"""
Run_pipeline plus the old Python-assembled PDE filter, kept for comparisons.

  PDEFilter_construct: the old filter, P1 (4-node tet) K and M assembled in Python
                       from the element rows of the stress export, solved with spsolve.
  PDEFilter_export2:   the current filter, Run_pipeline.PDEFilterAndGradient (mesh from
                       model.cdb, midside fill, one MAPDL run on pde_filter_in.txt), with
                       K and M exported by that same run.

Both return (filtered rows, K, M, xyz), with K, M and xyz in node-number order.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from Run_pipeline import (SHARED_DIR, PDE_DECK, PDEFilterAndGradient, fail, log,
                          read_cdb_mesh, read_mmf_vector)

PDE_DECK_KM = "pde_filter_km_in.txt"
KM_EXPORT_FILES = ("K_pde.mtx", "M_pde.mtx")


def cont3D4kmt(ex,ey,ez):
    C = np.array([
        [1 , ex[0], ey[0], ez[0]],
        [1 , ex[1], ey[1], ez[1]],
        [1 , ex[2], ey[2], ez[2]],
        [1 , ex[3], ey[3], ez[3]]
    ])
    Cinv = np.linalg.inv(C)
    V = np.abs(np.linalg.det(C)/6)

    B = Cinv[1:4,:]

    Ke = B.T @ B * V
    Me = V / 20 * (np.ones((4, 4)) + np.eye(4))

    return Ke,Me, V

def constructPDEMatrix(stress_file,r,nnod,nelm,npe):

    rows, cols, dataK, dataM = [], [], [], []

    #Construct matrices
    vols = []

    for element in range(nelm):
        index = element*npe

        nodes = [int(row[1])-1 for row in stress_file[index:index+npe]]

        ex = [row[2] for row in stress_file[index:index+npe]]
        ey = [row[3] for row in stress_file[index:index+npe]]
        ez = [row[4] for row in stress_file[index:index+npe]]

        Ke,Me,Ve = cont3D4kmt(ex,ey,ez)

        vols.append(Ve)

        for a in range(npe):
            for b in range(npe):
                rows.append(nodes[a]); cols.append(nodes[b])
                dataK.append(Ke[a, b]); dataM.append(Me[a, b])

    K = sparse.coo_matrix((dataK, (rows, cols)), shape=(nnod, nnod)).tocsr()
    M = sparse.coo_matrix((dataM, (rows, cols)), shape=(nnod, nnod)).tocsr()

    #Determine l_0
    l_0 = r / (2*np.sqrt(2))

    h = np.cbrt(6*np.sqrt(2) * np.array(vols))
    h10, h50, h90 = np.percentile(h, [10, 50, 90])
    kept = lambda he: 1 / (1 + (np.pi * l_0 / he)**2)
    log(f"  PDE filter: l_0={l_0:.4g}, edge p10/p50/p90 = {h10:.4g}/{h50:.4g}/{h90:.4g}, "
        f"node-to-node noise kept {kept(h10):.0%}/{kept(h50):.0%}/{kept(h90):.0%}")
    if kept(h10) > 0.9:
        log(f"  Warning: filter has almost no effect even in the finest 10% of elements; refine the mesh or raise r")
    # kept() only catches "the filter does nothing"; it says nothing about whether the
    # mesh resolves l_0. Below l_0 ~ h the solve is dominated by the discretisation, so
    # the P1 assembly here and the P2 matrices from MAPDL answer different questions.
    if l_0 < h50:
        log(f"  Warning: l_0={l_0:.4g} is below the median edge h={h50:.4g} (l_0/h={l_0/h50:.2f}); "
            f"the filtered field is mesh-dependent and the two implementations will not agree. "
            f"Raise r to >= {2*np.sqrt(2)*h50:.4g} or refine the mesh")

    return K,M,l_0

def PDEFilter_construct(stress_file_nodal,stress_file_elemental,r):
    #Find geometry variables
    npe = int(sum(1 for row in stress_file_elemental if row[0] == stress_file_elemental[0][0]))
    assert npe == 4, f"expected linear tets (4 rows per element), got {npe}"
    nelm = len(stress_file_elemental) // npe
    assert nelm * npe == len(stress_file_elemental), f"row count {len(stress_file_elemental)} not divisible into {nelm}×{npe}"
    eids = np.array([row[0] for row in stress_file_elemental]).reshape(nelm, npe)
    assert np.all(eids == eids[:, :1]), "element rows are not grouped in contiguous blocks"
    nnod = int(max(row[1] for row in stress_file_elemental))

    nids = np.array([row[1] for row in stress_file_nodal])
    assert np.array_equal(nids, np.arange(1, nnod + 1)), "nodal file is not one row per node, 1..nnod in order"
    assert len({row[1] for row in stress_file_elemental}) == nnod, "some node ids are unused by elements (A would be singular)"

    K,M,l_0 = constructPDEMatrix(stress_file_elemental,r,nnod,nelm,npe)

    stress_nodal = np.array([r[5:11] for r in stress_file_nodal])

    # solve once (nodal)
    A_pde = (l_0**2 * K + M).tocsc()
    stress_filt_nodal = spsolve(A_pde, M @ stress_nodal)

    out = []
    for row, s in zip(stress_file_nodal, stress_filt_nodal):
        out.append(row[:5] + tuple(s.tolist()))

    xyz = np.array([row[2:5] for row in stress_file_nodal])

    return out, K, M, xyz


def read_mmf_triplets(path):
    with open(path) as f:
        header = f.readline()
        symmetric = "symmetric" in header.lower()
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        nrow, ncol, nnz = (int(v) for v in line.split())
        data = np.loadtxt(f, max_rows=nnz)
    i = data[:, 0].astype(np.int64) - 1
    j = data[:, 1].astype(np.int64) - 1
    v = data[:, 2]
    if symmetric:
        off = i != j
        i, j, v = (np.concatenate([i, j[off]]),
                   np.concatenate([j, i[off]]),
                   np.concatenate([v, v[off]]))
    return i, j, v, nrow


def write_km_deck():
    # The pipeline deck plus two *EXPORTs, placed after M = Keff - Kt and before Kt is
    # overwritten by l_0^2 K + M, so K and M are exactly the matrices the pipeline solves with.
    if not (SHARED_DIR / PDE_DECK).exists():
        fail(f"{PDE_DECK} not found in {SHARED_DIR}")
    deck = (SHARED_DIR / PDE_DECK).read_text()
    anchor = "*AXPY,-1.0,,Kt,1.0,,Keff\n"
    if deck.count(anchor) != 1:
        fail(f"'{anchor.strip()}' not found exactly once in {PDE_DECK}")
    deck = deck.replace(anchor, anchor + "*EXPORT,Kt,MMF,K_pde.mtx\n*EXPORT,Keff,MMF,M_pde.mtx\n")
    (SHARED_DIR / PDE_DECK_KM).write_text(deck)


def PDEFilter_export2(stress_file_nodal,r):
    write_km_deck()
    for name in KM_EXPORT_FILES:
        p = SHARED_DIR / name
        if p.exists():
            p.unlink()

    out, _ = PDEFilterAndGradient(stress_file_nodal, r, deck=PDE_DECK_KM)
    for name in KM_EXPORT_FILES:
        if not (SHARED_DIR / name).exists():
            fail(f"{name} not written, see the MAPDL output of {PDE_DECK_KM}")

    # K and M come out in equation order; run_mapdl has checked that mapb.mtx is a
    # permutation of the mesh nodes, so this moves them to mesh_nids (node-number) order.
    mesh_nids, xyz, _ = read_cdb_mesh(SHARED_DIR / "model.cdb")
    n = len(mesh_nids)
    pos = np.searchsorted(mesh_nids, read_mmf_vector(SHARED_DIR / "mapb.mtx").astype(np.int64))
    mats = []
    for name in KM_EXPORT_FILES:
        i, j, v, nrow = read_mmf_triplets(SHARED_DIR / name)
        if nrow != n:
            fail(f"{name}: {nrow} rows, mesh has {n} nodes")
        mats.append(sparse.coo_matrix((v, (pos[i], pos[j])), shape=(n, n)).tocsr())
    K, M = mats

    return out, K, M, xyz
