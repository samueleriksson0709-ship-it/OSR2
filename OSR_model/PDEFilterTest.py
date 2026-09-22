from Run_pipeline import SHARED_DIR, log, fail, read_stress_file
from pathlib import Path
import subprocess
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


MAPDL_EXE_DIR = r"C:\Program Files\ANSYS Inc\v261\ansys\bin\winx64\ANSYS261.exe"


def run_mapdl_extract():
    export_files = ("K_pde.mtx", "M_pde.mtx", "mapb.mtx", "mapb_t.mtx",
             "pde_dims.txt", "econn.txt", "nxyz.txt")
    pde_deck = "pde_matrices_extract_in.txt"
    for name in export_files:
        p = SHARED_DIR / name
        if p.exists():
            p.unlink()

    if not (SHARED_DIR / pde_deck).exists():
        fail(f"{pde_deck} not found in {SHARED_DIR}")
    if not (SHARED_DIR / "model.cdb").exists():
        fail("model.cdb not found; add CDWRITE to the Mechanical Commands object")
    if not Path(MAPDL_EXE_DIR).exists():
        fail(f"MAPDL executable not found: {MAPDL_EXE_DIR}")

    res = subprocess.run(
        [MAPDL_EXE_DIR, "-b", "-np", "1", "-i", pde_deck, "-o", "pde_filter.out"],
        cwd=str(SHARED_DIR), timeout=7200,
    )
    for name in export_files:
        if not (SHARED_DIR / name).exists():
            fail(f"{name} not written (code {res.returncode}), see pde_filter.out")

    return int(float(np.loadtxt(SHARED_DIR / "pde_dims.txt")))


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


def build_csc(i, j, v, n):
    return coo_matrix((v, (i, j)), shape=(n, n)).tocsc()


def read_mmf_vector(path):
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("%")]
    data = np.loadtxt(lines[1:])
    return data if data.ndim == 1 else data[:, -1]
    
def fill_midside(rhs, known, econn, lut, xyz, tol=0.25):
    # Fill midside nodes by averaging the values of the two end nodes of each edge.
    acc = np.zeros_like(rhs)
    cnt = np.zeros(len(rhs))
    for conn in econn:
        nodes = np.unique(lut[conn[conn > 0]])
        kn = nodes[known[nodes]]
        un = nodes[~known[nodes]]
        if len(un) == 0 or len(kn) < 2:
            continue
        ia, ib = np.triu_indices(len(kn), 1)
        a, b = kn[ia], kn[ib]
        mid = 0.5 * (xyz[a] + xyz[b])
        L = np.linalg.norm(xyz[a] - xyz[b], axis=1)
        d = np.linalg.norm(xyz[un, None, :] - mid[None], axis=2) / L
        k = d.argmin(axis=1)
        ok = d[np.arange(len(un)), k] < tol
        for m, p in zip(un[ok], k[ok]):
            acc[m] += 0.5 * (rhs[a[p]] + rhs[b[p]])
            cnt[m] += 1
    filled = cnt > 0
    out = rhs.copy()
    out[filled] = acc[filled] / cnt[filled, None]
    return out, known | filled

if __name__ == "__main__":
    nn = run_mapdl_extract()
    ki, kj, kv, n = read_mmf_triplets(SHARED_DIR / "K_pde.mtx")
    mi, mj, mv, nm = read_mmf_triplets(SHARED_DIR / "M_pde.mtx")
    if n != nn or nm != nn:
        fail(f"matrix size mismatch: K={n} M={nm} mesh nodes={nn}")
    log(f"K: {n} rows, {len(kv)} nonzeros")
    log(f"M: {n} rows, {len(mv)} nonzeros")

    krow = np.zeros(n)
    np.add.at(krow, ki, kv)
    log(f"nullspace check |K@1|_max / |K|_max = {np.abs(krow).max() / np.abs(kv).max():.3e}")
    log(f"M.sum() = {mv.sum():.6g}  (should equal mesh volume 4.63e-4)")

    back = read_mmf_vector(SHARED_DIR / "mapb.mtx").astype(np.int64)
    back_t = read_mmf_vector(SHARED_DIR / "mapb_t.mtx").astype(np.int64)
    if not np.array_equal(back, back_t):
        fail("equation ordering differs between pdesteady.full and pdetrans.full")
    if len(back) != nn:
        fail(f"mapping length {len(back)} != node count {nn}")

    lut = np.full(back.max() + 1, -1, dtype=np.int64)
    lut[back] = np.arange(n)

    # PDE Filter test

    path = "C:/Users/samue/OSR/OSR_model/sigma_export_nodes_testBracket.txt"
    stress_file_nodal = read_stress_file(path)
    l_0 = 0.003

    nids = np.array([int(r[1]) for r in stress_file_nodal])
    stress_nodal = np.array([r[5:11] for r in stress_file_nodal], dtype=float)
    if nids.max() >= len(lut) or (lut[nids] < 0).any():
        fail("stress-file nodes not in PDE mesh")
    rows = lut[nids]

    rhs = np.zeros((n, 6))
    rhs[rows] = stress_nodal
    known = np.zeros(n, dtype=bool)
    known[rows] = True
    log(f"stress file covers {known.sum()} of {n} nodes before midside fill")

    econn = np.loadtxt(SHARED_DIR / "econn.txt").astype(np.int64).reshape(-1, 20)
    nxyz = np.loadtxt(SHARED_DIR / "nxyz.txt").reshape(-1, 3)
    xyz = nxyz[back - 1]

    rhs, known = fill_midside(rhs, known, econn, lut, xyz)
    if not known.all():
        fail(f"{(~known).sum()} nodes still without stress after midside fill")

    K = build_csc(ki, kj, kv, n)
    M = build_csc(mi, mj, mv, n)
    lu = splu((l_0**2 * K + M).tocsc())
    sol = lu.solve(M @ rhs)
    stress_filt_nodal = sol[rows]
    stress_filt_all = sol

    out_path = SHARED_DIR / "testStress.txt"
    meta = np.array([r[0:5] for r in stress_file_nodal], dtype=float)
    out = np.column_stack([meta, stress_filt_nodal])
    np.savetxt(
        out_path, out,
        fmt=["%8d", "%8d"] + ["%14.6E"] * 9,
        header="eid nid X Y Z SXX SYY SZZ SXY SYZ SXZ",
        comments="# ",
    )
    log(f"filtered stress written to {out_path}")
