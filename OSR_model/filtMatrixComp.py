"""Compare the K and M matrices of the old and the current PDE filter.

  construct (Run_pipelineWithOldFunc.PDEFilter_construct): P1 K, M assembled in Python from
            the element rows of the stress export, on the corner nodes 1..nnod
  export2   (Run_pipelineWithOldFunc.PDEFilter_export2):   K, M exported by the same MAPDL run as
            Run_pipeline.PDEFilterAndGradient, P2 (10-node tets of model.cdb), all nodes

The P2 matrices are restricted to the P1 space with the prolongation P (corner node -> itself,
midside node -> mean of its edge end nodes, the rule fill_midside uses). For straight-sided
elements P^T K_P2 P = K_P1 and P^T M_P2 P = M_P1 exactly, so the entries are compared one to
one. What is left comes from curved boundary elements (midside node off the edge midpoint)
and from the rounded coordinates of the stress export, so the difference is also split into
entries that only straight elements contribute to and entries a curved element contributes to.

K and M do not depend on r; r only enters the filtered stresses compared at the end.

Usage:  python filtMatrixComp.py [--stress-file FILE] [--r 0.003]
"""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import norm as sparse_norm

import Run_pipelineWithOldFunc as old
from Run_pipeline import SHARED_DIR, read_cdb_mesh, read_stress_file
from filtExtractComp import compare, load_nodal, pct, print_timing, run_timed

# SOLID187 / SOLID87 node order: corners I J K L, then the midside nodes of these edges.
TET10_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
CURVED_TOL = 1e-6   # midside offset from the edge midpoint, relative to the edge length


def prolongation(n1):
    # P (n2 x n1): P1 field on nodes 1..n1 (column = node id - 1) -> P2 field on the mesh of
    # model.cdb (row = position in mesh_nids, the order of export2's K and M). Also returns the
    # P2 row of each P1 node, C (n1 x n1, nonzero where a curved element contains both nodes)
    # and the midside offsets (n_elements x 6).
    mesh_nids, xyz, econn = read_cdb_mesh(SHARED_DIR / "model.cdb")
    if (econn[:, :10] <= 0).any() or (econn[:, 10:] != 0).any():
        raise SystemExit("model.cdb: only 10-node tets are supported")
    corners, mids = econn[:, :4], econn[:, 4:10]
    corner_ids = np.unique(corners)
    if not np.array_equal(corner_ids, np.arange(1, n1 + 1)):
        raise SystemExit(f"corner nodes of model.cdb ({len(corner_ids)}, max id {corner_ids.max()}) "
                         f"are not the nodes 1..{n1} of the stress export; different mesh?")

    lut = np.full(mesh_nids.max() + 1, -1, dtype=np.int64)
    lut[mesh_nids] = np.arange(len(mesh_nids))
    a = corners[:, [i for i, _ in TET10_EDGES]]
    b = corners[:, [j for _, j in TET10_EDGES]]
    rows = np.concatenate([lut[corner_ids], lut[mids].ravel(), lut[mids].ravel()])
    cols = np.concatenate([corner_ids - 1, a.ravel() - 1, b.ravel() - 1])
    vals = np.ones(len(rows))
    P = sparse.coo_matrix((vals, (rows, cols)), shape=(len(mesh_nids), n1)).tocsr()
    # A midside node shared by k elements got its two end nodes k times each.
    P = (sparse.diags(1.0 / np.asarray(P.sum(axis=1)).ravel()) @ P).tocsr()
    nnz_row = np.diff(P.indptr)
    if (nnz_row[lut[corner_ids]] != 1).any() or nnz_row.max() > 2:
        raise SystemExit("model.cdb: a node is both a corner and a midside node, or the "
                         "elements sharing a midside node disagree on its edge")

    xm, xa, xb = xyz[lut[mids]], xyz[lut[a]], xyz[lut[b]]
    offset = np.linalg.norm(xm - 0.5 * (xa + xb), axis=2) / np.linalg.norm(xa - xb, axis=2)
    curved = offset.max(axis=1) > CURVED_TOL
    cn = corners[curved] - 1
    S = sparse.coo_matrix((np.ones(cn.size), (cn.ravel(), np.repeat(np.arange(len(cn)), 4))),
                          shape=(n1, len(cn))).tocsr()
    return P, lut[corner_ids], (S @ S.T).tocsr(), offset


def matrix_checks(K, M, xyz):
    # From the old Run_pipeline.filter_comp: quantities both discretisations should share.
    ones = np.ones(K.shape[0])
    Mx = M @ ones
    vol = ones @ Mx
    out = {"n_nodes": K.shape[0], "volume": vol,
           "max|K 1| / max|K|": np.abs(K @ ones).max() / abs(K).max()}
    for d, name in enumerate("xyz"):
        u = xyz[:, d] - (Mx @ xyz[:, d]) / vol
        out[f"u'Ku/V ({name})"] = u @ (K @ u) / vol
        out[f"u'Mu ({name})"] = u @ (M @ u)
    return out


def pattern_diff(A, B):
    # Entries nonzero in B but not in A, and the other way round.
    d = ((B != 0).astype(np.int8) - (A != 0).astype(np.int8)).tocsr()
    d.eliminate_zeros()
    return int((d.data > 0).sum()), int((d.data < 0).sum())


def compare_matrix(name, A, B, C):
    # A: construct (P1), B: export2 restricted to P1, C: entries a curved element touches.
    extra, missing = pattern_diff(A, B)
    print(f"\n{name}: P^T {name}_P2 P (export2) vs {name}_P1 (construct)")
    print(f"  nonzeros {A.nnz} vs {B.nnz}, {extra} only in export2, {missing} only in construct")
    curved = (C != 0).astype(float)
    D = (B - A).tocsr()
    parts = (("all", A, D),
             ("straight elements only", A - A.multiply(curved), D - D.multiply(curved)),
             ("with a curved element", A.multiply(curved), D.multiply(curved)))
    print(f"  {'entries':<24s} {'count':>8s} {'rel L2 %':>10s} {'max|d| % of max|A|':>20s}")
    for label, a, d in parts:
        print(f"  {label:<24s} {np.count_nonzero(a.data):8d} {pct(sparse_norm(d), sparse_norm(a)):10.4f} "
              f"{pct(abs(d).max(), abs(a).max()):20.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stress-file", default=str(SHARED_DIR / "sigma_export_nodes_testBracket.txt"))
    ap.add_argument("--r", type=float, default=0.003, help="filter radius; l_0 = r/(2*sqrt(2))")
    args = ap.parse_args()

    stress_file_elemental = read_stress_file(args.stress_file)
    stress_file_nodal = load_nodal(args.stress_file)
    nids = np.array([row[1] for row in stress_file_nodal])
    print(f"Stress file: {args.stress_file}  ({len(stress_file_elemental)} rows, {len(nids)} unique nodes)")
    print(f"r = {args.r:g}, l_0 = {args.r / (2 * np.sqrt(2)):.6g}")

    methods = (
        ("construct", "P1, Python assembly", old, ("constructPDEMatrix", "spsolve"),
         lambda: old.PDEFilter_construct(stress_file_nodal, stress_file_elemental, args.r)),
        ("export2", "P2, MAPDL", old, ("read_cdb_mesh", "run_mapdl_onerun"),
         lambda: old.PDEFilter_export2(stress_file_nodal, args.r)),
    )
    results = []
    for name, desc, module, stages, run in methods:
        out, runs = run_timed(run, module, stages, 1)
        print_timing(f"{name} ({desc})", runs, stages)
        results.append(out)
    (out1, K1, M1, xyz1), (out2, K2, M2, xyz2) = results

    P, corner_rows, C, offset = prolongation(K1.shape[0])
    if P.shape[0] != K2.shape[0]:
        raise SystemExit(f"export2 matrices are {K2.shape[0]} x {K2.shape[0]}, model.cdb has {P.shape[0]} nodes")
    n_curved = (offset.max(axis=1) > CURVED_TOL).sum()
    print(f"\nMesh: {K1.shape[0]} corner nodes (P1), {K2.shape[0]} nodes (P2), {len(offset)} elements")
    print(f"  curved elements (midside offset > {CURVED_TOL:g} of the edge length): {n_curved}, "
          f"largest offset {offset.max():.4f}")
    print(f"  max |xyz stress export - xyz model.cdb| on the corner nodes: "
          f"{np.abs(xyz1 - xyz2[corner_rows]).max():.3e}")

    compare_matrix("K", K1, (P.T @ K2 @ P).tocsr(), C)
    compare_matrix("M", M1, (P.T @ M2 @ P).tocsr(), C)

    a, b = matrix_checks(K1, M1, xyz1), matrix_checks(K2, M2, xyz2)
    print("\nFull matrices, no restriction (old filter_comp checks)")
    print(f"  {'':<20s} {'construct (P1)':>16s} {'export2 (P2)':>16s} {'rel diff %':>11s}")
    for k in a:
        # Node counts differ by design and K 1 is ~0 for both, so no % for those two.
        rel = "" if k in ("n_nodes", "max|K 1| / max|K|") else f"{pct(abs(b[k] - a[k]), abs(a[k])):11.4f}"
        print(f"  {k:<20s} {a[k]:16.6g} {b[k]:16.6g} {rel:>11s}")

    s1 = np.array([row[5:11] for row in out1], dtype=float)
    s2 = np.array([row[5:11] for row in out2], dtype=float)
    compare(s1, s2, nids, "construct", "export2")


if __name__ == "__main__":
    main()
