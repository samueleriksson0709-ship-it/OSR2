"""Validate the stress gradient that pde_filter_onerun_in.txt writes (grad_onerun.txt).

Run after a PDE filter run, from the OSR_model folder:
    python validate_pde_gradient.py [folder]       (folder defaults to SHARED_DIR)

It reads what that run leaves in the folder (model.cdb, sol_onerun.txt, grad_onerun.txt,
pde_onerun_params.txt) and checks the MAPDL gradient against the filtered stress two ways:

1. Recompute (10-node tets only): differentiate the filtered field the way MAPDL does
   with ERESX,YES -- SOLID87 gradient at its 4 Gauss points, extrapolated linearly to the
   corner nodes, averaged over the elements at each node. MAPDL's nodal TG is that same
   quantity, so the two should agree to round-off; anything above ~1e-4 % means a
   different extrapolation setting (ERESX,NO copies Gauss point values instead) or wrongly
   mapped files. For information it also shows how far that extrapolated gradient is from
   the gradient evaluated exactly at the node: on straight-edged elements the two are equal,
   on curved ones (fillets) they differ, which is a property of the method, not an error.
2. Edge gradient (any quadratic element, no shape functions): along an element edge the
   field and the coordinates are quadratic in the edge parameter s, so at a corner
   df/ds = -3 f_a + 4 f_m - f_b and dx/ds likewise follow exactly from the three edge
   nodes. The three edges meeting at a corner give dx/ds_k . grad f = df/ds_k, k = 1..3,
   which fixes the element gradient at that corner exactly; averaged over the elements at
   the node it must match MAPDL on nodes that touch only straight-edged elements (where
   extrapolation is exact). Anything above ~1e-4 % there means wrong scaling, units, axis
   order or node mapping.
"""
import sys
from pathlib import Path

import numpy as np

from Run_pipeline import SHARED_DIR
from PDEFilterTest3 import read_cdb_mesh

COMPONENTS = ("SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ")

# (corner a, corner b, midside m) per edge, in Ansys node order
EDGES = {
    10: [(0, 1, 4), (1, 2, 5), (2, 0, 6), (0, 3, 7), (1, 3, 8), (2, 3, 9)],
    20: [(0, 1, 8), (1, 2, 9), (2, 3, 10), (3, 0, 11), (4, 5, 12), (5, 6, 13),
         (6, 7, 14), (7, 4, 15), (0, 4, 16), (1, 5, 17), (2, 6, 18), (3, 7, 19)],
}
N_CORNERS = {4: 4, 8: 8, 10: 4, 20: 8}

# SOLID187 natural coordinates (xi, eta, zeta) of the 10 nodes
TET10_NAT = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                      [.5, 0, 0], [.5, .5, 0], [0, .5, 0], [0, 0, .5], [.5, 0, .5], [0, .5, .5]])


def dN_tet10(xi, eta, zeta):
    L = np.array([1 - xi - eta - zeta, xi, eta, zeta])
    dL = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    dN = np.zeros((10, 3))
    for i in range(4):
        dN[i] = (4 * L[i] - 1) * dL[i]
    for k, (i, j, _) in enumerate(EDGES[10]):
        dN[4 + k] = 4 * (L[i] * dL[j] + L[j] * dL[i])
    return dN.T                                                     # 3 x 10


def read_params(path):
    params = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            params[key.strip().lower()] = float(value)
    return params


def load_run(folder):
    """Mesh, filtered stress and MAPDL gradient, all in mesh_nids order."""
    mesh_nids, xyz, econn = read_cdb_mesh(folder / "model.cdb")
    n = len(mesh_nids)
    lut = np.full(mesh_nids.max() + 1, -1, dtype=np.int64)
    lut[mesh_nids] = np.arange(n)

    data = np.loadtxt(folder / "sol_onerun.txt").reshape(-1, 7)
    sol = np.zeros((n, 6))
    sol[lut[data[:, 0].astype(np.int64)]] = data[:, 1:]

    # grad_onerun.txt holds the gradient of sol * gsc_ (MAPDL refuses D values above 1e6)
    gsc = read_params(folder / "pde_onerun_params.txt").get("gsc_", 1.0)
    g = np.loadtxt(folder / "grad_onerun.txt").reshape(-1, 19)
    rows = g[:, 0].astype(np.int64)                                 # node numbers 1..nrow_
    in_mesh = rows <= mesh_nids.max()
    in_mesh[in_mesh] = lut[rows[in_mesh]] >= 0
    grad = np.zeros((n, 3, 6))                                      # [mesh row, x/y/z, component]
    grad[lut[rows[in_mesh]]] = g[in_mesh, 1:].reshape(-1, 6, 3).transpose(0, 2, 1) / gsc
    return mesh_nids, xyz, econn, lut, sol, grad, gsc


def corner_mask(econn, lut, n):
    corner = np.zeros(n, dtype=bool)
    for nn, nc in N_CORNERS.items():
        rows = econn[(econn > 0).sum(1) == nn]
        if len(rows):
            corner[lut[rows[:, :nc]].ravel()] = True
    return corner


def curved_elements(econn, lut, xyz, tol=1e-6):
    """True for quadratic elements with a midside node off the chord of its edge."""
    curved = np.zeros(len(econn), dtype=bool)
    nnode = (econn > 0).sum(1)
    for nn, edges in EDGES.items():
        sel = np.flatnonzero(nnode == nn)
        if not len(sel):
            continue
        X = xyz[lut[econn[sel, :nn]]]
        for a, b, m in edges:
            chord = np.linalg.norm(X[:, a] - X[:, b], axis=1)
            off = np.linalg.norm(X[:, m] - 0.5 * (X[:, a] + X[:, b]), axis=1)
            curved[sel] |= off > tol * np.maximum(chord, 1e-300)
    return curved


def print_diff_table(ref, grad, nodes):
    print("  differences in % of the largest reference gradient of each component")
    print("  component  max|g| (ref)    max diff [%]   99th pct diff [%]")
    for c, name in enumerate(COMPONENTS):
        r, m = ref[nodes, :, c], grad[nodes, :, c]
        scale = np.abs(r).max() or 1.0
        rel = np.linalg.norm(m - r, axis=1) / scale
        print(f"  {name:9s}  {scale:12.4e}   {100 * rel.max():12.3g}   {100 * np.percentile(rel, 99):12.3g}")


def check_recompute(econn, lut, xyz, sol, grad, corner, curved):
    nnode = (econn > 0).sum(1)
    if not np.all(nnode == 10):
        print("  skipped: the mesh is not all 10-node tets "
              f"(node counts {sorted(set(nnode.tolist()))})")
        return
    n = len(sol)
    idx = lut[econn[:, :10]]                                        # (ne, 10)
    cidx = idx[:, :4].ravel()
    cnt = np.maximum(np.bincount(cidx, minlength=n), 1)[:, None, None]

    def average(g_corners):                                         # (ne, 4, 3, 6) -> (n, 3, 6)
        acc = np.zeros((n, 3, 6))
        np.add.at(acc, cidx, g_corners.reshape(-1, 3, 6))
        return acc / cnt

    def gradient_at(points):                                        # natural coords -> (ne, p, 3, 6)
        dn = np.stack([dN_tet10(*p) for p in points])
        J = np.einsum('pdk,ekx->epdx', dn, xyz[idx])
        return np.linalg.solve(J, np.einsum('pdk,ekm->epdm', dn, sol[idx]))

    # MAPDL with ERESX,YES: 4-point Gauss rule, linear extrapolation to the corners
    a, b = 0.5854101966249685, 0.1381966011250105
    P = np.full((4, 4), b) + (a - b) * np.eye(4)                    # P[j, i] = L_i at Gauss point j
    ref = average(np.einsum('ij,ejdm->eidm', np.linalg.inv(P), gradient_at(P[:, 1:])))
    at_node = average(gradient_at(TET10_NAT[:4]))

    print("  a) MAPDL gradient vs the same calculation in Python (should be round-off)")
    print_diff_table(ref, grad, corner)
    worst = np.unravel_index(np.argmax(np.linalg.norm(grad - ref, axis=1)[corner]), (corner.sum(), 6))
    print(f"  largest difference at mesh row {np.flatnonzero(corner)[worst[0]]}, component {COMPONENTS[worst[1]]}")

    on_curved = np.zeros(n, dtype=bool)
    on_curved[lut[econn[curved, :4]].ravel()] = True
    print(f"\n  b) for information: extrapolated vs evaluated at the node, on the {(corner & on_curved).sum()} "
          f"of {corner.sum()} corner nodes touching curved elements")
    print("     (zero on straight-edged elements; a method difference, not an error)")
    if (corner & on_curved).any():
        print_diff_table(at_node, ref, corner & on_curved)


def edge_gradient(econn, lut, xyz, sol):
    """Nodal gradient from the edges at each element corner, averaged over the elements."""
    n = len(sol)
    acc, cnt = np.zeros((n, 3, 6)), np.zeros(n)
    nnode = (econn > 0).sum(1)
    for nn, edges in EDGES.items():
        rows = econn[nnode == nn]
        if not len(rows):
            continue
        idx = lut[rows[:, :nn]]
        for c in range(N_CORNERS[nn]):
            # edges at corner c as (this corner, far corner, midside)
            inc = [(a, b, m) if a == c else (b, a, m) for a, b, m in edges if c in (a, b)]
            T = np.stack([-3 * xyz[idx[:, a]] + 4 * xyz[idx[:, m]] - xyz[idx[:, b]] for a, b, m in inc], 1)
            D = np.stack([-3 * sol[idx[:, a]] + 4 * sol[idx[:, m]] - sol[idx[:, b]] for a, b, m in inc], 1)
            # collapsed corners of degenerate elements have a zero-length edge: skip them
            lengths = np.linalg.norm(T, axis=2)
            ok = (lengths.min(1) > 1e-8 * lengths.max(1)) & \
                 (np.abs(np.linalg.det(T)) > 1e-8 * np.prod(lengths, axis=1))
            np.add.at(acc, idx[ok, c], np.linalg.solve(T[ok], D[ok]))   # T grad = D
            np.add.at(cnt, idx[ok, c], 1)
    return acc / np.maximum(cnt, 1)[:, None, None], cnt > 0


def check_edges(econn, lut, xyz, sol, grad, corner, curved):
    if not np.isin((econn > 0).sum(1), list(EDGES)).any():
        print("  skipped: no quadratic elements (linear elements have no midside nodes)")
        return
    ref, has = edge_gradient(econn, lut, xyz, sol)
    on_curved = np.zeros(len(sol), dtype=bool)
    on_curved[lut[econn[curved][econn[curved] > 0]]] = True
    nodes = corner & has & ~on_curved
    print(f"  {nodes.sum()} corner nodes touching only straight-edged elements "
          f"({(corner & on_curved).sum()} on curved elements left out)")
    if nodes.any():
        print_diff_table(ref, grad, nodes)


def main(folder):
    mesh_nids, xyz, econn, lut, sol, grad, gsc = load_run(folder)
    corner = corner_mask(econn, lut, len(mesh_nids))
    zero = np.all(grad == 0, axis=(1, 2))
    print(f"{len(mesh_nids)} mesh nodes, {corner.sum()} corner nodes, gsc_ = {gsc:.6g}")
    print(f"zero gradient rows: {zero[corner].sum()} corner nodes, {zero[~corner].sum()} midside nodes")
    curved = curved_elements(econn, lut, xyz)
    print("\n1) recompute from sol_onerun.txt with SOLID87 Gauss points + extrapolation (corner nodes)")
    check_recompute(econn, lut, xyz, sol, grad, corner, curved)
    print("\n2) gradient from exact edge derivatives at each element corner (no shape functions)")
    check_edges(econn, lut, xyz, sol, grad, corner, curved)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else SHARED_DIR)
