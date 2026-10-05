"""Validate the stress gradient that pde_filter_onerun_in.txt writes (grad_onerun.txt).

Run after a PDE filter run, from the OSR_model folder:
    python validate_pde_gradient.py [folder]       (folder defaults to SHARED_DIR)

It reads what that run leaves in the folder (model.cdb, sol_onerun.txt, grad_onerun.txt,
pde_onerun_params.txt) and checks the MAPDL gradient against the filtered stress two ways:

1. Recompute (10-node tets only): differentiate the filtered field with the SOLID187
   shape functions and average over the elements at each corner node. MAPDL's nodal TG is
   the same quantity, so the two should agree to round-off; anything above ~1e-4 % means
   MAPDL averages or extrapolates differently, or the files are mapped wrongly.
2. Edge derivatives (any quadratic element): along each element edge, the derivative of
   the field at a corner follows exactly from the three edge nodes,
   df/ds = -3 f_a + 4 f_m - f_b, and must equal grad f . dx/ds. The MAPDL gradient is an
   average over elements, so the match is close rather than exact; a slope far from 1 or
   a low R^2 means wrong scaling, units, axis order or node mapping.
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


def check_recompute(econn, lut, xyz, sol, grad, corner):
    nnode = (econn > 0).sum(1)
    if not np.all(nnode == 10):
        print("  skipped: the mesh is not all 10-node tets "
              f"(node counts {sorted(set(nnode.tolist()))})")
        return
    n = len(sol)
    idx = lut[econn[:, :10]]                                        # (ne, 10)
    dn = np.stack([dN_tet10(*p) for p in TET10_NAT])                # (10 points, 3, 10)
    J = np.einsum('pdk,ekx->epdx', dn, xyz[idx])                    # (ne, 10, 3, 3)
    dfdxi = np.einsum('pdk,ekm->epdm', dn, sol[idx])                # (ne, 10, 3, 6)
    g_el = np.linalg.solve(J, dfdxi)                                # element gradient at its nodes
    acc = np.zeros((n, 3, 6))
    np.add.at(acc, idx.ravel(), g_el.reshape(-1, 3, 6))
    cnt = np.bincount(idx.ravel(), minlength=n)
    ref = acc / np.maximum(cnt, 1)[:, None, None]

    print("  differences in % of the largest reference gradient of each component")
    print("  component  max|g| (ref)    max diff [%]   99th pct diff [%]")
    for c, name in enumerate(COMPONENTS):
        r, m = ref[corner, :, c], grad[corner, :, c]
        scale = np.abs(r).max() or 1.0
        rel = np.linalg.norm(m - r, axis=1) / scale
        print(f"  {name:9s}  {scale:12.4e}   {100 * rel.max():12.3g}   {100 * np.percentile(rel, 99):12.3g}")
    worst = np.unravel_index(np.argmax(np.linalg.norm(grad - ref, axis=1)[corner]), (corner.sum(), 6))
    print(f"  largest difference at mesh row {np.flatnonzero(corner)[worst[0]]}, component {COMPONENTS[worst[1]]}")


def check_edges(econn, lut, xyz, sol, grad, corner):
    d_exact, d_mapdl = [], []
    for nn, edges in EDGES.items():
        rows = econn[(econn > 0).sum(1) == nn]
        if not len(rows):
            continue
        idx = lut[rows[:, :nn]]
        for a, b, m in edges:
            ia, ib, im = idx[:, a], idx[:, b], idx[:, m]
            for end, sign in ((ia, 1.0), (ib, -1.0)):
                # derivative w.r.t. the edge parameter s in [0, 1] at this end of the edge,
                # exact for the quadratic (isoparametric) field even on curved edges
                if sign > 0:
                    dxds = -3 * xyz[ia] + 4 * xyz[im] - xyz[ib]
                    dfds = -3 * sol[ia] + 4 * sol[im] - sol[ib]
                else:
                    dxds = xyz[ia] - 4 * xyz[im] + 3 * xyz[ib]
                    dfds = sol[ia] - 4 * sol[im] + 3 * sol[ib]
                length = np.linalg.norm(dxds, axis=1)
                ok = (length > 1e-12 * np.abs(xyz).max()) & corner[end]   # skip collapsed edges
                t = dxds[ok] / length[ok, None]
                d_exact.append(dfds[ok] / length[ok, None])                # (k, 6)
                d_mapdl.append(np.einsum('kd,kdc->kc', t, grad[end[ok]]))
    if not d_exact:
        print("  skipped: no quadratic elements (linear elements have no midside nodes)")
        return
    x, y = np.vstack(d_exact), np.vstack(d_mapdl)
    print(f"  {len(x)} edge ends;  ideal: slope deviation 0 %, R^2 1 (nodal averaging keeps them slightly off)")
    print("  differences in % of the largest exact edge derivative of each component")
    print("  component  slope dev [%]     R^2    median diff [%]  95th pct diff [%]")
    for c, name in enumerate(COMPONENTS):
        xc, yc = x[:, c], y[:, c]
        slope = (xc @ yc) / (xc @ xc) if xc @ xc > 0 else np.nan
        r2 = 1 - ((yc - xc) ** 2).sum() / max(((xc - xc.mean()) ** 2).sum(), 1e-300)
        diff = 100 * np.abs(yc - xc) / (np.abs(xc).max() or 1.0)
        print(f"  {name:9s}  {100 * (slope - 1):12.3g}  {r2:9.5f}  {np.median(diff):14.3g}"
              f"  {np.percentile(diff, 95):16.3g}")


def main(folder):
    mesh_nids, xyz, econn, lut, sol, grad, gsc = load_run(folder)
    corner = corner_mask(econn, lut, len(mesh_nids))
    zero = np.all(grad == 0, axis=(1, 2))
    print(f"{len(mesh_nids)} mesh nodes, {corner.sum()} corner nodes, gsc_ = {gsc:.6g}")
    print(f"zero gradient rows: {zero[corner].sum()} corner nodes, {zero[~corner].sum()} midside nodes")
    print("\n1) recompute from sol_onerun.txt with SOLID187 shape functions (corner nodes)")
    check_recompute(econn, lut, xyz, sol, grad, corner)
    print("\n2) exact edge derivatives vs MAPDL gradient projected on the edge")
    check_edges(econn, lut, xyz, sol, grad, corner)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else SHARED_DIR)
