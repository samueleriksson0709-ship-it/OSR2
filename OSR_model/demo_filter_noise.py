"""Visual check of the PDE filter and its gradient: a smooth bump plus a checkerboard.

The stresses are made up on the mesh of model.cdb, so the right answer is known:
    SXX = A exp(-d^2 / (2 w^2))                                  smooth bump, d = |x - centre|
        + a sin(k x) sin(k y) sin(k z) * patch(d),  k = 2 pi / lam   3D checkerboard
The filter treats the six components independently, so SYY and SZZ carry the bump alone
and the checkerboard alone through the same MAPDL run; this separates the two effects
exactly (the filter is linear). The filter must
  * remove the checkerboard: in the interior its amplitude is multiplied by
    1 / (1 + 3 (l0 k)^2), l0 = r / (2 sqrt 2), the transfer of (l0^2 K + M) u = M f
    for a field with one wavenumber k along each of the three axes;
  * leave the bump almost unchanged when r << w;
and the gradient of the bump is zero at the centre and peaks on a ring at d = w with
|grad SXX| = A / (w sqrt e).

Run from the OSR_model folder after model.cdb has been written by the APDL block:
    python demo_filter_noise.py --center X Y Z [--w 15 --lam 5 --amp 20 --r 7]
X Y Z is a point on the face where the bump goes, in the length unit of the model
(Mechanical: Hit Point Coordinate). w, lam, r and the patch radius are in mm.
The mesh must be several times finer than lam and r around the centre (see the warnings).

Writes SHARED_DIR/demo_noise.csv for an External Data block:
  X,Y,Z        node coordinates, model unit
  Raw          bump + checkerboard (what goes into the filter)
  Filtered     output of the PDE filter
  Bump         the bump alone (the ideal result)
  GradFilt     |grad Filtered| from MAPDL, stress / model length unit
  GradBump     |grad Bump|, exact
"""
import argparse

import numpy as np

from Run_pipeline import SHARED_DIR, PDEFilterAndGradient, read_cdb_mesh, fail, log

OUTPUT_CSV = "demo_noise.csv"


def local_element_size(xyz, econn, lut, centre, radius):
    corners = lut[econn[:, :4]]
    near = np.linalg.norm(xyz[corners].mean(axis=1) - centre, axis=1) < radius
    if not near.any():
        return np.nan
    c = xyz[corners[near]]
    edges = [np.linalg.norm(c[:, i] - c[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)]
    return float(np.median(np.mean(edges, axis=0)))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--center", nargs=3, type=float, required=True, metavar=("X", "Y", "Z"),
                   help="centre of the bump, model length unit")
    p.add_argument("--amp-bump", type=float, default=100.0, help="bump height A")
    p.add_argument("--w", type=float, default=15.0, help="bump width (Gaussian sigma), mm")
    p.add_argument("--amp", type=float, default=20.0, help="checkerboard amplitude a")
    p.add_argument("--lam", type=float, default=5.0, help="checkerboard wavelength, mm")
    p.add_argument("--patch", type=float, default=None, help="checkerboard patch radius, mm (default 2 w)")
    p.add_argument("--r", type=float, default=7.0, help="filter radius, mm")
    args = p.parse_args()

    mesh_nids, xyz, econn = read_cdb_mesh(SHARED_DIR / "model.cdb")
    lut = np.full(mesh_nids.max() + 1, -1, dtype=np.int64)
    lut[mesh_nids] = np.arange(len(mesh_nids))

    # model.cdb is in the unit the model was solved in: m for MKS, mm for the mm-kg-N system
    extent = float(np.ptp(xyz, axis=0).max())
    mm = 1e-3 if extent < 5.0 else 1.0
    unit = "m" if mm == 1e-3 else "mm"
    w, lam, r = args.w * mm, args.lam * mm, args.r * mm
    patch = (args.patch if args.patch is not None else 2 * args.w) * mm

    centre = np.array(args.center)
    k_near = int(np.argmin(np.linalg.norm(xyz - centre, axis=1)))
    if np.linalg.norm(xyz[k_near] - centre) > 0.05 * extent:
        fail(f"--center is {np.linalg.norm(xyz[k_near] - centre):.4g} {unit} from the nearest node; "
             f"is it in the model unit ({unit})?")
    centre = xyz[k_near]

    h = local_element_size(xyz, econn, lut, centre, patch)
    log(f"Demo: model unit {unit}, bump centre at node {mesh_nids[k_near]} {np.round(centre, 6)}")
    log(f"  element size within {patch / mm:g} mm of the centre: {h / mm:.2f} mm (median edge)")
    if not h <= lam / 4:
        log(f"  WARNING: checkerboard wavelength {args.lam:g} mm < 4 elements; refine the mesh "
            f"to about {args.lam / 4:.2g} mm or increase --lam")
    if not h <= r / 3:
        log(f"  WARNING: filter radius {args.r:g} mm < 3 elements; refine the mesh "
            f"to about {args.r / 3:.2g} mm or increase --r")

    d = np.linalg.norm(xyz - centre, axis=1)
    k = 2 * np.pi / lam
    bump = args.amp_bump * np.exp(-d**2 / (2 * w**2))
    window = 0.5 * (1 - np.tanh((d - patch) / (0.15 * patch)))      # soft edge, no jump
    noise = args.amp * np.sin(k * xyz[:, 0]) * np.sin(k * xyz[:, 1]) * np.sin(k * xyz[:, 2]) * window
    raw = bump + noise

    # Every mesh node gets its exact value, so fill_midside has nothing to fill.
    # SXX = bump + checkerboard, SYY = bump, SZZ = checkerboard, shear components unused.
    rows = [(0, int(nid), *xyz[i], raw[i], bump[i], noise[i], 0.0, 0.0, 0.0)
            for i, nid in enumerate(mesh_nids)]
    out, grad_by_nid = PDEFilterAndGradient(rows, args.r * mm)
    filt, filt_bump, filt_noise = np.array([row[5:8] for row in out]).T
    grad_filt = np.linalg.norm(grad_by_nid[mesh_nids - 1, 0, :], axis=1)        # component SXX
    grad_bump = np.linalg.norm((xyz - centre) / w**2 * bump[:, None], axis=1)

    with open(SHARED_DIR / OUTPUT_CSV, "w") as f:
        f.write("X,Y,Z,Raw,Filtered,Bump,GradFilt,GradBump\n")
        for (x, y, z), *v in zip(xyz, raw, filt, bump, grad_filt, grad_bump):
            f.write(f"{x:15.7E},{y:15.7E},{z:15.7E}," + ",".join(f"{a:15.7E}" for a in v) + "\n")

    l0 = r / (2 * np.sqrt(2))
    inner = d < 0.7 * patch
    left = np.abs(filt_noise[inner]).max() / np.abs(noise[inner]).max()
    ring = d < 3 * w
    grad_filt_bump = np.linalg.norm(grad_by_nid[mesh_nids - 1, 1, :], axis=1)  # SYY = bump alone
    k_max = int(np.argmax(np.where(ring, grad_filt_bump, 0.0)))
    k_tot = int(np.argmax(np.where(inner, grad_filt, 0.0)))
    k_noise = int(np.argmax(np.where(inner, np.abs(grad_filt - grad_filt_bump), 0.0)))
    log(f"  bump peak:          {args.amp_bump:g} -> {filt_bump.max():.4g} after filtering")
    log(f"  checkerboard left:  {100 * left:.3g} % of its amplitude "
        f"(interior theory {100 / (1 + 3 * (l0 * k)**2):.3g} %)")
    log(f"  gradient ring:      max {grad_filt_bump[k_max] * mm:.4g} per mm at {d[k_max] / mm:.3g} mm from the centre "
        f"(exact for the raw bump: {args.amp_bump / (args.w * np.sqrt(np.e)):.4g} per mm at {args.w:g} mm)")
    log(f"  GradFilt column:    max {grad_filt[k_tot] * mm:.4g} per mm; leftover checkerboard adds up to "
        f"{abs(grad_filt[k_noise] - grad_filt_bump[k_noise]) * mm:.3g} per mm")
    log(f"  written {SHARED_DIR / OUTPUT_CSV}")


if __name__ == "__main__":
    main()
