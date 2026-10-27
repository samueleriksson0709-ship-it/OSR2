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

With --stress the bump is replaced by a real load case: the von Mises stress of an Ansys
export (element contributions averaged per node, converted to MPa with the unit declared
in loadcases.txt), e.g. the bracket solved with a downward force on its two top lugs.
The checkerboard is added on top exactly as above, centred on --center, which then does
not have to be a node; --amp is in MPa:
    python demo_filter_noise.py --stress [FILE] --center X Y Z [--patch 40 ...]
FILE is in SHARED_DIR, by default the first stress file of loadcases.txt. It must come
from the same solve as model.cdb.

Writes SHARED_DIR/demo_noise.csv for an External Data block:
  X,Y,Z        node coordinates, model unit
  Raw          bump + checkerboard (what goes into the filter)
  Filtered     output of the PDE filter
  Bump         the bump alone (the ideal result)                 Clean with --stress
  GradFilt     |grad Filtered| from MAPDL, stress / model length unit
  GradBump     |grad Bump|, exact                                GradClean with --stress:
               there is no exact gradient, so it is |grad| of the filtered clean field
"""
import argparse

import numpy as np

from Run_pipeline import (SHARED_DIR, INPUT_LOADCASES, PDEFilterAndGradient, read_cdb_mesh,
                          read_stress_file, parse_loadcases, von_mises_gradient, fail, log)
from validate_pde_gradient import EDGES

OUTPUT_CSV = "demo_noise.csv"


def fill_midside_mean(values, econn, lut):
    # MAPDL writes nodal TG, and the APDL block the stresses, at corner nodes only; give
    # each midside node the mean of the two corner nodes at the ends of its edge.
    nnode = (econn > 0).sum(1)
    for nn, edges in EDGES.items():
        rows = lut[econn[nnode == nn][:, :nn]]
        for a, b, m in edges:
            values[rows[:, m]] = 0.5 * (values[rows[:, a]] + values[rows[:, b]])
    return values


def clean_von_mises(fname, xyz, econn, lut, extent):
    # Von Mises stress in MPa of an Ansys export, one value per node of model.cdb.
    loadcases, *_, stress_unit, stress_scale, _ = parse_loadcases(SHARED_DIR / INPUT_LOADCASES)
    path = SHARED_DIR / (fname or loadcases[0][0])
    if not path.exists():
        fail(f"stress file not found: {path}")
    rows = read_stress_file(path)
    if not rows:
        fail(f"no stress rows in {path}")
    nid = np.array([r[1] for r in rows])
    if nid.max() >= len(lut) or (lut[nid] < 0).any() or \
            np.abs(np.array([r[2:5] for r in rows]) - xyz[lut[nid]]).max() > 1e-4 * extent:
        fail(f"{path.name} was not exported from the mesh in model.cdb; solve again so the "
             f"APDL block writes both")

    # NaN marks nodes without a value, so a node the export misses cannot pass unnoticed.
    k = lut[nid]
    acc = np.zeros((len(xyz), 6))
    np.add.at(acc, k, np.array([r[5:11] for r in rows]) * stress_scale)
    cnt = np.bincount(k, minlength=len(xyz))
    stress = np.full((len(xyz), 6), np.nan)
    stress[cnt > 0] = acc[cnt > 0] / cnt[cnt > 0, None]
    stress = fill_midside_mean(stress, econn, lut)
    if np.isnan(stress).any():
        fail(f"{np.isnan(stress[:, 0]).sum()} nodes of model.cdb have no stress in {path.name}")
    vm, _ = von_mises_gradient(stress, np.zeros(stress.shape + (3,)))
    log(f"  clean field: von Mises of {path.name} ({stress_unit} -> MPa), max {vm.max():.4g} MPa")
    return vm


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
                   help="centre of the bump (of the checkerboard with --stress), model length unit")
    p.add_argument("--amp-bump", type=float, default=100.0, help="bump height A")
    p.add_argument("--w", type=float, default=15.0, help="bump width (Gaussian sigma), mm")
    p.add_argument("--amp", type=float, default=20.0, help="checkerboard amplitude a (MPa with --stress)")
    p.add_argument("--lam", type=float, default=5.0, help="checkerboard wavelength, mm")
    p.add_argument("--patch", type=float, default=None, help="checkerboard patch radius, mm (default 2 w)")
    p.add_argument("--r", type=float, default=7.0, help="filter radius, mm")
    p.add_argument("--stress", nargs="?", const="", default=None, metavar="FILE",
                   help="von Mises stress of this Ansys export instead of the bump "
                        "(default FILE: first stress file of loadcases.txt)")
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
    if args.stress is None:
        k_near = int(np.argmin(np.linalg.norm(xyz - centre, axis=1)))
        if np.linalg.norm(xyz[k_near] - centre) > 0.05 * extent:
            fail(f"--center is {np.linalg.norm(xyz[k_near] - centre):.4g} {unit} from the nearest node; "
                 f"is it in the model unit ({unit})?")
        centre = xyz[k_near]
        log(f"Demo: model unit {unit}, bump centre at node {mesh_nids[k_near]} {np.round(centre, 6)}")
    else:
        # only the checkerboard is centred here, so a point in the air (e.g. on a pin axis) is fine
        if ((centre < xyz.min(axis=0) - 0.05 * extent) | (centre > xyz.max(axis=0) + 0.05 * extent)).any():
            fail(f"--center is outside the model; is it in the model unit ({unit})?")
        log(f"Demo: model unit {unit}, checkerboard centred at {np.round(centre, 6)}")

    h = local_element_size(xyz, econn, lut, centre, patch)
    log(f"  element size within {patch / mm:g} mm of the centre: {h / mm:.2f} mm (median edge)")
    if not h <= lam / 4:
        log(f"  WARNING: checkerboard wavelength {args.lam:g} mm < 4 elements; refine the mesh "
            f"to about {args.lam / 4:.2g} mm or increase --lam")
    if not h <= r / 3:
        log(f"  WARNING: filter radius {args.r:g} mm < 3 elements; refine the mesh "
            f"to about {args.r / 3:.2g} mm or increase --r")

    d = np.linalg.norm(xyz - centre, axis=1)
    k = 2 * np.pi / lam
    if args.stress is None:
        clean = args.amp_bump * np.exp(-d**2 / (2 * w**2))
    else:
        clean = clean_von_mises(args.stress, xyz, econn, lut, extent)
    window = 0.5 * (1 - np.tanh((d - patch) / (0.15 * patch)))      # soft edge, no jump
    noise = args.amp * np.sin(k * xyz[:, 0]) * np.sin(k * xyz[:, 1]) * np.sin(k * xyz[:, 2]) * window
    raw = clean + noise

    # Every mesh node gets its exact value, so fill_midside has nothing to fill.
    # SXX = clean + checkerboard, SYY = clean, SZZ = checkerboard, shear components unused.
    rows = [(0, int(nid), *xyz[i], raw[i], clean[i], noise[i], 0.0, 0.0, 0.0)
            for i, nid in enumerate(mesh_nids)]
    out, grad_by_nid = PDEFilterAndGradient(rows, args.r * mm)
    filt, filt_clean, filt_noise = np.array([row[5:8] for row in out]).T
    grad = fill_midside_mean(grad_by_nid[mesh_nids - 1], econn, lut)            # (n, 6, 3)
    grad_filt = np.linalg.norm(grad[:, 0, :], axis=1)                            # component SXX
    grad_filt_clean = np.linalg.norm(grad[:, 1, :], axis=1)                      # SYY = clean alone
    if args.stress is None:
        ref = "Bump"
        grad_ref = np.linalg.norm((xyz - centre) / w**2 * clean[:, None], axis=1)
    else:
        ref = "Clean"
        grad_ref = grad_filt_clean

    with open(SHARED_DIR / OUTPUT_CSV, "w") as f:
        f.write(f"X,Y,Z,Raw,Filtered,{ref},GradFilt,Grad{ref}\n")
        for (x, y, z), *v in zip(xyz, raw, filt, clean, grad_filt, grad_ref):
            f.write(f"{x:15.7E},{y:15.7E},{z:15.7E}," + ",".join(f"{a:15.7E}" for a in v) + "\n")

    l0 = r / (2 * np.sqrt(2))
    inner = d < 0.7 * patch
    left = np.abs(filt_noise[inner]).max() / np.abs(noise[inner]).max()
    k_tot = int(np.argmax(np.where(inner, grad_filt, 0.0)))
    k_noise = int(np.argmax(np.where(inner, np.abs(grad_filt - grad_filt_clean), 0.0)))
    if args.stress is None:
        log(f"  bump peak:          {args.amp_bump:g} -> {filt_clean.max():.4g} after filtering")
    else:
        peak = clean[inner].max()
        log(f"  clean peak:         {peak:.4g} MPa -> {filt_clean[inner].max():.4g} MPa after filtering, "
            f"within {0.7 * patch / mm:g} mm of the centre (checkerboard amplitude = {100 * args.amp / peak:.3g} % of it)")
    log(f"  checkerboard left:  {100 * left:.3g} % of its amplitude "
        f"(interior theory {100 / (1 + 3 * (l0 * k)**2):.3g} %)")
    if args.stress is None:
        ring = d < 3 * w
        k_max = int(np.argmax(np.where(ring, grad_filt_clean, 0.0)))
        log(f"  gradient ring:      max {grad_filt_clean[k_max] * mm:.4g} per mm at {d[k_max] / mm:.3g} mm from the centre "
            f"(exact for the raw bump: {args.amp_bump / (args.w * np.sqrt(np.e)):.4g} per mm at {args.w:g} mm)")
    log(f"  GradFilt column:    max {grad_filt[k_tot] * mm:.4g} per mm; leftover checkerboard adds up to "
        f"{abs(grad_filt[k_noise] - grad_filt_clean[k_noise]) * mm:.3g} per mm")
    log(f"  written {SHARED_DIR / OUTPUT_CSV}")


if __name__ == "__main__":
    main()
