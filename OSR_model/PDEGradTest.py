"""PDE filter in a single MAPDL run.

Same maths as PDEFilterTest1/2, but the mesh (coordinates and connectivity for the
midside fill) is read straight from model.cdb, so Python can prepare the stresses
before MAPDL starts. They are written by node number; pde_filter_onerun_in.txt
builds the matrices, reorders the stresses into equation order, solves and writes
the result with node numbers.
"""
import numpy as np
from Run_pipeline import read_stress_file, SHARED_DIR, Path, MAPDL_EXE_DIR, subprocess, fail, log, \
    parse_loadcases, write_gradient_csv, gradient_csv_path

# Element types the decks delete before building the PDE matrices (ESEL,S,ENAME,,169,177).
EXCLUDED_ENAMES = range(169, 178)


def run_mapdl(pde_deck, job, export_files):
    for name in export_files + (f"{job}.lock",):
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
        [MAPDL_EXE_DIR, "-b", "-np", "1", "-j", job,
         "-i", pde_deck, "-o", f"{job}.out"],
        cwd=str(SHARED_DIR), timeout=7200,
    )
    for name in export_files:
        if not (SHARED_DIR / name).exists():
            fail(f"{name} not written (code {res.returncode}), see {job}.out")


def read_mmf_vector(path):
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("%")]
    data = np.loadtxt(lines[1:])
    return data if data.ndim == 1 else data[:, -1]


def run_mapdl_onerun(rhs, mesh_nids, l_0):
    # Builds K and M and solves (l_0^2 K + M) sol = M rhs in one MAPDL run;
    # rhs is n x 6 with row k belonging to node mesh_nids[k], sol is returned the same way.
    n = rhs.shape[0]

    # Row k of rhs_nodes.txt holds node number k (zeros for numbers not in the mesh).
    nrow = int(mesh_nids.max())
    rhs_by_nid = np.zeros((nrow, 6))
    rhs_by_nid[mesh_nids - 1] = rhs
    gscale = float(np.abs(rhs).max()) or 1.0
    with open(SHARED_DIR / "pde_onerun_params.txt", "w") as f:
        f.write(f"nrow_={nrow}\nl0sq_={l_0**2:.16e}\ngsc_={1.0/gscale:.16e}\n")

    np.savetxt(SHARED_DIR / "rhs_nodes.txt", rhs_by_nid, fmt="%25.15E", delimiter="")
    run_mapdl("pde_filter_onerun_in.txt", "pdeonerun",
              ("mapb.mtx", "mapb_t.mtx", "sol_onerun.txt", "grad_onerun.txt"))

    back = read_mmf_vector(SHARED_DIR / "mapb.mtx").astype(np.int64)
    back_t = read_mmf_vector(SHARED_DIR / "mapb_t.mtx").astype(np.int64)
    if not np.array_equal(back, back_t):
        fail("equation ordering differs between pdesteady.full and pdetrans.full")
    data = np.loadtxt(SHARED_DIR / "sol_onerun.txt").reshape(-1, 7)
    if not np.array_equal(data[:, 0].astype(np.int64), back):
        fail("node numbers in sol_onerun.txt differ from mapb.mtx (mapping read wrongly by *VREAD)")
    if len(back) != n or not np.array_equal(np.sort(back), mesh_nids):
        fail(f"MAPDL mesh ({len(back)} nodes) differs from model.cdb mesh ({n} nodes)")

    # Equation order -> mesh_nids order (mesh_nids is sorted, so searchsorted gives the row).
    sol = np.empty_like(rhs)
    sol[np.searchsorted(mesh_nids, back)] = data[:, 1:]

    # Row k of grad_onerun.txt is node number k+1: node, then (dx, dy, dz) for each of the 6 components.
    g = np.loadtxt(SHARED_DIR / "grad_onerun.txt").reshape(-1, 19)
    if not np.array_equal(g[:, 0].astype(np.int64), np.arange(1, nrow + 1)):
        fail("grad_onerun.txt is not one row per node number 1..nrow_")
    grad_by_nid = g[:, 1:].reshape(-1, 6, 3) * gscale      # [node number - 1, component, x/y/z]

    return sol, grad_by_nid


def cdb_format_fields(fmt):
    # Finds every <count><i|e><width> in a format line, e.g. "3i9" -> ("3", "i", "9");
    # the count may be empty ("i9"), the width may not.
    fields, pos = [], 0
    while pos < len(fmt):
        j = pos
        while j < len(fmt) and fmt[j].isdecimal():
            j += 1
        k = j + 1
        while j < len(fmt) and fmt[j] in "ie" and k < len(fmt) and fmt[k].isdecimal():
            k += 1
        if k > j + 1:
            fields.append((fmt[pos:j], fmt[j], fmt[j + 1:k]))
            pos = k
        else:
            pos += 1
    return fields



def cdb_field_widths(fmt_line, kind):
    # "(3i9,6e21.13e3)" -> ([9, 9, 9, 21, ...], 3);  "(19i10)" -> ([10]*19, 19)
    widths, nint = [], 0
    for count, letter, width in cdb_format_fields(fmt_line.lower()):
        k = int(count) if count else 1
        widths += [int(width)] * k
        nint += k if letter == "i" else 0
    if not widths:
        fail(f"cannot parse {kind} format line: {fmt_line.strip()}")
    return widths, nint


def split_fixed(line, widths):
    out, pos = [], 0
    for w in widths:
        s = line[pos:pos + w].strip()
        if not s:
            break
        out.append(s)
        pos += w
    return out


def read_cdb_mesh(path):
    # Returns (node ids, xyz per node id row, connectivity (ne, 20) as node ids) for
    # the elements the PDE decks keep, i.e. everything except contact/target (169-177).
    ename = {}
    coords = {}
    elems = []
    with open(path) as f:
        lines = iter(f)
        for line in lines:
            key = line[:8].upper()
            if key.startswith("ET,"):
                p = line.split(",")
                ename[int(p[1])] = int(float(p[2]))
            elif key.startswith("ETBLOCK"):
                next(lines)                                  # format line
                for ln in lines:
                    p = ln.split()
                    if p[0] == "-1":
                        break
                    ename[int(p[0])] = int(p[1])
            elif key.startswith("NBLOCK"):
                widths, nint = cdb_field_widths(next(lines), "NBLOCK")
                for ln in lines:
                    s = ln.strip()
                    if s == "-1" or s.upper().startswith("N,"):
                        break
                    p = split_fixed(ln, widths)
                    xyz = [float(v) for v in p[nint:nint + 3]]
                    coords[int(p[0])] = xyz + [0.0] * (3 - len(xyz))  # trailing zeros are omitted
            elif key.startswith("EBLOCK"):
                if "SOLID" not in line.upper():
                    fail(f"only the SOLID EBLOCK format is supported: {line.strip()}")
                widths, _ = cdb_field_widths(next(lines), "EBLOCK")
                for ln in lines:
                    p = split_fixed(ln, widths)
                    if p[0] == "-1":
                        break
                    etype, nnode = int(p[1]), int(p[8])
                    nodes = [int(v) for v in p[11:]]
                    while len(nodes) < nnode:
                        nodes += [int(v) for v in split_fixed(next(lines), widths)]
                    if ename.get(etype) not in EXCLUDED_ENAMES:
                        elems.append((int(p[10]), nodes[:nnode]))

    if not elems:
        fail(f"no elements found in {path}")
    elems.sort()
    econn = np.zeros((len(elems), 20), dtype=np.int64)
    for k, (_, nodes) in enumerate(elems):
        econn[k, :len(nodes)] = nodes
    nids = np.unique(econn[econn > 0])
    missing = [n for n in nids if n not in coords]
    if missing:
        fail(f"{len(missing)} element nodes missing from NBLOCK, e.g. {missing[:5]}")
    xyz = np.array([coords[n] for n in nids])
    return nids, xyz, econn


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


def PDEFilter_onerun(stress_file_nodal, r):
    mesh_nids, xyz, econn = read_cdb_mesh(SHARED_DIR / "model.cdb")
    n = len(mesh_nids)

    lut = np.full(mesh_nids.max() + 1, -1, dtype=np.int64)
    lut[mesh_nids] = np.arange(n)

    nids = np.array([int(r[1]) for r in stress_file_nodal])
    if nids.max() >= len(lut) or (lut[nids] < 0).any():
        fail("stress-file nodes not in PDE mesh")
    rows = lut[nids]

    stress_nodal = np.array([r[5:11] for r in stress_file_nodal], dtype=float)
    rhs = np.zeros((n, 6))
    rhs[rows] = stress_nodal
    known = np.zeros(n, dtype=bool)
    known[rows] = True

    rhs, known = fill_midside(rhs, known, econn, lut, xyz)
    if not known.all():
        fail(f"{(~known).sum()} nodes still without stress after midside fill")

    l_0 = r / (2*np.sqrt(2))
    sol, grad_by_nid = run_mapdl_onerun(rhs, mesh_nids, l_0)
    stress_filt_nodal = sol[rows]

    out = []
    for row, s in zip(stress_file_nodal, stress_filt_nodal):
        out.append(row[:5] + tuple(s.tolist()))
    return out, grad_by_nid


if __name__ == "__main__":
    rows = read_stress_file(SHARED_DIR / "sigma_export_nodes_testBracket.txt")

    nodal = {}
    for row in rows:
        nodal.setdefault(row[1], row)
    stress_file_nodal = sorted(nodal.values(), key=lambda r: r[1])

    r = 0.003
    out, grad_by_nid = PDEFilter_onerun(stress_file_nodal, r)

    # grad_by_nid[int(row[1]) - 1] is (6, 3): component SXX..SXZ, then d/dx, d/dy, d/dz
    stress_scale = parse_loadcases(SHARED_DIR / "loadcases.txt")[9]
    write_gradient_csv(out, grad_by_nid, gradient_csv_path("sigma_export_nodes_testBracket.txt"), stress_scale)
