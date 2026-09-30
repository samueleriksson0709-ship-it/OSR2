"""PDE filter in a single MAPDL run.

Same maths as PDEFilterTest1/2, but the mesh (coordinates and connectivity for the
midside fill) is read straight from model.cdb, so Python can prepare the stresses
before MAPDL starts. They are written by node number; pde_filter_onerun_in.txt
builds the matrices, reorders the stresses into equation order, solves and writes
the result with node numbers.
"""
import re
import numpy as np
from Run_pipeline import read_stress_file, SHARED_DIR, fail, log
from PDEFilterTest2 import run_mapdl, fill_midside

# Element types the decks delete before building the PDE matrices (ESEL,S,ENAME,,169,177).
EXCLUDED_ENAMES = range(169, 178)


def _field_widths(fmt_line, kind):
    # "(3i9,6e21.13e3)" -> ([9, 9, 9, 21, ...], 3);  "(19i10)" -> ([10]*19, 19)
    widths, nint = [], 0
    for count, letter, width in re.findall(r"(\d*)([ie])(\d+)", fmt_line.lower()):
        k = int(count) if count else 1
        widths += [int(width)] * k
        nint += k if letter == "i" else 0
    if not widths:
        fail(f"cannot parse {kind} format line: {fmt_line.strip()}")
    return widths, nint


def _split_fixed(line, widths):
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
                widths, nint = _field_widths(next(lines), "NBLOCK")
                for ln in lines:
                    s = ln.strip()
                    if s == "-1" or s.upper().startswith("N,"):
                        break
                    p = _split_fixed(ln, widths)
                    xyz = [float(v) for v in p[nint:nint + 3]]
                    coords[int(p[0])] = xyz + [0.0] * (3 - len(xyz))  # trailing zeros are omitted
            elif key.startswith("EBLOCK"):
                if "SOLID" not in line.upper():
                    fail(f"only the SOLID EBLOCK format is supported: {line.strip()}")
                widths, _ = _field_widths(next(lines), "EBLOCK")
                for ln in lines:
                    p = _split_fixed(ln, widths)
                    if p[0] == "-1":
                        break
                    etype, nnode = int(p[1]), int(p[8])
                    nodes = [int(v) for v in p[11:]]
                    while len(nodes) < nnode:
                        nodes += [int(v) for v in _split_fixed(next(lines), widths)]
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

    # Row k of rhs_nodes.txt holds node number k (zeros for numbers not in the mesh).
    nrow = int(mesh_nids.max())
    rhs_by_nid = np.zeros((nrow, 6))
    rhs_by_nid[mesh_nids - 1] = rhs
    l_0 = r / (2*np.sqrt(2))
    with open(SHARED_DIR / "pde_onerun_params.txt", "w") as f:
        f.write(f"nrow_={nrow}\nl0sq_={l_0**2:.16e}\n")
    np.savetxt(SHARED_DIR / "rhs_nodes.txt", rhs_by_nid, fmt="%25.15E", delimiter="")

    run_mapdl("pde_filter_onerun_in.txt", "pdeonerun", ("sol_onerun.txt",))

    data = np.loadtxt(SHARED_DIR / "sol_onerun.txt").reshape(-1, 8)
    back, back_t, sol = data[:, 0].astype(np.int64), data[:, 1].astype(np.int64), data[:, 2:]
    if not np.array_equal(back, back_t):
        fail("equation ordering differs between pdesteady.full and pdetrans.full")
    if len(back) != n or not np.array_equal(np.sort(back), mesh_nids):
        fail(f"MAPDL mesh ({len(back)} nodes) differs from model.cdb mesh ({n} nodes)")

    sol_by_row = np.empty_like(sol)
    sol_by_row[lut[back]] = sol
    stress_filt_nodal = sol_by_row[rows]

    out = []
    for row, s in zip(stress_file_nodal, stress_filt_nodal):
        out.append(row[:5] + tuple(s.tolist()))
    return out


if __name__ == "__main__":
    rows = read_stress_file(SHARED_DIR / "sigma_export_nodes_testBracket.txt")

    nodal = {}
    for row in rows:
        nodal.setdefault(row[1], row)
    stress_file_nodal = sorted(nodal.values(), key=lambda r: r[1])

    r = 0.003
    out = PDEFilter_onerun(stress_file_nodal, r)
