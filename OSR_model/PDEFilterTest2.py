import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import os
from Run_pipeline import read_stress_file, SHARED_DIR, Path, MAPDL_EXE_DIR, subprocess, fail, log, coo_matrix

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


def run_mapdl_extract():
    run_mapdl("pde_matrices_extract_in.txt", "pdeextract",
              ("mapb.mtx", "mapb_t.mtx", "pde_dims.txt", "econn.txt", "nxyz.txt"))
    return int(float(np.loadtxt(SHARED_DIR / "pde_dims.txt")))


def run_mapdl_solve(rhs, l_0):
    # Solves (l_0^2 K + M) sol = M rhs in MAPDL; rhs is n x 6 in equation order.
    n = rhs.shape[0]
    with open(SHARED_DIR / "pde_params.txt", "w") as f:
        f.write(f"n_={n}\nl0sq_={l_0**2:.16e}\n")
    np.savetxt(SHARED_DIR / "rhs_pde.txt", rhs, fmt="%25.15E", delimiter="")
    run_mapdl("pde_filter_solve_in.txt", "pdesolve", ("sol_pde.txt",))
    sol = np.loadtxt(SHARED_DIR / "sol_pde.txt").reshape(-1, 6)
    if sol.shape != rhs.shape:
        fail(f"sol_pde.txt has shape {sol.shape}, expected {rhs.shape}")
    return sol


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

def PDEFilter_export(stress_file_nodal,r):
    nn = run_mapdl_extract()
    n = nn

    back = read_mmf_vector(SHARED_DIR / "mapb.mtx").astype(np.int64)
    back_t = read_mmf_vector(SHARED_DIR / "mapb_t.mtx").astype(np.int64)
    if not np.array_equal(back, back_t):
        fail("equation ordering differs between pdesteady.full and pdetrans.full")
    if len(back) != nn:
        fail(f"mapping length {len(back)} != node count {nn}")

    lut = np.full(back.max() + 1, -1, dtype=np.int64)
    lut[back] = np.arange(n)

    nids = np.array([int(r[1]) for r in stress_file_nodal])
    if nids.max() >= len(lut) or (lut[nids] < 0).any():
        fail("stress-file nodes not in PDE mesh")
    rows = lut[nids]

    stress_nodal = np.array([r[5:11] for r in stress_file_nodal], dtype=float)
    rhs = np.zeros((n, 6))
    rhs[rows] = stress_nodal
    known = np.zeros(n, dtype=bool)
    known[rows] = True

    econn = np.loadtxt(SHARED_DIR / "econn.txt").astype(np.int64).reshape(-1, 20)
    nxyz = np.loadtxt(SHARED_DIR / "nxyz.txt").reshape(-1, 3)
    xyz = nxyz[back - 1]

    rhs, known = fill_midside(rhs, known, econn, lut, xyz)
    if not known.all():
        fail(f"{(~known).sum()} nodes still without stress after midside fill")

    l_0 = r / (2*np.sqrt(2))
    sol = run_mapdl_solve(rhs, l_0)
    stress_filt_nodal = sol[rows]

    out = []
    for row, s in zip(stress_file_nodal, stress_filt_nodal):
        out.append(row[:5] + tuple(s.tolist()))

    #return out, K, M, xyz
    return out

rows = read_stress_file(SHARED_DIR / "sigma_export_nodes_testBracket.txt")

nodal = {}
for row in rows:
    nodal.setdefault(row[1], row)
stress_file_nodal = sorted(nodal.values(), key=lambda r: r[1])

r = 0.003
out = PDEFilter_export(stress_file_nodal, r)