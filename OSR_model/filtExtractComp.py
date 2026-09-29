"""Compare the two PDE filter implementations on the same stress field.

  Test1 (PDEFilterTest1.PDEFilter_splu):  K, M exported from MAPDL, solve in Python (splu)
  Test2 (PDEFilterTest2.PDEFilter_export): solve done inside MAPDL (*LSENGINE / *LSBAC)

Both get the same de-duplicated nodal stress rows and the same filter length
l_0 = r / (2*sqrt(2)) (the Run_pipeline convention). Test1 is used as the reference
for the percentage differences.

Usage:  python filtExtractComp.py [--stress-file FILE] [--r 0.003] [--repeats 1]
"""
import argparse
import time
from contextlib import contextmanager

import numpy as np

import PDEFilterTest1 as t1
import PDEFilterTest2 as t2
from Run_pipeline import SHARED_DIR, read_stress_file

COMPONENTS = ("SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ")


def load_nodal(path):
    # Same de-duplication as PDEFilterTest2: first row per node id, sorted by node id.
    nodal = {}
    for row in read_stress_file(path):
        nodal.setdefault(row[1], row)
    return sorted(nodal.values(), key=lambda r: r[1])


@contextmanager
def stage_timers(module, names):
    # Temporarily wrap module-level functions so the time spent in each is recorded.
    spent = {name: 0.0 for name in names}
    originals = {name: getattr(module, name) for name in names}

    def wrap(name, fn):
        def timed(*args, **kwargs):
            t0 = time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                spent[name] += time.perf_counter() - t0
        return timed

    for name, fn in originals.items():
        setattr(module, name, wrap(name, fn))
    try:
        yield spent
    finally:
        for name, fn in originals.items():
            setattr(module, name, fn)


def run_timed(run, module, stages, repeats):
    runs = []
    for _ in range(repeats):
        with stage_timers(module, stages) as spent:
            t0 = time.perf_counter()
            out = run()
            total = time.perf_counter() - t0
        runs.append({"total": total, **spent})
    return out, runs


def von_mises(s):
    sxx, syy, szz, sxy, syz, sxz = s.T
    return np.sqrt(0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2)
                   + 3.0 * (sxy**2 + syz**2 + sxz**2))


def pct(num, den):
    return 100.0 * num / den if den != 0 else float("nan")


def print_timing(label, runs, stages):
    print(f"\n{label}")
    for key in ("total",) + tuple(stages):
        vals = np.array([r[key] for r in runs])
        extra = f"   (min {vals.min():.2f} s, max {vals.max():.2f} s)" if len(vals) > 1 else ""
        print(f"  {key:<20s} {vals.mean():10.2f} s{extra}")
    other = np.mean([r["total"] - sum(r[k] for k in stages) for r in runs])
    print(f"  {'other (Python)':<20s} {other:10.2f} s")


def compare(s1, s2, nids):
    d = s2 - s1
    print("\nPer-component difference (Test2 vs Test1 reference)")
    print(f"  {'comp':<5s} {'rel L2 %':>10s} {'max|d| % of max|s1|':>21s} {'mean|d| % of mean|s1|':>23s}")
    for c, name in enumerate(COMPONENTS):
        a, dc = s1[:, c], d[:, c]
        print(f"  {name:<5s} {pct(np.linalg.norm(dc), np.linalg.norm(a)):10.4f} "
              f"{pct(np.abs(dc).max(), np.abs(a).max()):21.4f} "
              f"{pct(np.abs(dc).mean(), np.abs(a).mean()):23.4f}")
    print(f"  {'all':<5s} {pct(np.linalg.norm(d), np.linalg.norm(s1)):10.4f}")

    vm1, vm2 = von_mises(s1), von_mises(s2)
    dvm = vm2 - vm1
    i_peak1, i_peak2 = vm1.argmax(), vm2.argmax()
    i_worst = np.abs(dvm).argmax()
    # Per-node % differences are only meaningful where the stress is not ~0.
    sig = vm1 > 0.01 * vm1.max()
    node_pct = np.abs(dvm[sig]) / vm1[sig] * 100.0

    print("\nVon Mises")
    print(f"  peak Test1               {vm1[i_peak1]:.6e}  (node {nids[i_peak1]})")
    print(f"  peak Test2               {vm2[i_peak2]:.6e}  (node {nids[i_peak2]})")
    print(f"  peak difference          {pct(vm2[i_peak2] - vm1[i_peak1], vm1[i_peak1]):+.4f} %")
    print(f"  rel L2 difference        {pct(np.linalg.norm(dvm), np.linalg.norm(vm1)):.4f} %")
    print(f"  max |d| / peak           {pct(np.abs(dvm).max(), vm1.max()):.4f} %  (node {nids[i_worst]})")
    print(f"  per-node |d|/vm1 (nodes with vm1 > 1% of peak, {sig.sum()} of {len(vm1)}):")
    print(f"    mean {node_pct.mean():.4f} %   median {np.median(node_pct):.4f} %   "
          f"p95 {np.percentile(node_pct, 95):.4f} %   max {node_pct.max():.4f} %")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stress-file", default=str(SHARED_DIR / "sigma_export_nodes_testBracket.txt"))
    ap.add_argument("--r", type=float, default=0.003, help="filter radius; l_0 = r/(2*sqrt(2))")
    ap.add_argument("--repeats", type=int, default=1, help="runs per method for timing")
    args = ap.parse_args()

    stress_file_nodal = load_nodal(args.stress_file)
    nids = np.array([row[1] for row in stress_file_nodal])
    l_0 = args.r / (2 * np.sqrt(2))
    print(f"Stress file: {args.stress_file}  ({len(stress_file_nodal)} unique nodes)")
    print(f"r = {args.r:g}, l_0 = {l_0:.6g}, repeats = {args.repeats}")

    stages1 = ("run_mapdl_extract", "splu")
    stages2 = ("run_mapdl_extract", "run_mapdl_solve")

    (s1, _), runs1 = run_timed(lambda: t1.PDEFilter_splu(stress_file_nodal, l_0),
                               t1, stages1, args.repeats)
    out2, runs2 = run_timed(lambda: t2.PDEFilter_export(stress_file_nodal, args.r),
                            t2, stages2, args.repeats)
    s2 = np.array([row[5:11] for row in out2], dtype=float)

    if s1.shape != s2.shape:
        raise SystemExit(f"output shapes differ: Test1 {s1.shape}, Test2 {s2.shape}")
    if not np.array_equal(nids, [row[1] for row in out2]):
        raise SystemExit("Test2 output node order differs from input")

    print("\nRun time (mean over repeats)")
    print_timing("Test1 (Python splu solve)", runs1, stages1)
    print_timing("Test2 (MAPDL solve)", runs2, stages2)
    t1_mean = np.mean([r["total"] for r in runs1])
    t2_mean = np.mean([r["total"] for r in runs2])
    print(f"\n  Test2 / Test1 total time = {t2_mean / t1_mean:.3f}  "
          f"({pct(t2_mean - t1_mean, t1_mean):+.1f} %)")

    compare(s1, s2, nids)


if __name__ == "__main__":
    main()
