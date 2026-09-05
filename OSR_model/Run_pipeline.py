"""
Automated FEM--OSR pipeline.

The script prepares Ansys stress exports for the Fortran OSR solver, launches the
solver, and copies the damage outputs back to the shared Workbench directory.

Deduplication principle:
for each shared physical node, several element--node stress contributions may
exist. The script keeps the element contribution that maximises the OSR
activation measure max_t beta(t, alpha=0) for the prescribed loading definition.
For multiaxial cases, the selected (eid, nid) pairs are chosen from the first
stress file and kept consistently across all loadcase files, ensuring coherent
tensor superposition.
"""

import math
import os
import subprocess
import shutil
import sys
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import numpy as np 

# PATH CONFIGURATION
# Nothing here is tied to a particular machine or user account. The folders are
# derived from the location of this file, so the pipeline works from any
# checkout, and can be overridden with environment variables (the Workbench ACT
# extension sets them when it launches this script).
#
#   <OSR root>/OSR_model/Run_pipeline.py   -> FORTRAN_DIR = <OSR root>/OSR_model
#                                          -> SHARED_DIR  = <OSR root>/shared
SCRIPT_DIR = Path(__file__).resolve().parent
OSR_ROOT   = Path(os.environ.get("OSR_ROOT") or SCRIPT_DIR.parent)

FORTRAN_DIR = Path(os.environ.get("OSR_FORTRAN_DIR") or SCRIPT_DIR)
SHARED_DIR  = Path(os.environ.get("OSR_SHARED_DIR") or (OSR_ROOT / "shared"))
EXE_NAME    = "osr_pipeline_multi.exe" if os.name == "nt" else "osr_pipeline_multi"

# Optional: reuse the existing Workbench ACT configuration when available.
# The current script mainly uses it to retrieve omp_threads, while keeping the
# main paths explicit above for clarity during handover.
ACT_CONFIG_FILE = SHARED_DIR / "act_config.json"

INPUT_LOADCASES = "loadcases.txt"
OUTPUT_CSV      = "damage.csv"
OUTPUT_FIELD    = "damage_field.txt"
LOG_FILE        = SHARED_DIR / "pipeline_log.txt"
OUTPUT_CRITICAL = "critical_points.txt"

# Coordinates of the fixed supports, in the length unit of the export (these ones
# come from a model meshed in mm). Update them for a model solved in MKS, where
# the exported coordinates are in m, before turning USE_CLAMP_FILTER on.
CLAMP_CENTERS = [ (-1.4912e-11, -25.936 , 28.185),
    (-1.5118e-11, 2.6002, 2.9392),
    (-148.06, -22.569, 25.206),
    (-148.06, 16.431, -9.2965)
]

CLAMP_RADIUS = 15.0

A_OSR = 0.225
SIGMA_OE = 490.0
NSAMP_BETA = 128

# The OSR material parameters above are calibrated in MPa, but Ansys exports in the
# unit system the model was solved in: MKS gives Pa, the mm-kg-N system gives MPa.
# The unit is declared as an optional 4th token on the first data line of
# loadcases.txt ("n_cases omega mode [stress_unit]") and this table turns it into
# the factor that brings the export to MPa. The same table exists in
# osr_pipeline_multi.f90; the Fortran reader is what actually rescales the field,
# so the .txt and .bin files written here stay in the unit of the export.
STRESS_UNITS = {
    "MPA": 1.0, "N/MM2": 1.0, "NMM": 1.0,
    "PA": 1.0e-6, "N/M2": 1.0e-6, "MKS": 1.0e-6, "SI": 1.0e-6,
    "KPA": 1.0e-3,
    "GPA": 1.0e3,
    "PSI": 6.894757e-3,
    "KSI": 6.894757,
}
DEFAULT_STRESS_UNIT = "MPA"

# enable this only for geometries where fixed-support singularities must be excluded from fatigue post-processing
USE_CLAMP_FILTER  = False 
USE_DEDUPLICATION = True

# Change this string if the binary-cache format or the deduplication criterion
# changes. This prevents silently reusing old .bin files after code changes.
CACHE_VERSION = "beta_dedup_clamp__toggle_v1"

def load_optional_act_config():
    """Read the Workbench ACT config if it exists; otherwise return an empty dict."""
    if not ACT_CONFIG_FILE.exists():
        return {}
    try:
        with open(ACT_CONFIG_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        print(f"[WARNING] Could not read {ACT_CONFIG_FILE.name}: {e}")
        return {}

def get_omp_threads(default=8):
    """Return the OpenMP thread count, preferably from act_config.json."""
    cfg = load_optional_act_config()
    try:
        return str(int(cfg.get("omp_threads", default)))
    except Exception:
        return str(default)

def log(message, also_print=True):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {message}"
    try:
        with open(LOG_FILE, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except Exception:
        pass
    if also_print:
        print(line)

def fail(message):
    """Log an error and exit; always write to the log before quitting."""
    log(f"[ERROR] {message}")
    log("=" * 60, also_print=False)
    sys.exit(1)

def stress_tensor_from_row(row):
    _,_,_,_,_, sxx, syy, szz, sxy, syz, sxz = row 
    return [
        [sxx, sxy, sxz],
        [sxy, syy, syz],
        [sxz, syz, szz],
    ]

def add_scaled_tensor(out, scale, s):
    for i in range(3):
        for j in range(3):
            out[i][j] += scale*s[i][j] 

def stress_scale_to_mpa(token):
    """Conversion factor to MPa for a unit declared in loadcases.txt.

    A bare number is accepted as an explicit factor, for a unit not in the table.
    """
    key = token.strip().upper()
    if key in STRESS_UNITS:
        return STRESS_UNITS[key]
    try:
        factor = float(key)
    except ValueError:
        fail(f"unknown stress unit '{token}' in loadcases.txt. "
             f"Expected one of {', '.join(sorted(STRESS_UNITS))}, "
             f"or a positive conversion factor to MPa.")
    if factor <= 0.0:
        fail(f"stress unit factor must be positive, got '{token}'")
    return factor

def beta_osr_alpha0(sigma, scale=1.0):
    """compute beta with alpha = 0, for a stress tensor in the export's own unit"""
    tr = sigma[0][0] + sigma[1][1] + sigma[2][2]

    mean = tr/3.0
    dev = [
        [sigma[0][0]-mean, sigma[0][1], sigma[0][2]],
        [sigma[1][0], sigma[1][1]-mean, sigma[1][2]],
        [sigma[2][0], sigma[2][1], sigma[2][2]-mean],
    ]

    dev_contract = 0.0
    for i in range(3):
        for j in range(3):
            dev_contract += dev[i][j]*dev[i][j]

    sigma_eff = math.sqrt(1.5*dev_contract)*scale
    beta = (sigma_eff+A_OSR*tr*scale-SIGMA_OE)/SIGMA_OE

    return beta 

def is_near_clamp(x,y,z):  
    for cx, cy, cz in CLAMP_CENTERS:
        dist = math.sqrt((x-cx)**2+(y-cy)**2+(z-cz)**2)
        if dist < CLAMP_RADIUS:
            return True 
    return False 

def filter_clamps(rows):
    if not USE_CLAMP_FILTER :
        return rows, 0
    
    kept=[]
    remove = 0 

    for row in rows : 
        x, y, z = row[2], row[3], row[4]

        if is_near_clamp(x,y,z):
            remove+=1
            continue 
        kept.append(row)

    return kept, remove 

def write_rows_file(rows, dst_path):
    with open(dst_path, "w") as f:
        f.write("# eid nid X Y Z SXX SYY SZZ SXY SYZ SXZ\n")
        for row in rows:
            eid, nid, x, y, z, sxx, syy, szz, sxy, syz, sxz = row
            f.write(
                f"{float(eid):8.0f} {float(nid):8.0f}"
                f"{x:14.6E} {y:14.6E} {z:14.6E}"
                f"{sxx:14.6E} {syy:14.6E} {szz:14.6E}"
                f"{sxy:14.6E} {syz:14.6E} {sxz:14.6E}\n"
            )
    
def parse_line(line):
    line = line.strip()
    if not line:
        return None
    if line[0] in "#abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ":
        return None
    parts = line.split()
    if len(parts) < 11:
        return None
    try:
        eid = int(float(parts[0]))
        nid = int(float(parts[1]))
        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        sxx = float(parts[5])
        syy = float(parts[6])
        szz = float(parts[7])
        sxy = float(parts[8])
        syz = float(parts[9])
        sxz = float(parts[10])
        return (eid, nid, x, y, z, sxx, syy, szz, sxy, syz, sxz)
    except (ValueError, IndexError):
        return None

def read_stress_file(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parsed = parse_line(line)
            if parsed is not None:
                rows.append(parsed)
    return rows

def write_critical_points_file(damage_field_path, critical_path, n_top=20):
    
    #read damage_field.txt and write the n_top points with the smallest finite N_fail.
    
    rows = []

    with open(damage_field_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 7:
                continue

            try:
                eid = int(float(parts[0]))
                nid = int(float(parts[1]))
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                D = float(parts[5])
                N_fail = float(parts[6])

                # Keep only actual finite failures
                if math.isfinite(N_fail) and N_fail < 1.0e100:
                    rows.append((N_fail, eid, nid, x, y, z, D))

            except ValueError:
                continue

    rows.sort(key=lambda r: r[0])

    with open(critical_path, "w") as f:
        f.write("# rank eid nid X Y Z D N_fail\n")

        for rank, row in enumerate(rows[:n_top], start=1):
            N_fail, eid, nid, x, y, z, D = row

            f.write(
                f"{rank:4d} "
                f"{eid:8d} {nid:8d} "
                f"{x:14.6E} {y:14.6E} {z:14.6E} "
                f"{D:14.6E} {N_fail:14.6E}\n"
            )

    return len(rows), min(n_top, len(rows))

def write_deduplicated_file(rows, chosen_eids, dst_path):
    n_in = len(rows)
    kept_rows = []
    seen_nids = set() 

    for row in rows:
        eid, nid = row[0], row[1]

        if nid in chosen_eids and chosen_eids[nid] == eid and nid not in seen_nids:
            kept_rows.append(row)
            seen_nids.add(nid)

    kept_rows.sort(key=lambda r: r[1])
    write_rows_file(kept_rows, dst_path)

    return n_in, len(kept_rows)

#  Binary cache helper functions 
def get_kept_rows(rows, chosen_eids):
    """Return deduplicated rows without writing them to disk."""
    kept, seen = [], set()
    for row in rows:
        eid, nid = row[0], row[1]
        if nid in chosen_eids and chosen_eids[nid] == eid and nid not in seen:
            kept.append(row)
            seen.add(nid)
    kept.sort(key=lambda r: r[1])
    return kept

def bin_path_for(fname, fortran_dir):
    """Return the .bin path associated with a stress .txt file."""
    return fortran_dir / (Path(fname).stem + ".bin")

def cache_metadata_path(bin_path):
    """Sidecar JSON file used to validate the binary cache."""
    return bin_path.with_suffix(bin_path.suffix + ".meta.json")

def current_cache_metadata(mode, history_filename, stress_unit=DEFAULT_STRESS_UNIT):
    """Metadata that affects clamp filtering and beta-based deduplication."""
    return {
        "cache_version": CACHE_VERSION,
        "A_OSR": A_OSR,
        "SIGMA_OE": SIGMA_OE,
        "NSAMP_BETA": NSAMP_BETA,
        "USE_DEDUPLICATION": USE_DEDUPLICATION,
        "USE_CLAMP_FILTER": USE_CLAMP_FILTER,
        "CLAMP_RADIUS": CLAMP_RADIUS,
        "CLAMP_CENTERS": [list(c) for c in CLAMP_CENTERS],
        "mode": mode,
        "history_filename": history_filename if mode == "HISTORY" else "-",
        "stress_unit": stress_unit.upper(),
    }

def write_cache_metadata(bin_path, metadata):
    """Write cache metadata next to the .bin file."""
    with open(cache_metadata_path(bin_path), "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, sort_keys=True)

def binary_cache_valid(bin_dst_path, dependencies, expected_metadata):
    """
    True if the .bin exists, is newer than all relevant input files, and its
    sidecar metadata matches the current deduplication settings.
    """
    if not bin_dst_path.exists():
        return False

    bin_mtime = bin_dst_path.stat().st_mtime

    for dep in dependencies:
        if dep is None:
            continue
        if not dep.exists():
            return False
        if bin_mtime <= dep.stat().st_mtime:
            return False

    meta_path = cache_metadata_path(bin_dst_path)
    if not meta_path.exists():
        return False

    try:
        with open(meta_path, "r", encoding="utf-8") as f:
            cached_metadata = json.load(f)
    except Exception:
        return False

    return cached_metadata == expected_metadata

def write_binary_cache(kept_rows, bin_path):
    """
    Write deduplicated rows as a NumPy-compatible binary array (float64, shape N x 11).
    Format: eid nid X Y Z SXX SYY SZZ SXY SYZ SXZ, without a header.
    read by read_stress_field_binary / read_stress_field_auto (pipeline_aux.f90).
    """
    arr = np.array(kept_rows, dtype=np.float64)
    arr.tofile(bin_path)

def parse_loadcases_for_filenames(loadcases_path):
    filenames = []
    with open(loadcases_path, "r") as f:
        lines = [l.strip() for l in f if l.strip() and l.strip()[0] not in "#!"]
    if len(lines) < 2:
        return filenames
    try:
        n_cases = int(lines[0].split()[0])
    except (ValueError, IndexError):
        return filenames
    for i in range(1, min(1 + n_cases, len(lines))):
        parts = lines[i].split()
        if parts:
            filenames.append(parts[0])
    return filenames

def parse_loadcases(loadcases_path): # read the different values needed for the construction of the stress path

    # new supported format : 
    # SINUS or HISTORY mode : 
    # n_cases omega mode [stress_unit]
    # file amp phase mean 

    loadcases = []

    with open(loadcases_path, "r") as f:
        lines = [l.strip() for l in f if l.strip() and l.strip()[0] not in "#!"]

    header = lines[0].split()
    n_cases = int(header[0])
    omega = float(header[1])

    if len(header)>=3 :
        mode = header[2].upper()
    else : 
        mode = "SINUS" 

    # the stress unit is optional : a loadcases.txt written before it existed
    # declares nothing, which means MPa
    stress_unit = header[3] if len(header) >= 4 else DEFAULT_STRESS_UNIT
    stress_scale = stress_scale_to_mpa(stress_unit)

    if mode not in ("SINUS", "HISTORY"):
        fail(f"invalid loading mode '{mode}' in loadcases.txt"
             "Expected SINUS or HISTORY")

    for i in range(1, 1 + n_cases):
        parts = lines[i].split()

        fname = parts[0]
        mean = float(parts[3])
        amp = float(parts[1])
        phase = float(parts[2])

        loadcases.append((fname, mean, amp, phase))
    
    #debug 
    control_parts = lines[1+n_cases].split()
    n_cycles_max = int(float(control_parts[0]))
    n_steps_per_cycle = int(float(control_parts[1]))
    n_substeps_history = int(float(control_parts[2]))

    # history file line 
    history_filename = lines[2+n_cases].split()[0]

    history = None 

    if mode == "HISTORY" : 
        history_path = loadcases_path.parent/history_filename
        if not history_path.exists():
            fail(f"history file not found : {history_path}")
        history = read_history_file(history_path,n_cases)

    return (loadcases, omega, mode, history, n_cycles_max, n_steps_per_cycle,
            n_substeps_history, history_filename, stress_unit, stress_scale)

def select_reference_eids_by_max_beta(all_rows_by_case, loadcases, mode = "SINUS", history = None,
                                      stress_scale = 1.0):

    #Select one eid per nid by maximizing max_t beta(t, alpha=0).
    #all_rows_by_case: list of rows lists, one per load case.
    #loadcases: list of (fname, mean, amp, phase)
    #mode SINUS : lambda_k(theta) = mean_k + amp_k*sin(theta + phase_k) 
    #mode HISTORY : lambda_k values are read from history.txt

    # Build one lookup dictionary per loadcase for fast access by (eid, nid).
    case_maps = []
    for rows in all_rows_by_case:
        d = {}
        for row in rows:
            eid, nid = row[0], row[1]
            d[(eid, nid)] = row
        case_maps.append(d)

    # Candidate element contributions are taken from the first stress file.
    candidates_by_nid = defaultdict(list)
    for row in all_rows_by_case[0]:
        eid, nid = row[0], row[1]
        candidates_by_nid[nid].append(eid)

    chosen = {}

    for nid, candidate_eids in candidates_by_nid.items():
        best_eid = None
        best_beta = -1.0e300

        for eid in candidate_eids:
            # Check that the selected (eid, nid) pair exists in every loadcase.
            rows_for_candidate = []
            valid = True

            for case_map in case_maps:
                key = (eid, nid)
                if key not in case_map:
                    valid = False
                    break
                rows_for_candidate.append(case_map[key])

            if not valid:
                continue

            # Convert unit constraints in tensors
            tensors = [stress_tensor_from_row(row) for row in rows_for_candidate]

            beta_max = -1.0e300

            if mode == "HISTORY" : 
                if history is None : 
                    raise RuntimeError ("HISTORY mode selected but no history data provided ")
                
                for _, lambdas in history :
                    sigma = [
                        [0.0,0.0,0.0],
                        [0.0,0.0,0.0],
                        [0.0,0.0,0.0],
                    ]

                    for k, lam in enumerate(lambdas) : 
                        add_scaled_tensor(sigma, lam, tensors[k])
                    
                    beta = beta_osr_alpha0(sigma, stress_scale)

                    if beta > beta_max : 
                        beta_max = beta 
            
            else : 
                
                for m in range(NSAMP_BETA):
                    theta = 2.0 * math.pi * float(m) / float(NSAMP_BETA)
                    
                    sigma = [
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                    ]

                    for k, (_, mean, amp, phase) in enumerate(loadcases):
                        lam = mean + amp * math.sin(theta + phase)
                        add_scaled_tensor(sigma, lam, tensors[k])

                    beta = beta_osr_alpha0(sigma, stress_scale)

                    if beta > beta_max:
                        beta_max = beta

            if beta_max > best_beta:
                best_beta = beta_max
                best_eid = eid

        if best_eid is not None:
            chosen[nid] = best_eid

    return chosen

def read_history_file(history_path,n_cases):

    # expected format : 
    # time  lambda_1  lambda_2 ... lambda_n_cases 

    history = []

    with open(history_path, "r") as f : 
        for line in f:
            line=line.strip()

            if not line or line.startswith("#"):
                continue 

            parts = line.split()
            if len(parts) < n_cases + 1 : 
                continue 

            t = float(parts[0])
            lambdas = [float(x) for x in parts[1:1+n_cases]]

            history.append((t, lambdas))
    
    if not history : 
        raise RuntimeError("Empty or invalid history file : " + str(history_path))
    
    return history 

def main():


    SHARED_DIR.mkdir(parents=True, exist_ok=True)
    log("Python script started")
    log(f"  SHARED_DIR  = {SHARED_DIR}")
    log(f"  FORTRAN_DIR = {FORTRAN_DIR}")

    # Initial checks
    loadcases_src = SHARED_DIR / INPUT_LOADCASES
    if not loadcases_src.exists():
        fail(f"{loadcases_src} not found")

    stress_filenames = parse_loadcases_for_filenames(loadcases_src)
    if not stress_filenames:
        fail(f"No stress files found in {loadcases_src}")
    log(f"  {len(stress_filenames)} stress file(s) declared in loadcases.txt")

    for fname in stress_filenames:
        p = SHARED_DIR / fname
        if not p.exists():
            fail(f"Stress file not found: {p}")
        log(f"  [OK] {fname} found")

    exe_path = FORTRAN_DIR / EXE_NAME
    if not exe_path.exists():
        fail(
            f"Fortran solver not built on this machine: {exe_path}\n"
            f"         Open a terminal in {FORTRAN_DIR} and run 'make' "
            f"(gfortran with OpenMP support is required)."
        )
    log(f"  [OK] {EXE_NAME} found")

    # Deduplication with binary cache
    # parse_loadcases is called here so that mode/history_filename are available
    # even when the binary cache is used.
    loadcases, omega, mode, history, n_cycles_max, n_steps_per_cycle, \
        n_substeps_history, history_filename, stress_unit, stress_scale = \
        parse_loadcases(loadcases_src)

    log(f"  stress unit of the export = {stress_unit} "
        f"(x{stress_scale:g} -> MPa, the unit of the OSR material parameters)")

    if USE_DEDUPLICATION:

        # Check whether all binary caches are valid.
        # A cache is valid only if it is newer than the stress .txt file,
        # loadcases.txt, and history.txt in HISTORY mode. Its metadata must
        # also match the current clamp and beta-deduplication settings.
        history_src = SHARED_DIR / history_filename if mode == "HISTORY" else None
        cache_metadata = current_cache_metadata(mode, history_filename, stress_unit)

        all_cache_valid = all(
            binary_cache_valid(
                bin_path_for(fname, FORTRAN_DIR),
                dependencies=[SHARED_DIR / fname, loadcases_src, history_src],
                expected_metadata=cache_metadata,
            )
            for fname in stress_filenames
        )

        if all_cache_valid:
            # Cache hit: skip the full deduplication step.
            log("Binary cache valid: skipping deduplication (cache hit)")
            for fname in stress_filenames:
                bp = bin_path_for(fname, FORTRAN_DIR)
                n_pts = bp.stat().st_size // (11 * 8)
                log(f"  Using cached {bp.name} ({n_pts} points)")

        else:
            # Cache miss: read, filter, deduplicate, and write a new cache.
            log("Binary cache missing or stale: running full deduplication")
            all_rows_by_case = []
            raw_counts_by_case = []

            for fname in stress_filenames:
                rows = read_stress_file(SHARED_DIR / fname)
                n_raw_before_filter = len(rows)
                rows, n_removed = filter_clamps(rows)
                all_rows_by_case.append(rows)
                raw_counts_by_case.append(n_raw_before_filter)

                if USE_CLAMP_FILTER : 
                    log(f" {fname}: {len(rows)} rows after clamp filtering ({n_removed} removed)")
                else : 
                    log(f" {fname}: {len(rows)} rows read (clamp filter disabled)")

            chosen_eids = select_reference_eids_by_max_beta(
                all_rows_by_case, loadcases, mode=mode, history=history,
                stress_scale=stress_scale
            )
            log(f" {len(chosen_eids)} unique nodes selected using max beta criterion")

            for k, fname in enumerate(stress_filenames):
                try:
                    kept_rows = get_kept_rows(all_rows_by_case[k], chosen_eids)
                    n_raw = len(all_rows_by_case[k])
                    n_out = len(kept_rows)

                    # Text file written for human inspection.
                    write_rows_file(kept_rows, FORTRAN_DIR / fname)

                    # Binary cache read by Fortran through read_stress_field_auto
                    bp = bin_path_for(fname, FORTRAN_DIR)
                    log(bp)
                    log(f"Binary path = {bp}")
                    write_binary_cache(kept_rows, bp)
                    write_cache_metadata(bp, cache_metadata)

                    log(f"  {fname}: {n_raw} raw -> {n_out} points "
                        f"[.txt + .bin written to {FORTRAN_DIR.name}/]")

                except Exception as e:
                    fail(f"Error processing {fname}: {e}")
    else:
        if USE_CLAMP_FILTER :
            log("Mode: clamp filtering without deduplication")
        else : 
            log("Mode: no deduplication and clamp filtering disabled")

        
        for fname in stress_filenames:
            src = SHARED_DIR / fname
            dst = FORTRAN_DIR / fname
            
            try:
                rows = read_stress_file(src)
                n_raw = len(rows)
                rows, n_clamp_removed = filter_clamps(rows)
                #rows = rows[:20000]
                write_rows_file(rows, dst)

                if USE_CLAMP_FILTER : 
                    log(
                    f"  {fname}: {n_raw} raw -> {len(rows)} after clamp filter "
                    f"({n_clamp_removed} rows removed near clamps, no dedup)"
                )
                
                else : 
                    log(
                    f"  {fname}: {n_raw} raw rows copied "
                    )

            except Exception as e:
                fail(f"Error processing {fname}: {e}")

    shutil.copy2(loadcases_src, FORTRAN_DIR / INPUT_LOADCASES)
    log("[OK] loadcases.txt copied")

    if mode == "HISTORY":
        history_src = SHARED_DIR / history_filename

        if not history_src.exists() : 
            fail(f"{history_filename} not found in shared directory")
        
        shutil.copy2(history_src, FORTRAN_DIR / history_filename)
        log("[OK] history file copied")

    # Launch Fortran solver
    log("Launching osr_pipeline_multi.exe")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = get_omp_threads(default=8)
    log(f"  OMP_NUM_THREADS = {env['OMP_NUM_THREADS']}")

    try:
        with open(LOG_FILE, "a", encoding="utf-8") as log_f:
            log_f.write(" Fortran output \n")
            log_f.flush()
            result = subprocess.run(
                [str(exe_path)],
                cwd=str(FORTRAN_DIR),
                env=env,
                stdout=log_f,
                stderr=log_f,
                timeout=7200,
            )
            log_f.write(" End Fortran output \n")
        if result.returncode != 0:
            fail(f"Fortran exe returned code {result.returncode}")
    except subprocess.TimeoutExpired:
        fail("Fortran exe timed out (>2h)")
    except Exception as e:
        fail(f"Unexpected error launching exe: {e}")

    log("[OK] Fortran pipeline finished")

    # Copy output files
    damage_csv = FORTRAN_DIR / OUTPUT_CSV
    if not damage_csv.exists():
        fail(f"{OUTPUT_CSV} not produced by Fortran exe")
    shutil.copy2(damage_csv, SHARED_DIR / OUTPUT_CSV)
    log(f"[OK] {OUTPUT_CSV} copied to shared")

    damage_field = FORTRAN_DIR / OUTPUT_FIELD
    if damage_field.exists():
        shutil.copy2(damage_field, SHARED_DIR / OUTPUT_FIELD)
        log(f"[OK] {OUTPUT_FIELD} copied to shared")

        critical_path = SHARED_DIR/OUTPUT_CRITICAL

        try:
            n_failed, n_written = write_critical_points_file(damage_field, critical_path, n_top=20)

            log (f"[OK] {OUTPUT_CRITICAL} written to shared")
        
        except Exception as e : 
            log (f"[WARNING] could not write {OUTPUT_CRITICAL} : {e}")

    log("Pipeline completed successfully")
    log("=" * 60, also_print=False)
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        log(f"[FATAL ERROR] {type(e).__name__}: {e}")
        sys.exit(1)