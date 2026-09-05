#!/usr/bin/env python3
"""One-shot setup of the OSR toolchain on a new machine.

The project exchanges files between three programs that each need an absolute
path (Workbench ACT, the APDL command block, the Fortran solver). Those paths
are machine specific, so run this script once after cloning the repository:

    python configure_osr.py

It:
  * writes shared/act_config.json with paths valid on this machine,
  * picks a Python 3 interpreter that really runs and has numpy,
  * rewrites the *CFOPEN export paths in shared/APDL.txt,
  * reports anything still missing (Fortran build, file name mismatches).

Use --check to inspect without writing anything.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

OSR_ROOT = Path(__file__).resolve().parent
SHARED_DIR = OSR_ROOT / "shared"
FORTRAN_DIR = OSR_ROOT / "OSR_model"
CONFIG_FILE = SHARED_DIR / "act_config.json"
APDL_FILE = SHARED_DIR / "APDL.txt"
EXE_NAME = "osr_pipeline_multi.exe" if os.name == "nt" else "osr_pipeline_multi"

PLACEHOLDER = "__OSR_SHARED_DIR__"

# The OSR material parameters are calibrated in MPa. Ansys exports in whatever
# unit system the model was solved in, so loadcases.txt declares the unit of the
# export as an optional 4th token of its first data line. Keep this table in step
# with STRESS_UNITS in OSR_model/Run_pipeline.py and set_stress_unit in
# OSR_model/osr_pipeline_multi.f90.
STRESS_UNITS = {
    "MPA": 1.0, "N/MM2": 1.0, "NMM": 1.0,
    "PA": 1.0e-6, "N/M2": 1.0e-6, "MKS": 1.0e-6, "SI": 1.0e-6,
    "KPA": 1.0e-3,
    "GPA": 1.0e3,
    "PSI": 6.894757e-3,
    "KSI": 6.894757,
}
DEFAULT_STRESS_UNIT = "MPA"
SIGMA_OE = 490.0          # endurance limit of the material, in MPa
MAX_PLAUSIBLE_MPA = 100.0 * SIGMA_OE   # ~49 GPa: past any real fracture stress
STRESS_SAMPLE_ROWS = 20000

problems = []
notes = []


def ok(msg):
    print("  [OK]      " + msg)


def warn(msg):
    print("  [WARNING] " + msg)
    notes.append(msg)


def error(msg):
    print("  [PROBLEM] " + msg)
    problems.append(msg)


# ---------------------------------------------------------------------------
# Python interpreter
# ---------------------------------------------------------------------------

def probe(cmd, code):
    try:
        result = subprocess.run(cmd + ["-c", code],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        return result.returncode, result.stdout.decode("utf-8", "replace").strip()
    except Exception as exc:
        return -1, str(exc)


def is_store_alias(path):
    """Microsoft Store app-execution aliases are zero-byte stubs.

    They exist as far as os.path.exists is concerned but exit with code 9009
    without starting an interpreter, which is exactly how this project failed
    when it was moved to a second machine.
    """
    path = Path(path)
    if "windowsapps" not in str(path).lower():
        return False
    try:
        return path.stat().st_size == 0
    except OSError:
        return True


def candidate_interpreters():
    seen = []

    def add(cmd):
        if cmd and cmd not in seen:
            seen.append(cmd)

    add([sys.executable])

    env_python = os.environ.get("OSR_PYTHON")
    if env_python:
        add([env_python])

    if os.name == "nt":
        roots = []
        local_app = os.environ.get("LOCALAPPDATA")
        if local_app:
            roots.append(Path(local_app) / "Programs" / "Python")
        for var in ("ProgramFiles", "ProgramFiles(x86)"):
            if os.environ.get(var):
                roots.append(Path(os.environ[var]))
        roots.append(Path("C:/"))
        for root in roots:
            if not root.is_dir():
                continue
            try:
                entries = sorted(root.iterdir(), reverse=True)
            except OSError:
                continue
            for entry in entries:
                if entry.name.lower().startswith("python") and (entry / "python.exe").is_file():
                    add([str(entry / "python.exe")])

    import shutil
    for name in ("python", "python3"):
        found = shutil.which(name)
        if found:
            add([found])
    if os.name == "nt":
        launcher = shutil.which("py")
        if launcher:
            add([launcher, "-3"])

    # Store aliases last: they may work, but a real install is always better.
    real = [c for c in seen if not (len(c) == 1 and is_store_alias(c[0]))]
    alias = [c for c in seen if c not in real]
    return real + alias


def resolve_python():
    """Return (interpreter argv, numpy_ok).

    An interpreter is only accepted once it has actually been executed: a path
    that exists is not necessarily runnable.
    """
    fallback = None
    for cmd in candidate_interpreters():
        rc, out = probe(cmd, "import sys; sys.stdout.write(sys.executable)")
        if rc != 0:
            continue
        rc_np, _ = probe(cmd, "import numpy")
        if rc_np == 0:
            return cmd, True
        if fallback is None:
            fallback = cmd
    if fallback is not None:
        return fallback, False
    return None, False


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------

def step_python(check_only):
    print("\nPython interpreter")
    cmd, numpy_ok = resolve_python()
    if cmd is None:
        error("No working Python 3 interpreter found. Install one from "
              "python.org (not the Microsoft Store) and tick "
              "'Add python.exe to PATH'.")
        return None
    ok(" ".join(cmd))
    if not numpy_ok:
        error('numpy is not installed for this interpreter. Run:\n'
              '              "%s" -m pip install numpy' % cmd[0])
    if len(cmd) == 1 and is_store_alias(cmd[0]):
        warn("This is a Microsoft Store alias. It works right now, but Workbench "
             "may fail to launch it (return code 9009). A python.org install is "
             "more reliable.")
    return cmd


def step_config(python_cmd, check_only):
    print("\nshared/act_config.json")
    if python_cmd is None:
        warn("skipped: no usable Python interpreter was found.")
        return

    cfg = {
        "python_exe": python_cmd[0],
        "pipeline_script": str(FORTRAN_DIR / "Run_pipeline.py"),
        "shared_dir": str(SHARED_DIR),
        "omp_threads": os.cpu_count() or 8,
    }
    if len(python_cmd) > 1:
        cfg["python_args"] = python_cmd[1:]

    if check_only:
        print(json.dumps(cfg, indent=4, sort_keys=True))
        return

    SHARED_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_FILE, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=4, sort_keys=True)
        f.write("\n")
    ok("written for this machine (%d OpenMP threads)" % cfg["omp_threads"])


def apdl_export_names(text):
    """File names produced by the *CFOPEN lines of the APDL block."""
    names = []
    for match in re.finditer(r"^\s*\*CFOPEN\s*,\s*([^,\s]+)\s*,\s*([^,\s]+)",
                             text, re.MULTILINE | re.IGNORECASE):
        names.append("%s.%s" % (match.group(1), match.group(2)))
    return names


def step_apdl(check_only):
    print("\nshared/APDL.txt")
    if not APDL_FILE.exists():
        error("not found: %s" % APDL_FILE)
        return

    text = APDL_FILE.read_text(encoding="utf-8", errors="replace")
    target = str(SHARED_DIR)

    # Replace the placeholder and any path left over from another machine.
    patched = text.replace(PLACEHOLDER, target)
    patched = re.sub(
        r"(^\s*\*CFOPEN\s*,[^,\n]+,[^,\n]+,\s*)([^\r\n]+)",
        lambda m: m.group(1) + target,
        patched,
        flags=re.MULTILINE | re.IGNORECASE,
    )

    if patched == text:
        ok("export path already points at %s" % target)
    elif check_only:
        warn("export path would be rewritten to %s" % target)
    else:
        APDL_FILE.write_text(patched, encoding="utf-8", newline="\r\n")
        ok("export path set to %s" % target)

    notes.append(
        "Copy the contents of shared/APDL.txt into the Commands (APDL) object "
        "under Solution in Mechanical - editing the file alone does not update "
        "the model."
    )
    return apdl_export_names(patched)


def loadcases_lines():
    loadcases = SHARED_DIR / "loadcases.txt"
    if not loadcases.exists():
        return None
    return [l.strip() for l in loadcases.read_text(encoding="utf-8",
                                                   errors="replace").splitlines()
            if l.strip() and l.strip()[0] not in "#!"]


def declared_stress_files():
    lines = loadcases_lines()
    if lines is None:
        return None
    if len(lines) < 2:
        return []
    try:
        n_cases = int(lines[0].split()[0])
    except (ValueError, IndexError):
        return []
    return [lines[i].split()[0] for i in range(1, min(1 + n_cases, len(lines)))]


def declared_stress_unit():
    """(token, factor to MPa) declared in loadcases.txt, or (None, None)."""
    lines = loadcases_lines()
    if not lines:
        return None, None
    header = lines[0].split()
    token = header[3] if len(header) >= 4 else DEFAULT_STRESS_UNIT
    factor = STRESS_UNITS.get(token.upper())
    if factor is None:
        try:
            factor = float(token)
        except ValueError:
            factor = None
    return token, factor


def peak_stress_in_file(path):
    """Largest stress component of an export, or None if nothing could be read."""
    peak = None
    with open(path, "r", errors="replace") as f:
        for count, line in enumerate(f):
            if count >= STRESS_SAMPLE_ROWS:
                break
            parts = line.split()
            if len(parts) < 11 or line.lstrip()[:1] in ("#", ""):
                continue
            try:
                values = [abs(float(v)) for v in parts[5:11]]
            except ValueError:
                continue
            top = max(values)
            if peak is None or top > peak:
                peak = top
    return peak


def step_stress_unit():
    """The unit mismatch that makes every point fail at zero cycle.

    The damage rate is K*exp(L*beta), so a field exported in Pa but read as MPa
    puts every point far past the endurance surface and the solver reports
    D = 1, N_failure = 0 over the whole model.
    """
    print("\nStress unit declared in loadcases.txt")
    token, factor = declared_stress_unit()
    if token is None:
        error("shared/loadcases.txt not found")
        return
    if factor is None:
        error("unknown stress unit '%s'. Use one of %s, or an explicit factor to MPa."
              % (token, ", ".join(sorted(STRESS_UNITS))))
        return
    ok("%s (x%g -> MPa)" % (token, factor))

    declared = declared_stress_files() or []
    for name in declared:
        path = SHARED_DIR / name
        if not path.is_file():
            warn("%s not exported yet, unit not verified against real data." % name)
            continue
        peak = peak_stress_in_file(path)
        if peak is None:
            warn("%s holds no readable stress row." % name)
            continue
        peak_mpa = peak * factor
        if peak_mpa > MAX_PLAUSIBLE_MPA:
            header = (loadcases_lines() or [""])[0].split()[:3]
            suggestion = " ".join((header or ["1", "6.283185307", "SINUS"]) + ["PA"])
            error("%s peaks at %.3g MPa once converted, which no material sustains.\n"
                  "              The export is almost certainly in Pa (an Ansys model "
                  "solved in MKS) while loadcases.txt declares %s.\n"
                  "              Every point would report D = 1 and N_failure = 0. "
                  "Put the unit on the first data line of loadcases.txt:\n"
                  "                  %s"
                  % (name, peak_mpa, token, suggestion))
        else:
            ok("%s peaks at %.3g MPa, consistent with sigma_oe = %g MPa"
               % (name, peak_mpa, SIGMA_OE))


def step_names(exported):
    print("\nExported file names vs loadcases.txt")
    declared = declared_stress_files()
    if declared is None:
        error("shared/loadcases.txt not found")
        return
    if not exported:
        warn("no *CFOPEN found in APDL.txt")
        return

    exported_stress = [n for n in exported if not n.startswith("apdl_marker")]
    missing = [n for n in declared if n not in exported_stress]
    if missing:
        error("loadcases.txt declares files the APDL block never writes:\n"
              "              declared: %s\n"
              "              exported: %s\n"
              "              Super Solve waits 60 s for these files and then fails."
              % (", ".join(declared), ", ".join(exported_stress) or "<none>"))
    else:
        ok("%s" % ", ".join(declared))


def step_solver():
    print("\nFortran solver")
    exe = FORTRAN_DIR / EXE_NAME
    if exe.exists():
        ok(str(exe))
    else:
        error("%s not built. Open a terminal in %s and run 'make' "
              "(gfortran with OpenMP required)." % (EXE_NAME, FORTRAN_DIR))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true",
                        help="report only, do not modify any file")
    args = parser.parse_args()

    print("OSR root: %s" % OSR_ROOT)

    python_cmd = step_python(args.check)
    step_config(python_cmd, args.check)
    exported = step_apdl(args.check)
    step_names(exported or [])
    step_stress_unit()
    step_solver()

    print("\n" + "-" * 70)
    if problems:
        print("%d problem(s) still to fix:" % len(problems))
        for item in problems:
            print("  - " + item)
    else:
        print("Configuration complete.")
    if notes:
        print("\nRemember:")
        for item in notes:
            print("  - " + item)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
