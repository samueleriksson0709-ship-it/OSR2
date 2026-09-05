# -*- coding: utf-8 -*-
"""OSR "Super Solve" ACT action.

This module is executed by Workbench in IronPython 2.7, so it stays Python 2/3
compatible: no f-strings, no pathlib.

Nothing here is tied to a particular machine or user account. Paths are
resolved in this order:

  1. environment variables (OSR_ROOT, OSR_SHARED_DIR, OSR_PYTHON)
  2. the location of this file  (<root>/ansys_extension/SuperSolve/)
  3. the ACT extension install directory
  4. values stored in shared/act_config.json, when they still exist locally

The resolved values are written back to shared/act_config.json, so the config
file repairs itself the first time the extension runs on a new machine.
"""

import os
import subprocess
import json
import traceback
import time


IS_WINDOWS = (os.name == "nt")

# Cached resolution results (per Workbench session).
_SHARED_DIR = None
_PYTHON_CMD = None


# ----------------------------------------------------------------------------
# Path discovery
# ----------------------------------------------------------------------------

def _script_dir():
    """Directory containing this file, or None when it cannot be determined."""
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        pass
    try:
        # ACT injects ExtAPI as a global.
        return str(ExtAPI.ExtensionManager.CurrentExtension.InstallDir)
    except Exception:
        return None


def _looks_like_osr_root(path):
    if not path or not os.path.isdir(path):
        return False
    return (os.path.isdir(os.path.join(path, "shared")) or
            os.path.isdir(os.path.join(path, "OSR_model")))


def detect_osr_root():
    """Locate the OSR checkout without relying on a hard coded user name."""
    env_root = os.environ.get("OSR_ROOT")
    if env_root and _looks_like_osr_root(env_root):
        return os.path.abspath(env_root)

    here = _script_dir()
    if here:
        # <root>/ansys_extension/SuperSolve -> <root>
        # <root>/ansys_extension           -> <root>
        # <root>                           -> <root>
        for up in ("..%s.." % os.sep, "..", "."):
            candidate = os.path.abspath(os.path.join(here, up))
            if _looks_like_osr_root(candidate):
                return candidate

    return None


def shared_dir():
    """Absolute path of the shared folder used to exchange files with Ansys."""
    global _SHARED_DIR
    if _SHARED_DIR:
        return _SHARED_DIR

    env_shared = os.environ.get("OSR_SHARED_DIR")
    if env_shared and os.path.isdir(env_shared):
        _SHARED_DIR = os.path.abspath(env_shared)
        return _SHARED_DIR

    root = detect_osr_root()
    if root:
        candidate = os.path.join(root, "shared")
        if not os.path.isdir(candidate):
            try:
                os.makedirs(candidate)
            except Exception:
                pass
        if os.path.isdir(candidate):
            _SHARED_DIR = candidate
            return _SHARED_DIR

    raise Exception(
        "Could not locate the OSR shared folder.\n"
        "Set the OSR_ROOT environment variable to the OSR checkout "
        "(the folder containing 'shared' and 'OSR_model'), or reinstall the "
        "extension from <OSR root>/ansys_extension."
    )


def _shared_dir_quiet():
    """shared_dir() that returns None instead of raising (used by the loggers)."""
    try:
        return shared_dir()
    except Exception:
        return None


def config_file():
    return os.path.join(shared_dir(), "act_config.json")


# ----------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------

def _write_log(filename, msg):
    base = _shared_dir_quiet()
    if not base:
        return
    try:
        f = open(os.path.join(base, filename), "a")
        try:
            f.write(msg + "\n")
        finally:
            f.close()
    except Exception:
        pass


def log_debug(msg):
    _write_log("act_debug.txt", msg)


def log_error(msg):
    _write_log("act_error.txt", msg)


# ----------------------------------------------------------------------------
# Python interpreter discovery
# ----------------------------------------------------------------------------

def _probe_capture_dir():
    base = _shared_dir_quiet()
    if base:
        return base
    import tempfile
    return tempfile.gettempdir()


def _read_text(path):
    try:
        f = open(path, "r")
        try:
            return f.read()
        finally:
            f.close()
    except Exception:
        return ""


def _run(cmd):
    """Run a command, return (returncode, stdout, stderr) as text.

    This is used only to probe candidate Python interpreters. It deliberately
    avoids subprocess.Popen(...).communicate(): this file is executed by
    Workbench under IronPython 2.7, whose subprocess module does not
    implement communicate() (it raises NotImplementedError("Popen.communicate")
    for every candidate, good or bad). Capturing through plain files plus
    wait() is the pattern that already works for the pipeline subprocess
    below, so probing reuses it.
    """
    directory = _probe_capture_dir()
    tag = str(os.getpid()) + "_" + str(int(time.time() * 1000))
    out_path = os.path.join(directory, "_osr_probe_out_" + tag + ".tmp")
    err_path = os.path.join(directory, "_osr_probe_err_" + tag + ".tmp")

    try:
        out_f = open(out_path, "w")
        err_f = open(err_path, "w")
        try:
            process = subprocess.Popen(cmd, stdout=out_f, stderr=err_f)
            returncode = process.wait()
        finally:
            out_f.close()
            err_f.close()
        return returncode, _read_text(out_path), _read_text(err_path)
    except Exception as e:
        # A Microsoft Store alias with no app behind it, a missing file, ...
        return -1, "", str(e)
    finally:
        for p in (out_path, err_path):
            try:
                os.remove(p)
            except Exception:
                pass


def _is_store_alias(path):
    """True for a Microsoft Store 'app execution alias' stub.

    Those are zero byte reparse points under ...\\AppData\\Local\\Microsoft\\
    WindowsApps. os.path.exists() says True, but executing one when the Store
    app is not installed just returns 9009 ('command not recognized') without
    running anything.
    """
    if "windowsapps" not in path.replace("/", "\\").lower():
        return False
    try:
        return os.path.getsize(path) == 0
    except Exception:
        return True


def _which_all(name):
    """Every match for `name` on PATH (no shutil.which in IronPython 2.7)."""
    found = []
    exts = [""]
    if IS_WINDOWS:
        exts += os.environ.get("PATHEXT", ".EXE").split(os.pathsep)
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        if not directory:
            continue
        for ext in exts:
            candidate = os.path.join(directory.strip('"'), name + ext)
            if os.path.isfile(candidate) and candidate not in found:
                found.append(candidate)
    return found


def _python_candidates(configured):
    """Ordered list of candidate interpreter commands (each one an argv list)."""
    candidates = []
    deprioritized = []

    def add(cmd):
        if not cmd:
            return
        if cmd in candidates or cmd in deprioritized:
            return
        if len(cmd) == 1 and os.path.isabs(cmd[0]):
            if not os.path.isfile(cmd[0]):
                return
            if _is_store_alias(cmd[0]):
                deprioritized.append(cmd)
                return
        candidates.append(cmd)

    env_python = os.environ.get("OSR_PYTHON")
    if env_python:
        add([env_python])
    if configured:
        add([configured])

    # Typical per-user and system-wide CPython installs.
    roots = []
    local_app = os.environ.get("LOCALAPPDATA")
    if local_app:
        roots.append(os.path.join(local_app, "Programs", "Python"))
    for var in ("ProgramFiles", "ProgramFiles(x86)"):
        base = os.environ.get(var)
        if base:
            roots.append(base)
    roots.append("C:\\")
    for base in roots:
        if not os.path.isdir(base):
            continue
        try:
            entries = sorted(os.listdir(base), reverse=True)
        except Exception:
            continue
        for entry in entries:
            if not entry.lower().startswith("python"):
                continue
            add([os.path.join(base, entry, "python.exe")])

    for name in ("python", "python3"):
        for path in _which_all(name):
            add([path])

    # The Windows launcher resolves a real install even when PATH does not.
    if IS_WINDOWS:
        for launcher in _which_all("py"):
            add([launcher, "-3"])

    return candidates + deprioritized


def resolve_python(configured=None):
    """Return an argv prefix for a Python 3 interpreter that actually runs.

    Every candidate is probed by executing it: this is what separates a working
    interpreter from a Microsoft Store alias stub, which exists on disk but
    exits with 9009 without running anything.
    """
    global _PYTHON_CMD
    if _PYTHON_CMD:
        return _PYTHON_CMD

    tried = []
    numpy_missing = []
    candidates = _python_candidates(configured)

    # True once any candidate returns a real exit code. While this stays False
    # every failure came from the probe mechanism itself, not the interpreters.
    probe_works = False

    for cmd in candidates:
        rc, out, err = _run(cmd + ["-c", "import sys; sys.stdout.write(sys.executable)"])
        if rc >= 0:
            probe_works = True
        if rc != 0:
            tried.append("%s -> rc=%s %s" % (" ".join(cmd), rc, (err or "").strip()[:120]))
            continue

        rc_np, _, err_np = _run(cmd + ["-c", "import numpy"])
        if rc_np >= 0:
            probe_works = True
        if rc_np != 0:
            numpy_missing.append(out.strip() or " ".join(cmd))
            tried.append("%s -> numpy not installed" % " ".join(cmd))
            continue

        log_debug("Resolved python interpreter: " + " ".join(cmd))
        _PYTHON_CMD = cmd
        return _PYTHON_CMD

    if not probe_works:
        # Not one candidate produced an exit code, so launching probes is
        # broken here rather than every interpreter being unusable. Refusing to
        # run in that situation hides a working install behind a tooling
        # problem, so fall back to the best candidate that exists on disk and
        # let the pipeline itself report a genuine failure.
        for cmd in candidates:
            if len(cmd) > 1 or not os.path.isabs(cmd[0]):
                continue
            if _is_store_alias(cmd[0]) or not os.path.isfile(cmd[0]):
                continue
            log_error(
                "Could not probe any interpreter (subprocess capture is "
                "unavailable here). Falling back to: " + " ".join(cmd)
            )
            _PYTHON_CMD = cmd
            return _PYTHON_CMD

    message = ["No usable Python 3 interpreter with numpy was found."]
    if numpy_missing:
        message.append(
            "These interpreters work but have no numpy (the OSR pipeline needs it):\n  " +
            "\n  ".join(numpy_missing) +
            "\nInstall it with:  \"%s\" -m pip install numpy" % numpy_missing[0]
        )
    else:
        message.append(
            "Install Python 3 from python.org (not the Microsoft Store) and tick "
            "'Add python.exe to PATH', then run:  python -m pip install numpy"
        )
    message.append(
        "You can also point the extension at a specific interpreter by setting "
        "the OSR_PYTHON environment variable, or 'python_exe' in act_config.json."
    )
    if tried:
        message.append("Candidates tried:\n  " + "\n  ".join(tried))
    raise Exception("\n".join(message))


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

def _default_omp_threads():
    try:
        return int(os.environ.get("NUMBER_OF_PROCESSORS", "8"))
    except Exception:
        return 8


def load_config():
    """Load act_config.json, repairing any entry that does not exist locally.

    The config file is user and machine specific: a checkout copied from
    another computer still carries the previous user's absolute paths. Rather
    than failing, every stale entry is replaced by a locally valid one and the
    file is rewritten.
    """
    base = shared_dir()
    path = os.path.join(base, "act_config.json")

    cfg = {}
    if os.path.exists(path):
        try:
            f = open(path, "r")
            try:
                cfg = json.load(f)
            finally:
                f.close()
        except Exception as e:
            log_error("Could not parse " + path + ": " + str(e) + " (using defaults)")
            cfg = {}

    original = dict(cfg)
    root = detect_osr_root()

    # shared_dir: always the folder this config was loaded from.
    cfg["shared_dir"] = base

    # pipeline_script
    script = cfg.get("pipeline_script")
    if not script or not os.path.exists(script):
        if script:
            log_debug("Configured pipeline_script does not exist here: " + str(script))
        if root:
            cfg["pipeline_script"] = os.path.join(root, "OSR_model", "Run_pipeline.py")
        else:
            cfg["pipeline_script"] = ""

    # omp_threads
    try:
        cfg["omp_threads"] = int(cfg.get("omp_threads", _default_omp_threads()))
    except Exception:
        cfg["omp_threads"] = _default_omp_threads()

    # python_exe: probed, because a path that exists is not necessarily runnable.
    configured_python = cfg.get("python_exe")
    python_cmd = resolve_python(configured_python)
    cfg["python_exe"] = python_cmd[0]
    if len(python_cmd) > 1:
        cfg["python_args"] = python_cmd[1:]
    elif "python_args" in cfg:
        del cfg["python_args"]

    if cfg != original:
        try:
            f = open(path, "w")
            try:
                json.dump(cfg, f, indent=4, sort_keys=True)
                f.write("\n")
            finally:
                f.close()
            log_debug("act_config.json updated for this machine")
        except Exception as e:
            log_error("Could not write " + path + ": " + str(e))

    return cfg


def python_command(cfg):
    cmd = [cfg["python_exe"]]
    cmd.extend(cfg.get("python_args", []))
    return cmd


# ----------------------------------------------------------------------------
# Checks on the Ansys side
# ----------------------------------------------------------------------------

def _list_shared_txt(base):
    try:
        names = [n for n in sorted(os.listdir(base)) if n.lower().endswith(".txt")]
    except Exception:
        return "<could not list the shared folder>"
    if not names:
        return "<no .txt file in the shared folder>"
    return "\n  ".join(names)


def verify_stress_files_declared_exist():
    cfg = load_config()
    base = cfg["shared_dir"]
    loadcases_path = os.path.join(base, "loadcases.txt")

    if not os.path.exists(loadcases_path):
        raise Exception("loadcases.txt not found in shared folder: " + base)

    f = open(loadcases_path, "r")
    try:
        lines = [l.strip() for l in f
                 if l.strip() and not l.strip().startswith("#") and not l.strip().startswith("!")]
    finally:
        f.close()

    if len(lines) < 2:
        raise Exception("Invalid loadcases.txt: no stress file declared")

    n_cases = int(lines[0].split()[0])

    marker_path = os.path.join(base, "apdl_marker.txt")
    deadline = time.time() + 60.0

    while True:
        missing = []
        for i in range(1, 1 + n_cases):
            fname = lines[i].split()[0]
            path = os.path.join(base, fname)
            if not os.path.exists(path):
                missing.append(path)

        if not missing:
            return

        if time.time() > deadline:
            detail = [
                "Ansys stress export missing after Solution update:",
                "\n".join(missing),
            ]
            if not os.path.exists(marker_path):
                detail.append(
                    "apdl_marker.txt is missing too, so the Commands (APDL) block either "
                    "did not run or wrote its output somewhere else.\n"
                    "Open shared/APDL.txt and check that both *CFOPEN lines end with this "
                    "machine's shared folder:\n  " + base + "\n"
                    "Running configure_osr.py from the OSR root patches those paths "
                    "automatically."
                )
            else:
                detail.append(
                    "apdl_marker.txt is present, so the APDL block ran and wrote to the right "
                    "folder. The exported file name therefore does not match the name declared "
                    "in loadcases.txt.\n"
                    "*CFOPEN,<name>,txt,... in APDL.txt must produce exactly the file names "
                    "listed in loadcases.txt."
                )
            detail.append("Files currently in " + base + ":\n  " + _list_shared_txt(base))
            raise Exception("\n".join(detail))

        log_debug("Waiting for Ansys stress export files")
        time.sleep(2.0)


# ----------------------------------------------------------------------------
# ACT entry points
# ----------------------------------------------------------------------------

def ExecuteSuperSolve(ext):
    try:
        log_debug("ExecuteSuperSolve called")
        ExtAPI.Log.WriteMessage(" OSR Super Solve started")

        log_debug("Before update_solution_systems")
        update_solution_systems()
        log_debug("After update_solution_systems")

        log_debug("Before verify_stress_file")
        verify_stress_files_declared_exist()
        log_debug("After verify_stress_file")

        log_debug("Before run_osr_pipeline")
        run_osr_pipeline()
        log_debug("After run_osr_pipeline")

        log_debug("Before update_external_data_systems")
        update_external_data_systems()
        log_debug("After update_external_data_systems")

        ExtAPI.Log.WriteMessage(" Super Solve finished successfully ")
        log_debug("ExecuteSuperSolve finished successfully")

    except Exception as e:
        ExtAPI.Log.WriteError("Super Solve failed : " + str(e))
        log_error("ERROR in ExecuteSuperSolve")
        log_error(str(e))
        log_error(traceback.format_exc())


def RunOSROnly(ext):
    try:
        log_debug("RunOSROnly called")
        ExtAPI.Log.WriteMessage("Run OSR only started")

        log_debug("Before run_osr_pipeline")
        run_osr_pipeline()
        log_debug("After run_osr_pipeline")

        log_debug("Before update_external_data_systems")
        update_external_data_systems()
        log_debug("After update_external_data_systems")

        ExtAPI.Log.WriteMessage("Run OSR only finished successfully")
        log_debug("OSR only finished successfully")

    except Exception as e:
        ExtAPI.Log.WriteMessage("Run OSR only failed:" + str(e))
        log_error("ERROR in RunOSROnly")
        log_error(str(e))
        log_error(traceback.format_exc())


def RefreshDamageOnly(ext):

    try:
        log_debug("RefreshDamageOnly called")
        ExtAPI.Log.WriteMessage("Refresh damage only started")

        update_external_data_systems()

        ExtAPI.Log.WriteMessage("Refresh damage only finished successfully")
        log_debug("RefreshDamageOnly finished successfully")

    except Exception as e:
        ExtAPI.Log.WriteMessage("Refresh damage only failed:" + str(e))
        log_error("ERROR in RefreshDamageOnly")
        log_error(str(e))
        log_error(traceback.format_exc())


def _tail(path, max_chars=2000):
    try:
        f = open(path, "r")
        try:
            text = f.read().strip()
        finally:
            f.close()
    except Exception:
        return ""
    if len(text) > max_chars:
        text = "..." + text[-max_chars:]
    return text


def run_osr_pipeline():

    cfg = load_config()
    python_cmd = python_command(cfg)
    script = cfg["pipeline_script"]
    base = cfg["shared_dir"]
    omp_threads = str(cfg["omp_threads"])

    log_debug("Entering run_osr_pipeline")
    log_debug("python_cmd = " + " ".join(python_cmd))
    log_debug("script = " + script)
    log_debug("shared_dir = " + base)
    log_debug("omp_threads = " + omp_threads)

    if not script or not os.path.exists(script):
        raise Exception(
            "Run_pipeline.py not found: " + str(script) + "\n"
            "Set OSR_ROOT to the OSR checkout, or fix 'pipeline_script' in act_config.json."
        )

    if not os.path.exists(base):
        raise Exception("shared directory not found : " + base)

    fortran_dir = os.path.dirname(script)
    exe = os.path.join(fortran_dir, "osr_pipeline_multi.exe")
    if not os.path.exists(exe) and not os.path.exists(exe[:-4]):
        raise Exception(
            "The Fortran solver has not been built on this machine:\n  " + exe + "\n"
            "Open a terminal in " + fortran_dir + " and run 'make' (gfortran with OpenMP required)."
        )

    stdout_path = os.path.join(base, "osr_stdout.txt")
    stderr_path = os.path.join(base, "osr_stderr.txt")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = omp_threads
    # Let Run_pipeline.py find the same folders without hard coded paths.
    env["OSR_SHARED_DIR"] = base
    env["OSR_FORTRAN_DIR"] = fortran_dir
    root = detect_osr_root()
    if root:
        env["OSR_ROOT"] = root

    stdout_file = open(stdout_path, "w")
    stderr_file = open(stderr_path, "w")

    try:

        log_debug("launching OSR pipeline subprocess")

        process = subprocess.Popen(
            python_cmd + [script],
            cwd=fortran_dir,
            env=env,
            stdout=stdout_file,
            stderr=stderr_file
        )

        returncode = process.wait()

    finally:
        stdout_file.close()
        stderr_file.close()

    log_debug("Subprocess finished")
    log_debug("returncode = " + str(returncode))

    if returncode != 0:
        message = ["OSR pipeline failed with return code " + str(returncode) + "."]
        if returncode == 9009:
            message.append(
                "9009 means Windows could not start the interpreter at all - usually a "
                "Microsoft Store 'app execution alias' (…\\AppData\\Local\\Microsoft\\"
                "WindowsApps\\python.exe) rather than a real Python install.\n"
                "Install Python 3 from python.org and delete 'python_exe' from "
                "act_config.json so it is detected again."
            )
        err = _tail(stderr_path)
        if err:
            message.append("osr_stderr.txt:\n" + err)
        else:
            message.append("osr_stderr.txt is empty. See pipeline_log.txt in " + base + ".")
        raise Exception("\n".join(message))

    log_debug("Leaving run_osr_pipeline successfully")


def update_solution_systems():
    log_debug("Entering update_solution_systems")

    for system in GetAllSystems():
        try:
            name = str(system.DisplayText)
        except:
            name = "<unknown>"

        if "External Data" in name or "External" in name:
            log_debug("Skipping external system: " + name)
            continue

        # Try to get the solution component 
        try:
            solution = system.GetComponent(Name="Solution")
            log_debug("Found solution component in system: " + name)
        except Exception as e : 
            log_debug("skipping system without solution component: " + name)
            log_debug("Reason: " + str(e))
            continue
        
        # try to update the solution component directly 
        try : 
            log_debug("Updating solution component directly : " + name)
            solution.Update()
            log_debug("Update solution component : " + name)
            continue
        except Exception as e:
            log_debug("Direct Solution update failed for : " + name)
            log_error(str(e))

        # fallback : update whole system 
        try : 
            log_debug("fallback : updating the whole system " + name)
            system.Update()
            log_debug("Update whole system " + name)
        except Exception as e : 
            log_error("System update failed for: " + name)
            log_error(str(e))
            raise
    
    log_debug("Leaving update_solution_systems")


def update_external_data_systems():
    log_debug("Entering update_external_data_systems")

    for system in GetAllSystems():
        try:
            name = str(system.DisplayText)
        except:
            name = "<unknown>"

        if "External Data" in name or "External" in name:
            log_debug("Trying to update External Data system: " + name)

            # 1) Nettoyer le cache
            try:
                setup = system.GetComponent(Name="Setup")
                setup.Clean()
                log_debug("Cleaned setup for: " + name)
            except Exception as e:
                log_error("Could not clean External Data setup: " + name)
                log_error(str(e))

            # 2) Essayer update du Setup
            try:
                setup = system.GetComponent(Name="Setup")
                setup.Update()
                log_debug("Updated External Data setup: " + name)
                continue
            except Exception as e:
                log_error("Could not update External Data setup directly: " + name)
                log_error(str(e))

            # 3) Fallback : update du système External Data
            try:
                system.Update()
                log_debug("Updated External Data system: " + name)
            except Exception as e:
                log_error("Could not update External Data system: " + name)
                log_error(str(e))

    log_debug("Leaving update_external_data_systems")
