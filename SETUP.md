# Running OSR on a new machine

The project glues together three programs that each need an absolute path:
Workbench (ACT extension), the APDL command block inside Mechanical, and the
Fortran solver. Those paths used to be written into the sources for one
specific user account, which is why a fresh clone failed on a second computer.

Everything except the APDL block is now derived from the location of the
checkout, so the only per-machine step is running the configuration script.

## 1. Prerequisites

* **Python 3** installed **from python.org**, with `numpy`:

  ```
  python -m pip install numpy
  ```

  Do not rely on the Microsoft Store version. `...\AppData\Local\Microsoft\
  WindowsApps\python.exe` is an *app execution alias*: a zero byte stub that
  exists on disk but exits with code **9009** without ever starting an
  interpreter when Workbench launches it.

* **gfortran** with OpenMP (MSYS2 / MinGW-w64 on Windows).

## 2. Build the solver

```
cd OSR_model
make
```

This produces `osr_pipeline_multi.exe`. The binary is not in the repository, so
it has to be rebuilt on every machine.

## 3. Configure the paths

From the root of the checkout:

```
python configure_osr.py
```

It writes `shared/act_config.json`, picks a Python interpreter that actually
runs and has numpy, rewrites the export path inside `shared/APDL.txt`, and
reports anything still missing. Use `python configure_osr.py --check` to see
what it would do without touching any file.

## 4. Update the APDL block in Mechanical

`configure_osr.py` fixes the file, but Mechanical holds its own copy: open the
**Commands (APDL)** object under *Solution* and paste in the contents of
`shared/APDL.txt` again. Both `*CFOPEN` lines must end with this machine's
`shared` folder.

The exported file names must match the names declared in
`shared/loadcases.txt`; `configure_osr.py` checks this and reports a mismatch.

## 4b. Declare the unit of the export

The OSR material parameters (`sigma_oe = 490 MPa`) are calibrated in **MPa**, but
APDL exports the stresses in whatever unit system the model was *solved* in. A
model solved in MKS - the Workbench default, `/units,MKS` at the top of `ds.dat`,
lengths in m and `MP,EX` in Pa - exports **Pa**, so its numbers are a million
times too large for the solver. Under such a field every point sits far past the
endurance surface and the damage rate `K*exp(L*beta)` saturates, so the whole
model reports `D = 1` and `N_failure = 0` no matter how small the applied load.

The unit is therefore declared as an optional 4th token on the first data line of
`shared/loadcases.txt`:

```
# n_cases omega mode stress_unit
1 6.283185307   SINUS   PA
```

`PA`, `KPA`, `MPA`, `GPA`, `PSI`, `KSI` or an explicit factor to MPa are
accepted. Omitting the token means `MPA`, which is what a model solved in the
mm-kg-N system exports. `python configure_osr.py --check` compares the declared
unit against the magnitude of the exported file and reports a mismatch; the
solver refuses to integrate a field that cannot be in the declared unit.

## 5. Install the extension

In Workbench: *Extensions > Install Extension*, pick
`ansys_extension/SuperSolve.xml`, then load it. The toolbar
"Automatisation OSR" appears in the Project page.

## Overriding the automatic paths

| Variable | Meaning |
|---|---|
| `OSR_ROOT` | root of the checkout (the folder holding `shared` and `OSR_model`) |
| `OSR_SHARED_DIR` | shared exchange folder |
| `OSR_FORTRAN_DIR` | folder holding the solver executable |
| `OSR_PYTHON` | interpreter used to run `Run_pipeline.py` |

The ACT extension sets the first four automatically for the pipeline
subprocess; set them yourself only for a non-standard layout.

## Troubleshooting

Logs are written to the shared folder: `act_debug.txt`, `act_error.txt`,
`osr_stdout.txt`, `osr_stderr.txt`, `pipeline_log.txt`.

| Symptom | Cause |
|---|---|
| `return code 9009` | Python could not be started - Store alias, or a `python_exe` from another machine. Delete `python_exe` from `act_config.json` and re-run `configure_osr.py`. |
| `Ansys stress export missing`, `apdl_marker.txt` absent | The APDL block did not run, or its `*CFOPEN` path points at another machine's folder. |
| `Ansys stress export missing`, `apdl_marker.txt` present | The APDL block ran, but the exported file names differ from those in `loadcases.txt`. |
| `Fortran solver not built` | Run `make` in `OSR_model`. |
| `D = 1` and `N_failure = 0` at every point, whatever the load | The stress unit declared in `loadcases.txt` does not match the export. A model solved in MKS exports Pa: declare `PA` (see section 4b). The solver now stops with "the stress field is not consistent with the material parameters" instead of returning this field. |
| `Error: Line truncated` while running `make` | An old checkout without `-ffree-line-length-none` in `OSR_model/Makefile`. |
