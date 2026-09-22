""" import subprocess
from scipy.io import mmread
import numpy as np

def export_pde_matrices(mapdl_exe, work_dir):
    subprocess.run(
        [mapdl_exe, "-b", "-i", "pde_matrices_extract_in.txt", "-o", "pde_matrices_extract_out.txt"],
        cwd=work_dir, check=True)


K = mmread("K_pde.mtx").tocsr()
M = mmread("M_pde.mtx").tocsr()

print(abs(K @ np.ones(K.shape[0])).max() / abs(K).max())
print(M.sum())
 """

# We want to run:
# "C:/Program Files/ANSYS Inc/v261/ansys/bin/winx64/ANSYS261.exe" -b -i pde_matrices_extract_in.txt -o pde_matrices_extract_out.txt

#To find the path:
from pathlib import Path
hits = sorted(Path("C:/Program Files/ANSYS Inc").glob("v*/ansys/bin/winx64/ANSYS*.exe"))
print(hits)
