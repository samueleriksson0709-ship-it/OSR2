import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import os

def parse_line(line):
    """One row of the export, or None for headers, blanks and malformed rows.
    Kept byte-compatible with Run_pipeline.parse_line so both readers agree on
    what counts as data.
    """
    line = line.strip()
    if not line or line[0].isalpha() or line[0] == "#":
        return None
    parts = line.split()
    if len(parts) < 11:
        return None
    try:
        return (int(float(parts[0])), int(float(parts[1]))) + tuple(
            float(p) for p in parts[2:11])
    except (ValueError, IndexError):
        return None

def read_stress_file(path):
    with open(path, "r") as f:
        return [r for r in (parse_line(line) for line in f) if r is not None]


def cont3D4kmt(ex,ey,ez):
    C = np.array([
        [1 , ex[0], ey[0], ez[0]],
        [1 , ex[1], ey[1], ez[1]],
        [1 , ex[2], ey[2], ez[2]],
        [1 , ex[3], ey[3], ez[3]]
    ])
    Cinv = np.linalg.inv(C)
    V = np.abs(np.linalg.det(C)/6)

    B = Cinv[1:4,:]

    Ke = B.T @ B * V
    Me = V / 20 * (np.ones((4, 4)) + np.eye(4))
    Te = V / 4 * np.ones(4)

    return Ke,Me,Te, V

def constructPDEMatrix(stress_file,filter_element_frac,nnod,nelm,npe):

    rows, cols, dataK, dataM = [], [], [], []

    #Construct matrices
    K = np.zeros([nnod,nnod])
    M = np.zeros([nnod,nnod])
    T = np.zeros([nnod,nelm])

    V_tot = 0.0

    for element in range(nelm):
        index = element*npe

        nodes = [int(row[1])-1 for row in stress_file[index:index+npe]]

        ex = [row[2] for row in stress_file[index:index+npe]]
        ey = [row[3] for row in stress_file[index:index+npe]]
        ez = [row[4] for row in stress_file[index:index+npe]]
        
        Ke,Me,Te,Ve = cont3D4kmt(ex,ey,ez)

        V_tot += Ve

        for a in range(npe):
            for b in range(npe):
                rows.append(nodes[a]); cols.append(nodes[b])
                dataK.append(Ke[a, b]); dataM.append(Me[a, b])

    K = sparse.coo_matrix((dataK, (rows, cols)), shape=(nnod, nnod)).tocsr()
    M = sparse.coo_matrix((dataM, (rows, cols)), shape=(nnod, nnod)).tocsr()

    #Determine l_0
    L_charac = np.cbrt(V_tot)
    r = filter_element_frac * L_charac
    l_0 = r / (2*np.sqrt(2))

    h = np.cbrt(V_tot / nelm)
    if l_0 < 2 * h:
        print(f"Warning: l_0={l_0:.4g} under-resolved (h={h:.4g}); refine mesh or raise feature_frac")

    return K,M,l_0

path = "C:/Users/samue/OSR/shared/sigma_export_nodes_testBracket.txt"
stress_file = read_stress_file(path)
path = "C:/Users/samue/OSR/OSR_model/sigma_export_nodes_testBracket.txt"
stress_nodal_file = read_stress_file(path)

stress = np.array([row[5:11] for row in stress_file])
nodes = np.array([row[1] for row in stress_file])

""" ids = np.array([r[1] for r in stress_file], dtype=int)
print(ids.min(), ids.max(), len(np.unique(ids))) """

filter_element_frac = 0.5 #Choose approprietly

#Find geometry variables
npe = int(sum(1 for row in stress_file if row[0] == stress_file[0][0]))
nelm = int(len(stress_file) / npe)
assert nelm * npe == len(stress_file), f"row count {len(stress_file)} not divisible into {nelm}×{npe}"
nnod = int(max([row[1] for row in stress_file]))


K,M,l_0 = constructPDEMatrix(stress_file,filter_element_frac,nnod,nelm,npe)

# collapse duplicated per-element rows into one value per node
stress_nodal = np.zeros((nnod, 6))
count = np.zeros(nnod)
for r in stress_file:
    nid = int(r[1]) - 1
    stress_nodal[nid] += r[5:11]


# solve once (nodal)
A_pde = (l_0**2 * K + M).tocsc()
stress_filt_nodal = spsolve(A_pde, M @ stress_nodal)



# scatter back to original per-row structure, same order as stress_file
out = np.empty((len(stress_file), 6))
for k, r in enumerate(stress_file):
    nid = int(r[1]) - 1
    print(nid)
    out[k] = stress_filt_nodal[nid]

# write with identical columns: eid nid X Y Z SXX..SXZ
base, ext = os.path.splitext(path)
out_path = base + "_filtered" + ext
with open(out_path, "w") as f:
    f.write("# eid nid X Y Z SXX SYY SZZ SXY SYZ SXZ\n")
    for k, r in enumerate(stress_file):
        eid, nid = r[0], r[1]
        x, y, z = r[2], r[3], r[4]
        sxx, syy, szz, sxy, syz, sxz = out[k]
        f.write(f"{eid:>8.0f}. {nid:>7.0f}. "
                f"{x: .6E} {y: .6E} {z: .6E} "
                f"{sxx: .6E} {syy: .6E} {szz: .6E} "
                f"{sxy: .6E} {syz: .6E} {sxz: .6E}\n")





