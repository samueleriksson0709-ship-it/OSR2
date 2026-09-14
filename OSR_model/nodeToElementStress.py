
# OBS: Only for four node tetraed

# Using:
# sigma_e = sum_i w_i * sigma_i
#    w_i = 1 / V * intergral(N_i dV)

import numpy as np

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


def n_to_e(stress_file,npe,nelm,nnod):
    
    for element in range(nelm):
        index = element*npe

        
        ex = [row[2] for row in stress_file[index:index+npe]]
        ey = [row[3] for row in stress_file[index:index+npe]]
        ez = [row[4] for row in stress_file[index:index+npe]]
        

        sigma_e += sigma_n
     


"""     for element in elements:
        nodes = 
        sigma_n = stress(nodes)
        w_n = calc_weight(nodes)
        w_list(element) = w_n
        sigma_e = sum (w_n * sigma_n)

    return sigma_e, w_list """

path = "C:/Users/samue/OSR/shared/sigma_export_nodes_testBracket.txt"
stress_file = read_stress_file(path)

npe = int(sum(1 for row in stress_file if row[0] == stress_file[0][0]))
nelm = int(len(stress_file) / npe)
assert nelm * npe == len(stress_file), f"row count {len(stress_file)} not divisible into {nelm}×{npe}"
nnod = int(max([row[1] for row in stress_file]))

n_to_e(stress_file,npe,nelm,nnod)

