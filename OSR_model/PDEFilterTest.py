
import numpy as np
from itertools import permutations
from Run_pipeline import (read_stress_file, parse_loadcases_for_filenames, parse_loadcases, PDEFilter,
                          constructPDEMatrix, SHARED_DIR, INPUT_LOADCASES)


def stresses(rows):
    return np.array([row[5:11] for row in rows])


def with_stress(rows, S):
    return [row[:5] + tuple(s) for row, s in zip(rows, S.tolist())]


def nodal_from_elemental(elem):
    first = {}
    for row in elem:
        first.setdefault(row[1], row)
    return [first[n] for n in sorted(first)]


def rel_err(a, b):
    return np.abs(a - b).max() / np.abs(b).max()


def test_real_mesh(r):
    fname = parse_loadcases_for_filenames(SHARED_DIR / INPUT_LOADCASES)[0]
    elem = read_stress_file(SHARED_DIR / fname)
    nodal = nodal_from_elemental(elem)
    nnod = len(nodal)

    C = np.tile(np.arange(1.0, 7.0), (nnod, 1))
    e_const = rel_err(stresses(PDEFilter(with_stress(nodal, C), elem, r)), C)

    S = stresses(nodal)
    S_f = stresses(PDEFilter(nodal, elem, r))
    _, M, _ = constructPDEMatrix(elem, r, nnod, len(elem) // 4, 4)
    w = M @ np.ones(nnod)
    e_int = rel_err(w @ S_f, w @ S)

    e_id = rel_err(stresses(PDEFilter(nodal, elem, 1e-6 * r)), S)

    print(f"real mesh: constant {e_const:.1e}, integral {e_int:.1e}, identity {e_id:.1e}")
    assert e_const < 1e-10 and e_int < 1e-10 and e_id < 1e-6


def box_mesh(n, L=1.0):
    g = np.linspace(0.0, L, n + 1)
    nid = lambda i, j, k: 1 + i + (n + 1) * (j + (n + 1) * k)
    X = {nid(i, j, k): (g[i], g[j], g[k])
         for i in range(n + 1) for j in range(n + 1) for k in range(n + 1)}
    elem, eid = [], 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for perm in permutations(range(3)):
                    p = [i, j, k]
                    path = [nid(*p)]
                    for ax in perm:
                        p[ax] += 1
                        path.append(nid(*p))
                    eid += 1
                    elem += [(eid, m) + X[m] + (0.0,) * 6 for m in path]
    nodal = [(0, m) + X[m] + (0.0,) * 6 for m in sorted(X)]
    return elem, nodal, np.array([X[m] for m in sorted(X)])


def test_box(n, L=1.0, l_0=0.1):
    elem, nodal, X = box_mesh(n, L)
    x, y, z = X.T
    k1, k2 = 3 * np.pi / L, 2 * np.pi / L

    S = np.zeros((len(nodal), 6))
    S[:, 0] = np.cos(k1 * x)
    S[:, 1] = np.cos(k2 * y) * np.cos(k2 * z)
    S[:, 2] = x

    exact = np.zeros_like(S)
    exact[:, 0] = S[:, 0] / (1 + l_0**2 * k1**2)
    exact[:, 1] = S[:, 1] / (1 + 2 * l_0**2 * k2**2)
    exact[:, 2] = x + l_0 * (np.cosh((L - x) / l_0) - np.cosh(x / l_0)) / np.sinh(L / l_0)

    S_f = stresses(PDEFilter(with_stress(nodal, S), elem, 2 * np.sqrt(2) * l_0))
    err = np.abs(S_f - exact).max(axis=0)[:3]
    top = np.isclose(x, L)
    print(f"box n={n:3d}: cos {err[0]:.2e}, cos*cos {err[1]:.2e}, linear {err[2]:.2e}, "
          f"x at x=L after filter {S_f[top, 2].mean():.4f} (exact {exact[top, 2].mean():.4f})")
    return err

if __name__ == "__main__":
    pde_radius = parse_loadcases(SHARED_DIR / INPUT_LOADCASES)[-1]
    assert pde_radius is not None, "no filtRad on the header line of loadcases.txt"
    print(f"filter radius from loadcases.txt: {pde_radius}")

    test_real_mesh(r=pde_radius)
    errs = [test_box(n) for n in (8, 16, 24)]
    print("rates:", np.log(errs[0] / errs[1]) / np.log(2), np.log(errs[1] / errs[2]) / np.log(1.5))
