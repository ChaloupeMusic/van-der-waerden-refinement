#!/usr/bin/env python3
"""Optimistic ceilings of candidate templates outside the current one.

Base commit: 63bef9ca996786a41169de4a0c3d6b3b47b0fe96. Standard library only.
Run:  python3 -B verify.py          (assertions are load-bearing; never -O)
      python3 -B verify.py --write  (regenerates RESULT.json)

Every template here is "base coloring c0 + rare perturbation of rate p = k^-b".
To first order in the exponent (units of k log k) a template is a system of
linear inequalities in exponents; sup c is computed exactly by vertex
enumeration in Fractions. The derivations are in README.md; this file checks
the arithmetic of each optimistic system, plus negative controls (systems with
a constraint deliberately dropped must show the gain, so the solver sees gains).
"""
from fractions import Fraction as Q
from itertools import combinations
import json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

def sup(objective, cons, nvar, box=Q(1)):
    """max <objective, x> over {x : <alpha, x> <= rhs for (alpha, rhs) in cons,
    0 <= x_i <= box}, by enumerating vertices (exact)."""
    cons = list(cons)
    for i in range(nvar):
        e = [Q(0)] * nvar; e[i] = Q(-1); cons.append((e, Q(0)))
        e = [Q(0)] * nvar; e[i] = Q(1); cons.append((e, box))
    best = None
    for rows in combinations(cons, nvar):
        M = [list(r[0]) for r in rows]; rhs = [r[1] for r in rows]
        # Gaussian elimination in Fractions
        A = [M[i] + [rhs[i]] for i in range(nvar)]
        ok = True
        for c in range(nvar):
            p = next((i for i in range(c, nvar) if A[i][c] != 0), None)
            if p is None: ok = False; break
            A[c], A[p] = A[p], A[c]
            for i in range(nvar):
                if i != c and A[i][c] != 0:
                    f = A[i][c] / A[c][c]
                    A[i] = [x - f * y for x, y in zip(A[i], A[c])]
        if not ok: continue
        x = [A[i][nvar] / A[i][i] for i in range(nvar)]
        if all(sum(al * xi for al, xi in zip(alpha, x)) <= r for alpha, r in cons):
            v = sum(o * xi for o, xi in zip(objective, x))
            best = v if best is None or v > best else best
    return best

H = Q(1, 2)
# variables x = (a, b, c) unless stated; gamma = 1/2 fixed (balance optimum)
def current(drop=()):
    cons = []
    if "B" not in drop: cons.append(([Q(-1), Q(1), Q(0)], Q(0)))          # b <= a      (bands, P1)
    if "R" not in drop: cons.append(([Q(1), Q(1), Q(0)], Q(1)))           # a + b <= 1  (boxes/runs, P3)
    if "F" not in drop: cons.append(([Q(1), Q(0), Q(0)], H))              # a <= 1/2    (labels, F)
    cons.append(([Q(0), Q(-H), Q(2)], Q(0)))                              # 2c <= gamma b
    return sup([0, 0, 1], cons, 3)

def asymmetric(theta, structured_ones=True):
    """x = (a, b0, b1, c). Base with a fraction theta of ones on rich P.
    Rich mono-0 flips the ones (cost theta*b1), mono-1 the zeros ((1-theta)*b0)."""
    cons = []
    for bi in (1, 2):
        if bi == 2 and not structured_ones:
            continue                                                      # b1 unconstrained
        e = [Q(0)] * 4; e[bi] = Q(1); e[0] = Q(-1); cons.append((e, Q(0)))      # b_i <= a
        e = [Q(0)] * 4; e[bi] = Q(1); e[0] = Q(1); cons.append((e, Q(1)))       # a + b_i <= 1
    cons.append(([Q(1), 0, 0, 0], H))                                     # a <= 1/2
    if structured_ones:                                                   # else p1 = 0: no flip on ones,
        cons.append(([0, Q(0), -theta, Q(2)], Q(0)))                      # mono-0 impossible; 2c <= theta b1
    cons.append(([0, -(1 - theta), Q(0), Q(2)], Q(0)))                    # 2c <= (1-theta) b0
    return sup([0, 0, 0, 1], cons, 4, box=Q(1))

def fine_coarse():
    """x = (a, bc, s, c): flips of rate k^-bc on a set Sigma of density k^-s.
    Structured progressions see rate k^-(bc+s) (density averaging), rich ones too."""
    cons = [([Q(-1), Q(1), Q(1), Q(0)], Q(0)),            # bc + s <= a
            ([Q(1), Q(1), Q(1), Q(0)], Q(1)),             # a + bc + s <= 1
            ([Q(1), 0, 0, 0], H),                         # a <= 1/2
            ([Q(0), -H, -H, Q(2)], Q(0))]                 # 2c <= gamma (bc + s)
    return sup([0, 0, 0, 1], cons, 4)

def uncertainty():
    """Structured-event entropy exponent max(1-a, a) (bands log N/D vs boxes D log k):
    its minimum over a is the least possible exponent of H, whence b <= 1/2."""
    grid = [Q(i, 1000) for i in range(1001)]
    m = min(max(1 - a, a) for a in grid)
    return m

def group_design(a):
    """One flip per group needs #groups >= #distinct slow words / O(1):
    exponents (units of k log k * 1/k^{...}): log #keys ~ 2 log q + 3D log k,
    log #slow words ~ 6 log q. Feasible iff 6 log q <= 2 log q + 3 D log k,
    i.e. (in powers of k) 1 - a <= a. Returns feasibility for the given a."""
    return (1 - a) <= a

def dense_radii_c(logk):
    """Dense random perturbation: generic P is a fair coin per position,
    Pr[mono] >= 2^{1-k}, so 2 c k log k <= k log 2: c <= log 2 / (2 log k)."""
    import math
    return math.log(2) / (2 * logk)

def minkowski_cells(q=101, D=2, k=5, seed=3, nbox=6):
    """Exact check of the lemma used in the no-go: a box cell of the y-torus with
    more than (4k)^D points of {y~(n)} contains a k-term progression of Z/q^D.
    Brute force on small instances (q prime, lambda = 1)."""
    import random
    rng = random.Random(seed)
    N = q ** D
    def yt(n):
        return tuple(((Q(n, q ** i) + Q(1, 2)) % 1) - Q(1, 2) for i in range(1, D + 1))
    pts = [yt(n) for n in range(N)]
    out = []
    for _ in range(nbox):
        while True:
            lo = [Q(rng.randrange(-500, 300), 1000) for _ in range(D)]
            wd = [Q(rng.randrange(150, 450), 1000) for _ in range(D)]
            cell = [n for n in range(N) if all(l <= c < l + w for c, l, w in zip(pts[n], lo, wd))]
            if len(cell) > (4 * k) ** D:
                break
        S = set(cell)
        found = None
        for d in range(1, N):
            for n0 in cell:
                if all((n0 + j * d) % N in S for j in range(k)):
                    found = (n0, d); break
            if found: break
        assert found is not None
        out.append({"cell_points": len(cell), "bound_(4k)^D": (4 * k) ** D, "k_AP_found": list(found)})
    return out

def main():
    res = {}
    res["current"] = str(current())
    assert current() == Q(1, 8)
    # negative controls: the solver sees gains when a constraint is dropped
    res["control_drop_B"] = str(current(drop=("B",)))
    res["control_drop_F_R"] = str(current(drop=("F", "R")))
    assert current(drop=("B",)) == Q(1, 4) and current(drop=("F", "R")) == Q(1, 4)
    assert current(drop=("R",)) == Q(1, 8) and current(drop=("F",)) == Q(1, 8)
    # asymmetric flips
    asym = {str(t): str(asymmetric(t)) for t in (Q(1, 10), Q(1, 4), Q(1, 2), Q(3, 4))}
    res["asymmetric_with_structured_ones"] = asym
    assert max(asymmetric(t) for t in (Q(i, 20) for i in range(1, 20))) == Q(1, 8)
    res["asymmetric_if_no_structured_ones_theta=1/20"] = str(asymmetric(Q(1, 20), structured_ones=False))
    assert asymmetric(Q(1, 20), structured_ones=False) == Q(19, 80)   # -> 1/4 as theta -> 0
    res["fine_coarse"] = str(fine_coarse())
    assert fine_coarse() == Q(1, 8)
    res["uncertainty_min_exponent_of_H"] = str(uncertainty())
    assert uncertainty() == H
    res["group_design_feasible_at_a=49/100"] = group_design(Q(49, 100))
    res["group_design_feasible_at_a=1/2"] = group_design(Q(1, 2))
    assert not group_design(Q(49, 100))
    res["dense_radii_c_at_k=2^1024"] = f"{dense_radii_c(1024 * 0.6931471805599453):.3e}"
    res["minkowski_cells_contain_k_APs"] = minkowski_cells()
    text = json.dumps(res, indent=2, sort_keys=True) + "\n"
    if "--write" in sys.argv:
        (HERE / "RESULT.json").write_text(text); print("RESULT.json written"); return 0
    assert (HERE / "RESULT.json").read_text() == text, "output differs from RESULT.json"
    print("OK: all checks pass and match RESULT.json")
    return 0

if __name__ == "__main__":
    sys.exit(main())
