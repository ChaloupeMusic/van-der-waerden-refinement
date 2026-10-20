#!/usr/bin/env python3
"""Checks for the composed endpoint bound (PROOF.md). Not a formal proof.

Base: N1b head 8c2b34f (seven-eighths) on PR #12 head 029431bb29f2bc64bf3f6dd150ddd15bbdf7c7d5; lemmas A-E
of N1a (endpoint-constants, a8a9edc). Standard library only, no network. Run from anywhere:

    python3 -B verify.py

Exit code 0 and standard output byte-identical to RESULT.json. Assertions are load-bearing; never -O.

1. Monomial certificates A t^j e^{-dt} <= B for every real t >= 512 (sign of j - dt; value at 512 with
   e >= 2); bounds with sqrt(t) are squared first.
2. Exact budgets for Sections 4-9 of PROOF.md with certified rational bounds for logarithms, for the main
   version (Lemma D', 19^k) and for the fallback (N1a's Lemma D, 37^k).
3. Exhaustive finite checks: Lemma D' (and N1a's 37^k for comparison), the canonical partition and the
   counting core of the short-block regularity lemma.
4. Negative controls, a 90-digit decimal replay on a grid of t, and file integrity.
"""
import hashlib
import json
import math
import sys
from decimal import Decimal as Dec, getcontext, MAX_EMAX, MIN_EMIN, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
T0 = 512
AD, BH, AS, AP = 22, 5, 4096, 36               # u = sqrt(k)/(AD t^{5/2}), h0 = BH D t^3, Ms = AS D t, p = AP t^{7/2}/sqrt(k)
E2 = Q(1, 64)                                   # eps_reg = E2/t, eta = eps_reg/(24 D)
G = 2 + E2 / 2                                  # gamma = 1/2 - G/t (single source)
C = Q(403, 100)                                 # F = k (t/8 - (7/8) log t - C), with Lemma D' (19^k)
C_FALLBACK = Q(35, 8)                           # with N1a's Lemma D (37^k)
CLASS, CLASS_FALLBACK = 19, 37
RHO = Q(1, 2 ** 200)                            # rho = 1/u <= 2^-200


def require(cond, msg):
    if not cond:
        raise AssertionError(msg)


def pow2(n):
    return Q(1 << n) if n >= 0 else Q(1, 1 << (-n))


def dec_floor(x, digits=6):
    n = math.floor(x * 10 ** digits)
    return f"{n // 10 ** digits}.{n % 10 ** digits:0{digits}d}" if n >= 0 else str(float(x))


def dec_ceil(x, digits=6):
    n = math.ceil(x * 10 ** digits)
    return f"{n // 10 ** digits}.{n % 10 ** digits:0{digits}d}"


# ---------------------------------------------------------------- logarithms
def ln2_bounds(N=60):
    s = sum(Q(1, n * 2 ** n) for n in range(1, N + 1))
    return s, s + Q(1, (N + 1) * 2 ** N)


def ln_bounds(y, N=60):
    y = Q(y)
    if y == 1:
        return Q(0), Q(0)
    z = (y - 1) / (y + 1)
    s = 2 * sum(z ** (2 * m + 1) / (2 * m + 1) for m in range(N))
    return s, s + 2 * z ** (2 * N + 1) / ((2 * N + 1) * (1 - z * z))


def ln_up(y):
    """Upper bound for log y, y >= 1 rational: y = 2^e m with 1 <= m < 2."""
    y = Q(y)
    e = math.floor(math.log2(y.numerator) - math.log2(y.denominator))
    while pow2(e) > y:
        e -= 1
    while pow2(e + 1) <= y:
        e += 1
    return e * ln2_bounds()[1] + ln_bounds(y / pow2(e))[1]


def ln_lo(y):
    y = Q(y)
    e = math.floor(math.log2(y.numerator) - math.log2(y.denominator))
    while pow2(e) > y:
        e -= 1
    while pow2(e + 1) <= y:
        e += 1
    return e * ln2_bounds()[0] + ln_bounds(y / pow2(e))[0]


# ---------------------------------------------------------------- 1. monomials
def monomial(name, A, j, d, B, t0=T0):
    A, d, B = Q(A), Q(d), Q(B)
    require(A >= 0 and d > 0 and B > 0, f"{name}: bad parameters")
    require(j <= 0 or d * t0 - j >= 0, f"{name}: not monotone from t0")
    bound = A * Q(t0) ** j * pow2(-math.floor(d * t0))
    require(bound <= B, f"{name}: endpoint bound fails")
    return {"name": name, "A": str(A), "j": j, "d": str(d), "B": str(B), "t0": t0,
            "endpoint_upper": f"2^{-math.floor(d * t0)}*{A * Q(t0) ** j}"}


def monomials():
    L4 = ln_up(4)
    one = 1 + RHO
    specs = [
        ("u_at_least_8", (8 * AD) ** 2, 5, 1, 1),
        ("rho_below_2^-200", AD ** 2, 5, 1, RHO ** 2),
        ("M_at_least_k_over_33", 1056, 0, 1, 1),
        ("h0_plus_one_term", 1, -1, 1, (Q(1, 4) - Q(BH, AD) * one) ** 2),       # e^{-t/2} t^{-1/2}
        ("h0_below_2M", Q(33, 8) ** 2, 1, 1, Q(1, 2)),
        ("Ms_below_k", (Q(AS + 2, AD)) ** 2, -3, 1, 1),
        ("p_at_most_half", 4 * AP * AP, 7, 1, 1),
        ("log_t_below_t_over_32", 1, 1, Q(1, 32), 1),
        ("mesh_leading", 2, 1, 3, Q(1, 2)),
        ("mesh_rounding", 2, 0, 3, Q(1, 2)),
        ("Q_over_k", (6 * Q(AD, 8)) ** 2, 7, 1, Q(1, 10 ** 4)),                 # 6Q/k leading
        # first family small terms (units of 16 t^2/k), each <= 1/200
        ("F1_h0_plus_one", (64 * L4 * one / AD) ** 2, -1, 1, Q(1, 200) ** 2),
        ("F1_cube_root", (128 * L4 * one / (AD * Q(5, 2))) ** 2, 1, Q(1, 3), Q(1, 200) ** 2),
        ("F1_32Dt", (512 * one / AD) ** 2, 1, 1, Q(1, 200) ** 2),
        ("F1_2Dlog441", (32 * ln_up(441) * one / AD) ** 2, -1, 1, Q(1, 200) ** 2),
        ("F1_t3", 16, 3, 1, Q(1, 600)),
        ("F1_t2", 16 * ln_up(96) + 2, 2, 1, Q(1, 600)),
        ("F1_const", 32, 0, 1, Q(1, 600)),
        # rich: (5Q + D(23.25t + 10.81) + 16.5t + 18.1)/k, each term <= 1/1000
        ("rich_5Q", (5 * Q(AD, 8)) ** 2, 7, 1, Q(1, 1000) ** 2),
        ("rich_15t", 15, 1, 1, Q(1, 1000)),
        ("rich_D", ((Q(2325, 100) + Q(1081, 100 * T0)) * one / AD) ** 2, -3, 1, Q(1, 1000) ** 2),
        ("rich_lin", Q(165, 10) + Q(181, 10 * T0), 1, 1, Q(1, 1000)),
        # case (i): (log T_K + log 16)/k small; case (ii) residual terms
        ("caseI_6Q", (Q(33, 2)) ** 2, 7, 1, Q(1, 100) ** 2),
        ("caseI_lin", Q(225, 10) + 9, 1, 1, Q(1, 100)),
        ("caseI_D", (17 * one / AD) ** 2, -3, 1, Q(1, 100) ** 2),
        ("caseII_lin", (Q(245, 10) + Q(9, T0)) ** 2, -5, 1, Q(1, 100) ** 2),    # (24.5t+9)/(sqrt k t^{7/2})
    ]
    return [monomial(*s) for s in specs]


# ---------------------------------------------------------------- 2. exact budgets
def budgets(Cv, cls, label):
    L2lo, L2up = ln2_bounds()
    L4, L9, L36 = ln_up(4), ln_up(9), ln_up(AP)
    one = 1 + RHO
    r = {}
    # scales
    require(Q(BH, AD) * one < Q(1, 4), "h0 <= sqrt(kt)/4 coefficient")
    require(Q(9 * T0, 64) >= L36, "L >= t/4")
    require(Q(17 * T0, 256) >= Cv, "F >= kt/32")
    require(math.ceil(24576 * T0 + 4) <= 24577 * T0, "B0 <= 24577 t")
    # geometry: D(12 eta + 3/h0) = 12 E2/(24 t) + 3/(BH t^3) <= E2/t, and 4MH/k <= eta (192 D t <= k^2: trivial)
    require(12 * E2 / 24 + Q(3, BH * T0 ** 2) <= E2, "short-block regularity budget")
    # dichotomy: A0 - 1/t - gamma = E2/t^2 exactly (rational identity in 1/t, checked at 3 points)
    for tt in (Q(T0), Q(1000), Q(10 ** 6)):
        A0 = (Q(1, 2) - 1 / tt) * (1 - E2 / tt)
        require(A0 - 1 / tt - (Q(1, 2) - G / tt) == E2 / tt ** 2, "dichotomy identity")
        require(A0 - 1 / tt > Q(1, 2) - G / tt, "dichotomy margin")
    # first family: main terms + five small terms (<= 1/200 each) <= 1; sqrt(t) >= 22.6 for t >= 512
    require(Q(226, 10) ** 2 <= T0, "sqrt(t0) >= 22.6")
    f1 = 64 * BH * L4 * (1 + 3 * RHO) / AD ** 2 + 16 * L4 * one / (AD * Q(226, 10)) + Q(5, 200)
    require(f1 <= 1, "first-family budget")
    r["first_family_total_upper"] = dec_ceil(f1, 5)
    # second family: phi(h0) <= -D((BH - 9/2) t + (7/2) log t - c2) + t + log t, where the entropy gives
    # (1/2 + 4) t per unit of D, the tail -(1 - E2/t) BH t, and c2 = 2 log 9 + 4 log 280 + BH E2;
    # B0 D h0 + 1 <= 280 k/t uses h0 <= sqrt(kt)/4 and B0 <= 24577 t.
    require(Q(24577, 4 * AD) * one <= 280 - Q(1, 100), "B0 D h0 coefficient")
    c2 = 2 * L9 + 4 * ln_up(280) + BH * E2
    require(BH > Q(9, 2), "second-family tail below its entropy")
    lt0_lo, lt0_up = ln_lo(T0), ln_up(T0)
    lhs2 = 8 * ((BH - Q(9, 2)) * T0 + Q(7, 2) * lt0_lo - c2)                # D >= 8; increasing in t
    require(lhs2 >= 3 * T0 + lt0_up + ln_up(192) + Q(1, 8), "second-family margin at t0")
    require(8 * (BH - Q(9, 2)) >= 3, "second-family margin is not increasing")
    require(5 * 8 + 2 <= (1 - E2 / T0) * BH * 8 * T0, "phi decreasing from h0")
    # third family: M_s/128 >= (AS/128) D t, log A_s <= D(18.5 t + (1/2) log t + 6.6) + 1.5 t + (1/2) log t;
    # need D((AS/128 - 37/2) t - (1/2) log t - 6.6) >= 2.5 t + (1/2) log t + log 96 + 1/8 with D >= 1
    require(3 * L9 <= Q(66, 10) + Q(1, 1000), "third-family constant")
    a3 = Q(AS, 128) - Q(37, 2)
    require(a3 >= 3 and a3 * T0 - lt0_up / 2 - Q(66, 10) >= Q(5, 2) * T0 + lt0_up / 2 + ln_up(96) + Q(1, 8),
            "third-family budget")
    # rich: gamma L / 2 identity on an interpolation grid (T, A = log t, B = log AP)
    for T in (Q(1), Q(2), Q(5)):
        for A in (Q(0), Q(1)):
            for B in (Q(0), Q(1)):
                gamma = Q(1, 2) - G / T
                Lx = T / 2 - Q(7, 2) * A - B
                rhs = T / 8 - Q(7, 8) * A - B / 4 - G / 4 + Q(7, 4) * G * A / T + G / 2 * B / T
                require(gamma * Lx / 2 == rhs, "gamma L / 2 expansion")
    # errors from PR #9's displayed factors (D <= sqrt k, D + 4 <= 2 sqrt k, J <= k^6)
    Dc_c = L9 + 2 * ln_up(7) + ln_up(3) + ln_up(6) + ln_up(19)            # D-constant in E_c
    Dc_p = L9 + 2 * ln_up(7) + L9                                         # D-constant in E_p
    c_c = ln_up(16) + 3 * ln_up(3) + 2 * 2 * L2up + 3 * ln_up(6) + 3 * ln_up(19)
    c_p = ln_up(100)
    require(Dc_c <= Q(1193, 100) and Dc_p <= Q(829, 100) and c_c <= Q(231, 10) and c_p <= Q(461, 100),
            "error constants")
    k_ebar = (ln_up(2 * cls) + ln_up(5)) / 2                              # (1/2) log(2*cls*5)
    kappa = Cv - L36 / 4 - G / 4 - k_ebar
    require(kappa > 6 * Q(1, 1000) + Q(1, 1000), f"rich budget ({label})")   # six monomials of <= 1/1000 above
    cut = Cv - L36 / 4 - G / 4 - ln_up(2 * cls) / 2
    require(cut >= Q(8, 10), f"rich cut ({label})")
    r["kappa_lower"] = dec_floor(kappa)
    r["cut_margin_lower"] = dec_floor(cut, 4)
    # cases (i) and (ii)
    require(1 >= L2up + 3 * Q(1, 100), "case (i) budget")                 # t/512 >= 1 at t0
    brun = Q(AS + 1, AD) * one * (3 + Q(1, T0))
    require(brun <= 560, "run coefficient")
    # residual terms over sqrt(k) t^{7/2}: 560 t^{-4}; D(16t+5) <= 17 (1+rho) u t gives 17 (1+rho)/(AD t^5);
    # (24.5 t + 9) t^{-7/2} e^{-t/2} is the monomial caseII_lin (<= 1/100)
    caseII_D = 17 * one / (AD * Q(T0) ** 5)
    require(6 * Q(AD, 8) + 560 * Q(1, T0 ** 4) + caseII_D + Q(1, 100) < Q(15, 32) * AP, "case (ii) budget")
    # Lemma C tail: 4(D+1) e^{-log 64 - t/2} <= 8 sqrt(k) e^{-t/2}/64 = 1/8 and 2 e^{-2 Gamma} <= 1/8 (Gamma >= 1)
    require(Q(8, 64) + Q(1, 8) <= Q(1, 4), "Lemma C tail")
    # error coefficients in t (PR #9 factors, see PROOF.md Section 7)
    require(1 + 16 + Q(3, 2) + 3 + 4 == Q(51, 2) and 1 + 16 == 17 and 1 + 6 + 2 + 9 + 12 == 30, "error t-coefficients")
    # thresholds
    require(739 * L2lo >= T0, "k >= 2^739 gives t >= 512")
    return r


def thresholds():
    L2lo, L2up = ln2_bounds()
    out = {}
    for Cv, t9, n9 in ((C, 704, 1016), (C_FALLBACK, 731, 1055)):
        g9 = Q(t9, 72) - Q(7, 8) * ln_up(t9) - Cv
        require(g9 > 0 and t9 > 63 and n9 * L2lo >= t9, "1/9 threshold")
        require(Q(t9 - 1, 72) - Q(7, 8) * ln_lo(t9 - 1) - Cv < 0, "1/9 threshold is not minimal at this integer t")
        out[f"one_ninth_C_{Cv}"] = {"t": t9, "k": f"2^{n9}", "margin_lower": dec_floor(g9)}
    # coefficient at k = 2^1024 (increasing in t), main version
    t1024 = 1024 * L2lo
    coeff = Q(1, 8) - (Q(7, 8) * ln_up(t1024) + C) / t1024
    require(coeff >= Q(1112, 10000), "coefficient at 2^1024")
    out["coefficient_at_2^1024_lower"] = dec_floor(coeff, 5)
    return out


# ---------------------------------------------------------------- 3. finite checks
def check_lemma_D_prime(kmax=300):
    cases = 0
    for K in range(1, 2001):
        require(4 * K * K + 14 * K + 1 <= CLASS ** K, "4K^2+14K+1 <= 19^K")
    for k in range(2, kmax + 1):
        for h in range(1, k):
            Kh = -(-k // h)
            prod = 1
            for a in range(h):
                Ka = -(-(k - a) // h)
                require(Ka >= 1 and Kh <= Ka + 1, "K_h <= K_a + 1")
                prod *= 2 * Ka * (2 * Kh + 5) + 1
            require(prod <= CLASS ** k, "class product above 19^k")
            require(prod <= CLASS_FALLBACK ** k, "class product above 37^k")
            cases += 1
    return cases


def canonical_partition(k, h):
    n = -(-k // h)
    base, rho = divmod(k, n)
    return [base + 1] * rho + [base] * (n - rho)


def check_partition(kmax=300):
    cases = 0
    for k in range(2, kmax + 1):
        for h in range(1, k // 2 + 1):
            lengths = canonical_partition(k, h)
            require(sum(lengths) == k and all(-(-h // 2) <= ell <= h for ell in lengths), "partition")
            cases += 1
    return cases


def check_regularity_core(smax=14):
    cases = 0
    for s in range(1, smax + 1):
        for a in [x for x in range(1, s + 1) if math.gcd(x, s) == 1]:
            for ell in range(1, 2 * s + 3):
                visits = [0] * s
                for z in range(ell):
                    visits[(z * a) % s] += 1
                require(max(visits) <= -(-ell // s), "visit count")
                doubled = visits + visits
                for g in range(1, s + 1):
                    window = sum(doubled[:g])
                    worst = window
                    for c in range(1, s):
                        window += doubled[c + g - 1] - doubled[c - 1]
                        worst = max(worst, window)
                    require(worst <= g * -(-ell // s), "arc count")
                    cases += 1
    return cases


# ---------------------------------------------------------------- 4. negative controls
def negative_controls():
    rej = {}
    L36 = ln_up(AP)
    # (a) AD = 21: the first-family main term alone exceeds 1
    rej["AD_21"] = 64 * BH * ln_lo(4) / 21 ** 2 > 1
    # (b) BH = 4: the second-family tail 4 D t no longer beats the entropy 4.5 D t (D-coefficient positive)
    rej["BH_4"] = Q(1, 2) + 4 - (1 - E2 / T0) * 4 > 0
    # (c) AP = 35: (15/32) * 35 < 6 * AD/8 = 16.5
    rej["AP_35"] = Q(15, 32) * 35 < 6 * Q(AD, 8)
    # (d) C = 4.02 with D', C = 4.35 with D: kappa < 0
    rej["C_4.02"] = Q(402, 100) - L36 / 4 - G / 4 - (ln_lo(2 * CLASS) + ln_lo(5)) / 2 < 0
    rej["C_4.35_with_D"] = Q(435, 100) - L36 / 4 - G / 4 - (ln_lo(2 * CLASS_FALLBACK) + ln_lo(5)) / 2 < 0
    # (e) 37^k in place of 19^k at C = 4.03: kappa < 0, so Lemma D' is load-bearing for 403/100
    rej["37k_at_C_4.03"] = C - L36 / 4 - G / 4 - (ln_lo(2 * CLASS_FALLBACK) + ln_lo(5)) / 2 < 0
    # (f) N1b's eligibility threshold l/2 (tail h/(2t^2)): the tail 2.5 D t loses to the entropy 4.5 D t
    rej["threshold_half"] = Q(BH, 2) < Q(9, 2)
    # (g) G = 2 (eps_reg ignored in gamma): A0 - 1/t - gamma < 0 at t0
    tt = Q(T0)
    rej["G_2"] = (Q(1, 2) - 1 / tt) * (1 - E2 / tt) - 1 / tt - (Q(1, 2) - 2 / tt) < 0
    # (h) a monomial claimed from t0 = 30 is refused
    try:
        monomial("control", (8 * AD) ** 2, 5, 1, 1, t0=30)
        rej["monomial_t0_30"] = False
    except AssertionError:
        rej["monomial_t0_30"] = True
    require(all(rej.values()), f"negative control not rejected: {rej}")
    return rej


# ---------------------------------------------------------------- 5. decimal replay
def decimal_replay():
    getcontext().prec = 90
    getcontext().Emax = MAX_EMAX
    getcontext().Emin = MIN_EMIN
    ONE = Dec(1)
    TOL = ONE - Dec(10) ** -60
    LN2, LN4, LN5, LN9 = (Dec(x).ln() for x in (2, 4, 5, 9))

    def ceil(x):
        return x.to_integral_value(rounding=ROUND_CEILING)

    def floor(x):
        return x.to_integral_value(rounding=ROUND_FLOOR)

    def dq(x):
        return Dec(x.numerator) / Dec(x.denominator)

    fails, checks, grid = 0, 0, list(range(512, 2001, 37)) + [5000, 10 ** 4, 10 ** 5]
    for tt in grid:
        for Cv, cls in ((C, CLASS), (C_FALLBACK, CLASS_FALLBACK)):
            t = Dec(tt); k = t.exp(); lt = t.ln()
            u = (t / 2).exp() / (AD * t ** Dec("2.5")); D = ceil(u); M = floor(k / 32)
            h0, Ms = ceil(BH * D * t ** 3), ceil(AS * D * t)
            eps_reg = dq(E2) / t; eta = eps_reg / (24 * D); B0 = ceil(24576 * t + 4)
            gamma = ONE / 2 - dq(G) / t
            logp = Dec(AP).ln() + Dec("3.5") * lt - t / 2; p = logp.exp(); L = -logp
            F = k * (t / 8 - Dec("0.875") * lt - dq(Cv))
            logP = LN2 + 2 * t; Qh = F / D + logP
            A0 = (ONE / 2 - ONE / t) * (1 - eps_reg)
            Lam = LN4 * (h0 + (2 * M).sqrt() + 2 * t * (2 * M) ** (ONE / 3))
            logTglob = 4 * D * Lam + 32 * D * t + 2 * D * Dec(441).ln()
            phi = (D + 2) * h0.ln() + 2 * D * LN9 + 4 * D * (B0 * D * h0 + 1).ln() - (1 - eps_reg) * h0 / t ** 2
            logAs = t + (D + 1) * h0.ln() + D * (LN9 + 2 * (k * (4 * M + 4) + 1).ln()) + D * (LN9 + 2 * (3 * k ** 7).ln())
            Gamma = gamma * k * L / 2 - (F + D * logP)
            Ec = 6 * Qh + k * (Dec(2 * cls).ln() + LN5 / 2) + D * (Dec("25.5") * t + Dec("11.93")) + 30 * t + Dec("23.1")
            Ep = 4 * Qh + k * LN5 / 2 + D * (17 * t + Dec("8.29")) + 2 * t + Dec("4.61")
            logTK = 6 * Qh + D * (16 * t + 5) + Dec("4.5") * t + 6
            logBrun = h0.ln() + (Ms + 1).ln() + Ms * (LN2 + 3 * t)
            c = [u >= 8, p <= ONE / 2, L >= t / 4, h0 < 2 * M, 2 * M / k * TOL <= ONE / 16,
                 4 * M / k ** 3 <= eta, D * (12 * eta + 3 / h0) <= eps_reg, A0 - ONE / t > gamma,
                 logTglob + Dec(96).ln() + ONE / 8 + t <= 2 * M / t ** 2,
                 phi <= -2 * t - Dec(192).ln() - ONE / 8,
                 logAs + Dec(96).ln() + ONE / 8 + t <= Ms / 128,
                 Gamma - (Ec + Ep) / 2 >= Dec(64).ln() + t / 2, Gamma > (Ec - Ep) / 2,
                 Dec(16).ln() + logTK + k * LN2 <= k * L / 128,
                 logTK + logBrun + Dec(16).ln() <= Dec(15) / 32 * p * k]
            checks += len(c)
            fails += sum(1 for ok in c if not ok)
    require(fails == 0, "decimal replay found a failing inequality")
    return {"t_values": len(grid), "versions": 2, "checks": checks, "failures": fails}


# ---------------------------------------------------------------- 6. integrity
def integrity():
    manifest = json.loads((HERE / "MANIFEST.json").read_text())
    for name, digest in sorted(manifest.items()):
        if name == "RESULT.json":
            continue
        require(hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest, f"hash mismatch: {name}")
    return sorted(name for name in manifest if name != "RESULT.json")


def main():
    report = {
        "base": {"n1b_head": "8c2b34f", "n1a_head": "a8a9edc",
                 "pr12_head": "029431bb29f2bc64bf3f6dd150ddd15bbdf7c7d5"},
        "parameters": {"u": f"sqrt(k)/({AD} t^(5/2))", "M": "floor(k/32)", "h0": f"ceil({BH} D t^3)",
                       "Ms": f"ceil({AS} D t)", "eps_reg": f"{E2}/t", "eta": "eps_reg/(24 D)",
                       "B0": "ceil(24576 t + 4)", "delta_s": "1/16", "gamma": f"1/2 - ({G})/t",
                       "p": f"{AP} t^(7/2)/sqrt(k)", "C": str(C), "C_fallback": str(C_FALLBACK), "t0": T0},
        "monomial_certificates": monomials(),
        "budgets_main_19k": budgets(C, CLASS, "main"),
        "budgets_fallback_37k": budgets(C_FALLBACK, CLASS_FALLBACK, "fallback"),
        "thresholds": thresholds(),
        "lemma_D_prime_cases": check_lemma_D_prime(),
        "partition_cases": check_partition(),
        "regularity_core_cases": check_regularity_core(),
        "negative_controls_rejected": negative_controls(),
        "decimal_replay": decimal_replay(),
        "integrity_files": integrity(),
        "conclusion": "W_r(k) > (k^(1/8) / (e^(403/100) (log k)^(7/8)))^(k floor(log2 r)) for k >= 2^739 "
                      "(e^(35/8) with N1a's Lemma D); c = 1/9 for k >= 2^1016 (2^1055); "
                      "conditional on the interfaces of PR #12, Section 2, and PR #9's displayed factors",
    }
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    sys.stdout.write(text)
    stored = HERE / "RESULT.json"
    if "--write" in sys.argv:
        stored.write_text(text)
    else:
        require(stored.read_text() == text, "output differs from RESULT.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
