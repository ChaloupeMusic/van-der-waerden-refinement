"""Exact checks for: algebraic band rules leave whole single-label blocks unsplit.

Standard library only.  Run:  python3 -B verify.py   (exit 0 = all checks pass;
the last line is "Test passé : all checks pass" and the JSON summary printed
before it must equal RESULT.json; --write regenerates RESULT.json).

The construction (README, Theorem 1).  Fix e >= 1 and L0 | lam e / 2.  Choose
the low digits r so that f = y_{D-1}(r) lies in [f*, f* + 2/q) with
f* = 1/4 + L0/(2 lam e), and put phi = 2f, b_f = lam e/2 + L0.  On the
single-label progression n_j = r + q^(D-1)(m + j e) with q y~_D(n_0) = sigma + f,
sigma = sigma_0 + 2 tau inside a "stretch", every band is EXACTLY

    Q_tau(j) = (sigma + j lam e)^2 + b_f j + I_0 + tau ,

and the discriminant of Q_tau in j is K_0 + 8 lam e L0 tau: an arithmetic
progression in tau.  Hence Q_tau takes on Z/l (l prime, l not | 2 lam e L0) the
values  squares + kappa(tau)  with kappa running over all of Z/l, and the
Chinese remainder theorem gives a block avoiding every per-prime condition.

What is checked:
  1. Algebra.  On random instances (several q = P^w, D = 3, 4, lam, e, L0, low
     digits), for every tau of the stretch and every j < k: the band computed
     from the source's coordinates equals Q_tau(j), and the discriminant equals
     K_0 + 8 lam e L0 tau.  All in Fractions.
  2. U-phase lemma.  phase(sigma) = frac(x_D(n_0) k^2 lam) moves by exactly
     k^2/q (mod 1) when sigma moves by 1, so U-spoiled starts form intervals.
  3. Toy instance in the right regime (q = 1009^6, D = 3, k = 24,
     lam = 720720, e = 256): for four rules of the class (a quadratic-character
     rule mod a prime above the P2 range, a Boolean mix with a periodic part
     mod 77, a polynomial condition combined with a character, and the CRT
     union), and for several box words: the constructed block keeps one full
     label (source U mesh and adaptive V mesh, 200-digit decimals), its bands
     are Q_tau(j), every per-prime indicator is 0 and t mod L0 is constant on
     it, so the rule puts all of its keys on the same side.  Distinct blocks
     have disjoint key sets.
  4. Beyond P2: the moduli of rules 1 and 3 exceed the largest step of any
     single-label block, so no step makes their band constant (P2 and its
     extension L | lam e do not apply).
  5. Negative controls, all rejected: L0 = 0 (f = 1/4) freezes the
     discriminant; r = 0 (f = 0) makes every block meet 0 mod l; a start
     sigma_0 outside the stretch breaks the band formula; a rule violating the
     per-prime hypothesis (T = squares) has no good residue; tau* + 1 meets the
     rule; b_f + 1 in the formula is detected; a step with L0 not | lam e / 2
     moves t mod L0.

Mesh, label and coordinate helpers are transcribed from the verifier of
radar-99's single-label-obstructions package (same base commit, same source
definitions).
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

getcontext().prec = 200
FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print('FAIL', msg)


# ------------------------------------------------------------ source objects
def centered(fr):
    """Representative of a circle point in [-1/2, 1/2)."""
    r = fr - math.floor(fr)
    if r >= Fraction(1, 2):
        r -= 1
    return r


def ytilde(n, q, D, lam):
    return [centered(Fraction(lam * n, q ** i)) for i in range(1, D + 1)]


def xhat(n, q, D):
    return [Fraction(n % q ** i, q ** i) for i in range(1, D + 1)]


def band(n, q, D, lam):
    return math.floor(sum((q * v) ** 2 for v in ytilde(n, q, D, lam)))


def mesh(q, D, k):
    alpha = Decimal(1) / (Decimal(q) ** (D + 1))
    R = (1 + 1 / (2 * alpha)).ln()
    m = int((Decimal(k * k) * R).to_integral_value(rounding='ROUND_CEILING'))
    return alpha, R, m


def v_index(Y, alpha, R, m):
    """Index of the source's adaptive V interval containing Y in [-1/2,1/2)."""
    Yd = Decimal(Y.numerator) / Decimal(Y.denominator)
    if Yd < 0:
        a = Yd + Decimal('0.5')
        side = 'L'
        l = int((m * (1 + a / alpha).ln() / R).to_integral_value(rounding='ROUND_FLOOR'))
    else:
        a = Decimal('0.5') - Yd
        side = 'R'
        l = int((m * (1 + a / alpha).ln() / R).to_integral_value(rounding='ROUND_CEILING')) - 1
    return side, l


def interval_bounds(side, l, alpha, R, m):
    a_l = alpha * ((Decimal(l) * R / m).exp() - 1)
    a_l1 = alpha * ((Decimal(l + 1) * R / m).exp() - 1)
    if side == 'L':
        return Decimal('-0.5') + a_l, Decimal('-0.5') + a_l1
    return Decimal('0.5') - a_l1, Decimal('0.5') - a_l


def u_index(X, k, lam):
    return math.floor(X * k * k * lam)


def label(n, q, D, lam, k, alpha, R, m):
    xs = xhat(n, q, D)
    ys = ytilde(n, q, D, lam)
    return (tuple(u_index(x, k, lam) for x in xs),
            tuple(v_index(y, alpha, R, m) for y in ys))


# ---------------------------------------------------------------- arithmetic
def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:                  # deterministic below 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime(n):
    while not is_prime(n):
        n += 1
    return n


def squares(l):
    return {x * x % l for x in range(l)}


def is_nonresidue(t, l):
    return t % l != 0 and pow(t % l, (l - 1) // 2, l) == l - 1


def crt(residues):
    x, M = 0, 1
    for l, a in residues:
        t = ((a - x) * pow(M, -1, l)) % l
        x, M = x + M * t, M * l
    return x % M, M


# ---------------------------------------------------------------- the rules
# A per-prime condition T_l is one of
#   ('single', rho)          T = {rho}
#   ('shiftNR', c)           T = c + (non-residues)        [quadratic character]
#   ('set', frozenset)       T given explicitly            [small l only]
#   ('poly', coeffs, S)      T = {t : P(t) mod l in S}     [small l only]
def in_T(l, spec, t):
    kind = spec[0]
    if kind == 'single':
        return t % l == spec[1] % l
    if kind == 'shiftNR':
        return is_nonresidue(t - spec[1], l)
    if kind == 'set':
        return t % l in spec[1]
    if kind == 'poly':
        v = 0
        for c in reversed(spec[1]):
            v = (v * t + c) % l
        return v in spec[2]
    raise ValueError(kind)


def T_set(l, spec):
    return {t for t in range(l) if in_T(l, spec, t)}


def hypothesis_G(l, spec, k):
    """(G)(i): T - squares != Z/l, or (G)(ii): 2k|T| < l.  Returns the label."""
    if spec[0] in ('single', 'shiftNR'):
        return 'G(i)'                # proved in the README for every l >= 5
    T = T_set(l, spec)
    if len({(t - s) % l for t in T for s in squares(l)}) < l:
        return 'G(i)'
    if 2 * k * len(T) < l:
        return 'G(ii)'
    return None


def rule_member(rule, t):
    bits = tuple(int(in_T(l, spec, t)) for l, spec in rule['primes'])
    return rule['F'](t % rule['L0'], bits), bits


# -------------------------------------------------------- the construction
class Instance:
    def __init__(self, P, w, D, k, lam, e, L0, force=False):
        self.P, self.q, self.D, self.k, self.lam, self.e, self.L0 = P, P ** w, D, k, lam, e, L0
        self.le = lam * e
        assert self.le % 2 == 0 and (force or (self.le // 2) % L0 == 0)
        self.bf = self.le // 2 + L0
        self.gamma = self.le ** 2
        self.Ts = self.le // (8 * L0)
        self.fstar = Fraction(1, 4) + Fraction(L0, 2 * self.le)
        self.lam_inv_low = pow(lam, -1, self.q ** (D - 1))
        self.lam_inv_q = pow(lam, -1, self.q)

    def word(self, low_digits):
        """r for prescribed digits u_0..u_{D-3}; digit u_{D-2} = ceil(f* q)."""
        q, D = self.q, self.D
        top = math.ceil(self.fstar * q)
        u = sum(d * q ** t for t, d in enumerate(low_digits)) + top * q ** (D - 2)
        r = (self.lam_inv_low * u) % q ** (D - 1)
        f = Fraction(self.lam * r, q ** (D - 1))
        f -= math.floor(f)
        yt = ytilde(r, q, D, self.lam)[:D - 1]
        C = sum((q * v) ** 2 for v in yt)
        return dict(r=r, f=f, phi=2 * f, alpha0=C + f * f,
                    floor_lr=(self.lam * r) // q ** (D - 1))

    def stretch_start(self, W, sigma_lo):
        """sigma_0 >= sigma_lo, same parity steps, with frac(alpha0 + phi sigma_0) < zeta."""
        zeta = 2 * W['phi'] - 1
        psi = W['alpha0'] + W['phi'] * sigma_lo
        psi -= math.floor(psi)
        # psi(sigma_lo + 2i) = frac(psi + i zeta); the first i landing in [0, zeta)
        i = 0 if psi < zeta else math.ceil((1 - psi) / zeta)
        s = sigma_lo + 2 * i
        val = W['alpha0'] + W['phi'] * s
        assert val - math.floor(val) < zeta
        return s, math.floor(val), zeta

    def K0(self, sigma0, I0):
        return self.bf ** 2 + 4 * self.le * self.bf * sigma0 - 4 * self.gamma * I0

    def Q(self, sigma0, I0, tau, j):
        return (sigma0 + 2 * tau + j * self.le) ** 2 + self.bf * j + I0 + tau

    def disc(self, sigma0, I0, tau):
        s = sigma0 + 2 * tau
        B = 2 * s * self.le + self.bf
        A = s * s + I0 + tau
        return B * B - 4 * self.gamma * A

    def n0(self, W, sigma):
        m = (self.lam_inv_q * (sigma - W['floor_lr'])) % self.q
        return W['r'] + self.q ** (self.D - 1) * m

    def good_residues(self, l, spec, sigma0, I0):
        """tau mod l such that Q_tau(j) mod l avoids T_l for all j < k."""
        if l <= 5000:
            out = []
            for tau in range(l):
                if not any(in_T(l, spec, self.Q(sigma0, I0, tau, j)) for j in range(self.k)):
                    out.append(tau)
            return out
        # large l: closed form through kappa(tau) = -Delta_0(tau)/(4 gamma)
        K0 = self.K0(sigma0, I0)
        slope = 8 * self.le * self.L0
        if spec[0] == 'shiftNR':
            kappa = spec[1] % l                     # values squares + c miss c + NR
        elif spec[0] == 'single':
            nr = next(x for x in range(2, l) if is_nonresidue(x, l))
            kappa = (spec[1] - nr) % l              # rho - kappa non-residue
        else:
            raise ValueError('large l needs a closed form')
        tau = ((-4 * self.gamma * kappa - K0) * pow(slope, -1, l)) % l
        return [tau]


def build_block(inst, W, sigma_lo, rule):
    """One block in the window starting at sigma_lo; returns a record."""
    sigma0, I0, zeta = inst.stretch_start(W, sigma_lo)
    res = []
    for l, spec in rule['primes']:
        g = inst.good_residues(l, spec, sigma0, I0)
        if not g:
            return dict(ok=False, why=f'no good residue mod {l}')
        res.append((l, g[0]))
    tau, M = crt(res) if res else (0, 1)
    if tau >= inst.Ts:
        return dict(ok=False, why='tau beyond the stretch')
    sigma = sigma0 + 2 * tau
    ns = [inst.n0(W, sigma) + inst.q ** (inst.D - 1) * j * inst.e for j in range(inst.k)]
    return dict(ok=True, sigma0=sigma0, I0=I0, tau=tau, sigma=sigma, ns=ns, M=M, zeta=zeta)


# ============================================================ 1. algebra
rng = random.Random(20261020)
alg_instances = alg_points = 0
for trial in range(60):
    P = rng.choice([101, 103, 1009, 1013])
    w = rng.randint(2, 3)
    D = rng.choice([3, 4])
    k = rng.randint(3, 9)
    lam = rng.choice([6, 12, 30, 60, 210])
    L0 = rng.choice([1, 2, 3, 5, 6])
    e = rng.randint(1, 6)
    if (lam * e // 2) % L0 or lam * e < 8 * L0 or 64 * k * lam * e > P ** w:
        continue
    inst = Instance(P, w, D, k, lam, e, L0)
    q = inst.q
    low = [rng.randrange(q // 8 + 1, 3 * q // 8) for _ in range(D - 2)]
    W = inst.word(low)
    check(W['f'] >= inst.fstar and W['f'] < inst.fstar + Fraction(2, q),
          f'alg: f outside its window, trial {trial}')
    check(math.floor(W['phi'] * inst.le) == inst.bf, f'alg: floor(phi lam e) != b_f, trial {trial}')
    sigma_lo = q // 8 + rng.randrange(q // 16)
    sigma0, I0, zeta = inst.stretch_start(W, sigma_lo)
    if sigma0 + 2 * inst.Ts + k * inst.le >= 3 * q // 8:
        continue
    alg_instances += 1
    K0 = inst.K0(sigma0, I0)
    for tau in range(inst.Ts):
        sigma = sigma0 + 2 * tau
        n0 = inst.n0(W, sigma)
        check(q * ytilde(n0, q, D, lam)[-1] == sigma + W['f'], f'alg: q y~_D(n_0) != sigma + f, trial {trial}')
        for j in range(k):
            n = n0 + q ** (D - 1) * j * e
            check(band(n, q, D, lam) == inst.Q(sigma0, I0, tau, j),
                  f'alg: band != Q_tau(j), trial {trial} tau {tau} j {j}')
            alg_points += 1
        check(inst.disc(sigma0, I0, tau) == K0 + 8 * inst.le * L0 * tau,
              f'alg: discriminant not K_0 + 8 lam e L0 tau, trial {trial}')
check(alg_instances >= 20, f'alg: only {alg_instances} instances')
print(f'1. algebra: {alg_points} bands on {alg_instances} random instances equal Q_tau(j); '
      f'discriminants form the predicted progression')


# ============================================================ 2. U-phase lemma
P, w, D, k, lam = 1009, 6, 3, 24, 720720
q = P ** w
lam_inv_q = pow(lam, -1, q)
r0 = rng.randrange(q ** (D - 1))
flr = (lam * r0) // q ** (D - 1)


def phase(sigma):
    m = (lam_inv_q * (sigma - flr)) % q
    X = Fraction(r0 + q ** (D - 1) * m, q ** D) * k * k * lam
    return X - math.floor(X)


for _ in range(200):
    s = rng.randrange(q)
    dlt = rng.randrange(1, 10 ** 6)
    diff = phase(s + dlt) - phase(s) - Fraction(dlt * k * k, q)
    check(diff == math.floor(diff), 'U-phase: increment is not k^2/q mod 1')
print('2. U-phase: phase(sigma + d) - phase(sigma) = d k^2/q mod 1 on 200 random pairs')


# ============================================================ 3. toy instance
e = 2 ** 18
alpha, R, mm = mesh(q, D, k)
H = Fraction(1, k * k)
check(lam <= Fraction(q, 400 * k ** 3), 'toy: lam above q/(400 k^3)')
check((6 * k + 12) * lam * e + 8 <= q * H / 32 - 2, 'toy: hypothesis (S2) fails')
check(32 * k * lam * e < q, 'toy: hypothesis (S3) fails')
W_toy = math.floor((q * H / 32 - 2) / ((k + 2) * lam * e)) - 4
check(W_toy >= 2, 'toy: W < 2 although (S2) holds')

# longest V interval (upper bound), hence the largest step of a single-label block
a_top = alpha * ((Decimal(mm) * R / mm).exp() - 1)
a_prev = alpha * ((Decimal(mm - 1) * R / mm).exp() - 1)
Vmax = a_top - a_prev
e_single_max = int(Decimal(q) * Vmax / ((k - 1) * lam)) + 1   # every single-label step is below

l_big = next_prime(10 ** 9)
l_mid = next_prime(10 ** 6)
rules = [
    dict(name='quadratic character mod a prime above the P2 range',
         L0=1, primes=[(l_big, ('shiftNR', 12345))], F=lambda t0, b: b[0] == 1),
    None,                                         # Boolean mix, built below
    dict(name='cubic condition mod 1201 OR character mod a prime near 10^6',
         L0=1, primes=[(1201, ('poly', [7, 2, 0, 1], frozenset({0, 1, 2}))),
                       (l_mid, ('shiftNR', 3))],
         F=lambda t0, b: max(b) == 1),
    dict(name='CRT union {t = 0 mod l}, l = 17, 19, 23',
         L0=1, primes=[(17, ('single', 0)), (19, ('single', 0)), (23, ('single', 0))],
         F=lambda t0, b: max(b) == 1),
]
# Boolean mix: an arbitrary truth table on (t mod 77, four indicators)
sq29 = squares(29)
t2 = next(t for t in range(4, 29)
          if len({(x - s) % 29 for x in (3, t) for s in sq29}) < 29)
table_rng = random.Random(77)
table = {(t0, bits): table_rng.random() < 0.5
         for t0 in range(77) for bits in
         [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]}
rules[1] = dict(name='arbitrary Boolean rule on (t mod 77, four per-prime indicators)',
                L0=77, primes=[(17, ('single', 5)), (19, ('shiftNR', 5)),
                               (23, ('single', 0)), (29, ('set', frozenset({3, t2})))],
                F=lambda t0, b: table[(t0, b)])

# box words: coordinate 1 (low digit) and coordinate D (moving) V intervals
targets = [Fraction(1, 8) + Fraction(t, 97) for t in range(1, 23)]


def digit_inside(lo, hi):
    u = int((lo * q).to_integral_value(rounding='ROUND_CEILING'))
    if Decimal(u + 1) / q > hi:
        return None
    return u % q


words1 = []
for t1 in targets[:8]:
    side, l = v_index(t1, alpha, R, mm)
    lo, hi = interval_bounds(side, l, alpha, R, mm)
    d = digit_inside(lo, hi)
    if d is not None and (side, l) not in [w_[0] for w_ in words1]:
        words1.append(((side, l), d))
check(len(words1) >= 6, 'toy: too few coordinate-1 words')
intervals_D = []
for tD in targets[12:22:3]:
    side, l = v_index(tD, alpha, R, mm)
    intervals_D.append(interval_bounds(side, l, alpha, R, mm))

summary_rules = []
n_pairs = n_blocks = 0
for ri, rule in enumerate(rules):
    keysets = []
    inst = Instance(P, w, D, k, lam, e, rule['L0'])
    Mrule = math.prod(l for l, _ in rule['primes'])
    hyps = [hypothesis_G(l, spec, k) for l, spec in rule['primes']]
    check(all(h is not None for h in hyps), f'toy rule {ri}: hypothesis (G) fails')
    check(all(l % 2 and (2 * lam * e * rule['L0']) % l for l, _ in rule['primes']),
          f'toy rule {ri}: a prime divides 2 lam e L0')
    check(rule['L0'] * (Mrule + 1) * 8 <= lam * e, f'toy rule {ri}: hypothesis (S1) fails')
    blocks = spoiled = 0
    for (vw1, d1) in words1:
        W = inst.word([d1])
        check(v_index(ytilde(W['r'], q, D, lam)[0], alpha, R, mm) == vw1,
              f'toy rule {ri}: coordinate-1 word not realized')
        for lo, hi in intervals_D[:2]:
            sigma_lo = math.ceil(Fraction(lo) * q - W['f'])
            sigma_hi = math.floor(Fraction(hi) * q - W['f'])
            win = (k + 2) * inst.le
            got = 0
            for wi in range(6):
                start = sigma_lo + wi * win
                if start + win >= sigma_hi:
                    break
                B = build_block(inst, W, start, rule)
                check(B['ok'], f'toy rule {ri}: construction failed ({B.get("why")})')
                if not B['ok']:
                    continue
                check(B['sigma'] + (k - 1) * inst.le < start + win, f'toy rule {ri}: block leaves window')
                labs = {label(n, q, D, lam, k, alpha, R, mm) for n in B['ns']}
                if len(labs) != 1:
                    spoiled += 1
                    continue
                bands = [band(n, q, D, lam) for n in B['ns']]
                check(bands == [inst.Q(B['sigma0'], B['I0'], B['tau'], j) for j in range(k)],
                      f'toy rule {ri}: bands differ from Q_tau')
                mem = [rule_member(rule, t) for t in bands]
                check(all(bits == (0,) * len(rule['primes']) for _, bits in mem),
                      f'toy rule {ri}: an indicator is 1 on the block')
                check(len({t % rule['L0'] for t in bands}) == 1, f'toy rule {ri}: t mod L0 moves')
                check(len({v for v, _ in mem}) == 1, f'toy rule {ri}: block is split')
                vword = next(iter(labs))[1]
                keysets.append({(vword, t) for t in bands})
                blocks += 1
                got += 1
                if got == 2:
                    break
    for i in range(len(keysets)):
        for j in range(i + 1, len(keysets)):
            n_pairs += 1
            check(not (keysets[i] & keysets[j]), f'toy rule {ri}: two blocks share a key')
    n_blocks += len(keysets)
    check(spoiled <= 4, f'toy rule {ri}: {spoiled} spoiled windows')
    check(blocks >= 16, f'toy rule {ri}: only {blocks} unsplit blocks')
    summary_rules.append(dict(rule=rule['name'], L0=rule['L0'],
                              primes=[l for l, _ in rule['primes']],
                              hypotheses=hyps, modulus=rule['L0'] * Mrule,
                              unsplit_blocks_checked=blocks))
    print(f'3.{ri + 1} {rule["name"]}: {blocks} single-label blocks, every key on one side')

print(f'3.  disjointness: {n_blocks} blocks, {n_pairs} pairs within a rule, no shared key')


# ============================================================ 4. beyond P2
for ri in (0, 2):
    L = rules[ri]['L0'] * math.prod(l for l, _ in rules[ri]['primes'])
    check(math.gcd(L, lam) == 1 and L > e_single_max,
          f'beyond P2: rule {ri} modulus within reach of a constant-band step')
check(e_single_max < l_big, 'beyond P2: l_big not above the largest single-label step')
print(f'4. beyond P2: largest single-label step < {e_single_max}; rule moduli '
      f'{l_big} and {1201 * l_mid} are coprime to lam and larger')


# ============================================================ 5. controls
ctrl = 0
# (a) L0 = 0: f = 1/4 + O(1/q), b_f = lam e/2, the discriminant does not move
inst = Instance(P, w, D, k, lam, e, 1)
W = inst.word([words1[0][1]])
sigma0, I0, _ = inst.stretch_start(W, q // 4)
frozen = set()
f0 = Fraction(math.ceil(q / 4), q)                  # the f* of L0 = 0
bf0 = inst.le // 2
for tau in range(20):
    s = sigma0 + 2 * tau
    A = s * s + math.floor(W['alpha0'] - W['f'] ** 2 + f0 ** 2 + 2 * f0 * s)
    B = 2 * s * inst.le + bf0
    frozen.add(B * B - 4 * inst.gamma * A)
check(len(frozen) <= 2, 'control (a): with L0 = 0 the discriminant should not move')
ctrl += 1
# (b) r = 0: bands are exact squares, every block with l not | e meets 0 mod l
for l in (17, 19, 23):
    for _ in range(50):
        s = rng.randrange(q // 8, 3 * q // 8)
        check(any(((s + j * inst.le) ** 2) % l == 0 for j in range(k)),
              'control (b): an r = 0 block avoids 0 mod l')
ctrl += 1
# (c) a start outside the stretch breaks the formula
inst = Instance(P, w, D, k, lam, e, 1)
W = inst.word([words1[1][1]])
s = q // 4
psi = W['alpha0'] + W['phi'] * s
psi -= math.floor(psi)
zeta = 2 * W['phi'] - 1
if psi <= 1 - Fraction(1, 64):
    s += 2 * (math.floor((1 - Fraction(1, 64) - psi) / zeta) + 1)
psi = W['alpha0'] + W['phi'] * s
check(psi - math.floor(psi) > 1 - Fraction(1, 64), 'control (c): bad start not found')
I0bad = math.floor(W['alpha0'] + W['phi'] * s)
bad = 0
for tau in range(0, inst.Ts, max(1, inst.Ts // 40)):
    n = inst.n0(W, s + 2 * tau)
    if any(band(n + q ** (D - 1) * j * e, q, D, lam) != inst.Q(s, I0bad, tau, j) for j in range(k)):
        bad += 1
check(bad > 0, 'control (c): a bad stretch start went undetected')
ctrl += 1
# (d) T = squares mod 29 violates (G); the construction finds no residue
check(hypothesis_G(29, ('set', frozenset(squares(29))), k) is None, 'control (d): (G) should fail')
check(inst.good_residues(29, ('set', frozenset(squares(29))), sigma0, I0) == [],
      'control (d): a good residue exists for T = squares')
ctrl += 1
# (e) tau* + 1 meets the Boolean rule's indicators
inst = Instance(P, w, D, k, lam, e, 77)
W = inst.word([words1[2][1]])
lo, hi = intervals_D[0]
B = build_block(inst, W, math.ceil(Fraction(lo) * q - W['f']), rules[1])
hit = any(any(rule_member(rules[1], inst.Q(B['sigma0'], B['I0'], B['tau'] + 1, j))[1])
          for j in range(k))
check(hit, 'control (e): tau* + 1 should meet an indicator')
ctrl += 1
# (f) b_f + 1 in the formula is detected
n = B['ns'][3]
check(band(n, q, D, lam) != inst.Q(B['sigma0'], B['I0'], B['tau'], 3) + 3,
      'control (f): b_f + 1 not detected')
ctrl += 1
# (g) the periodic part needs L0 | lam e/2: with L0 = 17 (prime to lam e) t mod 17 moves
inst17 = Instance(P, w, D, k, lam, e, 17, force=True)
W17 = inst17.word([words1[2][1]])
B17 = build_block(inst17, W17, math.ceil(Fraction(lo) * q - W17['f']),
                  dict(L0=17, primes=[], F=lambda t0, b: t0 == 0))
check(len({band(n, q, D, lam) % 17 for n in B17['ns']}) > 1,
      'control (g): t mod 17 should move when 17 does not divide lam e/2')
ctrl += 1
print(f'5. controls: {ctrl} negative controls rejected')


# ============================================================ summary
summary = {
    'base_commit': '63bef9ca996786a41169de4a0c3d6b3b47b0fe96',
    'theorem': {
        'band_polynomial': 'Q_tau(j) = (sigma_0 + 2 tau + j lam e)^2 + b_f j + I_0 + tau, b_f = lam e/2 + L0',
        'discriminant': 'K_0 + 8 lam e L0 tau',
        'class': 'X_s on word W = {t : F_W(t mod L0, (1[t mod l in T_l])_l) = 1}, F_W arbitrary',
        'per_prime_hypothesis': ['G(i): T_l - squares != Z/l', 'G(ii): 2k|T_l| < l'],
        'hypotheses': ['L0 | lam e/2', 'primes odd, not dividing lam e L0',
                       'S1: 8 L0 (M + 1) <= lam e', 'S2: (6k+12) lam e + 8 <= q H/32 - 2',
                       'S3: 32 k lam e < q'],
        'unsplit_blocks_at_least': '(k^2/2 - 4)^(D-1) * W, W = floor((qH/32 - 2)/((k+2) lam e)) - 4 >= 2',
        'consequence': {'inside': 'K_in blocks inside X_s: monochromatic for every draw; K_in >= 1 gives P(no mono) = 0',
                        'outside': 'K_out blocks outside X_s: sub-prescription union bound >= K_out (1-p)^k, P(no mono) <= exp(-K_out (1-p)^k)',
                        'split': 'K = K_in + K_out; K_in = 0 forces K_out = K'},
    },
    'toy_instance': {'q': '1009^6', 'D': D, 'k': k, 'lam': lam, 'e': e, 'W': W_toy,
                     'largest_single_label_step_below': e_single_max, 'l_big': l_big,
                     'rules': summary_rules},
    'finite_checks': {'algebra_instances': alg_instances, 'algebra_band_points': alg_points,
                      'u_phase_pairs': 200, 'disjoint_block_pairs': n_pairs,
                      'negative_controls_rejected': ctrl},
}
out = json.dumps(summary, indent=2, sort_keys=True) + '\n'
result = Path(__file__).with_name('RESULT.json')
if '--write' in sys.argv:
    result.write_text(out)
elif not result.exists() or result.read_text() != out:
    check(False, 'RESULT.json differs from the computed summary')
print(out)
if FAIL:
    print(f'{len(FAIL)} check(s) failed')
    sys.exit(1)
print('Test passé : all checks pass')
