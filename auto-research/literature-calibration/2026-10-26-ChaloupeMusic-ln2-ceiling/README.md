# Literature calibration: the `ln 2` ceiling of the uniform coefficient

**Status:** three elementary propositions and a map of the 2024–2026
literature. No new bound on `W_r(k)`. Proposition 2 says what an improvement
past `ln 2` would settle; it is not an impossibility result.
**Base commit:** `63bef9ca996786a41169de4a0c3d6b3b47b0fe96`.
**Depends on:** the source's digit product (Proposition 7.1 of the pinned
manuscript, `\label{transfer:product}`), used as stated.

## Normalization

All logarithms `ln` are natural, as in Fox–Hunter ("All logarithms are base
`e` unless otherwise specified"). For `r ≥ 2` and `α > 0`,

```
r^{α k ln k} = exp(α k ln k ln r) = k^{(α ln 2) k log₂ r},
```

so a bound `W_r(k) ≥ r^{α k ln k}` is the coefficient `c = α ln 2` in this
repository's form `W_r(k) > k^{c k ⌊log₂ r⌋}` (at least that, since
`⌊log₂ r⌋ ≤ log₂ r`).

## Where the literature stands, in this normalization

| Result | Statement | `c` here |
|---|---|---|
| Erdős–Lovász 1975; Szabó 1990; Kozik–Shabanov 2016 | `W_r(k) ≥ β r^{k−1}` | 0 |
| Berlekamp 1968; Campos–Fox–Schildkraut 2026 (arXiv 2608.20824) | `W_2(k) ≥ (1−o(1)) k 2^{k−1}`, all `k` | 0 |
| Blankenship–Cummings–Taranchuk 2018 | `w(r,k) > p(w(r−⌈r/p⌉,k) − 1)`, `p ≤ k` prime | 0 |
| Hunter 2025 (Israel J. Math., arXiv 2301.06212) | `(3^{r/3})^{(1−o(1))k}`, `r ≥ 6`, `3 ∣ r` | 0 |
| Kim 2026 (arXiv 2610.32000) | explicit blocks for `3 ≤ r ≤ 6`, `9 ≤ k ≤ 28`; proved not to change asymptotics | 0 |
| Fox–Hunter 2026, Thm 1 (arXiv 2606.02541) | `W_3(k) > 2^{k log* k / 4}` | 0 (superexponential, but `log* k = o(ln k)`) |
| **Fox–Hunter 2026, Thm 2** | `W_r(k) ≥ r^{(1−ε)k ln k}` for `k ≥ k₀(ε)`, `r ≥ (ln k)^{3/ε}` | **`(1−ε) ln 2`**, only for `r ≥ (ln k)^{3/ε}` |
| OpenAI 2026 (pinned source) | uniform in `r ≥ 2`, absolute threshold | `10⁻⁵` |
| This repository (conditional) | uniform in `r ≥ 2` | `1/9`; sup `1/8` of the current template |

The density side, relevant to a deterministic progression-sparse set
(Porte 1): Behrend 1946 and Rankin 1961; Elsholtz–Hunter–Proske–Sauermann
2024 (arXiv 2406.12290), the first improvement of Behrend's constant for
three-term progressions, `C = 2√(ln(24/7) ln 2) < 2√(ln 4)`; Kelley–Meka 2023
and Leng–Sah–Sawhney 2024 on the upper side. Green 2022 and Hunter 2022
treat the asymmetric `w(3,k)`. The pinned source cites most of these; this
repository did not.

## Proposition 1 (the uniform coefficient is a two-color exponent)

Let `N_cyc(k)` be the largest `N` such that some `C : ℤ/Nℤ → {0,1}` has no
monochromatic `a, a+d, …, a+(k−1)d` with `d ≢ 0 (mod N)`. Put

```
c_cyc  = liminf_k  ln N_cyc(k)      / (k ln k),
c_2    = liminf_k  ln (W_2(k) − 1)  / (k ln k),
c_unif = sup { c : ∃K ∀k ≥ K ∀r ≥ 2,  W_r(k) > k^{c k ⌊log₂ r⌋} }.
```

Then `c_cyc ≤ c_unif ≤ c_2`, and `N_cyc(k) ≤ W_2(k) − 1`.

*Proof.* If `c < c_cyc`, then `N_cyc(k) > k^{ck}` for all large `k`, and the
digit product gives `W_r(k) > N_cyc(k)^{⌊log₂ r⌋} > k^{ck⌊log₂ r⌋}` for every
`r ≥ 2` at once, so `c ≤ c_unif`. The case `r = 2` of the uniform statement
gives `c_unif ≤ c_2`. Finally, a progression in `{0, …, N−1}` with step
`1 ≤ d` and `(k−1)d ≤ N−1` has `d ≢ 0 (mod N)`, so the restriction of a cyclic
coloring is a valid interval coloring. ∎

So the uniform question **is** the two-color question, between its cyclic
and its interval forms. Fox–Hunter's Theorem 2 does not move `c_unif`: its
range `r ≥ (ln k)^{3/ε}` never contains `r = 2`. (PR #16, *Remark
(consolidation)*, notes the same thing in words: the uniform statement is
pinned by small `r`.)

## Proposition 2 (the `ln 2` ceiling)

Suppose `c_cyc > ln 2`, and write `c_cyc > (1+2δ) ln 2` with `δ > 0`. Then for
every large `k` and **every** `r = 2^m`, `m ≥ 1`, the interval
`[⌊r^{(1+δ)k ln k}⌋]` contains a set of density at least `1/r` with no
`k`-term arithmetic progression.

*Proof.* For large `k`, `N = N_cyc(k) ≥ k^{(1+δ)(ln 2)k}`. The digit product
colors `[0, N^m)` with `2^m = r` colors and no monochromatic `k`-term
progression, and `N^m ≥ exp((1+δ)(ln 2) k m ln k) = r^{(1+δ)k ln k}`.
Restrict to the first `M = ⌊r^{(1+δ)k ln k}⌋ ≤ N^m` integers: the restriction
still has no monochromatic progression, and one of its `r` classes has at
least `M/r` elements. ∎

Fox–Hunter (arXiv 2606.02541, §1, after Theorem 2): *"when `r` is not too
large with respect to `k` (i.e., when `r ≤ k^{O(1)}`), the authors are not
aware of any way to construct a set `S ⊂ [r^{(1+Ω(1))k log k}]` with density
`1/r` which does not contain a `k`-term arithmetic progression. This serves
as a bottleneck for improving the exponent in Theorem 2 by a factor of
`1 + Ω(1)`."* Proposition 2 produces that set for all `r = 2^m`, including
every `2^m ≤ k^{O(1)}`.

Reading. A cyclic two-color coefficient above `ln 2 ≈ 0.693` would settle a
construction Fox and Hunter state they do not know how to make. Between
`1/8` and `ln 2` no obstruction is known from the literature: the `1/8`
ceilings of PR #16, and the obstructions of PR #14 behind them, bound the
current template, not the problem. In these units the conditional `1/9` is `16.0 %` of `ln 2`, the
template's `1/8` is `18.0 %`, and the trust-reduced `1/16` is `9.0 %`.

The digit product needs a **cyclic** coloring; the source's two-color
theorem (`\label{perturb:cyclic}`) is cyclic, so any improvement of its
template is automatically covered.

## Proposition 3 (where Fox–Hunter already pass a given `c`)

For a fixed `0 < c < ln 2`, set `ε = 1 − c/ln 2`. For `k ≥ k₀(ε)` and
`r ≥ (ln k)^{3/ε}`, Fox–Hunter give `W_r(k) ≥ r^{(1−ε)k ln k} =
k^{c k log₂ r} ≥ k^{c k ⌊log₂ r⌋}`.

| `c` | `ε` | Fox–Hunter reach `c` once `r ≥` |
|---|---|---|
| `1/16` | `0.9098` | `(ln k)^{3.297}` |
| `1/9` | `0.8397` | `(ln k)^{3.573}` |
| `1/8` | `0.8197` | `(ln k)^{3.660}` |

`k₀(ε)` is not explicit in their paper, so this is an asymptotic statement
only. Against the explicit thresholds of PR #15 (`2^{960}`, `2^{739}`) no
comparison is possible from their text.

## A direction the literature closes

Fox–Hunter's random shifted product (their Lemma 4.2) enlarges a robust
coloring by a factor exponential in `k` per application, at a
double-exponential loss in the density parameter, and they iterate it about
`(log* k)/2` times. On that accounting it cannot move `c`, which needs a
`k^{Ω(k)}` gain: feeding it the source's two-color coloring would change
lower-order factors only. This is a reading of their stated accounting, not
a proof.

## Reproduce

```sh
python3 -B verify.py
```

Expected: `OK: all checks pass and match RESULT.json` (standard library,
under a second). It checks the normalization, then the mechanics of
Proposition 2 exhaustively on small cases: the largest `N` with a cyclic
two-coloring (`4` for `k = 3`, `11` for `k = 4`, `20` for `k = 5`, within
the searched ranges), its digit product on `[0, N^m)` for `m = 2, 3` without
monochromatic progressions, and a class of density at least `2^{−m}` without
a `k`-term progression. A coloring with a monochromatic progression, and its
product, are caught. Finally it recomputes the table of Proposition 3.
