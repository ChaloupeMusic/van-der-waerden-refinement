# Optimistic ceilings of templates outside the current one

**Status:** structural research, negative with one conditional exception. No
new bound. The numbers are first-order optimistic ceilings (exponents in units
of `k log k`), not theorems about `W_2(k)`.
**Base commit:** `63bef9ca996786a41169de4a0c3d6b3b47b0fe96`.

## Frame

Every known superexponential construction for `W_2(k) > k^{ck}` is a base
coloring `c₀` plus a rare perturbation of rate `p = k^{−b}`. Two first-moment
costs decide it:

* **rich** progressions (`c₀` balanced): monochromatic means changing the
  minority, `Pr ≥ p^{k/2}` when the positions carry distinct units of
  randomness, so `2c ≤ γb` with `γ ≤ 1/2`;
* **structured** progressions (`c₀` constant): they must receive a
  perturbation, `Pr[none] ≈ e^{−pk}`, so their events need entropy
  `H ≪ pk = k^{1−b}`.

Hence `c ≤ (1/4)(1 − log_k H)` for this shape, and **`c > 1/8` needs `H ≪
k^{1/2}`**.

## The torus uncertainty principle

On the `D`-dimensional torus (`N = q^D`, `D = k^a`) the structured events of
case (ii) pay for two words, each shown necessary for union bounds with iid
flips by the concurrent obstructions P1 and P3 (`short-period-sensitivity`
track):

* the band word, three real coefficients at resolution `q²`:
  `6 log q = 6 log N/D = k^{1−a+o(1)}` (P1, `b < a`);
* the box word, two parameters per coordinate: `(D−1)·log(k²/2) =
  k^{a+o(1)}` (P3, `a + b < 1`).

So `H ≥ log N/D + D log k`, whose exponent `max(1 − a, a)` is at least `1/2`
**for every `a`**: `b ≤ 1/2` is one phenomenon seen from both sides, not two
coincidences. A template beyond `1/8` must break this product, with keys
whose slow words cost `o(log N/D)`, or a box/run structure costing `o(D)`.
Separately, the label entropy `D^{2−o(1)}` (dilation `λ`) keeps `a < 1/2`
(constraint F); dropping only the box word gains nothing (`b < a ≤ 1/2`).

## Candidates

| template | what changes | optimistic system | ceiling |
|---|---|---|---|
| current (PR #7) | — | `b ≤ a`, `a + b ≤ 1`, `a ≤ 1/2`, `2c ≤ b/2` | **1/8** |
| asymmetric flips, rates `p₀ ≠ p₁`, base with a fraction θ of ones | rich progressions pay only the flips of one color | structured progressions of both colors exist (every label cell contains step-1 progressions), so both `b₀, b₁ ≤ 1/2`; `2c ≤ min(θb₁, (1−θ)b₀)` | **1/8** at θ = 1/2 (1/40 at θ = 1/10) |
| asymmetric, **if no structured progression has mostly ones** | `p₁ = 0`: a progression containing a one is never monochromatic 0 | `2c ≤ (1−θ)b₀`, `b₀ ≤ 1/2` | **(1−θ)/4 → 1/4** (19/80 at θ = 1/20) |
| fine / coarse split, flips only on a set Σ of density `k^{−s}` | structured progressions should see Σ densely | step-1 progressions tile `G`, so density(Σ) ≥ min over structured S of `|S∩Σ|/k`; rate seen by both sides `k^{−(b+s)}` | **1/8** |
| one flip per group (covering designs) | a deterministic transversal of the slow words | disjoint groups give ≤ `#keys/g = q^{2+o(1)}` groups for `q^{6−o(1)}` distinct slow words, each containing O(1) groups: needs `1 − a ≤ a` | infeasible for `a < 1/2` |
| random radii or offsets per box (Green 2022, Hunter 2022) | randomness in the geometry | dense: each position of a generic P is a fair coin, `Pr[mono] ≥ 2^{1−k}`, `c ≤ log 2/(2 log k) → 0`; sparse: the current template | 0 or 1/8 |
| LCM dilations (Fox–Hunter) | control of rational returns | already the `λ` of the current template | 1/8 |
| finite-field norm colorings (Shi–Dong) | coloring of `F_{q^D}` | progressions of `[N]` do not map to progressions of `(Z/q)^D` (carries); no transfer | — |
| iterated products, direct `r`-color templates | more colors | the uniform-in-`r` coefficient is the `r = 2` one; a direct 3-coloring gives `(1−1/r)b/2 → 1/6` for `W_3` only | 1/8 for `r = 2` |

## The one exception, stated exactly

The asymmetric template reaches `1/4` if there is a set `A = c₀^{−1}(1)` with

1. **AP-sparse:** every `k`-term progression has at most `θk` points in `A`
   (so monochromatic 1 needs `(1−θ)k` flips, and no progression lies in `A`),
2. **generic hitting:** every progression disjoint from `A` belongs to a
   family whose flip events have entropy `≪ k^{1−b₀}`.

No label-based `A` works: a label cell contains step-1 progressions of length
`k`, which would lie in `A`. A random `A` fails (2): a generic progression
misses a random set of density `δ` with probability `(1−δ)^k`, and `N²` of
them would need `δ ≳ 2c log k`, contradicting (1). The object needed is a
**sparse hitting set of the high-entropy progressions containing no dense
part of any progression**, the dual of the open "structured transversal"
(Porte 1: a sparse transversal of the low-entropy slow words). Both give
`1/4` optimistically, and both are open.

## Remark (consolidation)

Fox–Hunter (Theorem 2, cited in the upstream introduction): `W_r(k) ≥
r^{(1−ε)k log k}` for `r ≥ (log k)^{3/ε}`, i.e. a coefficient `(1−ε) ln 2` in
the normalization `k^{ck⌊log₂ r⌋}`, already `ln 2/2 ≈ 0.35` at `r ≥ (log k)^6`.
The uniform-in-`r` statement is pinned by small `r`.

## Reproduce

```sh
python3 -B verify.py
```

Expected: `OK: all checks pass and match RESULT.json` (standard library,
under a second). Each optimistic system is solved exactly by vertex
enumeration in Fractions. Negative controls: dropping `b ≤ a` or both
`a ≤ 1/2` and `a + b ≤ 1` must give `1/4`, and dropping only one of the
latter must give `1/8` (they do), so the solver sees gains when they exist.

## Scope

First-order optimistic costs only: lower-order terms, local-lemma weights and
the outer-coloring construction are not modeled. A ceiling `≤ 1/8` here means
that even the favorable budget does not pass `1/8`. "No candidate found" is
not an impossibility theorem, and the exception above is open, not refuted.
