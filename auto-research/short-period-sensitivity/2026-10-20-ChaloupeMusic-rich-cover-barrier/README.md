# Rich-side barrier for cover bounds: c ≤ b/4, conditional on a rank statement

**Status:** research note (conditional barrier, with exact certificates and a
scoped measurement). Not a new bound, and not a statement about van der
Waerden numbers: it is about what a family of proofs can show.
**Base commit:** `63bef9ca996786a41169de4a0c3d6b3b47b0fe96`, on top of the
`short-period-sensitivity` track (`2026-10-20-erwan-binding-lever`).

## Question

In the short-period region the coefficient enters only through the rich
condition `2c < γb` with `γ < 1/2`. Can the rich side give more, for the
keys of the construction (V-box, `⌊‖qỹ‖²⌋`), by counting fewer than `N²`
rich events?

Every rich estimate in this repository (source, #4, #6, #7) is a **cover
bound**: a family 𝒢 of key sets such that every rich `(P, β)` has some
`G ∈ 𝒢` inside its minority key set `F_{P,β}`, and
`Pr[some rich P mono] ≤ Σ_𝒢 p^{|G|}`.

## Result

1. **Duality (rigorous, trivial).** If a probability measure ν on rich
   `(P, β)` satisfies `ν(F ⊇ S) ≤ p^{|S|}` for every nonempty key set `S`,
   every cover costs at least 1.
2. **Small sets (rigorous).** With ν uniform on ≥ N²/4 non-short-period
   progressions and β the majority color (`|F| ≤ k/2`), the bound holds for
   `|S| = 1` (box fibre ≤ `Π(q|I_i| + 1)`) and for `2 ≤ |S| < s₀`,
   `s₀ = ⌊(4 − 2/log k)D/b⌋` (two keys, two positions, two fibres).
3. **Linear rigidity is false (exact certificate).** Treating
   `Θ = (|A|², ⟨A,B⟩, |B|²)` as free, the band at position `t` is linear in
   `z = (A, B, Θ)` with integer row `r_t = (−2w_t, −2t w_t, 1, 2t, t²)`. One
   would like any `O(D)` rows to have full rank. They do not: by pigeonhole on
   the `2^D` sums `⟨g′,B⟩`, `g′ ∈ {0,1}^D`, some `g ∈ {−1,0,1}^D` has
   `|⟨g,B⟩| ≤ D2^{−D}`; when `2^D > Dk`, `⟨g, w_t⟩` takes ≤ `D+2` values, and
   on a level set `J_m` (≥ `k/(D+2)` positions) the vectors
   `(ug, vg, 2um, vm, 0)` annihilate every row. `verify.py` certifies this
   exactly (`D = 12, k = 300`: 138 positions of rank ≤ 25). This does **not**
   refute spread: with the true constraint `Θ = Θ(A,B)`, the band moves at first
   order by `2(u + vt)⟨φ_t, g⟩` along that direction, and `⟨φ_t, g⟩` is
   linear in `t` and vanishes at most once (checked).
4. **Tube lemma (rigorous).** Fix `s₀` keys, their positions and wrap vectors
   in `P′`; let `d = 2D + 3 − rank(r_t)`. At most `(1 + 4q√D)^d` progressions
   fit. Proof: by Cauchy–Binet, `Π σ_i² = det(M_ρ M_ρᵀ) ≥ 1`, with
   `σ₁ ≤ 2nk²`, so `σ_ρ ≥ (2nk²)^{−n}`; the solutions lie within
   `τ ≤ 1/(4q)` of an affine subspace of dimension `d` (because `a < 1/2`);
   points of the lattice `Λ² = {(ỹ(n₀), ỹ(d))}` are `1/q`-separated, so a
   packing argument concludes.
5. **Proposition (rigorous reduction).** Let `a < 1/2`, `0 < b < 1`,
   `c > b/4`, `θ < 2 − b/(2c)`. Assume
   > **(G_θ)** there are ≥ N²/4 non-short-period progressions such that, for
   > each one and every set `J` of `s₀` positions,
   > `rank(r_t : t ∈ J) ≥ 2D + 3 − θD`.

   Then for large `k` every rich cover costs ≥ 1. No full rank is needed:
   a deficit of a fixed fraction of `D` is tolerated. For `|S| ≥ s₀` the
   count is `k^{s₀}(k+2)^{Ds₀}(1 + 4q√D)^{θD}`, whose logarithm is
   `O(k^{2a} log k) + θck log k`, below `(2c − b/2)k log k`.

So the rich-side barrier `c ≤ b/4 = γb/2` at `γ = 1/2` holds for every cover
bound with these keys, **conditional on (G_θ)**, a purely combinatorial
statement about Beatty vectors `w_t = round(A + tB)`. The only deficit
mechanism found is that of Dirichlet: `r` vectors `g` costing `2r` on a
common level set, which shrinks by ≈ 0.5–0.7 per added `g`. On `s₀` positions
that predicts `O(log k)`, far below `θD = θk^a`.

## Measurement of (G_θ) (scoped)

| progression | D | k | s₀ | deficit |
|---|---|---|---|---|
| affine, no coordinate wraps (positive control) | 8 | 400 | 40 | 16 = 2D |
| quasi-affine, 3 wrapping coordinates (positive control) | 8 | 400 | 40 | 10 = 2(D−3) |
| generic, seed 13: worst searched / 20 random sets | 8 | 400 | 40 | 5 / 0 |
| generic, seed 14 | 8 | 400 | 40 | 6 / 0 |

Domain: progressions with `A, B` drawn from dyadic rationals of denominator
`2^30`. Searched sets: intersections of the most frequent level sets of
`⟨g, w_t⟩` for the 40 best `g ∈ {−1,0,1}^D` (smallest `|⟨g,B⟩|`), up to 6
vectors, truncated to `s₀`. Stopping rule: 400 sets per level. Rank exact
modulo `2^61 − 1` (a lower bound for the rational rank, so an upper bound
for the deficit). Coverage: not exhaustive, and `D ≤ 8` is far from
`D = k^a`. "No worse deficit found" is not a proof. The positive controls
show the measurement does see clusters where they exist.

A brute force at small `q` would not test the right regime: with `q ≲ 400`
and `D ≤ 3` there are more keys than elements, whereas everything hinges on
`#keys = N^{o(1)}`.

## Reproduce

```sh
python3 -B verify.py
```

Expected: `OK: all checks pass and match RESULT.json` (≈ 7 s, standard
library only). It replays: the exact Dirichlet certificate (three instances),
the first-order check, the four deficit measurements, on small instances the
fibre bound, the `1/q` separation and the Gram determinant ≥ 1, and the
three regimes of the proposition at `a = 49/100, b = 97/200, c = 13/100,
θ = 13/100, k = 10^400` (margin `6·10⁻⁴` in units of `k log k`, which is
`2c − b/2 − θc`).

## What is not claimed

No lower or upper bound on `W_r(k)`. (G_θ) is open. The barrier only
concerns cover bounds with these keys. A failure of (G_θ) would be good
news: it would yield cheaper rich covers.
