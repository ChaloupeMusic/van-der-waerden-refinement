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

**No-go for the natural class (cell colorings chosen at random or by a local
lemma).** Suppose `A` is a union of cells of a partition of the `y`-torus into
boxes, the cells being colored by a union bound or a local lemma over events
indexed by cell words, as in the source.

1. *(Minkowski, sketch.)* There is an absolute `C` such that a box cell
   holding more than `(Ck)^D` points of `{ỹ(n)}` contains a `k`-term
   progression. Sketch: the representatives form a lattice of covolume `1/N`
   whose fundamental domain has diameter `≈ √D/q`. In coordinates where the
   cell is narrower than the lattice spacing, restrict to one layer, a coset
   of the sublattice stationary in those coordinates. In the others the
   cell's middle half-box contains a lattice point, and pigeonholing the
   points over `(4k)^{D′}` sub-boxes gives a nonzero lattice vector `v` with
   `|v_i| ≤ ℓ_i/(4k)`, so `n + jv` stays in the cell for `j < k`. The layer
   bookkeeping is not written out: a complete proof is an open obligation,
   and `verify.py` only checks the statement exactly on small instances
   (where the step `d = 1` already works). The constant plays no role below.
   So every 1-cell of a `k`-AP-free `A` has `≤ (Ck)^D` points.
2. *(Entropy.)* A cell of `≤ (Ck)^D` points determines its point up to
   `(Ck)^D`. A generic progression with two points in `A` is therefore fixed,
   up to `(Ck)^{2D}k²`, by its cell word: there are `≥ N^{2−o(1)}` distinct
   words, and `≥ N^{1−o(1)}` through any one 1-cell. Every event about the
   colors of `k` cells has probability `≥ e^{−O(k)}`, so both the union bound
   and the local-lemma condition (incident weight per cell `≪ 1`) fail, since
   `log N = ck log k ≫ k`.

What survives is a deterministic or algebraic choice of `A`, outside every
probabilistic selection used in this repository. (Bohr sets are one convex
cell and fall under 1. Band-level sets `{n : band(n) ∈ Λ}` have cell words of
entropy `≈ 2 log N` on generic progressions, so only a deterministic hitting
theorem for pseudo-random band sets could serve; none is known.)

Equivalently, the exception asks for an outer coloring whose monochromatic
progressions are **all of color 0** and all in a family of event entropy
`≪ k^{1/2}`: color 1 must be `k`-AP-free (indeed `θ`-sparse on every
progression) while still meeting every non-structured progression. The
current outer coloring has monochromatic structured progressions of both
colors (one per label cell), which is exactly what pins the two rates to the
same `k^{−1/2}`. Three further obstructions, all heuristic, for where such an
`A` could come from:

* **tiling:** step-1 progressions tile `G`, so some structured progression sees
  `A` at least at its global density; a set denser on generic progressions than
  on every structured one does not exist;
* **deterministic hitting** of all `N²` generic progressions is only available
  through a local lemma over low-entropy label words, i.e. coarse cells, which
  contain whole progressions; finer structure is random-like on generic
  progressions, and a random-like `A` misses `N²(1 − μ)^k ≫ 1` of them unless
  `μ ≳ 2c log k`;
* **geometric regions** `A = ỹ^{−1}(R)`: AP-freeness forbids `R` from
  containing an orbit that is slow in any subset of coordinates (cylinders,
  slabs, small boxes all fail), while `k ≪ 2^D` points cannot equidistribute
  in dimension `D = k^a`, so no discrepancy argument forces every generic orbit
  into a thin `R`.

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
