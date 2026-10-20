# Structured flips on single-label progressions: what a frozen key set must hit

**Status:** structural research. One partial positive result (a structured set
passes (I) on the families `P1` and `P3`), one negative result (an exact
family that defeats it and every rule of its kind), and one exact reduction.
No new bound is claimed, and no `X` passing (I) on all of case (ii) is known.

**Base commit:** `63bef9ca996786a41169de4a0c3d6b3b47b0fe96`.

**Sources:** the OpenAI manuscript *Quantitative Superexponential Bounds for
van der Waerden Numbers*, pinned in `source-manifest.json` (commit
`fd4aeeb`; Section 2 for the meshes and dilation, Section 6 for the keys).
Also `short-period-runs-proof.md` (case (ii)) and the single-label families
`P1`, `P2`, `P3` of the open contribution
`short-period-sensitivity/2026-10-20-ChaloupeMusic-single-label-obstructions`.

## Reproduce

```sh
python3 -B auto-research/structured-flips/2026-10-20-ChaloupeMusic-single-label-lenses/verify.py
```

Expected: exit code 0, last line `Test passé : all checks pass`, and a JSON
summary equal to `RESULT.json`. It takes about 20 seconds and uses only the
standard library.

## The question

Every proof in this repository colors `C(n) = c₀(L(n)) + E_{κ(n)}`, with iid
Bernoulli(`p`) flips `E` on the keys `κ(n) = (V-box of ỹ(n), ⌊‖qỹ(n)‖²⌋)`.
On case (ii) (short-period progressions) this costs `b < a` (band words) and
`a + b < 1` (heavy runs). The single-label families `P1` and `P3` show that
iid flips cannot avoid either cost.

Porte 1 freezes a set `X` of keys instead: `c₁ = c₀ ⊕ 1_X`, then iid flips on
the other keys. Two conditions suffice for every `c < 1/4`:

* **(I)** every case (ii) progression has a key of `X` at a position of its
  dominant color;
* **(II)** no rich progression is cheapened by `X`.

The precise form of (II) for the union bound is weaker than "few keys in `X`".
A rich `P` can be monochromatic in color `β` only if every key of `X` on `P`
sits at a position with `c₀ ≠ β`. If `X ∩ P` meets both colors of `c₀`, then
`P` is never monochromatic. The cost of target `β` is then
`p^{#(c₀⁻¹(¬β) ∩ P) − t}`, with `t = |X ∩ P|`.

## Results

### Lemma 1 (exact lens on a single-label block)

Take `n_j = r + q^{D−1}(m + je)` for `j < k`, with `0 ≤ r < q^{D−1}`, inside
one full label (so `qỹ_D` does not wrap). Put `s = qỹ_D(n₀)` and
`C = Σ_{i<D} (qỹ_i(r))²`. Then

```
band(n_j) = ⌊θ + v j⌋ + λ²e²j²,   θ = C + s²,   v = 2sλe,   frac(s) = y_{D−1}(r).
```

*Proof.* Lemma A of `P1/P3` gives `f(n_j) = C + (s + jλe)²`. The term
`λ²e²j²` is an integer and comes out of the floor. The fractional part of
`q y_D(n) = (λn mod q^D)/q^{D−1}` is `(λn mod q^{D−1})/q^{D−1} = y_{D−1}(r)`.
Centering changes `q y_D` by an integer. ∎

Consequently, on every single-label block, the whole real part of the band
word is a Beatty sequence `⌊θ + vj⌋`. The key does not see `θ` (through `C`)
nor `frac(v)` (through `y_{D−1}(r)`). Both range over essentially every real
value inside one box.

### Lemma 2 (exact shifted squares)

Let `D ≥ 3` and `0 ≤ u < q/2`, and let `r = u λ⁻¹ mod q^{D−1}`. Then on every
single-label block of this `r` with `qỹ_D ≥ 0`:

```
band(n_j) = x_j² + u²,   x_j = ⌊s⌋ + jλe ∈ ℤ.
```

*Proof.* `λr mod q^i = u` for `i ≤ D − 1`, so `qỹ_i(r) = u/q^{i−1}`. Hence
`C = u²(1 + q⁻² + …)` and `frac(s) = u/q^{D−1}`. The remainder
`u²(q⁻² + …) + 2 frac(s) x_j + frac(s)²` lies in `[0, 1)`, because
`u < q/2`, `x_j < q/2` and `D ≥ 3`. ∎

`P1` is the case `u = 0`. For `D = 2` the identity is false (the verifier
finds counterexamples).

### Lemma 3 (blocking the squares)

For an odd prime `ℓ`, let `SQ_ℓ` be the set of squares mod `ℓ`, 0 included.
It has `(ℓ+1)/2` elements. Some `R' ⊂ ℤ/ℓ` with
`|R'| ≤ ⌈log₂ ℓ⌉ + 1` meets every translate `c + SQ_ℓ`.

*Proof.* Take `s = ⌈log₂ ℓ⌉ + 1` and a uniformly random `s`-subset. It
misses `c + SQ_ℓ` with probability
`Π_{i<s} ((ℓ−1)/2 − i)/(ℓ − i) < 2^{−s} ≤ 1/(2ℓ)`. A union bound over the
`ℓ` translates finishes the proof. ∎

### Theorem A (a residue rule passes (I) whenever `ℓ ∤ e`)

Let `ℓ` be an odd prime with `ℓ ∤ λ` (for instance `ℓ > h₀`). Let `R'` be
as in Lemma 3, `R = R' + {−1, 0, 1}`, and `k ≥ (ℓ − 1)² + 1`. Put

```
X_ℓ = {(B, t) : t mod ℓ ∈ R}   (every box B).
```

Then every single-label block with `ℓ ∤ e` has a key in `X_ℓ`, whatever
`r`. In particular this holds on all of `P3` (`e = 1`). The key density is
`|R|/ℓ ≤ 3(⌈log₂ ℓ⌉ + 1)/ℓ`.

*Proof.* By Dirichlet, there are `1 ≤ β ≤ ℓ − 1` and an integer `a` with
`|βv − a| ≤ 1/ℓ`. For `0 ≤ i < ℓ`, the index `βi ≤ (ℓ−1)² < k` is
admissible, and `|(βv − a) i| < 1`. So
`band(n_{βi}) ∈ g(i) + {−1, 0, 1}`, where
`g(i) = ⌊θ⌋ + a i + (λeβ)² i²`.

Since `ℓ ∤ λeβ`, the leading coefficient of `g` is a nonzero square mod
`ℓ`. Completing the square shows that `g` takes exactly the values of a
translate of `SQ_ℓ` on `ℤ/ℓ`. So `g(i₀) ∈ R'` for some `i₀`, and then
`band(n_{βi₀}) mod ℓ ∈ R`. ∎

### Corollary A (`P1 ∪ P3` at `b_X → 1/2`)

Add a patch on exact squares: `P_H = {(B, x²) : x ≥ q/8, x ∈ H}`, where
`H = ⋃_{p∈𝒫} pℤ` and `𝒫` is a set of primes in `(h₀, k]` with
`Π_{p∈𝒫} p > E_max`. Then `X = X_ℓ ∪ P_H` has a key on every `P1` block
(`r = 0`, any `e ≤ E_max`) and on every `P3` block.

*Proof.* For `P1` with `ℓ | e`, some `p ∈ 𝒫` does not divide `e`, and
`x₀ + jλe ≡ 0 (mod p)` has a solution `j < p ≤ k`. Every other block is
covered by Theorem A.

Take `ℓ` prime in `(√k/2, √k]`. This prime exceeds `h₀ = D⌈log k⌉²` once `k`
is large, because `a < 1/2`. The patch is negligible: above `q²/64`, perfect
squares are at most a fraction `4/q` of the bands of any box. So the key density is `O(log k/√k)`, which
means `b_X = 1/2 − o(1) > a` on these two families.

**The iid obstructions `b < a` (from `P1`) and `a + b < 1` (from `P3`) do not
survive structured flips on those families.** The limit `√k` comes from the
Dirichlet step, not from the families.

*Measured, not proved.* With `k = 2ℓ`, an exhaustive enumeration covers every
real lens (all `θ`, all `v`, all nonzero squares `λ²e²`, all mechanical
words of length `k`). It finds that every lens covers at least `(ℓ+1)/2`
residues. It also finds that a set `R` of 4, 4 and 5 residues meets every lens
for `ℓ = 13, 17, 19` (and 5 residues for `ℓ = 23`, a run outside the
verifier). Four is minimal for `ℓ = 13`. If "at least `(ℓ+1)/2` residues per lens" holds for all `ℓ` with
`k ≥ 2ℓ`, the counting of Lemma 3 gives `|R| = O(log k)`, hence
`b_X → 1` on `P3` and on every family with `ℓ ∤ e`.

### Proposition B (the shifted family defeats every such rule)

Let `X` be any key set whose band rule in each box depends only on the band
mod `L`, with `L ≤ E_max`, plus any patch on perfect squares. Let
`1 ≤ u ≤ √q/2`. Then the family of Lemma 2 contains a single-label block with
`e = L` on which `X` is constant: either `X` misses all its keys, so (I)
fails, or `X` contains all of them, so the block is monochromatic for every
flip assignment.

*Proof.* `x_j ≡ x₀ (mod L)`, so `band(n_j) = x_j² + u² ≡ x₀² + u²`. A number
`x² + u²` with `x ≥ q/8 ≥ u²/2` is never a perfect square: `z² = x² + u²`
with `z > x` would give `u² = (z − x)(z + x) ≥ 2x + 1 > u²`. Good blocks
with `e = L` exist by the argument of Lemma C of `P1/P3`, as for `P2`; the
verifier exhibits one for `u = 1`, `L = 13`. ∎

This is `P2` shifted to every `u`, and immune to square patches. On the real
instance (`q = 2⁶¹ − 1`, `D = 3`, `k = 145`, `ℓ = 13`), the verifier exhibits a
single-label block with `u = 1` and `e = 13` that misses `X_ℓ ∪ P_H`
entirely.

The rule of Corollary A therefore fails (I) on the real single-label family
`{r arbitrary, ℓ | e}`. Unions of residue rules have to cover the
`e`-direction, where `e` runs up to `E_max = q^{1−o(1)}`. Each prime of the
union handles only `ℓ ∤ e` and costs `≈ log k/ℓ`. The total cost is then
`≈ log E_max / k = k^{−a+o(1)}`, which is `b_X ≤ a` again. This is a count,
not a theorem. A general version for residue rules at `r ≠ 0`, CRT unions
included, is in preparation in a companion contribution.

### Lemma R (the exact reduction behind (I))

Fix `r` and a maximal interval `J` of `m` on which `r + q^{D−1}m` keeps one
full label. Then (I) on the single-label blocks inside `J` says exactly that
`S_J = {m ∈ J : κ(r + q^{D−1}m) ∈ X}` meets every `k`-AP contained in `J`.
This is a restatement: those APs are exactly the blocks.

A box is covered by such intervals as `r` varies, of length up to
`qH/λ = q/(k²λ)`. Averaging therefore shows that the point density of `X`
in every box of the `P1` region is at least `ρ_AP(Y, k)` over the long
intervals. Here `ρ_AP(Y, k)` is the least density of a set meeting every
`k`-AP contained in `[Y]`, and `log Y = (1 − o(1)) log q =
(c − o(1)) k^{1−a} log k`.

Known constructions give `ρ_AP ≈ log Y/(k log k)`: unions of residue classes
modulo primes up to `k`, or "some base-`k` digit is 0". That is `k^{−a+o(1)}`,
the same exponent as random. So:

* **Porte 1 through a sparse `X` needs a new additive-combinatorics object:**
  sets of density `k^{−1/2−δ}` meeting every `k`-AP of
  `[exp(Θ(k^{1−a} log k))]`, a power of `k` below every known construction.
  On top of that, the pullback has to work for every lens at once.
* A lower bound `ρ_AP(Y, k) ≥ log Y / k^{1+o(1)}` would close this route in
  the `P1` region. The remaining route would then be a dense `X` (density
  `k^{−a}`) that is anti-random on rich progressions: `X ∩ P` must
  essentially never equal the `c₀`-minority of a rich `P`, while a random
  set of that density does so on about `N² k^{−aγk}` progressions, which is
  `≫ 1` as soon as `2c > aγ`.

The analogous families moving in coordinate `i < D` (`d = q^{i−1}e`) carry
the same lens up to a quadratic perturbation below `10⁻⁵`. They are not
checked here.

### On (II)

On the real instance, 2000 uniformly random progressions of `ℤ/N` (generic
rich ones) meet `X_ℓ` with the binomial mean and variance of density `7/13`.
The patch is never met. This is the typical behaviour only. The tail, which
is what (II) needs, cannot be tested at this size.

A rule that is random-like on rich progressions gives at best
`c < γ min(b_X, b)/2`. So the rigorous `b_X → 1/2` of Corollary A would sit
exactly at `1/8`, even if it passed (I) everywhere.

## Checks in `verify.py`

1. Lemma 1 on 534 positions (random `q`, `D`, `λ`, `r`, `m`, `e`, Fractions).
2. Lemma 2 on 979 positions (`q` up to `2⁶¹ − 1`, `D = 3, 4`). Control:
   `D = 2` breaks it.
3. Lemma 3: stored `R'` for 10 primes (exhaustive), and the counting
   inequality for the 665 odd primes below 5000. Control: no single residue
   blocks.
4. Theorem A on 1500 random lenses (`ℓ = 13`, `k = 145`), with the Dirichlet
   witness of the proof recomputed. Control: an `ℓ | e` lens misses `R`.
5. Measured table, `k = 2ℓ`, exhaustive over all real lenses, for
   `ℓ ∈ {13, 17, 19}`. Control: no 3-element set works for `ℓ = 13`.
6. The real construction (`q = 2⁶¹ − 1`, `D = 3`, `λ = 840`, `k = 145`, the
   source's `U` mesh and adaptive `V` mesh with 200-digit decimals):
   * 30 `P3` and 30 `P1` single-label blocks meet `X_ℓ ∪ P_H`;
   * the `P1` blocks with `13 | e` and `x₀²` outside `R` are met only through
     the patch;
   * the killer block of Proposition B (`u = 1`, `e = 13`) misses
     `X_ℓ ∪ P_H` entirely.
7. Typical rich behaviour (2000 progressions), as above.

The toy instance uses `λ = 840`, chosen coprime to `13`, and not the dilation
the source would compute for `k = 145`. The lemmas do not depend on the value
of `λ`. The finite checks illustrate the lemmas; the proofs are the text
above.

## What is not claimed

* No `X` passing (I) on all of case (ii) is known, at any density below
  random.
* No lower bound on `ρ_AP` is proved. "No construction known" is not an
  impossibility theorem.
* The measured `k = 2ℓ` table is a finite enumeration for three primes, not a
  proof for general `ℓ`.
* Nothing is claimed about (II) beyond typical behaviour.
* Nothing changes for the current results: `1/9` and every `c < 1/8` stand as
  they are.
