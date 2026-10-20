# Track: structured-flips

Goal: replace the iid flips of the current template, on part of the keys, by a
deterministic set `X` of frozen keys (`c₁ = c₀ ⊕ 1_X`, iid flips elsewhere),
so that the short-period case (ii) stops costing `b < a`. If `X` meets every
case (ii) key word on its dominant color (I) and never covers the minority of
a rich progression (II), every `c < 1/4` follows. This is the only known
route past `1/8` with the current keys (Porte 1).

Open questions:

* a set `X` that passes (I) on the whole single-label family at key density
  `k^{-1/2-δ}`, not only on the sub-families `P1` and `P3`;
* the minimal density of a set meeting every `k`-AP of `[Y]` when
  `log Y = k^{Θ(1)}`: a power of `k` below `log Y / k` would be needed;
* a dense `X` (density `k^{-a}`) whose large deviations on rich progressions
  beat those of a random set.

Related: `template-ceilings` (the dual, asymmetric flips), the
`short-period-sensitivity` contributions on `P1`, `P2`, `P3` (open PRs), and
`barrier-note.md`.

## Index

| Date | Handle | Result | Status | PR |
|---|---|---|---|---|
| 2026-10-20 | ChaloupeMusic | Exact lens identity on single-label blocks. A residue rule passes (I) on `P1 ∪ P3` at key density `O(log k/√k)` (`b_X → 1/2 > a`). The exact shifted family `x² + u²` defeats every rule of this kind. (I) reduces, fiber by fiber, to `k`-AP hitting sets in `[q^{1-o(1)}]` | structural research (partial positive, negative, reduction) | – |
