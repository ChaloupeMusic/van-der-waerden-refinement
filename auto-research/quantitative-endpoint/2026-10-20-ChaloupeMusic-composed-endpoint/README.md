# Composed endpoint bound: seven-eighths loss with retuned constants

**Status:** conditional bound with exact checks. It combines two contributions built on open PR #12 and
depends on them. Independent review is pending.

```math
W_r(k)>\left(\frac{k^{1/8}}{e^{403/100}(\log k)^{7/8}}\right)^{k\lfloor\log_2 r\rfloor}
\qquad(k\ge2^{739},\ r\ge2),
```

with `W_r(k) > k^{(k/9)⌊log_2 r⌋}` for `k >= 2^1016`, and a coefficient above `0.1112` at `k = 2^1024`.
Both are conditional on PR #12's Section 2 interfaces and PR #9's displayed finite factors.

| | bound on `W_2(k)^{1/k}` | from | `c = 1/9` from |
|---|---|---|---|
| PR #12 | `k^{1/8} / (e^{20} log k)` | `2^1024` | `2^4096` |
| `2026-10-20-ChaloupeMusic-endpoint-constants` (N1a) | `k^{1/8} / (e^{9/2} log k)` | `2^1024` | `2^1163` |
| `2026-10-20-ChaloupeMusic-seven-eighths` (N1b) | `k^{1/8} / (e^{35/2} (log k)^{7/8})` | `2^739` | `2^2504` |
| **this note** | `k^{1/8} / (e^{403/100} (log k)^{7/8})` | `2^739` | `2^1016` |
| this note, without Lemma D′ | `k^{1/8} / (e^{35/8} (log k)^{7/8})` | `2^739` | `2^1055` |

Each version beats every earlier bound on that bound's own range, even before the `(log k)^{1/8}` factor.
No coefficient above `1/8` is claimed.

## What is combined

- **N1b's architecture.** A remainder-free canonical partition of long-period progressions into blocks of
  lengths in `[h/2, h]`. With it, `M = ⌊k/32⌋`, `D ∝ sqrt(k)/t^{5/2}` and `p ∝ t^{7/2}/sqrt(k)`.
- **N1a's five scale-free lemmas**, restated in [PROOF.md](PROOF.md), Section 2: the dilation bound (A),
  local-lemma weights `Pbar·e^{1/8}` (B), the asymmetric rich split (C), PR #9's products in `5^{k/2}` and
  `37^k` (D), and the s-light threshold `k/16` (E).
- **One new arithmetic lemma (D′).** PR #9's class product is bounded by `19^k` instead of `37^k`, using
  each class's own length `K_a` and `K_h <= K_a + 1`. It lowers the constant from `35/8` to `403/100`.
  Without D′ the proof gives `e^{35/8}`, and `verify.py` checks that version too.
- **Eligibility threshold.** The anchored patterns of N1b's (I2) are asked to have at least
  `(1 − ε_reg)ℓ` regular positions instead of `ℓ/2`. Canonical blocks satisfy this by N1b's Lemma 5.2.
  It halves the second-family tail exponent, so `h_0 = ⌈5Dt^3⌉` suffices.

## Reproduce

From the repository root:

```sh
python3 -B auto-research/quantitative-endpoint/2026-10-20-ChaloupeMusic-composed-endpoint/verify.py
```

Expected: exit code 0, and standard output byte-identical to `RESULT.json` (under a second, standard
library only). The verifier checks:

- 26 monomial certificates, valid for every real `t >= 512`;
- every budget of PROOF.md in exact rationals, for both versions, with certified logarithm bounds;
- Lemma D′, against `19^k` and against `37^k`, for every `k <= 300` and `h < k`;
- the canonical partition and the counting core of the short-block lemma;
- a 90-digit replay of 15 inequalities at 44 values of `t`, for both versions;
- the hashes of the package files.

Negative controls rejected by `verify.py`:

- `AD = 21`: the first family fails;
- `BH = 4`: the second family fails;
- `AP = 35`: case (ii) fails;
- `C = 4.02` with D′, and `C = 4.35` with D;
- `37^k` at `C = 4.03` (so D′ is load-bearing);
- N1b's `ℓ/2` eligibility threshold with `BH = 5`;
- `G = 2` (`γ` ignoring `ε_reg`);
- a monomial claimed from `t_0 = 30`.

Mutation tests, run by hand on a copy whose `MANIFEST.json` was regenerated (with the delivered manifest,
any edit of `verify.py` fails first on its own hash):

| Mutation | Inequality that fails |
|---|---|
| `AD = 21` | first-family budget |
| `BH = 4` | second-family tail against its entropy |
| `AS = 1024` | third-family budget |
| `AP = 35` | case (ii) budget |
| `AP = 37` | rich budget |
| `C = 4.02` | rich budget (main) |
| `C_FALLBACK = 4.35` | rich budget (fallback) |
| `CLASS = 37` | rich budget (main) |
| `G = 2` | dichotomy identity |
| `E2 = 1/2` | rich budget |
| `T0 = 400` | `ρ <= 2^{−200}` |

## Track index

The track README is not edited, because PR #12's `MANIFEST.json` pins it. The proposed line is:

```
2026-10-20, `ChaloupeMusic`: [Composed endpoint bound](2026-10-20-ChaloupeMusic-composed-endpoint/README.md).
Status: conditional bound with exact checks, combining the two 2026-10-20 contributions. Independent review is pending.
```

## Scope

- Imported, as stated by their authors: N1b's Lemmas 5.1–5.3 and its interface extensions (I1)–(I2),
  N1a's Lemmas A–E, PR #12's Section 2 interfaces and Family 3, and PR #9's displayed factors.
- Lemma D′, like N1a's Lemma D, re-bounds a displayed product of PR #9. It does not re-derive the encoding
  that produces it.
- The local lemma, Hoeffding's inequality and PR #9's encodings are not formalized.
- Finite checks corroborate the proof and are not its premises.
- N1a and N1b are separate packages on separate branches. This one cites them and does not contain them.
