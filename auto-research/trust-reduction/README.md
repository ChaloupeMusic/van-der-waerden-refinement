# Trust reduction

Goal: shrink what each coefficient depends on. Every refinement in this
repository is conditional on correctness, not on a conjecture: it assumes the
OpenAI source's lemmas and adds arguments reviewed only by AI-assisted source
reviews. A contribution here proves an added link in full, or shows that a
coefficient does not need it.

Open questions:

- A complete write-up of the rich-event chain (PRs #4, #6, #7) at the level of
  the source, which would carry the same treatment from `1/16` up to `1/9`.
- A Lean statement of the parameter transfer (the side-condition table).

Related tracks, in open PRs #13–#18: `independent-audit`, `short-period-sensitivity`.

## Index

| Date | Handle | Result | Status | PR |
|---|---|---|---|---|
| 2026-10-26 | ChaloupeMusic | every c < 1/16 (explicit 1/17) from the source's lemmas alone; (O), (R), (P2) proved, (K) not used | research note | — |
