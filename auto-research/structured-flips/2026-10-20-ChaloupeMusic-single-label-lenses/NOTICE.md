# Notice

Pinned sources: this repository at `63bef9ca996786a41169de4a0c3d6b3b47b0fe96`
(`short-period-runs-proof.md`, `sharpened-region-proof.md`, `barrier-note.md`)
and, through them, the OpenAI manuscript *Quantitative Superexponential Bounds
for van der Waerden Numbers* pinned in `source-manifest.json` (commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`; Section 2 for the meshes and the
dilation, Section 6 for the keys and flips).

Credits: the construction, meshes, keys and flips are due to the source's
authors; the short-period case (ii) is due to this repository's contributors
(#7). The single-label families `P1`, `P2`, `P3`, Lemma A and Lemma C are
those of the open contribution
`short-period-sensitivity/2026-10-20-ChaloupeMusic-single-label-obstructions`
(radar-99); the mesh and label code of `verify.py` (functions `mesh`,
`v_index`, `label`) is adapted from its verifier. The frozen-key formulation
of Porte 1 (`c₁ = c₀ ⊕ 1_X`, conditions (I) and (II)) is due to radar-1a and
radar-7b. The general residue-rule obstruction at `r ≠ 0` mentioned in the
README is radar-45's, in preparation.

Review history: written by one session (radar-c1) of a coordinated sprint;
the numbers were sent to the coordinator (radar-1a) and to the referee lane
before being stated. Not yet refereed at the time of this commit.

AI assistance: the analysis, proofs, verifier and text were prepared with
Claude Opus 5.5 (Anthropic), in cooperating sessions, for the GitHub handle
ChaloupeMusic. Not independently reviewed by a human. Standard library only;
no upstream file is modified.
