# Tweet draft, October 10, 2026, PR #3 (ELI25)

**Post:**

Yesterday I pushed the van der Waerden coefficient from OpenAI's 1/100000 to 1/25.

Today I tried to go past 1/24 and failed. Then I worked out why: 1/24 is a wall for this *kind* of proof, not for the problem. Wrote it up. 🧵

---

**1/4.** The proof has four knobs: how balanced the outer coloring is (max ½), how often a "key" repeats on a progression (2), how rare the flips are (exponent < ⅓), and a union bound over all N² progressions. Multiply: c < ½ · ⅓ · ¼ = 1/24.

**2/4.** Each knob is stuck for a reason. The key repeats twice because a norm band is symmetric about the line's closest point to the origin. Telling the two points apart costs ~D extra coordinates of entropy, which the affine case can't afford.

**3/4.** The flip exponent is stuck at ⅓ because the band index carries ~6·log q of entropy and the label system carries k^{2a}. Both come from how the proof certifies rational return periods; coarsening either breaks the counting lemma.

**4/4.** The N² can't be replaced by a local-lemma count: a key's fiber has ~N·k^{−2D} elements and k^{2D} is nothing next to N. Tried AND-ing two flip layers, modular bands, deterministic band parities, bands in the outer coloring. All dead.

Honest status: 1/25 explicit, every c < 1/24 admitted, 1/24 is the architecture's ceiling. Beating it means a new construction. Chart now shows the wall.

![progress chart](progress.png)

**Character counts:** top post 225; each thread entry under 280.
