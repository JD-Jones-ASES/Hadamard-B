# Audit — port lane b, 2026-09-05 (the general branch, and the product theorem)

Skeptic lane `hb-audit-b`, on `desk/ports-turn54-b` at `94853e6`. Mandate:
refute it. Nothing enters `note/NOTE-B.md` on the author lane's say-so — every
theorem below was re-derived here from §1.0–§1.1 before the author's proof was
graded, every certificate was replayed from the worktree root, and one banked
number in each certificate was tampered in a scratch copy to confirm it fails
closed.

**Verdict: the mathematics survives.** Every theorem, lemma and corollary in
§§1.12–1.13 is true as stated; every count in §2.6 and in certs 30/31/32
reproduces under an independent third implementation written for this audit;
all three certificates pass and all three fail closed under tampering. Six
findings, none of them a false theorem: one unstated hypothesis in a proof
(F1, with a two-line repair that removes the hypothesis outright), one
scope sentence that is **refuted as worded** by §1.12's own scope paragraph
three paragraphs later (F2), one dangling prior-art pointer (F4), and three
editorial (F3, F5, F6). Recommendation: **REPAIR-FIRST** — F1, F2 and F4
before any push.

---

## 1. Verdict table

| # | item | label claimed | verdict |
| --- | --- | --- | --- |
| 1 | §1.12 Theorem F (general branch) clause **(a)** — `QᵀQ = 4i·I`, `col(Q) = ℝ⁴⊗V_S` exactly | PROVEN | **SOUND** — re-derived §2.1 |
| 2 | §1.12 clause **(b)** — `PPᵀ = 4i·I`, rows blockwise in `V_S`, `PᵀP = I₄⊗M` | PROVEN | **SOUND-WITH-REPAIR** — the printed proof ("the computations of §1.8 verbatim") silently imports **D3**, whose hypothesis `w > 2s` is not in §1.12's setting. **F1.** The conclusion is true with no `w` hypothesis at all; a trace argument is given in §2.2 |
| 3 | §1.12 clause **(c)** — `ĈᵀĈ = ĈĈᵀ = I₄⊗dev(F)`, `Ĉ` commutes with `I₄⊗P_S` | PROVEN | **SOUND** — re-derived §2.3; conjugation-closure of `S` is exactly what the commutation needs, and the note says so |
| 4 | §1.12 clause **(d)** — `E = −(1/4i)·PĈᵀQ` | PROVEN | **SOUND** — one line from (H4) and (a) |
| 5 | §1.12 clause **(e)** — `EEᵀ = 4s·I` automatically; (H3)+(H4) ⟺ `E ∈ {±1}^{4s×4s}` | PROVEN | **SOUND** — re-derived §2.2; the substantive step `F̂(χ) = N − 4iw = 4s` on `S` is correct and is correctly identified as the place the house form was used in §1.8 |
| 6 | §1.12 clause **(f)** — `E` sees `σ` only through the S-part | PROVEN | **SOUND**, and correctly *weakened*: §1.8(f)'s kernel-dimension clause (`4(i−s) = 4`) is `i = s+1`-only and was dropped |
| 7 | §1.12 clause **(g)** — transport across `w`, hence across orders | PROVEN | **SOUND**; §1.8(g)'s numeric row-sum clause `Σ_q r_q² = N` is house-only (in the general branch it is `4s` when `χ₀ ∈ S`) and was correctly compressed away |
| 8 | §1.12 the reread of "blockwise in `V_S`" when `χ₀ ∈ S`; the **dropping** of §1.8(b)'s forced-completion clause | — | **SOUND, and the load-bearing edit of the port.** §1.8's "each block sums to zero" and its `h_J ⊗ 1_i` completion are both false when `χ₀ ∈ S`; both are gone. This is where a careless port would have broken |
| 9 | §1.12 **row-factorisation lemma** `WWᵀ = 4i·4s·(I₄⊗P_S)`; `‖e‖² = 4s` for every admissible `p`; kit ⟺ `4s` mutually orthogonal flat rows | PROVEN | **SOUND** — re-derived §2.4, including the converse leg |
| 10 | §1.12 **dual form** `p` flat ⟺ `Ĉp ∈ Q·{±1}^{4s}`, and `P = −(1/4s)·E·Qᵀ·Ĉ` | PROVEN, "re-derived here" | **SOUND** — re-derived §2.5. The proof printed is complete and is this repository's; the provenance claim is honest |
| 11 | §1.12 the `(7,12)/S3` **border kill** at every order `N = 4(12w+7)`, and the `χ₆` twin with it | PROVEN + PROVEN-BY-CERTIFICATE | **REFUTED-AS-WORDED** (**F2**) — "**no `Q ∈ {±1}^{48×28}` satisfies (H1)**" is, on the plain reading of (H1) as §1.1 states it (which quantifies `M` existentially), contradicted by §1.12's own scope paragraph and by cert 31 [C], which exhibit an admissible `Q` at `(7,12)/S2`. The *tier-scoped* statement is **SOUND** and is independently reproduced here (§2.6): alphabet `= {±𝟙}`, `M_{S3}` off-zero values `{−8,−4,8}`, twin `= {1,4,5,6,7,8,11}` and its two alternating vectors. §5's row already states it correctly. Repair R2 |
| 12 | §1.12 **scope paragraph** — no kit known in the general branch at 2092; the `Q` quantifier open at S1/S2/S4, empty at S3; nothing about `H(2092)` | — | **SOUND**. The source lab's 40 494-point negative is reported as "tens of thousands", called a bounded negative search, and explicitly *not* certified here — no strengthening |
| 13 | §1.13 **Theorem G** `a_{2k} = −(b₀c₀d₀)·a_k b_k c_k d_k` | PROVEN | **SOUND** — the group-ring proof re-derived line by line in §2.7. Every step checks, including the two that carry the prefactor (`U* = U + J + e` for the skew seed; `(−1)^{|V|} = v₀ = b₀` for a symmetric one) |
| 14 | §1.13 **closure sign law** `Π_j π_j = closure·(−(b₀c₀d₀))^L` | PROVEN | **SOUND** — re-derived §2.8 |
| 15 | §1.13 the arithmetic at `n = 7` (`L = 3`, `2³ = 1 ∈ M₀`, `closure = +1`) and at `n = 523` (`87` orbits, `L = 87`, `2⁸⁷ = 463`, `−463 = 60 ∈ M₀`, `closure = −1`) | computed | **SOUND** — recomputed independently here (§2.9). `523` prime, `|ℤ*₅₂₃| = 522 = 2·3²·29`, `M₀ = {1, 60, 462}` |
| 16 | §1.13 census `528/288/1440/3456/4224 = 9 936`; each unprefixed variant on exactly half | COMPUTATIONAL-EVIDENCE (cert 32) | **SOUND** — `n = 7, 9, 11` re-enumerated here by brute force with no meet-in-the-middle: `528 / 288 / 1440`, prefixed form in **all**, each unprefixed variant in **exactly half** (264 / 144 / 720) |
| 17 | §1.13 the Williamson control — 14 784 quadruples, constant-product corollary in all, doubling relation in **none** | COMPUTATIONAL-EVIDENCE (cert 32) | **SOUND** — `n = 7, 9, 11` re-enumerated here: `960 / 2112 / 1920`, constant-product in all, prefixed doubling in **0**. The author's correction to the source decision is confirmed (§2.10) |
| 18 | §1.13 "no novelty claimed"; arXiv:1811.05094 REPORTED-UNVERIFIED "(§4)" | — | **SOUND in stance, broken in reference** (**F4**): §4 contains no good-matrices entry and no `1811.05094`, so both `(§4)` pointers in §1.13 and the one in §5 dangle, and an unread source is asserted about outside the section whose stated rule is that unread sources are declared there |
| 19 | §2.6 σ-layer: the six tables, the two support-completeness claims, the `(3,8)` coarsening | PROVEN-BY-CERTIFICATE (cert 30) | **SOUND** — all three enumerations reproduced independently (§2.11): **exactly one** Galois-stable size-11 support of `Ẑ₁₆` containing `χ₀` (`{0,1,3,4,5,7,9,11,12,13,15}`, = the banked one); **exactly four** such size-7 supports of `Ẑ₁₂` (= the four banked, with `S2`, `S3` as named); **all 56** size-3 subsets of `𝔽̂₂³` and **all 12** Galois-stable size-3 subsets of `Ẑ₄×Ẑ₂` are real twists of `W∖{1}`. The house `(3,4)` target `[4w+12, 4w, 4w, 4w] = [532, 520, 520, 520]` at `w = 130` recomputes from §1.8's profile |
| 20 | §2.6 `(7,12)` column layer: alphabet sizes `66, 32, 2, 24` of rank `7, 7, 1, 7`; `S2` table admissible; multiset forced | PROVEN-BY-CERTIFICATE (cert 31) | **SOUND** — the four alphabets and their ranks recomputed here by a third implementation (projector fixed point over all `2¹²`, exact `Fraction` rank): `66, 32, 2, 24` and `7, 7, 1, 7`, in the note's order |
| 21 | §2.6 "what the two certificates do not say" | — | **SOUND** — seeds, kits, counts and `H(2092)` all disclaimed, in that order |
| 22 | §5 label table — five rows appended | — | **SOUND**; §5's rows are more precise than §1.12's prose (see F2). Every label used is earned: PROVEN has a written proof, PROVEN-BY-CERTIFICATE has a cert that runs green here, COMPUTATIONAL-EVIDENCE is on the censuses only |
| 23 | §2.4 forward pointer | — | **SOUND** — a dated parenthetical, nothing else in §2.4 touched |
| 24 | cert 30 `certs/30-general-branch-sigma` | — | **SOUND** — replays, tampers closed |
| 25 | cert 31 `certs/31-s2-column-table` | — | **SOUND** — replays, tampers closed |
| 26 | cert 32 `certs/32-good-product-theorem` | — | **SOUND** — replays, tampers closed. On the desk question the lane raised: cert 32 **should stay**. DISCLOSURE line 98 says "Every computational claim in the note carries a certificate"; §1.13's counts are computational claims; without cert 32 that sentence goes false |
| 27 | README, PROVENANCE, DISCLOSURE, pre-existing text | — | **SOUND-WITH-REPAIR** — see F3/F5; `git diff main --stat` is **5 533 insertions, 0 deletions** across 12 files, so no pre-existing line was altered anywhere |
| 28 | non-vocabulary labels `DEAD-BY-CERTIFICATE`, `REPORTED-UNVERIFIED` | — | **F3** — neither exists on `main`; NOTE-B's preamble names a closed six-label vocabulary |
| 29 | section numbering `§1.8 → §1.12`, `§2.4 → §2.6` | — | **F5** — the gaps are lane a's reserved numbers. Harmless on merge, a visible defect if this branch ships alone |

Nothing checked here is **NOT-CHECKED** except the two things the port itself
declines to certify — the source laboratory's 40 494-point flat-row census
(reported, not claimed) and the `n = 13, 15` legs of cert 32's enumeration,
which I replayed inside the certificate but did not re-enumerate by a second
program (the `n = 7, 9, 11` legs I did, and they agree).

---

## 2. Re-derivations

Written from §1.0–§1.1 only. Notation as there: `Ḡ = G/K` of order `i`,
`w = |K|`, `N = 4(wi + s)`, `V_S = span{v_χ : χ ∈ S} ⊂ ℝ^i`,
`P_S = (1/i)Σ_{χ∈S} v_χ v_χ*`, `M = M_S = 4i·P_S`, `Π := I₄⊗P_S`,
`F = Σ_q PAF_{σ_q} = N·δ₀ − w·M`.

First, the arithmetic the whole section turns on. `M(c) = 4Σ_{χ∈S} χ(c)` and
`(4i·P_S)[c,c'] = 4Σ_{χ∈S} χ(c)χ̄(c') = 4Σ_{χ∈S} χ(c−c')`, so `M = 4i·P_S`
is the same object read two ways; it is real because `S` is conjugation-closed
and has rank `|S| = s` because distinct characters are orthogonal. Hence
`M̂(χ) = 4i` on `S` and `0` off it, and

> `F̂(χ) = N − w·M̂(χ) = N − 4iw = 4(wi+s) − 4iw = 4s` for every `χ ∈ S`.

That single value `4s` is what replaces the house computation in §1.8, and it
is the only place `M = 4i·P_S` is consumed. The author identifies exactly this
step; the identification is correct.

**2.1 (a).** `QQᵀ = I₄⊗M` has nonzero spectrum `{(4i)^{4s}}`. `QᵀQ` is
`4s × 4s` with the same nonzero spectrum, so all `4s` of its eigenvalues equal
`4i` and `QᵀQ = 4i·I_{4s}`. Then `Q/√(4i)` has orthonormal columns, so
`QQᵀ/4i = Π` is the orthogonal projector onto `col(Q)`; `im Π = ℝ⁴⊗V_S`;
equal projectors have equal images. So `col(Q) = ℝ⁴⊗V_S` **exactly**, and in
particular every `i`-block of every column lies in `V_S`. Note what is *not*
available here and is available in §1.8: `1ᵀM1 = i·M̂(χ₀) = 0` fails when
`χ₀ ∈ S`, so "each block sums to zero" is false in the general branch. The
port says so. ∎

**2.2 (b), (e) — and the repair for F1.** From (H4) `EQᵀ = −PĈᵀ`; right-
multiplying by `Q` and using (a) gives (d), `E = −(1/4i)PĈᵀQ`. Then

```
  EEᵀ = (1/16i²)·PĈᵀ(I₄⊗M)ĈPᵀ = (1/4i)·PĈᵀΠĈPᵀ
      = (1/4i)·PΠĈᵀĈPᵀ            [ (c): Ĉ commutes with Π ]
      = (1/4i)·PΠ(I₄⊗dev F)Pᵀ = (4s/4i)·PΠPᵀ = (s/i)·A,   A := PΠPᵀ .
```

§1.8 now invokes **D3** to get `PPᵀ = 4i·I`, and D3's hypothesis is `w > 2s`.
That hypothesis is nowhere in §1.12's setting — the setting names only `w > s`,
and only in passing, as the width at which Theorem E′ makes `M = 4i·P_S` the
*only* admissible Gram. **This is F1: the printed proof of (b) has an unstated
hypothesis.** It is repairable without adding one, because a trace suffices.
Put `B := PΠ^⊥Pᵀ ⪰ 0`, so `PPᵀ = A + B`. Taking traces in (H3)
`EEᵀ + w·PPᵀ = N·I_{4s}`, and using `tr(PPᵀ) = 16is` because `P ∈ {±1}^{4s×4i}`:

```
  (s/i)·tr A + w·16is = 4s·N = 16swi + 16s²   ⟹   (s/i)·tr A = 16s²
                                              ⟹   tr A = 16is = tr(PPᵀ) .
```

So `tr B = 0` with `B ⪰ 0`, hence `B = 0`: `P(I₄⊗P_{S^c}) = 0`, the rows of
`P` are blockwise in `V_S`, and `PPᵀ = A`. Then (H3) reads
`(s/i + w)·PPᵀ = N·I`, i.e. `PPᵀ = 4i(wi+s)/(s+wi)·I = 4i·I`. Since the `4s`
rows of `P` are then an orthogonal basis of `ℝ⁴⊗V_S` (dimension `4s`),
`PᵀP = 4i·Π = I₄⊗M`, and `Pᵀ` is an admissible column table. Conversely, given
(a) and (b), `EEᵀ = (s/i)·4i·I = 4s·I` automatically, (H3) becomes the identity
`4s + 4iw = N`, and (H4) holds for the `E` of (d) because
`EQᵀ = −(1/4i)PĈᵀ(I₄⊗M) = −PΠĈᵀ = −PĈᵀ`; so (H3)+(H4) reduce to
`E ∈ {±1}^{4s×4s}`. **No hypothesis on `w` is used anywhere above.** ∎

**2.3 (c).** `Ĉᵀ` is the standard array over `Ḡ` of `(σ₀∘(−), −σ₁, −σ₂, −σ₃)`
at the same `κ(ρ)` (the §1.0 block bookkeeping, unchanged by the port), and
each of those seeds has the same PAF as its original, so Lemma 2 gives
`ĈᵀĈ = ĈĈᵀ = I₄⊗dev(F)`. For the commutation: each developed block is
diagonal in the character basis, hence commutes with `P_S`; and
`R v_χ = χ(κρ)·v_{χ̄}`, so `R` preserves both `V_S` and `V_{S^c}` **because `S`
is conjugation-closed** — which Galois stability supplies. Products of blocks
that commute with `P_S` commute with `P_S`. ∎

**2.4 The row-factorisation lemma.** `W := ĈᵀQ` is `4i × 4s`, and

```
  WWᵀ = Ĉᵀ(QQᵀ)Ĉ = 4i·ĈᵀΠĈ = 4i·ΠĈᵀĈ = 4i·Π(I₄⊗dev F) = 4i·4s·Π .
```

An admissible row `p` is `±1` of length `4i` with every `i`-block in `V_S`, so
`‖p‖² = 4i` and `Πp = p`; hence `Σ_b ⟨p, W[:,b]⟩² = pᵀWWᵀp = 4i·4s·4i =
(4i)²·4s`. By (d) the row of `E` under `p` is `e_b = −(1/4i)⟨p, W[:,b]⟩`, so
`‖e‖² = 4s` for **every** admissible `p`, whatever `Q, κ(ρ), σ` are. Flatness —
`e ∈ {±1}^{4s}`, i.e. `⟨p, W[:,b]⟩ ∈ {±4i}` for every `b` — is therefore a
condition on `p` alone. For the "iff": `4s` mutually orthogonal admissible flat
rows give a `P` with `PPᵀ = 4i·I` and rows blockwise in `V_S`, i.e. (b), whence
(e) supplies (H3) and (H4); conversely a kit's `P` has exactly those rows. ∎
*(The note's aside "a vector of `4s` reals of squared norm `4s` lies in
`{±1}^{4s}` iff every entry is `±1`" is a tautology as written; the content is
the fixed norm, which is what the display before it establishes. Harmless.)*

**2.5 The dual form.** For admissible `p`: `ΠĈp = ĈΠp = Ĉp` by (c), so
`Ĉp ∈ im Π = col(Q)` by (a). The columns `q_b` of `Q` are orthogonal of squared
norm `4i`, so `Ĉp = Σ_b (⟨Ĉp, q_b⟩/4i)·q_b`, and
`⟨Ĉp, q_b⟩ = ⟨p, Ĉᵀq_b⟩ = ⟨p, W[:,b]⟩ = −4i·e_b`, i.e. `Ĉp = −Qe`. Flatness is
`e ∈ {±1}^{4s}`, and `−{±1}^{4s} = {±1}^{4s}`, so `p` flat ⟺
`Ĉp ∈ Q·{±1}^{4s}`. Row by row this is `PĈᵀ = −EQᵀ`, which is (H4); right-
multiplying by `Ĉ` and using `PĈᵀĈ = P(I₄⊗dev F) = 4s·P` (rows of `P` in `V_S`)
gives `4s·P = −EQᵀĈ`, i.e. `P = −(1/4s)·E·Qᵀ·Ĉ`. ∎ Two lines, as claimed.

**2.6 The `(7,12)/S3` kill, re-derived and recomputed.** `S3 = {0,1,2,5,7,10,11}`
is the union of the character-order `1`, `6`, `12` orbits of `Ẑ₁₂`
(`{0} ∪ {2,10} ∪ {1,5,7,11}`), Galois-stable, of size `7`. I computed
`M_{S3}` from `M(c) = 4Σ_{χ∈S3}χ(c)` and the alphabet
`B_{S3} = {v ∈ {±1}¹² : M_{S3}·v = 48·v}` by exhaustive search over all `2¹²`
sign vectors — a **third** implementation, sharing no code with either of the
certificate's two:

```
  M_S3 = [28, 8, 8, -4, -8, 8, -4, 8, -8, -4, 8, 8]     off-zero values {-8,-4,8}
  B_S3 = { +1, -1 }  (the two constants),  rank 1
```

By (a) every `12`-block of every column of `Q` is in `B_{S3}`, so each column is
`(ε₀𝟙, ε₁𝟙, ε₂𝟙, ε₃𝟙)`: at most `16` distinct vectors, spanning dimension `4`,
while `QᵀQ = 48·I₂₈` needs `28` mutually orthogonal ones. Independently, all
twelve rows of a superblock are then *equal*, so their pairwise inner products
are `28`, while (H1) demands `M_{S3}(c) ∈ {−8,−4,8}` for `c ≠ 0`. Either
obstruction alone kills it, at any `σ`, any `w`, any `κ(ρ)` — none of the three
enters — hence at every `N = 4(12w+7)`. *(The note's "equal or antipodal" is a
harmless superset: antipodal cannot occur.)* For the twin: `6 + S3 =
{1,4,5,6,7,8,11}`, **seven** elements, Galois-stable; I confirmed
`M_{6+S3}(c) = (−1)^c·M_{S3}(c)` at every `c`, that its alphabet is exactly the
two alternating vectors, and that `±28` occurs nowhere off zero. The author's
correction to the source decision's six-element transcription is right. ∎

The one thing wrong here is the *sentence*, not the theorem — see F2.

**2.7 Theorem G.** In `ℤ[ℤ_n]` with `J = Σ_k x^k`, `e = x⁰`, write a seed as
`X = J − 2U` for the `{0,1}` indicator `U` of its `−1` positions. Then
`XX* = nJ − 4|U|J + 4UU*`, and the good condition
`Σ_q X_qX_q* = 4n·e` gives the exact integer identity

```
  Σ_q U_qU_q* = n·e + (σ − n)·J ,      σ = Σ_q |U_q| .
```

Reduce mod 2. For a **symmetric** seed, `V* = V` and `VV* = V²`, whose
coefficient at `x^{2k}` is `v_k` (squaring is a ring homomorphism in
characteristic 2, and `2` is invertible since `n` is odd) — write `V^{[2]}`.
For the **skew-type** seed, `u₀ = 0` and `u_{−k} = 1 + u_k`, so
`U* = Σ_k u_{−k}x^k = (J − e) + U + 0 = U + J + e` and
`UU* = U² + |U|J + U`. With `T = U+V+W+Z` and `n ≡ 1 (mod 2)` the identity
collapses to

```
  T^{[2]} + U = e + ε·J ,     ε = 1 + |V| + |W| + |Z|  (mod 2).
```

At `x^{2k}` with `k ≠ 0` (so `2k ≠ 0`, `n` odd): `t_k + u_{2k} = ε`. In signs,
`(−1)^{t_k} = a_kb_kc_kd_k`, and for a symmetric seed
`(−1)^{|V|} = Π_k b_k = b₀` because the nonzero indices pair off; hence
`(−1)^ε = −b₀c₀d₀` and

```
  a_{2k} = −(b₀c₀d₀)·a_k b_k c_k d_k ,   k ≠ 0 .   ∎
```

The prefactor is exactly the `(−1)^{|V|+|W|+|Z|}` the reduction produces; the
source consult's normalised form is the `b₀ = c₀ = d₀ = +1` special case, and
the note says so. The Williamson variant falls out of the same identity with
`U* = U` throughout: `T^{[2]} = e + (σ+1)J`, so `t_k` is *constant* off zero and
`a_kb_kc_kd_k = −(a₀b₀c₀d₀)` — a different conclusion, which is exactly the
control the certificate runs. Skewness of `A` is load-bearing, and both the
note and the certificate say where.

**2.8 The closure sign law.** `a_{2^{j+1}} = ε·a_{2^j}·π_j` with
`ε = −(b₀c₀d₀)` and `π_j = b_{2^j}c_{2^j}d_{2^j}`, so
`a_{2^L} = ε^L·(Π_{j<L} π_j)·a₁`. If the seeds are constant on the orbits of
`⟨M₀,−1⟩` then `a_{2^L} = a₁` when `2^L ∈ M₀` and `a_{2^L} = a_{−m} = −a_m =
−a₁` when `2^L = −m ∈ −M₀`; dividing by `a₁ = ±1` and using `ε² = 1` gives
`Π_j π_j = closure·ε^L`. One parity condition, whatever `L` is. ∎

**2.9 The two arithmetics, recomputed.** `n = 7`, `M₀ = {1}`: `⟨M₀,−1⟩ = {1,6}`
has `3` orbits on `ℤ₇∖{0}`, doubling is a single `3`-cycle on them, `2³ = 1 ∈ M₀`,
`closure = +1`, and under the normalisation `Π_j π_j = (−1)³ = −1`. `n = 523`:
prime, `|ℤ*₅₂₃| = 522 = 2·3²·29`, the unique order-3 subgroup is
`M₀ = {1, 60, 462}`, `|⟨M₀,−1⟩| = 6`, `522/6 = 87` orbits all of size `6`, the
least `L` with `2^L ∈ ⟨M₀,−1⟩` is **87**, `2⁸⁷ = 463`, `−463 = 60 ∈ M₀` so
`2⁸⁷ ∈ −M₀`, `closure = −1`, and `Π_j π_j = (−1)·(−1)⁸⁷ = +1`. All of it
reproduces. The note's "**one** parity condition, not 87" is the right reading.

**2.10 The censuses, re-enumerated.** By brute force over all skew and all
symmetric seeds with `b₀` free — no meet-in-the-middle, no prefilter, sharing
no code with cert 32:

| `n` | good quadruples | prefixed form holds | unprefixed `−` | unprefixed `+` |
| --- | --- | --- | --- | --- |
| 7 | **528** | 528 | 264 | 264 |
| 9 | **288** | 288 | 144 | 144 |
| 11 | **1 440** | 1440 | 720 | 720 |

and, for the Williamson control (four symmetric seeds):

| `n` | Williamson quadruples | constant-product | prefixed doubling | unprefixed `−` doubling |
| --- | --- | --- | --- | --- |
| 7 | **960** | 960 | **0** | **48** |
| 9 | **2 112** | 2112 | 0 | 0 |
| 11 | **1 920** | 1920 | 0 | 0 |

The author's correction to the source decision is **confirmed independently**:
the theorem's own prefixed form holds in `0` of the 960 at `n = 7`, and the
"48 of 960" belongs to the *unprefixed* variant (48 to each sign). The note and
PROVENANCE record the difference; the control is strengthened by it, not
weakened. The closure prediction at `n = 7` holds in **528 / 528** here.

**2.11 The support enumerations of §2.6.** Independently reproduced:
Galois-stable size-11 subsets of `Ẑ₁₆` containing `χ₀` — **exactly one**,
`{0,1,3,4,5,7,9,11,12,13,15}`, which is the `S` of the banked table
`A_11_16_32`; Galois-stable size-7 subsets of `Ẑ₁₂` containing `χ₀` — **exactly
four**, and they are the four banked (`B1…B4`), with `B2` carrying
`S2 = {0,1,4,5,7,8,11}` and `B3` carrying `S3` as §1.12 names them; all **56**
size-3 subsets of `𝔽̂₂³` and all **12** Galois-stable size-3 subsets of
`Ẑ₄×Ẑ₂` are real-character twists of `W∖{1}` for an order-4 subgroup `W`. The
alphabets and ranks at the four `(7,12)` supports: `66, 32, 2, 24` of rank
`7, 7, 1, 7`, in the note's stated order. Every banked table satisfies
`Σ_q r_q² = 4s` — which is `F̂(χ₀)` when `χ₀ ∈ S`, as it is in all six — and
`N = 2092 = 4(32·16+11) = 4(43·12+7) = 4(65·8+3)`.

---

## 3. Replays

From the worktree root, `PYTHONDONTWRITEBYTECODE=1`, one process at a time.

| cert | command | checks | FAIL | exit | wall |
| --- | --- | --- | --- | --- | --- |
| 30 | `python certs/30-general-branch-sigma/run.py` | **59** | 0 | **0** | 0.16 s |
| 31 | `python certs/31-s2-column-table/run.py` | **50** | 0 | **0** | 1.15 s |
| 32 | `python certs/32-good-product-theorem/run.py` | **73** | 0 | **0** | 1.25 s |

Matches the author's report exactly. No `--full` exists in any of the three;
every claim each makes is on its default path.

**Tamper tests** (scratch copy of the worktree, one banked number each):

| cert | tamper | result |
| --- | --- | --- |
| 30 | `data/general-branch-sigma-tables.json`, one `σ` entry `4 → 2` | **exit 1**, 3 FAIL: the file digest, the aggregate-PAF/parity/box/norm path *independently*, and the "unperturbed table still passes" guard |
| 31 | `data/q-7_12-S2.json`, one `Q` block entry `1 → −1` | **exit 1**, 8 FAIL: digest, alphabet membership, `QᵀQ = 48I`, `(H1)`, the forced multiset, `WWᵀ`, and both negative controls |
| 32 | `EXPECTED_GOOD[7] 528 → 529` in `run.py` | **exit 1**, 1 FAIL: the enumeration reports 528 against the banked 529 |

All three **fail closed**, and in cert 30's and 31's case the mathematics
catches the tamper even without the digest — the digest is not the only guard.

**Self-containment.** `certs/30,31,32/run.py` import only
`argparse copy hashlib itertools json math os random shutil subprocess sys
time fractions`; every path is built from `__file__` and stays inside the
repository (`data/`, `verify/verify.py`, a scratch `out/` under cert 32 that is
deleted). Nothing is imported from the private lab, and cert 31 carries its own
`648`-representative kit engine rather than importing cert 17's. `git status`
is clean after all three runs and after cert 32's five `verify/verify.py`
assemblies at orders 28, 36, 44, 52, 60.

**Untouched.** `git diff main --stat`: **12 files, 5 533 insertions, 0
deletions.** Certs 01–25, `verify/`, `tools/`, `data/`'s pre-existing files,
`DISCLOSURE.md`, `CITATION.cff`, `LICENSE`, `PLAN.md` and `DECISIONS.md`
anywhere: untouched. The only pre-existing files that changed are
`note/NOTE-B.md` (three pure insertions: §§1.12–1.13 after §1.8, the §2.4
forward pointer plus §2.6, five §5 rows), `README.md` (three replay lines plus
a paragraph, both at the end of their blocks) and `PROVENANCE.md` (a dated
paragraph at the end). Every hunk read.

**Names.** `grep -riE "jones|jd-|jdj|\bJD\b"` over the whole worktree returns
only the pre-existing GitHub URLs in `README.md` and `PROVENANCE.md`; **no
added line contains a person's name**, and no mathematical text anywhere does.
Crediting in the new PROVENANCE paragraph is by station throughout — "the Sol
station, an external reviewer", "a desk lane of the source laboratory", "the
desk audit of that laboratory's Cursor lane `PR-0051`", "an unaudited lane".
The `data/` provenance strings cite lab paths, which is the established house
practice on `main`. The two banked digests in PROVENANCE's table match
`sha256sum` on the files.

**DISCLOSURE.** Untouched and still accurate: the "thirteen `--full` flags"
sentence names certs 06/08/11/13/14/15/19/20/21/22/23/24/25 and certs 30–32
have none; "Every computational claim in the note carries a certificate"
survives *because* cert 32 exists.

---

## 4. Findings

**F1 — §1.12(b): an unstated hypothesis in the proof.** *(substantive, repair
supplied)* The proof of clauses (b)/(e) says "the computations of §1.8
verbatim". §1.8's forward derivation of `PPᵀ = 4i·I` goes through **D3**
(§1.3), whose hypothesis is `w > 2s`. §1.12's setting names no such condition —
it names `w > s`, and only as the width at which Theorem E′ forces the Gram.
As written, §1.12 therefore either has a silent `w > 2s` or a gap. No claim in
the port is affected (every cell at 2092 has `w > 2s`, and §2.6 says so), but
the theorem is stated in the general branch and should be true there. The trace
argument of §2.2 above closes it and needs no `w` hypothesis at all — which is
strictly better than adding one.

**F2 — §1.12, the kill sentence is refuted as worded.** *(substantive)* Line
1155: "Then **no `Q ∈ {±1}^{48×28}` satisfies (H1)**". (H1) as §1.1 states it
quantifies `M` existentially ("there is a symmetric, `Ḡ`-invariant `M` …"), so
the plain reading is that the whole `(7,12)` cell has no column table — and
that is contradicted by this same section's scope paragraph ("An admissible
column table exists at `(7,12)/S2`") and by cert 31 [C]. The intended,
tier-scoped statement is true and is what §5's row already says. One clause
fixes it.

**F3 — two labels outside the note's declared vocabulary.** *(editorial)*
NOTE-B's preamble names six honesty labels. `DEAD-BY-CERTIFICATE` (§2.6, and
in cert 31's `NOTES.md`/`run.py` output) and `REPORTED-UNVERIFIED` (§1.13, §5)
are new label-shaped tokens with no precedent on `main`. A reader who takes the
preamble at its word cannot resolve them.

**F4 — the §1.13 prior-art pointer dangles.** *(publication-relevant)* §1.13
and §5 both send the reader to §4 for arXiv:1811.05094 and the good-matrices
product relationship. §4 contains neither, and §4 was deliberately not touched.
§4's own opening rule is that a source not read firsthand is declared *there*.
As it stands, the only unread-source declaration in the note lives outside the
prior-art section, in a theorem section, behind a pointer that does not resolve.

**F5 — section-number gaps.** *(merge-conditional)* Movement I runs
§1.8 → §1.12; Movement II runs §2.4 → §2.6. The gaps are lane a's reserved
numbers. If both lanes merge, nothing to do; if this branch ships alone, the
note publishes with three phantom sections in Movement I and one in Movement II.
Related and already flagged by the author: NOTE-B's opening Movement-I bullet
and README's "The three movements" do not mention the general branch or the
product theorem, and DISCLOSURE's opening summary does not either.

**F6 — one ambiguous antecedent.** *(editorial)* §1.12's MEASURED paragraph:
"…checked at `(2,4)` on all four banked records of **cert 18** … and at `(3,4)`
on kits found by **that certificate's** own engine". The `(3,4)` kits are cert
**31**'s, not cert 18's. §2.6's parallel sentence ("its own engine") is fine.

Two things I looked for and did **not** find: any strengthening of the source
laboratory's claims (the 40 494-point flat-row census is reported as a bounded
negative search and explicitly not certified; the `Q` quantifier is left open at
S1/S2/S4; no seed statement anywhere), and any statement about `H(2092)` or
about a cell (every scope paragraph disclaims both, in the note, in the README
and in all three certificate epilogues). On the lane's own desk question: keep
cert 32.

---

## 5. Repairs

| file | old → new | kind |
| --- | --- | --- |
| `note/NOTE-B.md` §1.12, proof of Theorem F (line ~1093) | `(d),(e),(b),(f),(g): the computations of §1.8 verbatim,` → `(d),(f),(g): the computations of §1.8 verbatim. (b) and (e): §1.8 reaches PPᵀ = 4i·I through D3, which assumes w > 2s; no such hypothesis is available here, and none is needed. Write A = PΠPᵀ and B = PΠ^⊥Pᵀ ⪰ 0 with Π = I₄⊗P_S. From (d), EEᵀ = (s/i)·A. Taking traces in (H3) and using tr(PPᵀ) = 16is gives (s/i)·tr A = 16s², so tr A = 16is = tr(PPᵀ), hence B = 0 and PΠ^⊥ = 0; then (H3) reads (s/i + w)·PPᵀ = N·I, i.e. PPᵀ = 4i·I, and the 4s rows of P are an orthogonal basis of ℝ⁴⊗V_S, so PᵀP = I₄⊗M. Conversely (a)+(b) give EEᵀ = 4s·I, (H3) becomes 4s + 4iw = N, and EQᵀ = −PΠĈᵀ = −PĈᵀ is (H4). The one substantive step is…` | text (a proof gap) |
| `note/NOTE-B.md` §1.12 line 1155 | `Then **no \`Q ∈ {±1}^{48×28}\` satisfies (H1)**:` → `Then **no \`Q ∈ {±1}^{48×28}\` satisfies (H1) with \`M = M_{S3}\`** — the S3 tier admits no column table:` | text (scope) |
| `note/NOTE-B.md` §2.6 line 1540 | `- **\`(7,12)/S3\` is DEAD-BY-CERTIFICATE at every order** (§1.12):` → `- **\`(7,12)/S3\` admits no column table at any order** — **PROVEN-BY-CERTIFICATE** (§1.12):` | label |
| `note/NOTE-B.md` §1.13 line ~1293 and §5's Theorem G row | `**REPORTED-UNVERIFIED** here` → `read only in the abstract here, and so **NOT CLAIMED** either way` — *and* add the source to §4 (one line: arXiv:1811.05094, good matrices, abstract only, no comparison made), so that the two `(§4)` pointers resolve | label + text |
| `note/NOTE-B.md` §1.12 MEASURED paragraph (line ~1138) | `on kits found by that certificate's own engine` → `on kits found by cert 31's own engine` | text |
| `certs/31-s2-column-table/{NOTES.md,run.py}` | `DEAD-BY-CERTIFICATE` (4 sites) → `dead at the column layer, PROVEN-BY-CERTIFICATE` | label |
| `note/NOTE-B.md` lines 3–34; `README.md` "The three movements"; `DISCLOSURE.md` opening | one clause each for the general branch and the product theorem | text (at merge) |
| `note/NOTE-B.md` §§1.12/1.13/2.6 | renumber **only if lane a does not merge** | label (merge-conditional) |

---

## 6. Recommendation

**REPAIR-FIRST.** The port is faithful, the proofs are right, the certificates
are real and fail closed, the scope discipline is the best in the note, and the
two corrections the lane makes to the source laboratory's own decision text are
both confirmed here by independent enumeration. But two of the repairs are not
cosmetic — F1 leaves a hypothesis unstated in a theorem the section advertises
as hypothesis-free, and F2's sentence is contradicted by its own section three
paragraphs later. Apply R1, R2 and R4 (with the §4 line), then push; F3, F5 and
F6 can ride the same commit or the merge.

*Audit lane `hb-audit-b`, on `desk/ports-turn54-b`. Independent re-derivations
and re-enumerations were written from `note/NOTE-B.md` §§1.0–1.1 and run
outside the repository; nothing in this audit was taken from the author lane's
report except the list of what to check.*
