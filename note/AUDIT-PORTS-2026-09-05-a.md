# Audit — the turn-54 port, lane a (skeptic)

Branch `desk/ports-turn54-a`, worktree `wt/hb-ports-a`, author commits
`7007603` (certs 26–29) and `0e72a58` (NOTE-B, README, PROVENANCE, port note).
Mandate: **refute it**. Nothing entering `note/NOTE-B.md` is taken on the
author's word; every theorem and lemma below was re-derived here from the
definitions of §1.0–§1.1, every certificate was replayed from the worktree
root, and every banked number was tampered with in a scratch copy to confirm
it fails closed.

**Result: nothing is refuted.** The four ported statements are sound as
mathematics and the four certificates pass and fail closed. Seven repairs are
recorded, all textual or symbolic; one of them (**F1**) is a wrong symbol in a
displayed inequality, one (**F2**) is a stale source-laboratory certificate
number printed by a public certificate, and the rest are scope qualifiers and
staleness. Recommendation: **REPAIR-FIRST**.

---

## 1. Verdict table

| item | verdict | note |
| --- | --- | --- |
| **Theorem T** (§1.9), the identity `Hᵀ = D̃·BGS(x₀∘(−), −x₁, x₂, x₃; ρ; Eᵀ, (S′Q)ᵀ, (PS′)ᵀ)·D̃` | **SOUND** | re-derived below (§2.1); reproduced **entrywise** by an independent assembler written here from §1.0 alone, on six arbitrary `(n,i,w,s,ρ)` shapes with random seeds, corners and tables, none Hadamard |
| **Theorem T's branch-freeness** — via Theorem A's *only if* half, not via `EEᵀ = 4s·I` | **SOUND** | the argument is exactly right and is the load-bearing repair D-068.8a asked for; no house-branch hypothesis (`i = s+1`, `w > 2s`, `n ≥ 3`) enters at any step |
| **the `(H″)ᵀ` identity** (§1.9), exact, no conjugation | **SOUND** | re-derived (§2.1, F1′); reproduced entrywise by the independent assembler on the same six shapes |
| **Remark R** (§1.9) — `H″` Hadamard iff `P` annihilates `Ĉ`'s off-diagonal part | **SOUND-WITH-REPAIR** | the derivation is correct and is in fact an *iff*, which the note conservatively states as "only when". **F4**: "all seven records" is the `--full` path; the default path exercises three |
| **`H ~ Hᵀ` ⟺ fixed point of `T`** (§1.9), PROVEN | **SOUND** | a restatement of Theorem T, correctly hedged ("an involution up to the family's gauge") |
| **the sporadicity conjecture** (§1.9), CONJECTURE | **SOUND** | labelled, and claimed nowhere else |
| **the small-order census** (§1.9), COMPUTATIONAL-EVIDENCE | **SOUND-WITH-REPAIR** | agreement-is-not-equivalence is stated three times and the four undecided instances are reported as undecided. **F5**: "eleven instances … seven proved" is the `--full` count; the default run reports **9 instances, 5 proved, 4 undecided** |
| **the source laboratory's 44-instance census**, cited and not replayed | **SOUND** | correctly quarantined: finder-side, negatives resting on a refinement's completeness, "nothing in this note rests on it" |
| **Sylow-2 rigidity** (§1.10), PROVEN | **SOUND** | every case re-derived below (§2.2); the `t = 523` conclusion independently brute-forced here |
| **the 12 groups / 24 `(E, e*)` pairs** (§1.10) | **SOUND** | orbit counts `1,2,1,2,1` re-derived by hand from `Aut` of the five groups of order 8; the twelve summands `1+1+3+1+3+7+3+1+1+1+1+1 = 24` reproduce term by term in cert 27's own census |
| **`t` an odd prime `> 7`** (§1.10) | **SOUND** | the hypothesis is exactly what `n_t ∈ {1,2,4,8}`, `n_t ≡ 1 (mod t)` needs, and it is used only there |
| **DLFH 2000 CITED, not re-derived** (§1.10) | **SOUND** | said three times, in the section, the label table and the certificate's verdict block |
| **the `t = 3` RDS census** (§1.10), a control | **SOUND** | explicitly outside the theorem's range and explicitly not offered as an instance; its own scope ("all 12 groups of order 24 **with a normal `ℤ₃`**") is stated, so `S₄` is not silently swept in |
| **Lemma F (i)–(iv)** (§1.11), PROVEN | **SOUND** | re-derived below (§2.3); the trace/eigenvalue argument, the quadratic `p f² − N(p+1)f + N² ≥ 0` and its roots `N/p`, `N` are all correct, and both branches (`f > m` and `f ≤ m`) close |
| **the "Admissible primes" paragraph** (§1.11) | **SOUND-WITH-REPAIR** | the prime arithmetic is right (independently reproduced). **F1: `\|(AAᵀ)_{rr′}\| ≤ 4 < m` is the wrong symbol — it must be `4 < p`**, and the conclusion then needs `p ≥ 5` |
| **22 admissible / 293 excluded primes at 2092; 131, 349 excluded; `523 ⟹ f = 0, m = 4`** | **SOUND** | reproduced here by an independent sieve: 315 odd primes `≤ 2092`, 22 admissible, 293 excluded, 131 and 349 among them, `2092 = 4·523` |
| **the 16-block circulant shape, "not the Williamson array"** (§1.11) | **SOUND** | the distinction is real and correctly drawn; the GS array's exclusion (back-circulants; `2a = 0` impossible for odd `\|G\|`) is right |
| **the records theorem** (§1.11), PROVEN + PROVEN-BY-CERTIFICATE | **SOUND** | the implication `Σ_C (\|C\| mod p) > N/p ⟹ no order-`p`` is re-derived below and is airtight; the label is earned — see §4 on the missing ADR and §3 on the audit surface |
| **`Aut` 2-part NOT COMPUTED**, `H′(716)`/`H″(716)` never run | **SOUND** | stated in §1.11, in §3.4's and §3.6's inserted paragraphs, in the label table and in the certificate's verdict |
| **the `(2,4)` existence theorem** (§2.5), PROVEN + PROVEN-BY-CERTIFICATE | **SOUND** | all six steps re-derived below (§2.4); every combinatorial input independently reproduced here |
| **step (3)'s "take `Q′` signed-Sylvester"** (§2.5) | **SOUND-WITH-REPAIR** | **F3**: column signed permutations alone do *not* reduce `H(8)` to Sylvester — there are 480 right-orbits, as §2.4 itself says. Nothing is lost, because the conclusion is stated for `A(Q′)` and step (6) is uniform over all 105 matchings; the sentence needs one clause |
| **"only the rank-1 case is used"**, rank-2/3 surplus (§2.5) | **SOUND** | and it is the right way round: the proof is uniform, the certificate settles the surplus, and no case distinction is load-bearing |
| **the compound label PROVEN + PROVEN-BY-CERTIFICATE** at §2.5 | **SOUND** | the two finite inputs the paper proof does not re-derive (the 448-type table; "`Π` is always a translation") are named explicitly, which is what makes the compound honest |
| **§2.5 as its own section rather than §2.4's second half** | **SOUND** | §2.4 already carries a three-part proposition, Theorem 3 and the general-branch remark; the split is the right call and both directions forward-point |
| **the §3.4 repair** ("source-laboratory intel" → Theorem T) | **SOUND** | accurate: Theorem T is proved here and cert 26's **default** path checks it entrywise on that very record. Non-strengthening: "Nothing in this section rests on it either way" |
| **the §3-preamble, §1.4, §2.3, §2.4(c), §3.4, §3.6 forward pointers** | **SOUND** | every one is a pointer or a non-strengthening remark; §3's preamble says outright "none of the separations uses it, and it strengthens none of them" |
| **§5 label table, nine new rows** | **SOUND** | one row per new statement, labels matching the section headers, `--full`/not-run boundaries carried into the table |
| **README's two movement paragraphs, abstract, replay block, `--full` paragraph** | **SOUND-WITH-REPAIR** | **F6**: "Lemma F … forces a 16-block circulant array at `p = 523`" reads unconditional; it needs "for a hypothetical `H(2092)` admitting an automorphism of order 523" |
| **PROVENANCE's new dated section** | **SOUND** | stations only — "a Cursor cloud lane", "the desk", "the source laboratory"; no person anywhere; the four file digests in its table all reproduce (`sha256sum`) |
| **DISCLOSURE.md untouched** | **SOUND-WITH-REPAIR** | the *stance* is genuinely unaffected. **F7**: two staleness spots — the dated results summary now stops at 2026-09-03, and the "audit vs replay" sentence's cert list omits 29 |
| **cert 26** replay + tamper | **SOUND** | 46 checks / 16.1 s default, 64 checks / 44.8 s `--full`, exit 0 both; tampering one T-image pin fails closed (exit 1) |
| **cert 27** replay + tamper | **SOUND** | 14 checks / 13.1 s default, 14 checks / 13.7 s `--full`, exit 0 both; tampering the 8384 count fails closed |
| **cert 28** replay + tamper | **SOUND-WITH-REPAIR** | 14 checks / 0.77 s, exit 0; tampering the 448-type count fails closed. **F2**: the run prints and comments the *laboratory's* `0022` where this repository's cell-24 certificate is **18** |
| **cert 29** replay + two tampers | **SOUND** | 40 checks / 1.28 s, exit 0; tampering a class size and (separately) a matrix digest each fail closed, the second proving the digest is recomputed in-process |
| **self-containment of all four certs** | **SOUND** | standard library only; reads confined to `verify/`, `data/` and the certificate's own directory; generated matrices written to `tempfile.mkdtemp` and removed; the worktree is clean after every replay |
| **no pre-existing cert or theorem text modified** | **SOUND** | `git diff main --stat` touches only `NOTE-B.md`, `README.md`, `PROVENANCE.md`, the four new certificate directories and the port note; every hunk in `NOTE-B.md`'s pre-existing text is a forward pointer, the §3.4 repair, or an abstract line |
| **owner-name rule** | **SOUND** | `git grep -i` for the surname across the worktree returns only `LICENSE`'s copyright line, `CITATION.cff`'s author field and the repository URLs — nothing in any mathematical text |

---

## 2. My re-derivations

These are written from `note/NOTE-B.md` §1.0–§1.1 and from the definitions in
the sections themselves. They are not paraphrases of the author's proofs.

### 2.1 Theorem T

Index the core's blocks by `(I,J) ∈ {0,1,2,3}²` as in §1.0.

*(F1) — the block facts.* `(XR)[g,h] = x(ρ−g−h)` and `(XᵀR)[g,h] = x(g+h−ρ)`
depend on `g+h` only, so **every off-diagonal block is a symmetric matrix**.
Reading the array off §1.0, `block(1,0) = −BR = −block(0,1)`,
`block(2,0) = −CR = −block(0,2)`, `block(3,0) = −DR = −block(0,3)`,
`block(2,1) = −DᵀR = −block(1,2)`, `block(3,1) = CᵀR = −block(1,3)`,
`block(3,2) = −BᵀR = −block(2,3)`: all six above-diagonal positions.
Hence `(Cᵀ)`'s `(I,J)` block is `block(J,I)ᵀ = block(J,I) = −block(I,J)` for
`I ≠ J`, and `Aᵀ` on the diagonal, i.e.

`Cᵀ = C^{sw}(x₀∘(−), x₁, x₂, x₃; ρ)`,

using `dev(x)ᵀ = dev(x∘(−))` — only the type-1 seed reverses, because the
off-diagonal blocks are already symmetric. ✔ as stated.

*(F2) — the switch is a seed negation up to conjugation.* Conjugating by
`S₀₁ = diag(−1,−1,+1,+1)⊗I_n` multiplies block `(I,J)` by `ε_Iε_J` with
`ε = (−1,−1,+1,+1)`; that is `−1` exactly when one index is in `{0,1}` and the
other is not — the **eight** blocks `(0,2),(0,3),(1,2),(1,3)` and their
transposes. Negating `x₁` sends `B ↦ −B`, flipping precisely the **four**
remaining off-diagonal blocks `(0,1),(1,0),(2,3),(3,2)`. Together: all twelve,
diagonal untouched. So `C^{sw}(y) = S₀₁·C(y₀,−y₁,y₂,y₃)·S₀₁`. ✔

*(F3) — the border.* `H = [E, P̃ ; Q̃, C]` gives `Hᵀ = [Eᵀ, Q̃ᵀ ; P̃ᵀ, Cᵀ]`, and
with `D̃ = diag(I_{4s}, S₀₁)`,
`D̃HᵀD̃ = [Eᵀ, Q̃ᵀS₀₁ ; S₀₁P̃ᵀ, S₀₁CᵀS₀₁]`. By (F1) and (F2),
`S₀₁CᵀS₀₁ = C(x₀∘(−), −x₁, x₂, x₃; ρ)` (same `ρ`, since `R` is untouched).
The new top-right strip has entry `ε_I·Q[iI+κ(g)][r]` at `(r,(I,g))`, which is
`((S′Q)ᵀ)[r][iI+κ(g)]` — a **coset-constant** row table, because
`S₀₁ = S′⊗I_w`. The new bottom-left strip is `((PS′)ᵀ)[iI+κ(g)][c]`,
likewise coset-constant. ∎

*Branch-freeness.* `D̃` is a `±1` diagonal, so `D̃HᵀD̃` is Hadamard exactly when
`H` is; the displayed data is a standard-orientation instance of the §1.0
ansatz; Theorem A's *only if* half returns (H1)–(H4) for it. **Nothing but
`HHᵀ = N·I` is used.** I checked specifically that no step touches `w > 2s`,
`i = s+1`, `n ≥ 3` or `EEᵀ = 4s·I`. The author's claim is exact.

*(F1′) — the `(H″)ᵀ` identity.* `C^{sw} = I₄⊗A − (C − I₄⊗A)`, so
`(C^{sw})ᵀ = I₄⊗Aᵀ + (C − I₄⊗A)`: diagonal reversed, off-diagonal **unchanged**
— that is `C(x₀∘(−), x₁, x₂, x₃; ρ)`. With the border transposed as in (F3)
and no `S₀₁` anywhere, `(H″)ᵀ = BGS(x₀∘(−), x₁, x₂, x₃; ρ; Eᵀ, Qᵀ, Pᵀ)`
exactly. ✔

*Remark R.* `H″` shares `E`, `P`, `Q` and the seeds' PAFs with `H`, so
(H1),(H2),(H3) are untouched and only (H4) can fail. Coset-summing commutes
with negating blocks, so `H″`'s compressed core is `Ĉ^{sw}`, and subtracting
the two (H4)s gives `P·(Ĉ − Ĉ^{sw})ᵀ = 2·P·(offdiag Ĉ)ᵀ = 0`. Given that `H`
is Hadamard this is an **iff**, not merely "only when"; the note's weaker
phrasing is safe. At `(1,1)` the quotient is trivial, `Ĉ` is the `4×4` array of
the row sums `(±2,0,0,0)`, hence `±2I₄` with no off-diagonal part — automatic,
as claimed, and that covers every `H″` §3 names.

*Independent entrywise check.* I wrote my own `dev`/`R`/core/border assembler
from §1.0 and evaluated both identities on six random shapes
`(n,i,w,s,ρ) = (6,2,3,2,5), (6,3,2,1,5), (5,1,5,2,0), (8,2,4,1,7), (8,8,1,3,3),
(6,3,2,3,0)` with random `x`, `E`, `P`, `Q` — none Hadamard, none satisfying
(H1)–(H4). **Both identities held exactly in all six**, confirming the author's
"table-free and seed-free" claim independently of cert 26.

### 2.2 Sylow-2 rigidity

`n_t | 8` and `n_t ≡ 1 (mod t)` force `n_t = 1` once `t > 7`, so
`E = ℤ_t ⋊ P`; at `t = 523`, `gcd(8, 522) = 2`.

With `f` the projection of `g = D − De*` along `ℤ_t` and `a,b,c,d` its values
on coset representatives of `⟨e*⟩`, `Σf(q)² = 8t` becomes
`a²+b²+c²+d² = 4t` with all four **odd**. I re-derived the two obstructions
the table asserts as impossible for **every** odd `t`:

* `ℤ₈`: `ab+bc+cd−da = b(a+c) + d(c−a)`. Write `a+c = 2u`, `c−a = 2v`; then
  `u+v = c` is odd, so `u,v` have opposite parity, and with `b,d` odd
  `bu+dv ≡ u+v ≡ 1 (mod 2)`. So the expression is `2·(odd) ≡ 2 (mod 4)` and can
  never vanish. ✔ (the note's "`≡ 2 mod 4`")
* `D₄`: `ac = bd` and `ad = −bc` multiply to `a²cd = −b²cd`, i.e.
  `cd(a²+b²) = 0`; `a,b` odd makes `a²+b² > 0`, so `cd = 0`, impossible for odd
  `c,d`. ✔
* `ℤ₂³`: the three conditions give `a² = b² = c² = d²`, so `4t = 4a²` and `t`
  is a square. ✔
* `ℤ₄×ℤ₂`, both `e*` cases: I did **not** take the note's word. I brute-forced
  all odd `(a,b,c,d)` with `a²+b²+c²+d² = 2092` against `ac+bd = 0` (the square
  case) and against `(a₀+a₂)(a₁+a₃) = 0`, `a₀a₂+a₁a₃ = 0` (the non-square
  case): **0 solutions each**. And `523` is not a square. So the `t = 523`
  conclusion does not depend on the "sum of two squares" step at all.
* `Q₈`: only `a²+b²+c²+d² = 4t`, and the ordered odd four-square
  representations of 2092 number **8384**, which I recomputed independently.

The 12 groups: `Aut`-orbits of index-2 subgroups are `1, 2, 1, 2, 1` for
`ℤ₈, ℤ₄×ℤ₂, ℤ₂³, D₄, Q₈` — I re-derived each (the two `ℤ₄`s in `ℤ₄×ℤ₂` are
swapped by `Aut`; `Aut(Q₈) ≅ S₄` is transitive on its three `ℤ₄`s) — giving
`5 + 7 = 12`. A central involution of `E` is a central involution of `P` lying
in the kernel of `P → Aut(ℤ_t)`, and cert 27's per-group tally prints exactly
the note's twelve summands. ✔

### 2.3 Lemma F

(i) On a `p`-cycle the product of the signs is `+1` because `P^p = I`, so a
diagonal unsigns it; at a fixed point `e^p = e` and `e^p = 1` give `e = 1` for
odd `p`. (ii) `PHQᵀ = H` gives `PH = HQ`, so `P = HQH⁻¹ = HQHᵀ/N` and
`tr P = tr Q`; with the signs removed those traces count fixed points.

(iii) Rows split into `f` fixed and `m = (N−f)/p` orbits, likewise columns. By
(iv), `AAᵀ + p·BBᵀ = N·I_f` with `A` the `f×f` fixed/fixed block and `B` the
`f×m` fixed-row block (constant on column orbits). If `f > m` then
`rank(p·BBᵀ) ≤ m < f`, so `N·I_f − AAᵀ` is singular and `AAᵀ` carries the
eigenvalue `N` with multiplicity `≥ f−m`; every eigenvalue of `AAᵀ` is `≥ 0`
and `tr(AAᵀ) = f²` (an `f×f` sign matrix), so

`f² ≥ N(f−m) = N·((p+1)f − N)/p`, i.e. `p f² − N(p+1)f + N² ≥ 0`,

whose roots are `N/p` and `N`; `f = N` forces `P = Q = I`. If instead `f ≤ m`
then `pf ≤ N−f`, i.e. `f ≤ N/(p+1) ≤ N/p`. Either way `f ≤ N/p` for a
nontrivial automorphism. ∎ ✔ exactly as printed.

At `N = 2092`, `p = 523`: `N − f = pm` gives `f ≡ N ≡ 0 (mod 523)`, and
`f ≤ N/p = 4`, so `f = 0`, `m = 4` — four row orbits and four column orbits,
each block a `523×523` circulant. ✔

**The records implication.** `π(i)` is preserved by signed row/column
permutations, so an order-`p` automorphism permutes each `π`-class; within a
class the action is `p`-cycles plus fixed points, so the class contributes at
least `|C| mod p` fixed rows. Summing and comparing with (iii):
`Σ_C (|C| mod p) > N/p ⟹` no automorphism of order `p`. With classes of size 1
and 2 only and `p ≥ 3`, the sum is `4·1 + 332·2 = 668 = N > N/p`. ✔ At 716,
`4 + 356·2 = 716 = N`. ✔ Since a group with an element of odd order has one of
odd prime order, `Aut` is a 2-group. ✔ The four singletons are fixed
individually and each `τ`-pair is preserved setwise, which is the last clause.

### 2.4 The `(2,4)` existence theorem

*(1)* With `P` and `Q` anti-periodic, `(PĈᵀ)[r,(I,c)]` collapses to
`(P′C*ᵀ)[r,(I,c)]` under the displayed definition of `C*`, and `EQᵀ` restricted
to `c ∈ {0,1}` is `EQ′ᵀ`; the `c ∈ {2,3}` components are the negatives, so
(H4) ⟺ `EQ′ᵀ + P′C*ᵀ = 0`. ✔

*(2)* `W := (1/4)P′C*ᵀ` has entries `(1/4)(±2 ± 2) ∈ {0, ±1}` given two `±2`
per row of `C*`, and `WWᵀ = (1/16)P′(C*ᵀC*)P′ᵀ = (1/16)·8·8·I = 4I` — using
that `C*C*ᵀ = 8I` and `C*ᵀC* = 8I` are equivalent for a square matrix. So `W`
is a `W(8,4)`; `E = −(1/2)WQ′` and `P′ = (1/2)WC*` invert the construction. ✔

*(3)* `wQ′ ∈ {±2}⁸ ⟺ w = (1/4)vQ′ᵀ` for `v ∈ {±1}⁸`, so admissibility is
"the transform of a sign function has a 4-point support of values `±4`". I
exhausted all 256 sign functions on `𝔽₂³` myself: **112** qualify, their
supports are **exactly the 14 affine planes**, and **all 112** have an odd
number of minus signs in the spectrum. `112 = 14 × 8` ✔.

*(4)* `C*`'s column supports are pairs; `(1/2)(±W[r,k₁] ± W[r,k₂]) ∈ {±1}` iff
exactly one of the two is nonzero, i.e. `supp(w)` meets each pair once. ✔

*(6)* For a perfect matching `Π` of `𝔽₂³` with differences `D = {a_i+b_i}`, the
four differences sum to the sum of all eight points, which is `0`; no `d ∈ D`
is `0`; so no **odd** subset of `D` sums to `0` (a 3-subset summing to zero
would force the fourth difference to vanish). A linear `ℓ` with `ℓ|_D ≡ 1`
therefore exists — the consistency condition for that inhomogeneous system is
exactly "every zero-summing subset has even size" — and the solution set is an
affine space of dimension `3 − rk D`. A plane `ℓ = ε` meets `{a,b}` once iff
`ℓ(a+b) = 1`, so **both** planes of the class are transversal. A plane and its
complement, each with its four sign origins, give eight rows: within a plane,
`Σ (1−2e_u)(1−2e_{u′}) = 4−2−2 = 0`; across, disjoint supports. That is a
`W(8,4)`. ∎ ✔

*Independent reproduction.* I recomputed, from scratch: **105** perfect
matchings; the rank distribution **7 / 42 / 56** for `rk D = 1 / 2 / 3`; the
number of transversal parallel classes **4 / 2 / 1**, matching `2^{3−rk D}` on
every one of the 105; and `|GL(3,2)| = 168`, `|AGL(3,2)| = 1344`. All agree
with §2.5 and with cert 28.

---

## 3. Replays and tampering

All from the worktree root, Python 3.14 on this desk, 2026-09-05. Two-process
limit and the twenty-minute cap respected; nothing ran longer than 45 s.

| run | checks | exit | wall | note claims |
| --- | --- | --- | --- | --- |
| `verify/verify.py --selftest` | 39 lines, `SELFTEST: PASS` | 0 | 0.6 s | — |
| `certs/26-theorem-t-structure/run.py` | 46 | 0 | **16.1 s** | 46 checks, 15.5–16.9 s ✔ |
| `certs/26-theorem-t-structure/run.py --full` | 64 | 0 | **44.8 s** | 64 checks, 45.4 s ✔ |
| `certs/27-sylow2-rigidity/run.py` | 14 | 0 | **13.1 s** | 14 checks, ≈ 12.6 s ✔ |
| `certs/27-sylow2-rigidity/run.py --full` | 14 | 0 | **13.7 s** | 13.9 s ✔ |
| `certs/28-ag32-transversal/run.py` | 14 | 0 | **0.77 s** | 14 checks, 0.7 s, no `--full` ✔ |
| `certs/29-lemma-f-records/run.py` | 40 | 0 | **1.28 s** | 40 checks, 1.3 s ✔ |

`--full` on cert 27 widens the box to `[−9,9]⁴` — the run prints `10000` odd
4-tuples against the default's `4096` in `[−7,7]⁴`, matching §1.10. Cert 29's
`--full` is present and priced (`--full`, `--matrix`), and was **not** run.

**Tampering (scratch copy under `/tmp`, never the worktree):**

| tamper | result |
| --- | --- |
| cert 26: one hex digit of the order-668 T-image pin | **exit 1**, `FAILED: order 668: the T-image is Hadamard … and its canonical digest matches the pin` |
| cert 27: `8384 → 8385` | **exit 1**, `SOME CHECK FAILED` |
| cert 28: `448 → 447` types | **exit 1**, `[FAIL] 448 types = 112 S-parts × 4 κ(ρ)` |
| cert 29: `pins.json` class size `332 → 331` | **exit 1**, three `[FAIL]`s on the partition totals |
| cert 29: one hex digit of `H(668)`'s digest in `pins.json` | **exit 1**, `[FAIL] H(668): assembled digest == pin` — proving the digest is **recomputed in-process** and not merely echoed |

**Corroboration of the banks.** All four `pins.json` digests reproduce against
certificates that predate this port (`bdeb5059…` in certs 01/02/06/08/13,
`600849b0…` in 02/06/08/13/15, `af1c285c…` in 13/15/19 and in `data/`,
`3adcb1bb…` in 01/02/11/14/15) — so the four matrices cert 29 reasons about are
the same objects §3 profiles. The four `sha256` file pins in PROVENANCE's table
and in cert 26's `FILE_PINS` all reproduce under `sha256sum`.

**Audit surface, stated plainly.** Certs 20–25 bank their heavy artefact — the
full profile — in `data/`, digest-pinned, from two implementations, and their
default paths bind it. Cert 29 banks only a **four-line summary** of the class
partition; the partition itself is in no file in this repository. The default
path therefore binds the *matrices* to their digests and checks the summary's
internal consistency, but cannot bind the summary to any computation. That is
exactly what the note and `NOTES.md` say, and the label is still earned (§4),
but it is a thinner surface than the rest of the repository offers, and
DISCLOSURE's description of what an audit path checks does not currently
describe it (**F7**).

---

## 4. Fidelity to the source, and the D-073 question

I read the source laboratory's `DECISIONS.md` (read-only) at **D-062** and
**D-068**, and its `intel/fleet-2026-09-01/skeptic-pass/NOTE-B-draft-2.md`.

* The port introduces **no strengthening**. D-068.2's census figures
  ("24/25 inequivalent at orders ≥ 44") appear in §1.9 unchanged and quarantined;
  D-068.4's `w > 4` repair survives verbatim in §2.4 and §2.5's scope; D-068.8a
  (branch-free Theorem T) and D-068.8e (`N` for the order) are applied as
  adopted; remark R is restored as D-068.8a required.
* **Cert 29 is the laboratory's code where it matters.** I diffed the ported
  `run.py` against `Hadamard-2060/certs/0028-lemma-f-records/run.py`: the
  `|T4|` route-B implementation, `classes_of`, `fmin_killed`,
  `implication_kills_all` and the entire `stage_full` are **byte-identical**.
  So §1.11's phrase "cert 29's own standard-library route B under `--full`" is
  faithful, not a flattering re-attribution — the `--full` path that ran at the
  laboratory is character-for-character the one offered here, and it gates on
  the same digest the port pins.
* **The `--full` run exists and its walls are real.** The four logs in
  `certs/0028-lemma-f-records/replay-2026-09-03/` end `ALL CHECKS PASS` at
  `9044.0`, `9062.2`, `9049.1` and `11273.6` s — the four numbers §1.11,
  `NOTES.md`, README and PROVENANCE quote, summing to `38 428.9 s ≈ 10.7`
  core-hours as tabled.
* **On the brief's D-073 citation.** The author is right that D-073 is
  PR-0049's route-B tier work and says nothing about Lemma F. But the label
  does not need a new ADR: **D-068 item 5 is self-executing** — "*`--full` …
  IS the missing second implementation; when it has run the label becomes
  PROVEN-BY-CERTIFICATE*". It has run, and I verified the logs. So
  PROVEN-BY-CERTIFICATE is licensed by the lab's own adopted text, and what is
  outstanding is D-068 item 10's bookkeeping (an OPEN to close), not a
  licensing gap. **This is not a blocker for the push.**
* **Cert 26's substitution is sound and correctly disclosed.** Porting
  `transpose_smallorders.py` would have put a finder-side
  individualisation–refinement search into the trust chain and banked lab-only
  data; refusing that is right. The replacement census is smaller and says so,
  reports four instances as *undecided* rather than equivalent, and runs no
  isomorphism search. The three eight-hex T-image prefixes the laboratory
  recorded reproduce here from this repository's own records — a real
  corroboration, and the `D̃MᵀD̃` vs `Mᵀ` distinction is flagged in three places,
  which it needed to be.
* **Lab clone, one observation.** `Hadamard-2060` is not byte-clean: `PLAN.md`
  carries an uncommitted turn-54 rewrite, mtime `2026-09-04 18:39:55`. That is
  **after** this lane's last commit (`0e72a58`, `18:26:56`) and its content is
  desk campaign state, so it is not this port's doing — but the author's
  "byte-clean" line is stale and the desk should know. Nothing else in the lab
  clone is modified.

---

## 5. Findings

**F1 — wrong symbol in a displayed inequality (§1.11, "Admissible primes").**
The note reads `|(AAᵀ)_{rr′}| ≤ 4 < m`. From Lemma F(iv),
`(AAᵀ)_{rr′} = −p·(BBᵀ)_{rr′}` off the diagonal, so the integer
`(BBᵀ)_{rr′}` vanishes because `4 < **p**`, not because `4 < m`. As written the
claim is false at the very configuration it describes: with `f = 4` and
`p = 523` at `N = 2092` one has `m = (N−4)/p`, and `4 < m` is not the operative
bound — indeed at `f = 0, m = 4` the stated inequality reads `4 < 4`. The
mathematics is fine; the symbol is not. The wording is inherited verbatim from
the laboratory's draft, so this is a defect the port should fix rather than
carry.
*Repair — `note/NOTE-B.md` §1.11:*
old `` `|(AAᵀ)_{rr′}| ≤ 4 < m` forces `(BBᵀ)_{rr′} = 0` ``
new `` `|(AAᵀ)_{rr′}| ≤ 4 < p` forces `(BBᵀ)_{rr′} = 0` (every odd prime but 3) ``
*and the same sentence in `certs/29-lemma-f-records/NOTES.md`.*

**F2 — a source-laboratory certificate number printed by a public
certificate.** `certs/28-ag32-transversal/run.py:526` prints
`"borders built by this construction for cert 0022's H(88) seed quadruples"`,
and line 527 comments `certs/0022-cell24-instances/run.py`. In *this*
repository the cell-24 certificate is **18**; `0022` is the laboratory's
private numbering and resolves to nothing here. The certificate's own docstring
(lines 43, 68) already says "cert 18", so the file contradicts itself, and the
port note's claim that cert 28 "cites cert 18 where the lab cited its own
`0022`" is only two-thirds true.
*Repair — `certs/28-ag32-transversal/run.py`:*
line 526: old `cert 0022's H(88) seed quadruples` → new `cert 18's H(88) seed quadruples`
line 527: old `certs/0022-cell24-instances/run.py` → new `certs/18-cell24-instances/run.py`

**F3 — an unearned WLOG (§2.5, step (3)).** "Column signed permutations of `Q′`
do not change `(1/2)wQ′ ∈ {±1}`, so take `Q′` signed-Sylvester." Column signed
permutations do **not** reduce every `H(8)` to the Walsh matrix: §2.4 states in
its own proof that there are **480** right-orbits under exactly that group, one
per `(AG(3,2)` structure, plane-sign class) pair. What is true is that every
`H(8)` is a signed **row** permutation of Sylvester up to column operations, and
that the conclusion — which the note then correctly states for the structure
`A(Q′)`, not for the standard one — is invariant under relabelling; step (6)
is uniform over all 105 matchings, so no generality is lost. But the reason
given does not carry the reduction.
*Repair — `note/NOTE-B.md` §2.5, step (3):*
old `` so take `Q′` signed-Sylvester. ``
new `` so take `Q′` signed-Sylvester up to a relabelling of the eight columns of `Qᵀ` — i.e. identify the label set with `𝔽₂³` through `Q′`'s own affine structure `A(Q′)`, which is one of the 30 (§2.4); step (6) is uniform over all 105 matchings, so the choice of structure is immaterial. ``

**F4 — remark R's seven records are the `--full` path (§1.9 and §5).** The
default run of cert 26 evaluates remark R on **three** records (668, 716, 1916);
1388, 1436, 1676 and 1772 appear only under `--full`. §1.9's remark and §5's
label row both say "all seven records" with no path qualifier, although the
neighbouring sentence about the Theorem-T identities does carry one.
*Repair — `note/NOTE-B.md` §1.9, remark R:* after
`(cert 26 [4], which evaluates the condition and the verdict on all seven records`
insert ` — 668, 716 and 1916 by default, the other four under --full —`;
*and the matching clause in §5's remark-R row.*

**F5 — the census count is the `--full` count (§1.9, COMPUTATIONAL-EVIDENCE).**
"eleven instances … Seven are proved inequivalent … four are undecided" is what
`--full` reports; the **default** run prints `census totals: 9 instances, 5
PROVEN inequivalent to their transpose, 4 undecided by the profile`, the two
`H(88)` being `--full`-only. `--full` has been run here, so nothing is
unsupported — but a reader replaying the default path meets a different number
in the one place the note does not warn them.
*Repair — `note/NOTE-B.md` §1.9:*
old `` and the four `(2,4)` matrices of §2.3 — ``
new `` and the four `(2,4)` matrices of §2.3, the two `H(88)` of them under `--full` (the default path reports nine instances, five separated) — ``
*and the same qualifier in `note/PORTS-2026-09-05-a.md`'s "11 of them, 7 separated".*

**F6 — README states Lemma F unconditionally.** "**Lemma F** (NOTE-B.md §1.11;
cert 29) forces a 16-block circulant array at `p = 523`" omits both the order
and the hypothesis. §1.11 itself is careful ("At `N = 2092`, `p = 523`"; "it is
a constraint on a differently-shaped `H(2092)`"), and the README paragraph is
the summary a reader meets first.
*Repair — `README.md`:*
old `` forces a 16-block circulant array at `p = 523` ``
new `` forces a 16-block circulant array on any `H(2092)` admitting an automorphism of order `p = 523` ``

**F7 — DISCLOSURE.md is now stale in two places.** The author left it untouched
and is right that no **stance** changed — nothing ported claims novelty,
priority or existence. But DISCLOSURE is also the results summary and the
description of the trust chain, and both have drifted:
(a) its dated additions stop at 2026-09-03, so three theorems and an existence
theorem that README and NOTE-B carry are absent;
(b) its "the word *replay* belongs to the optional `--full` paths of certs 06,
08, 11, 13, 14, 15, 19, 20, 21, 22, 23, 24 and 25" sentence defines which
certificates audit by default — and cert **29** is now exactly such a
certificate, with the added wrinkle that its audit surface is a summary rather
than a banked artefact (§3 above).
*Repair — `DISCLOSURE.md`:* add an "Added 2026-09-05:" sentence naming Theorem T
(with remark R), the Sylow-2 rigidity theorem, Lemma F and the records claim,
and the `(2,4)` existence theorem, each with its label and with the
`H(2092)` statements marked conditional and non-existential; and extend the
`--full` list to `… 24, 25 and 29`, noting that cert 29's default path audits a
**pinned summary** of a partition computed at the source laboratory rather than
a banked profile in `data/`.

**F8 (minor, no repair required).** Three wording slips I record but do not
insist on: §1.10's opening "two groups out of twenty-four" compares groups with
`(E, e*)` pairs (there are 12 groups and 24 pairs — the theorem statement and
README both get it right); §3.4's repaired passage now runs two consecutive
sentences beginning "Nothing"; and `certs/28-ag32-transversal/run.py`'s
`shutil.rmtree` at line 582 is not in a `finally`, so a failing run leaves its
private temporary directory behind.

---

## 6. Recommendation

**REPAIR-FIRST.** No statement is refuted, no label is unearned, no scope
sentence overreaches, and no separation, cell or existence claim was
strengthened in transit. Every certificate passes from the worktree root, fails
closed under tampering, is self-contained on the standard library, and leaves
the tree clean. The mathematics of §§1.9, 1.10, 1.11 and 2.5 stands after an
independent re-derivation of every step and an independent recomputation of
every finite input I could reach in seconds.

What stands between this branch and a push is seven text repairs: **F1** (a
wrong symbol in a displayed inequality), **F2** (a private cert number printed
by a public certificate), **F3** (a WLOG whose stated reason does not carry it),
**F4**/**F5** (two `--full`-only counts presented without the qualifier),
**F6** (an unconditional README sentence), and **F7** (DISCLOSURE's summary and
trust-chain description, now stale). F1 and F2 are the two I would not push
without.

Two items for the desk that are not repairs to this branch: the D-073 pointer
in the brief was wrong and **needs no replacement**, because D-068 item 5
licenses the label on its own terms and the run is on disk (D-068 item 10's
OPEN should simply be closed in the laboratory's record); and the laboratory
clone carries an uncommitted `PLAN.md` rewrite timestamped after this lane
finished, which is desk work, not this port's.
