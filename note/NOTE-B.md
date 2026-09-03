# NOTE-B — Bordered Goethals–Seidel arrays: theorems, instances, and separations at orders 668, 716, 1676, 1772 and 2060

This note is the mathematics of the repository, in three movements:

- **I. The theorems** (§1) — an exact characterisation of when a
  Goethals–Seidel array over a finite abelian group extends to a
  Hadamard matrix through a *coset border* of width `4s`; the forced
  parameter classification **within the house-Gram branch, under
  `w > 2s`** (with `n ≥ 3` the surviving cells are `(0,1)`, `(1,1)`
  and `i = s+1` with `s` odd); the **Gram rigidity theorem**, which
  determines every admissible `M` — at `i = s+1` exactly the
  real-character twists of the house form, and no others (§1.2.1),
  sharpened to `w > s` with the `w = s` boundary classified and its
  escape realised by seeds (Theorem E′, §1.7); the structure theorem
  for **border kits**, which reduces (H3) and (H4) to the single
  condition `E ∈ {±1}` and makes a kit an order-independent object,
  with the exhaustive `(3,4)` census that makes that cell
  **one-layer** — at 2092 a seed quadruple *is* the matrix
  (Theorem F, §1.8); and the complete resolution of the `s = 1`
  border system.
- **II. The instances** (§2) — twelve publicly posted matrices
  re-verified through the theorem's hypotheses, and twelve matrices
  constructed here, every one replayable from `certs/`. Four of them
  inhabit the **even-`s`** branch of the rigidity theorem, at the
  cell `(s,i) = (2,4)` outside the house form (§2.3), where the
  border layer is proved never to obstruct (§2.4).
- **III. Existence plus separation** (§3) — a proof that order 668
  carries at least **four** Hadamard equivalence classes and that
  orders 716, 1676, 1772 and 2060 each carry at least **three**,
  separated by an exactly computed invariant; and the invariant theory that makes the
  separations honest, including two documented traps and one
  sound-but-blind published statistic.

Every claim carries one of the honesty labels **PROVEN**,
**PROVEN-BY-CERTIFICATE** (an explicit artifact plus a green run of
`verify/verify.py`), **MEASURED** (a machine run on stated inputs),
**COMPUTATIONAL-EVIDENCE** (reproducible but not a proof),
**BOUNDED-NEGATIVE-SEARCH** (a search that closes only its stated
class), or **CONJECTURE**. §5 is the label table. §4 is prior art and
credit; it contains this repository's only hedged sentence.

---

## 1. Movement I — the theorems

### 1.0 Notation

`G` a finite abelian group, written additively, `|G| = n`.
`K ≤ G` a subgroup, `i = [G:K]`, `w = |K| = n/i`, `κ : G → Ḡ := G/K`
the quotient map, `s ≥ 0` an integer, `v = n + s`, `N = 4v`.
`ρ ∈ G` a fixed element (the "reflection shift").

For `x : G → {±1}`, `PAF_x(t) = Σ_{u∈G} x(u)x(u+t)`; note
`PAF_x(−t) = PAF_x(t)`. For a quadruple `x_0,…,x_3` write
`Σ PAF(t) = Σ_q PAF_{x_q}(t)`. Coset sums:
`σ_q(c) = Σ_{g ∈ κ^{-1}(c)} x_q(g)`, a function `Ḡ → ℤ`. Row sums:
`r_q = Σ_g x_q(g) = Σ_c σ_q(c)`.

**Type-1 development.** `X = dev(x)` is the `G×G` matrix
`X[g,h] = x(h−g)`. **Reflection.** `R = R_ρ` is the permutation matrix
`R[k,h] = [k+h = ρ]`; `R = Rᵀ`, `R² = I`. Then `(XR)[g,h] = x(ρ−g−h)`
and `(XᵀR)[g,h] = x(g+h−ρ)`. On `G = ℤ_n`, `ρ = −1` gives the
classical back-diagonal and `ρ = 0` is group inversion. A shift
`g ↦ g+α` sends `ρ ↦ ρ−2α`, so the invariant is the class
`[ρ] ∈ G/2G`: `ρ` can be normalised to `0` exactly when `ρ ∈ 2G` —
automatic for odd `|G|`, and not guaranteed when `|G|` is even.

**The core.** `C = GS(x_0,x_1,x_2,x_3; ρ)` is the `4G × 4G` array

```
       [  A     BR     CR     DR  ]
       [ -BR     A    DᵀR   -CᵀR  ]        A = dev(x_0), B = dev(x_1),
       [ -CR   -DᵀR    A     BᵀR  ]        C = dev(x_2), D = dev(x_3)
       [ -DR    CᵀR  -BᵀR     A   ]
```

We call this orientation *standard* — a convention of this
repository, not of the literature; negating the six transposed
blocks gives the other valid orientation (used, e.g., by SageMath).
The proofs below are orientation-symmetric, but **every theorem in
this note is stated, and every primary trust-chain assembly in this
repository uses, the standard orientation only** (cert 05
additionally carries an alternate-orientation control instance); the
analogues for the other orientation are routine and are not derived
here.

**The border ansatz.** Given a corner `E ∈ {±1}^{4s×4s}`, a row table
`P ∈ {±1}^{4s×4i}` and a column table `Q ∈ {±1}^{4i×4s}`, set

```
H = [ E   P̃ ]      P̃[r, (J,h)] = P[r][iJ + κ(h)]
    [ Q̃   C ]      Q̃[(I,g), c] = Q[iI + κ(g)][c]
```

i.e. the first `4s` rows and columns are a corner plus strips that are
**constant on each K-coset inside each superblock**. `i = 1` (border
strips constant on all of `G`) recovers the classical bordered arrays
of Wallis–Whiteman and Spence (§4); `i ≥ 2` — strips that see the
coset structure of a proper subgroup — is the object of study here.

### 1.1 Theorem A (exact characterisation) — PROVEN

**Theorem A.** With the data above, `H Hᵀ = N·I_N` **if and only if**

- **(H1)** there is a symmetric, `Ḡ`-invariant `i×i` integer matrix
  `M` (i.e. `M[c,c'] = M(c−c')` for `M : Ḡ → ℤ`) with
  `Q Qᵀ = I_4 ⊗ M`, and
- **(H2)** `Σ PAF(t) = −M(κ(t))` for all `t ≠ 0` (at `t = 0` the
  sum is `4n` automatically), and
- **(H3)** `E Eᵀ + w·P Pᵀ = N·I_{4s}`, and
- **(H4)** `E Qᵀ + P Ĉᵀ = 0`, where `Ĉ = GS(σ_0,σ_1,σ_2,σ_3; κ(ρ))`
  is the Goethals–Seidel array of the **coset sums over the quotient
  group Ḡ**.

`M(0) = 4s` is automatic (rows of `Q` are `±1` of length `4s`), so
the "`−4s` on `K∖{0}`" tier of the profile is **forced by the
ansatz**, not chosen. (H4) depends on the seeds only through `σ` —
the content of the compression lemma (Lemma 3), which also makes (H4)
a `4s×4i` condition rather than a `4s×4n` one.

*Proof.* Three lemmas, then the block computation.

**Lemma 1.** For `x, y : G → ℤ`, with `X = dev(x)`, `Y = dev(y)`:
(a) `XY = YX` (convolution over an abelian group);
(b) `X Yᵀ = Yᵀ X`, and `X Xᵀ = Xᵀ X = dev(PAF_x∘(−))`, i.e.
`(X Xᵀ)[g,h] = PAF_x(g−h)`;
(c) `X R Yᵀ = Y R Xᵀ` (both equal `f(ρ−g−h)` with
`f(t) = Σ_u x(u)y(t−u)`);
(d) `X R Y = Yᵀ R Xᵀ` (both equal `Σ_u x(ρ−g−h+u)y(u)`).
Each is a one-line change of summation variable; (c) and (d) are
where `R` being the *reflection* `k+h = ρ` (rather than an arbitrary
involution) is used. ∎

**Lemma 2 (core Gram).** `C Cᵀ = I_4 ⊗ Σ`, where
`Σ[g,h] = Σ PAF(g−h)`.

*Proof.* Diagonal blocks: `R Rᵀ = I` collapses every reversed block,
and `XᵀX = XXᵀ` (Lemma 1b) makes all four diagonal blocks equal
`Σ_q X_q X_qᵀ = Σ`. Off-diagonal blocks: each is a sum of four
products which cancel in pairs of the four displayed types —
`−A R Bᵀ + B R Aᵀ` (zero by 1c), `CD − DC` (zero by 1a),
`B Cᵀ − Cᵀ B` (zero by 1b), `−A R D + Dᵀ R Aᵀ` (zero by 1d) — two
such cancelling pairs for each of the six off-diagonal block
positions above the diagonal. Nothing but Lemma 1 is used, so **no
hypothesis on the seeds enters**: Lemma 2 holds for arbitrary
sequences. ∎

That last sentence is why the border is cheap: the 16-block sign
bookkeeping is seed-independent, and the aggregate `Σ` is the only
channel through which the seeds reach the Gram.

**Lemma 3 (compression).** For every `(I,J)` and every `g`,
`Σ_{h ∈ κ^{-1}(c)} C[(I,g),(J,h)]` depends on `g` only through
`κ(g)`, and the resulting `4i×4i` matrix is
`Ĉ = GS(σ_0,…,σ_3; κ(ρ))`, the Goethals–Seidel array over `Ḡ`.

*Proof.* Three cases, each a bijection of a coset onto a coset:
`Σ_{h∈κ^{-1}(c)} x(h−g) = σ_x(c − κg)`;
`Σ_{h∈κ^{-1}(c)} x(ρ−g−h) = σ_x(κρ − κg − c)`;
`Σ_{h∈κ^{-1}(c)} x(g+h−ρ) = σ_x(κg + c − κρ)`.
These are precisely the entries of the array over `Ḡ` at `(κg, c)`,
and the block signs are untouched by summation. ∎

The function-level content of Lemma 3 — coset-summing intertwines
the PAF with the quotient PAF — is the **Đoković–Kotsireas
compression device**, stated by them at full abelian generality
(§4); what Lemma 3 adds is only the matrix-level bookkeeping that
the compressed core is again a Goethals–Seidel array over `G/K`,
reflection included, which is what (H4) consumes.

**Corollary (Σ̄ law).** `Ĉ Ĉᵀ = I_4 ⊗ Σ̄` with
`Σ̄(ē) = Σ_{t ∈ κ^{-1}(ē)} Σ PAF(t)`; in the house form of Theorem B,
`Σ̄ = 4(w(i−s−1)+s)·I_i + 4w·J_i`. (This identity is
Đoković–Kotsireas compression applied to the aggregate PAF — see §4;
nothing in it is new here, and it is stated only for its coupling
role in (H4).)

*The four blocks of `H Hᵀ`.*

```
H Hᵀ = [ E Eᵀ + P̃ P̃ᵀ      E Q̃ᵀ + P̃ Cᵀ ]
       [      (…)ᵀ         Q̃ Q̃ᵀ + C Cᵀ ]
```

**Top-left.** `P̃[r,(J,h)] = P[r][iJ+κ(h)]` and each class `(J,c)`
has exactly `w` members, so `P̃ P̃ᵀ = w·P Pᵀ`. Top-left `= N·I_{4s}`
is exactly **(H3)**.

**Bottom-right.** `(Q̃ Q̃ᵀ)[(I,g),(I',g')] = (Q Qᵀ)[iI+κg, iI'+κg']`
and, by Lemma 2, `(C Cᵀ)[(I,g),(I',g')] = δ_{I,I'}·Σ PAF(g−g')`.
Requiring `N δ_{I,I'} δ_{g,g'}`:

- `I ≠ I'` ⟹ `Q Qᵀ` has zero cross-superblock blocks;
- `I = I'`, `g = g'` ⟹ `(Q Qᵀ)[iI+c,iI+c] = N − 4n = 4s`, automatic;
- `I = I'`, `g ≠ g'` ⟹ `(Q Qᵀ)[iI+κg, iI+κg'] = −Σ PAF(g−g')`.

The left side depends on `g,g'` only through their cosets and the
right side only through `g−g'`, and every `t` in the coset
`κ(g)−κ(g')` is realised; so the two conditions are jointly
equivalent to (H1)+(H2). This also shows what the ansatz forces on
the profile: the value `−4s` on `K∖{0}` (the diagonal of a `±1` Gram
is constant), and off the origin the aggregate PAF **factors through
`G/K`** (`Q Qᵀ` cannot see anything finer than a coset). In the
house branch all nonzero quotient classes coalesce to the common
value `+4`, giving the two-tier profile of Theorem B; a general
admissible `M` may distinguish nonzero classes. This is the sentence
the whole construction turns on.

**Top-right.** `(E Q̃ᵀ)[r,(I,g)] = (E Qᵀ)[r, iI+κg]`, and by Lemma 3
`(P̃ Cᵀ)[r,(I,g)] = (P Ĉᵀ)[r, iI+κg]`. Vanishing is exactly **(H4)**.
Bottom-left is the transpose of top-right. ∎

### 1.2 Theorem B (the house form) — PROVEN

**Theorem B.** Take `M = (4s+4)·I_i − 4·J_i`. Then (H1)+(H2) become
exactly the **two-tier PAF profile**

```
Σ_q PAF_q(g) = 4n·δ_0 − 4s·[g ∈ K∖{0}] + 4·[g ∉ K]
```

together with `Q Qᵀ = I_4 ⊗ ((4s+4)I_i − 4J_i)`, and these plus (H3)
and (H4) are sufficient *and* necessary within the ansatz.

Degenerations. `s = 0` (take `K = G`, `i = 1`): the border is empty,
the profile becomes `Σ PAF(t) = 4n δ_0`, and Theorem B **is the
classical Goethals–Seidel theorem over an abelian group**
(Wallis–Whiteman 1972, Theorem 11 — see §4). `s = 1, i = 1`: the
profile is `Σ PAF(t) = −4` for `t ≠ 0`, and Theorem B **is the
Wallis–Whiteman/Spence bordered construction** (§4).

**On the choice of `M`.** Theorem A admits in principle other Gram
matrices `M`, and §1.2.1 determines all of them: under `w > 2s`
every admissible `M` is `4i` times an orthogonal projector, and at
`i = s+1` it is a real-quotient-character twist of the house form,
the twist family being complete. So "any admissible `M` is house" is
**false as a literal statement** — a character twist `x_q ↦ ψ x_q`
with `ψ² = 1`, `ψ|_K = 1` produces valid instances whose Gram is
non-house, and an explicit `N = 1916` witness is banked (§3.3
explains why, when `ψ(ρ) = 1`, that witness is a diagonal
conjugation of the house instance and not a new matrix) — but it is
false in a completely described way. The twists are the whole story,
known in advance rather than quotiented out on faith. Any exhaustive
treatment of admissible profiles at `i ≥ 4` must still quotient by
the twist; what changes is that quotienting now provably loses
nothing, and that the `(s,i) = (1,2)` picture of Theorem D below —
the house Gram in the genuine branch, one degenerate alternative —
is the shape of every cell rather than an accident of the smallest.

#### 1.2.1 Theorem E (Gram rigidity) — PROVEN

Write `M̂(χ) = Σ_{c∈Ḡ} M(c)χ(c)` for `χ` in `Ĝ̄`, the character
group of `Ḡ` (of order `i`, with values in `μ_e`, `e = exp Ḡ`), and
`v_χ ∈ ℂ^i` for the vector `v_χ(c) = χ(c)`. For a `Ḡ`-invariant
`Y[c,c'] = Y(c−c')` with `Y(−c) = Y(c)`, substituting `u = c−c'`,

```
(Y v_χ)[c] = Σ_{c'} Y(c−c')χ(c') = χ(c)·Σ_u Y(u)χ(−u)
           = Ŷ(χ̄)·χ(c) = Ŷ(χ)·χ(c),
```

the last step by `Y(−u) = Y(u)`. So the `i` vectors `v_χ` are an
eigenbasis and the eigenvalue **multiset** of `Y` is
`{ Ŷ(χ) : χ ∈ Ĝ̄ }`, all real. Used below without comment.

**Theorem E.** Assume (H1) and (H2), with `s ≥ 1` and `w > 2s`. Put
`S = { χ ∈ Ĝ̄ : M̂(χ) ≠ 0 }`. Then `s ≤ i`, `|S| = s`, `S` is stable
under complex conjugation and under the Galois action `χ ↦ χ^k`
(`k` coprime to `e`), and

```
       M(c)  =  4 · Σ_{χ ∈ S} χ(c)          for every c ∈ Ḡ,
```

i.e. `M = 4i·P_S`, the `4i`-multiple of the orthogonal projector
onto `span{ v_χ : χ ∈ S }`. In particular the eigenvalue multiset of
every admissible `M` is `{ 0^(i−s), (4i)^s }` — the pinning holds at
every cell `(s,i)`, not only at `i = s+1`. Neither (H3) nor (H4) is
used, so this is a statement about the seed-and-column-table layer
alone.

*Proof.* Five steps: a Parseval transfer, the arithmetic of `M`, the
rank, the Galois stability of `S`, and an integrality lock.

**(L1) Parseval transfer.** For `x : G → {±1}` and a character `ψ`
of `G`, using `ψ(g)ψ̄(h) = ψ(g−h)` and substituting `g = h+t`,

```
|x̂(ψ)|² = Σ_{g,h} x(g)x(h) ψ(g)ψ̄(h) = Σ_{h,t} x(h+t)x(h) ψ(t)
        = Σ_{t∈G} ψ(t)·PAF_x(t).
```

Take `ψ = χ∘κ`, sum over the quadruple, and split off `t = 0`, where
`Σ_q PAF_q(0) = 4n`. By (H2) the rest is `−Σ_{t≠0} χ(κt)M(κt)`; add
back and subtract the `t = 0` term, and use that each fibre
`κ^{-1}(c)` has exactly `w` elements, so
`Σ_{t∈G} χ(κt)M(κt) = w·M̂(χ)`:

```
       Σ_q |x̂_q(χ∘κ)|²  =  4n + M(0) − w·M̂(χ)  =  N − w·M̂(χ),
```

using `M(0) = 4s` from (L2) and `N = 4n+4s`. The left side is a sum
of squared absolute values, hence `≥ 0`:

> **the Parseval window.** `M̂(χ) ≤ N/w = 4i + 4s/w` for **every**
> `χ ∈ Ĝ̄`, and under `w > 2s` this reads `M̂(χ) < 4i + 2`.

This is the single channel through which the seeds constrain the
Gram; everything below is the window plus integrality.

**(L2) Arithmetic of `M`.** `M(c−c') = ⟨Q[iI+c], Q[iI+c']⟩` is an
inner product of two `±1` rows of length `4s`. Hence `M(0) = 4s`,
`|M(c)| ≤ 4s`, and, with `k` the number of disagreeing coordinates,
`M(c) = 4s − 2k ∈ 2ℤ`. `M(−c) = M(c)` is symmetry read through
`Ḡ`-invariance; `M = Q_0 Q_0ᵀ` is PSD; and
`4·rank M = rank(I₄⊗M) = rank(QQᵀ) = rank Q ≤ 4s` because `Q` has
`4s` columns, so `rank M ≤ s`. The evenness is what makes the
argument close, and nothing finer is available: `M(c) mod 4` is not
determined, since the parity of `k` is free.

**(L3) The rank is exactly `s`, and `s ≤ i`.** Write `r = rank M`,
so `r ≤ min(s,i)`. `M` is PSD, so `trace M = i·M(0) = 4si` is a sum
of `r` positive eigenvalues, each `< 4i+2` by (L1). Hence
`4si < r·(4i+2)`, i.e. `r > 2si/(2i+1)`. Bootstrap twice. *First,*
`r ≤ i` forces `i > 2si/(2i+1)`, i.e. `2i+1 > 2s`, i.e. **`s ≤ i`**
— re-proved here without (H3). *Second,* `2si/(2i+1) = s − s/(2i+1)`
and `s ≤ i` makes `s/(2i+1) < 1`, so `r > s−1`; with `r ≤ s` and
`r ∈ ℤ`, **`r = s`**. (This is the "rank exactly `s`" of §1.2, and
"exactly" is not a slip: it needs the trace-versus-window argument,
not just `rank ≤ s` from the four orthogonal superblock spans.)

**(L4) `S` is Galois-stable.** `|S| = rank M = s` by the dictionary
and (L3). Let `σ_k ∈ Gal(ℚ(ζ_e)/ℚ)` be `ζ_e ↦ ζ_e^k`. Since
`M(c) ∈ ℤ`,
`σ_k(M̂(χ)) = Σ_c M(c)σ_k(χ(c)) = Σ_c M(c)χ(c)^k = M̂(χ^k)`, and
`σ_k` fixes `0`, so `M̂(χ) = 0 ⟺ M̂(χ^k) = 0`. Complex conjugation
is `σ_{−1}`.

**(L5) The integrality lock.** Put `M_S(c) := 4·Σ_{χ∈S}χ(c)`.
Because `σ_k` permutes `S` by (L4), `Σ_{χ∈S}χ(c)` is a
Galois-invariant algebraic integer, hence a rational integer, so
`M_S(c) ∈ 4ℤ`; conjugation-stability gives `M_S(−c) = M_S(c)`;
character orthogonality `Σ_c χ'(c)χ(c) = i·[χ' = χ̄]` gives
`M̂_S(χ) = 4i·[χ ∈ S]`, so `M_S = 4i·P_S` and `M_S(0) = 4s`. Set
`Y := M_S − M`, a `Ḡ`-invariant symmetric matrix. Three facts:

1. **`Y` has zero diagonal.** `Y(0) = 4s − 4s = 0`.
2. **`Y` has even entries.** `M_S(c) ∈ 4ℤ`, `M(c) ∈ 2ℤ`.
3. **`λ_min(Y) ≥ −4s/w`.** `Ŷ(χ) = M̂_S(χ) − M̂(χ)` vanishes off `S`
   (both terms do) and equals `4i − M̂(χ) ≥ 4i − (4i + 4s/w)` on
   `S`, by (L1).

So `C := Y/2` is symmetric with **integer** entries, **zero
diagonal**, and `λ_min(C) ≥ −2s/w > −1` — that last strict
inequality being exactly `w > 2s`. But an integer symmetric matrix
with zero diagonal and a nonzero entry `x = C[a,b]`, `a ≠ b`, has
`λ_min ≤ −1`: interlace against the principal submatrix
`[[0,x],[x,0]]`, whose eigenvalues are `±|x|` with `|x| ≥ 1`;
equivalently `uᵀCu = −|x|` for the unit vector
`u = (e_a − sgn(x)·e_b)/√2`. Hence `C = 0`, `M = M_S`, and the
spectrum is `M̂_S`'s. ∎

The proof is sharp in `w`: at `w = 2s` the bound reads
`λ_min(C) ≥ −1` and the interlacing step is silent. Two hypotheses
earn their place, and it is worth saying where. Galois-stability is
not decoration — with `S` merely conjugation-stable, `Σ_{χ∈S}χ(c)`
would be a real algebraic integer that need not be rational, `Y`
would not be an integer matrix, and (L5) would collapse. And the
zero diagonal is free, because `M_S(0) = 4s` and `M(0) = 4s` are
forced by *different* mechanisms — character orthogonality on one
side, `±1` rows of length `4s` on the other. A zero-trace integer
matrix cannot be spectrally small unless it is zero; that is the
whole theorem.

**Corollary E1 (the twist family is complete at `i = s+1`).**
Assume in addition `i = s+1`. Then `Ĝ̄ ∖ S = {χ₀}` is a singleton,
`χ₀² = 1`, and

```
       M  =  χ₀ · ( (4s+4)·I_i − 4·J_i ),   i.e. M(c) = χ₀(c)·M_house(c).
```

Conversely every real character of `Ḡ` occurs: the Lemma-T twist of
§1.4 by `ψ = χ₀∘κ`, with `Q' = D̄Q` and `D̄` diagonal on the row
index, `D̄[iI+c] = χ₀(c)`, carries any (H1)+(H2) pair to one with
Gram `χ₀·M`. **The admissible Grams at `i = s+1` are therefore
exactly the `|Ĝ̄[2]|` real-character twists of the house form, and
there are no others.**

*Proof.* `|Ĝ̄ ∖ S| = i − s = 1`, and the complement is
conjugation-stable (Theorem E), so `χ̄₀ = χ₀`: `χ₀` is `±1`-valued
and `χ₀² = 1`. Then
`M(c) = 4Σ_{χ≠χ₀}χ(c) = 4(i·[c=0] − χ₀(c)) = χ₀(c)·4(i·[c=0] − 1)`,
using `χ₀(0) = 1` and `χ₀(c)² = 1`; at `i = s+1` the `c = 0` value
is `4s` and the rest is `−4χ₀(c)`, which is `χ₀·M_house`. For the
converse, `PAF_{ψx}(t) = ψ(t)·PAF_x(t)` (Lemma T) turns (H2)'s
`−M(κt)` into `−χ₀(κt)M(κt)`, and

```
(Q'Q'ᵀ)[(I,c),(I',c')] = χ₀(c)χ₀(c')·δ_{I,I'}M(c−c')
                       = δ_{I,I'}·χ₀(c−c')·M(c−c')
```

because `χ₀` is real — so (H1) transports with the same twisted
Gram. ∎

A non-real `χ₀` never arises, and the reason is the cleanest way to
see why the answer had to be a real character: if `χ̄₀ ≠ χ₀` then
both would lie in the kernel of the real matrix `M`, forcing
`rank M ≤ i−2 = s−1` against (L3).

**What Corollary E1 does not settle.** It classifies the *Grams*,
i.e. the (H1)+(H2) layer. Which twists carry a complete bordered
instance, (H3) and (H4) included, is a different question. When
`ψ(ρ) = 1` the whole instance transports and manufactures no new
matrix (§1.4's proposition; at `N = 1916` all four twists occur —
§3.3 and cert 09, whose four Grams are now known to be the complete
admissible list at that cell, cert 12). When `ψ(ρ) = −1` nothing
here decides it. And Theorem E says nothing about assembled
matrices: not about Hadamard equivalence, not about orders 668 or
716.

**Corollary E2 (the seeds' aggregate spectrum is two-valued).**
Under the hypotheses of Theorem E, for every `χ ∈ Ĝ̄`,

```
   Σ_q |x̂_q(χ∘κ)|²  =  4s    if χ ∈ S      (s characters)
                    =  N     if χ ∉ S      (i − s characters),
```

immediately from (L1) with `M̂(χ) ∈ {4i, 0}` and `N − 4iw = 4s`. The
absolute values are not cosmetic — `Ĝ̄` can hold non-real characters,
and does at `(3,4)` when `Ḡ = ℤ₄`. For *real* `χ` the four
`x̂_q(χ∘κ)` are integers with `Σ_q x̂_q² = 4s` on `S`, and since
`x̂_q(χ∘κ) ≡ n (mod 2)` they are all even when `n` is even — so at
`(3,4)` exactly three of them are `±2` and one is `0`, a sharp
checkable constraint on every `(3,4)` record. At `(s,i) = (1,2)`
this is (D-d)'s `Σ_q δ_q² = 4` (§1.5); in general it pins the
quadruple's aggregate energy at every "live" lifted character to
`4s`, with no slack at all.

Label for Theorem E, Corollary E1 and Corollary E2: **PROVEN**
(paper-grade), with the `(3,4)` classification **MEASURED** and the
general-branch cells **PROVEN-BY-CERTIFICATE** in cert 12.

### 1.3 Theorem C (the parameters are forced) — PROVEN

Assume the house form and the non-degeneracy hypothesis `w > 2s`.

| | statement |
| --- | --- |
| **D3** | `w > 2s ⟹ E` is a `4s×4s` Hadamard matrix and `P Pᵀ = 4i·I_{4s}` |
| **D1** | `i ≤ s+1` |
| **D2** | for `s ≥ 1`: `i = s+1 ⟹ i` even `⟹ s` odd |
| **D4** | `i ≥ 2 ⟹ w(s+1−i) ≤ s` |
| **D6** | `i = 1 ⟹ n(s−1) ≤ s`, so `i = 1` dies for `s ≥ 2` (once `n ≥ 3`) |
| **D5** | `Σ_q r_q² = 8n − 4w(s+1) + 4s`; at `i = s+1` this equals `N` |

`D2` is stated for `s ≥ 1` only: at `s = 0` the (empty-border) case
`i = 1 = s+1` is realised with `i` odd, so the parity clause needs
the border to exist.

`D2` does not, however, need the house branch: it is a corollary of
Theorem E (§1.2.1). Under (H1)+(H2) and `w > 2s`, at `i = s+1`
Corollary E1 gives `M v_{χ₀} = 0` for a *real* character `χ₀`, so
`|Q_Iᵀ v_{χ₀}|² = v_{χ₀}ᵀ M v_{χ₀} = 0` legitimately over `ℝ`, hence
`Σ_{c∈Ḡ} χ₀(c)·Q[iI+c][r] = 0` — a signed sum of `i` units, which
vanishes only for `i` even, so `s = i−1` is odd. The house proof
below is this argument with `χ₀` the trivial character and `v_{χ₀}`
the all-ones vector.

*Proof of D3.* Off-diagonal of (H3):
`(E Eᵀ)[r,r'] = −w·(P Pᵀ)[r,r']`. Entries of `P Pᵀ` are sums of `4i`
signs, hence even; entries of `E Eᵀ` are bounded by `4s`. If
`w > 2s` then `|w·(even)| ≤ 4s < 2w` forces that even number to be
`0`; hence `E Eᵀ` is diagonal with diagonal `4s` (`E` is Hadamard)
and `P Pᵀ = 4i·I_{4s}`. A byproduct: `4s ≤ 4i`, so `w > 2s ⟹ s ≤ i`. ∎

*Proof of D1, D2.* `M = (4s+4)I_i − 4J_i` has eigenvalues
`4s+4−4i` (once) and `4s+4` (`i−1` times). `M` is a Gram matrix,
hence PSD: `i ≤ s+1` (D1). `Q` is `4i×4s`; the four superblock row
groups have mutually orthogonal spans in `ℝ^{4s}`, each of rank
`rank(M)`, so `4·rank(M) ≤ 4s`, i.e. `rank(M) ≤ s`. If `i = s+1`
(with `s ≥ 1`, so the border is nonempty) the all-ones vector is in
the kernel of `M`, so the `i` rows of each group sum to `0`; a sum of
`i` signs is `0` only if `i` is even, so `s` is odd (D2). ∎

*Proof of D4, D5, D6.* Summing the profile over `G` gives
`Σ_q r_q² = 4n − 4s(w−1) + 4(n−w) = 8n − 4w(s+1) + 4s` (D5); at
`i = s+1`, `wi = n` turns this into `4n+4s = N`. Non-negativity at
`i = 1` (`w = n`) gives `n(s−1) ≤ s` (D6). For D4, Fourier-transform
the coset sums: for a character `χ` of `Ḡ`,
`Σ_q |σ̂_q(χ)|² = Σ_ē χ(ē) Σ_{t∈κ^{-1}(ē)} Σ PAF(t)`. In the house
form the `ē = 0` term is `4n − 4s(w−1)` and each `ē ≠ 0` term is
`4w`, so for `χ ≠ 1` the total telescopes to
`4n − 4sw + 4s − 4w ≥ 0`, i.e. `w(s+1−i) ≤ s` — which needs `i ≥ 2`
for a nontrivial `χ` to exist. ∎

**Corollary (classification).** Under `w > 2s` and `n ≥ 3`, exactly
three cells survive:

```
(s,i) = (0,1)      plain Goethals–Seidel
(s,i) = (1,1)      the classical bordered array (Wallis–Whiteman / Spence)
 i = s+1, s odd    the coset-border family
```

*Proof.* If `i = 1`: D6 gives `s ≤ 1`. If `i ≥ 2`: D4 gives
`w(s+1−i) ≤ s < w/2`, so `s+1−i < 1/2`, so `i ≥ s+1`; with D1,
`i = s+1`; with D2, `s` odd. ∎

**The house-branch qualifier is necessary.** `D1` (`i ≤ s+1`) is not
a consequence of (H1)–(H4) alone, and the counterexamples are in the
repository's own bank. The decoded `(s,i) = (1,1)` record at order
1676 (`G = ℤ₄₁₈`), re-read on the subgroup `K ≤ G` of index
`i = 11` — the column table repeated along each coset, the row table
likewise — satisfies all four hypotheses with `M ≡ 4` on
`Ḡ = ℤ₁₁`, `w = 38 > 2s`, and `11 > s+1 = 2`. Its Gram is `4J₁₁`,
which is `4i·P_S` for `S = {1}`, exactly Theorem E's prediction: the
trivial character is the only Galois-stable singleton of `ℤ̂₁₁`, so
no other Gram is available there. The indices `11, 19, 22, 38` all
work at 1676 and `13, 17, 26, 34` at 1772; cert 12 machine-checks
all eight. What they do not do is produce anything: each re-reading
leaves `P̃` and `Q̃`, hence the assembled matrix, unchanged — it is
the `i = 1` instance in `i`-fold coordinates, checked as such in
cert 12. The *general* statement, that every (H1)–(H4) instance
with `i > s+1` collapses this way, is proved only at
`(s,i) = (1,2)`, where it is (D-a′) of §1.5; at general `i` it is
open. What is proved here is the negative, and it is why Theorem C's
classification is stated inside the house branch rather than for
(H1)–(H4).

So `i = s+1` is a theorem for `s ≥ 2`, and at `s = 1` both `i = 1`
and `i = 2` are legal — and §1.4–§1.5 show they are *the same
problem*. The twelve records of §2 land as `(0,1)×5`, `(1,1)×4`,
`(3,4)`, `(5,6)`, `(7,8)`: the classification is exactly saturated.

### 1.4 Lemma T and the twist — PROVEN

**Lemma T.** Let `ψ` be a character of `G` with `ψ² = 1` and kernel
`K` of index `2`. Then `PAF_{ψx}(t) = ψ(t)·PAF_x(t)`. Hence a
quadruple satisfying the `s=1, i=1` profile (`Σ PAF(t) = −4`,
`t ≠ 0`) twists into one satisfying the `s=1, i=2` profile (`−4` on
`K∖0`, `+4` off `K`), and back. The map is exact only at `s = 1`
(it turns `−4s` into `+4s`, and the target off-`K` value is `+4`).

**Proposition (a twist with `ψ(ρ) = 1` is a diagonal conjugation).**
Let `ψ` be a character of `G` with `ψ² = 1`, `ψ|_K = 1` and
`ψ(ρ) = 1`. Then the twisted instance `x'_q = ψ x_q`, `P' = P D̄`,
`Q' = D̄ Q`, `E' = E` assembles to `S H S` with
`S = diag(I_{4s}, I_4 ⊗ diag(ψ(g)))`.

*Proof.* `A'[g,h] = ψ(h−g)x_0(h−g) = ψ(g)ψ(h)A[g,h]`;
`(BR)'[g,h] = ψ(ρ−g−h)x_1(ρ−g−h) = ψ(ρ)ψ(g)ψ(h)(BR)[g,h]`;
`(XᵀR)'[g,h] = ψ(g+h−ρ)x(g+h−ρ) = ψ(g)ψ(h)(XᵀR)[g,h]` — using
`ψ² = 1` and `ψ(ρ) = 1`. Hence `C' = D C D`. The border strips are
constant on `K`-cosets and `ψ` factors through `Ḡ`, so `P̃' = P̃ D`
and `Q̃' = D Q̃`, while `E' = E`. Therefore `H' = S H S`. ∎

(The hypothesis `ψ² = 1` is not a convenience: at `i = s+1` it is
**forced by rank** — a non-real kernel character would put its
conjugate in the kernel of the real matrix `M` too, dropping
`rank M` below `s` against Theorem E's step (L3). §1.2.1.)

So the twist **never manufactures a new matrix when `ψ(ρ) = 1`**,
however non-house its Gram. When `ψ(ρ) = −1` the argument breaks:
the conjugation produces the array with all twelve off-diagonal
blocks negated, and that sign pattern is not of rank one
(`s_I s_J = −1` for all `I ≠ J` is impossible since
`(s₁s₂)(s₁s₃)(s₂s₃) = (s₁s₂s₃)² = +1`). At orders 668, 716,
1676 and 1772 — where `G = ℤ₁₆₆`, `ℤ₁₇₈`, `ℤ₄₁₈` and `ℤ₄₄₂`, `ρ = −1` in
each, and the unique `ψ` has `ψ(ρ) = −1` — the twisted matrix is in fact
**provably a new equivalence class** (§3.4, §3.6, §3.7, §3.8). Whether a
twist with `ψ(ρ) = −1` always leaves the equivalence class is not
claimed in general; 668, 716, 1676 and 1772 are the **four** proven
instances (§3.4, §3.6, §3.7, §3.8), which is every decoded `(1,1)`
order in this repository and still not a theorem.

**Remark (the twist and the orientation switch).** With `ψ² = 1`,
`D = diag(ψ(g))` and `D₄ = I₄⊗D`, the conjugation above reads
`D₄·C(x;ρ)·D₄ = C(ψx;ρ)` when `ψ(ρ) = 1`, and when `ψ(ρ) = −1` it
produces instead `C^{sw}(ψx;ρ)` — the array with its twelve
off-diagonal blocks negated. That array is exactly the object certs
13, 14, 20 and 23 study. At 668, 716, 1676 and 1772 they build `H″`, the
decoded record with its twelve off-diagonal core blocks negated and
the border untouched, verify the sign-pattern identity
`S·H″·S = H_alt` (the alternate Goethals–Seidel orientation with the
border strips signed by superblock) cell by cell, and prove `H″`
inequivalent to the decoded record, to the cert-02 twist `H′`, and —
at 668 — to the Hall switch as well. So the orientation switch and the `ψ(ρ) = −1` twist
share a core and differ in how the border is transported, and that
difference is enough to change the equivalence class: one twisted
core, two inequivalent bordered completions. (The identity
`D̃·H_t·D̃ = H″` for the transport instance `H_t` — twisted seeds with
the decoded border carried along `ψ` — is recorded in
`certs/13-668-orientation/NOTES.md` as a source-laboratory
computation. It is not re-derived in this repository, and nothing
stated here rests on it.)

### 1.5 Theorem D (the `s = 1, i = 2` border system, resolved) — PROVEN

Setting: `G` abelian of even order `n ≥ 4` (equivalently `w ≥ 2`;
the vacuous `n = 2` case is excluded); `K ≤ G` of index `2`;
`w = n/2`; `s = 1`; `v = n+1`; `N = 4v`; `κ : G → ℤ₂`;
`ε = +1` if `κ(ρ)=0` else `−1`;
`r_q = σ_q(0) + σ_q(1)`, `δ_q = σ_q(0) − σ_q(1)`,
`d = (δ_0, εδ_1, εδ_2, εδ_3)`; and

```
          [  y0   y1   y2   y3 ]
Λ(y)  =   [ −y1   y0   y3  −y2 ]        Λ(y)Λ(y)ᵀ = (Σ y_q²)·I₄
          [ −y2  −y3   y0   y1 ]
          [ −y3   y2  −y1   y0 ]
```

**(D-a) The Gram is forced.** `QQᵀ = I₄⊗M` gives
`4·rank(M) ≤ rank(Q) ≤ 4`, so `rank(M) ≤ 1`; with `M(0) = 4` forced
and `M` PSD, `det M = 0` gives `M(1) = ±4`. The only admissible
Grams are `M = 8I₂ − 4J₂` (the house form — the genuine `i = 2`
branch) and `M = 4J₂` (the degenerate branch: the `i = 1`
construction written with `i = 2` bookkeeping). This settles
`(s,i) = (1,2)` by rank alone, for every `w ≥ 2`; Theorem E
(§1.2.1) reaches the same two Grams — they are the twists of the
house form by the two real characters of `ℤ₂` — but needs
`w > 2s = 2`. The two agree wherever both apply, which is a useful
independent check on each.

**(D-a′) The degenerate branch collapses — and this is exactly what
that means.** In the branch `M = 4J₂`, (H2) reads `Σ PAF(t) = −4` for
**every** `t ≠ 0`: both tiers coalesce, and the profile is the `i = 1`
profile. The claim is stronger than that, and is about the border, not
only the seeds. Write `P₁[r][J] = P[r][2J]` and `Q₁[I] = Q[2I]`.

- ***`Q` pairs up.*** `(QQᵀ)[2I, 2I+1] = M(1) = +4` and the rows of `Q`
  are `±1` of length `4`, so `Q[2I+1] = Q[2I] = Q₁[I]`; the
  cross-superblock zeros then make `Q₁` a 4×4 Hadamard matrix. (This is
  (D-b) with the sign reversed: `+4` forces equality where `−4` forces
  negation.)
- ***The compressed blocks.*** Read Lemma 3's three cases on `Ḡ = ℤ₂`:
  a developed block of `Ĉ` contributes `δ_q`, a reflected or
  transposed-reflected one contributes `εδ_q` — which is the `ε` in
  `d`. So at every block position `(I,J)`, for all `c, c′`,

  ```
  Ĉ[(I,1),(J,c′)] − Ĉ[(I,0),(J,c′)] = −(−1)^{c′}·Λ(d)[I][J]
  Ĉ[(I,c),(J,0)] + Ĉ[(I,c),(J,1)]  =              Λ(r)[I][J]
  ```

  — the (D-d) table differenced, and the same table summed. Both sides
  are linear in the eight coset sums, so both are identities of `ℤ[σ]`
  (cert 10 checks them on an interpolation grid, for `ε = ±1`).
- ***`P` pairs up.*** By the first item `E Qᵀ` has `col(2I+1) = col(2I)`,
  so by (H4) `P Ĉᵀ` must too; by the second item that difference is
  `−Λ(d)·b_r` with `b_{rJ} = P[r][2J] − P[r][2J+1]`. The Fourier
  identity of D4's proof, at the character trivial on `K` and
  nontrivial on `G`, reads `Σ_q δ_q² = Σ̄(0) − Σ̄(1)`; the degenerate
  profile gives `Σ̄(0) = 4n − 4(w−1)` and `Σ̄(1) = −4w`, so
  `Σ_q δ_q² = 4n + 4 = N`, and `ε² = 1` gives `Σ_q d_q² = N` too. Hence
  `Λ(d)Λ(d)ᵀ = N·I₄` is nonsingular, `b_r = 0`, and
  **`P[r][2J+1] = P[r][2J]`** — again a theorem, not a choice. (Exact
  mirror of (D-c): there `Σ_q δ_q² = 4` while `Σ_q r_q² = N`, and it is
  `Λ(r)` that is inverted, killing the *sum* instead of the
  difference.)
- ***So the strips are constant on all of `G`.***
  `P̃[r,(J,h)] = P₁[r][J]` and `Q̃[(I,g),c] = Q₁[I][c]`, with no
  dependence on `h` or `g` — which *is* the `i = 1` ansatz, on
  `(G, K′ = G)`. The hypotheses transport termwise: (H2) is the `i = 1`
  profile; `Q₁Q₁ᵀ = 4I₄` is (H1) with `M₁ = (4)`; `P Pᵀ = 2·P₁P₁ᵀ`
  turns (H3) into `EEᵀ + n·P₁P₁ᵀ = N·I₄`; and (H4) read at `c = 0`,
  with the row-sum identity above, becomes `E Q₁ᵀ + P₁Λ(r)ᵀ = 0`, i.e.
  `E = −¼·P₁·Λ(r)ᵀ·Q₁` — the `i = 1` border equation, with the row-sum
  vector `r` where the genuine branch has `d`. (The `c = 1` component
  is the same equation: the two rows of a block of `Ĉ` have equal row
  sums.) Conversely, doubling any `i = 1` border gives (H1)–(H4) with
  `M = 4J₂`.

Doubling is therefore a **bijection** between the `(s,i) = (1,1)`
borders on `(G, ρ)` and the degenerate-branch `(s,i) = (1,2)` borders
on `(G, K, ρ)` for the same seeds, and it does not move the assembled
matrix: `H` and `H₁` are equal entry for entry. The degenerate branch
contains no seed quadruple, no border and no matrix that the classical
`s = 1, i = 1` construction does not already contain; `K` is invisible
in its data. ∎

**(D-b)** In the genuine branch, `M(1) = −4` forces
`Q[2I+1] = −Q[2I]`; writing `u_I = Q[2I]`, cross-superblock
orthogonality makes `U = (u_0;…;u_3)` a 4×4 Hadamard matrix — `Q` is
its ±-row-doubling.

**(D-c)** `EQᵀ` then has `col(2I+1) = −col(2I)`, so `PĈᵀ` must too;
summing the two rows of a developed block over `Ḡ = ℤ₂` gives the
seed's row sum, and the condition reads `Λ(r)·a_r = 0` with
`a_{rJ} = P[r][2J] + P[r][2J+1]`. D5 gives `Σ_q r_q² = N ≠ 0`, so
`Λ(r)` is nonsingular and `a_r = 0`: **`P[r][2J+1] = −P[r][2J]` is a
theorem, not a choice.** Writing `p[r][J] = P[r][2J]`, (H3) becomes
`EEᵀ + 2w·ppᵀ = N·I₄`; entries of `ppᵀ` are even and
`|EEᵀ|_off ≤ 4`, so for `w ≥ 2` the split is forced: `p` is 4×4
Hadamard and `E` is Hadamard. (Sharper than D3, which needs
`w > 2s`: the doubling turned `w` into `2w`.)

**(D-d)** With (D-b) and (D-c), (H4) collapses to a single 4×4
equation: writing `D_I[J] = Ĉ[2I][2J] − Ĉ[2I][2J+1]`, the block
table gives `D = Λ(d)` exactly, so

> **(H4) ⟺ `E Uᵀ + p Λ(d)ᵀ = 0` ⟺ `E = −¼ · p · Λ(d)ᵀ · U`.**

Parseval on the character non-trivial on `G` but trivial on `K`
forces `Σ_q δ_q² = 4`; given that, `EEᵀ = 4I` is automatic. The
whole border reduces to: *choose 4×4 Hadamard `p, U` such that
`−¼·p·Λ(d)ᵀ·U` has entries `±1`.*

**(D-e) Transport.** A fixed `(p, U, E)` pins `d` uniquely (`Λ` is
injective) — the exact analogue of the `i = 1` template pinning the
row-sum vector. The integral vectors `y ∈ ℤ⁴` with `Σ y_q² = 4`
comprise eight even vectors `±2e_j` and sixteen odd vectors
`(±1,±1,±1,±1)`, and the border equation can be evaluated on all
24; but **seed-derived data is always even** — each coset sum
satisfies `σ_q(c) ≡ w (mod 2)`, so `δ_q = σ_q(0) − σ_q(1) ∈ 2ℤ` —
so the realizable vectors are exactly `d = ±2e_j`. Every one of the
24 abstract arguments admits a border — **PROVEN-BY-CERTIFICATE**, by
the complete census of cert 10 over all `768²` ordered pairs of 4×4
Hadamard `(p,U)`: every `±2e_j` admits exactly `½` of the pairs
(`294 912`) and every odd vector exactly `¾` (`442 368`).
*Mechanism.* The four rows of a 4×4 Hadamard matrix lie in one of the
two weight-parity cosets of `{±1}⁴ ≅ F₂⁴` (the even-weight subgroup
`F₂³` and its complement), and so do its four columns; write
`W = p·Λ(d)ᵀ/2`, integral with `W Wᵀ = 4I₄` (from `Λ(d)ᵀΛ(d) = 4I₄`
and `p pᵀ = 4I₄`), so its rows have square norm 4. Either every row of
`W` is a sign vector, all of one class, and then `E` is `±1` iff that
class differs from the column class of `U` — exactly half the `U`; or
every row is a spike `±2e_k`, and then every `U` is admitted. For
`d = ±2e_j` only the first case occurs (`½`); for odd `d` the 768 `p`
split `384/384` between the two cases, giving `½·½ + ½·1 = ¾`. (The
spike case is the second half of the mechanism, and it is what
separates `¾` from `½`.) For genuine `i = 2` data only the `±2e_j`
part of that census is in play. Borders transport across the
doubling: with `P[r][2J+c] = (−1)^c·P₁[r][J]`, `Q[2I] = Q₁[I]`,
`Q[2I+1] = −Q₁[I]`, `E = E₁`, an `i = 1` border for seeds with
row-sum vector `r` is a valid `i = 2` border exactly when `d = r`,
since the transported system reads `P₁Λ(d)ᵀ = P₁Λ(r)ᵀ` and `Λ` is
injective with `P₁` invertible. The four `(s,i) = (1,1)` records
here all have `r = (2,0,0,0)`, so their transports land at
`d = (2,0,0,0)`.

**One sentence:** *at `s = 1` the `i = 1` and `i = 2` border systems
are the same 4×4 system `E Uᵀ + p Λ(·)ᵀ = 0`; only the argument
changes, from the row-sum vector `r` to the twisted vector `d`.*

Theorem D is stated for the standard orientation only, and its
sufficiency direction inherits from Theorem A. Clauses (D-b)/(D-c)
are machine-validated on the four `i = 2` instances of §2 (cert 02),
and all four clauses (D-a)–(D-d) on the `H(52)` gate instance
(cert 03); all four §2 instances have the house Gram, so none of them
exercises (D-a′) — whose content is precisely that its branch holds
nothing new to exercise. The theorem is *gated*: a new
`s = 1, i = 2` instance on the non-cyclic group `G = ℤ₂×ℤ₂×ℤ₃` (the
`ε = +1` branch, which the four instances do not exercise) was built
from scratch and verifies green — `H(52)`, cert 03. Label: PROVEN
(paper-grade) + MEASURED ((D-b)/(D-c) on the four `i = 2` instances,
cert 02; all four clauses on the gate instance, cert 03) +
PROVEN-BY-CERTIFICATE (the gate, cert 03; the `768²` census of (D-e)
and the compressed-block identities of (D-a′)/(D-d), cert 10).

### 1.6 Corollary (the index-2 collapse) — PROVEN

**Corollary.** Let `G` be abelian of even order with a **unique**
subgroup `K` of index 2 (equivalently, 2-rank 1), and `ψ` the
character with kernel `K`. Then `K` is characteristic, so
`ψ∘α = ψ` for every `α ∈ Aut(G)`; hence `x ↦ ψx` is a bijection
between the `s=1, i=1` and `s=1, i=2` seed problems that commutes
with every automorphism of `G` — in particular it preserves
invariance under any multiplier subgroup of `Aut(G)`. The two
problems are one problem, and any exhaustive statement proved on one
side transports to the other.

(With Theorem D this closes the `s = 1` layer completely: the
`i = 2` **seed problem** is a character-twist reparametrization of
the `i = 1` seed problem whenever the index-2 subgroup is unique,
and its border system is the `i = 1` system in different
coordinates. The bijection is of seed problems, not of matrices —
it does **not** imply the assembled Hadamard matrices are
equivalent, and at order 668 §3.4 proves they are not. The underlying Parseval mechanism is that the
trivial character's constraint is tight at `i = 1` while at `i = 2`
the binding constraint is the character trivial on `K` but not on
`G`, which contributes `Σ_q δ_q² = 4n − 4(s+1)w + 4s`.)

### 1.7 Theorem E′ (Gram rigidity under `w > s`) and the boundary — PROVEN

Theorem E (§1.2.1) pins the admissible Gram to `M = 4i·P_S` under
`w > 2s`. The lock closes twice as early once one more integrality
fact is used.

**Lemma (mod 4).** For every `x : G → {±1}` and every `t ≠ 0`,
`PAF_x(t) ≡ n (mod 4)`. Consequently, under (H2), `M(c) ∈ 4ℤ` for
every `c ∈ Ḡ`.

*Proof.* `PAF_x(t) = n − 2d(t)` with `d(t) = #{u : x(u) ≠ x(u+t)}`.
The permutation `u ↦ u+t` of `G` splits into cycles, and along a
cycle the number of sign changes of `x` is even — the product of the
consecutive products `x(u)x(u+t)` around the cycle is `1`. So `d(t)`
is even. Summing over the quadruple, `Σ_q PAF_q(t) ≡ 4n ≡ 0 (mod 4)`,
and (H2) gives `M(κt) ∈ 4ℤ`; every `c ≠ 0` is some `κ(t)`, and
`M(0) = 4s`. ∎

**Theorem E′.** Assume (H1), (H2), `s ≥ 1` and **`w > s`**. Then
`S := {χ ∈ Ĝ̄ : M̂(χ) ≠ 0}` has `|S| = s`, is closed under
conjugation and under the Galois action `χ ↦ χ^k`, and

```
       M(c)  =  4 · Σ_{χ ∈ S} χ(c)          (c ∈ Ḡ),
```

i.e. `M = 4i·P_S`, with spectrum `{0^(i−s), (4i)^s}`. In particular
every conclusion of Theorem E and of Corollaries E1 and E2 holds
under `w > s`.

*Proof.* Steps (L1), (L2) and (L4) of §1.2.1 are unchanged; the
Parseval window `M̂(χ) ≤ N/w = 4i + 4s/w` of (L1) is what carries
the weight. Two steps used `w > 2s` and are re-run with `w ≥ s+1`.
(i) `s ≤ i`: otherwise `rank M ≤ i` and the trace bound
`4si = trace M ≤ i·(4i + 4s/w)` gives `s(1 − 1/w) ≤ i`, while
`s(1 − 1/w) ≥ s²/(s+1) > s − 1`; so `s − 1 < i`. (ii) `rank M = s`:
with `r = rank M`, the same bound `4si ≤ r·(4i + 4s/w)` gives
`r ≥ sn/(n+s) > s − 1`, because `n = wi ≥ (s+1)s > s(s−1)`; so
`r = s`, `|S| = s`, and (L4) makes `S` Galois-stable. Now set
`C := (M_S − M)/4` — note the `4`, where (L5) could only afford a
`2`. `M_S(c) ∈ 4ℤ` exactly as in (L5), and `M(c) ∈ 4ℤ` by the mod-4
lemma, so `C` is an **integer** matrix; it has zero diagonal, since
`M_S(0) = 4s = M(0)`; and its eigenvalues are `0` off `S` and
`(4i − M̂(χ))/4 ≥ −s/w > −1` on `S`, by the Parseval window. The
interlacing step that closes (L5) — an integer symmetric matrix with
zero diagonal and a nonzero off-diagonal entry has `λ_min ≤ −1` —
now gives `C = 0`, so `M = M_S = 4i·P_S`. ∎

*(That is the whole gain: the mod-4 lemma lets the integrality lock
divide by 4 instead of 2, which halves the width the strict
inequality `λ_min(C) > −1` needs.)*

**Theorem E′ (boundary).** Let `w = s`, and let `M` satisfy the
conditions that are necessary for any admissible Gram — `M(0) = 4s`,
`M ∈ 4ℤ`, `|M(c)| ≤ 4s`, `M` PSD, `M̂(χ) ≤ N/w`, `rank M ≤ s`. Then
exactly one of:

- **(a)** `s < i`, and `M = 4i·P_S` — rigid, as above;
- **(b)** `s = i`, and `M = 4i·I − 4C` where
  `C(c) = ξ(c)·[c ∈ Ḡ₀∖{0}]` for a subgroup `Ḡ₀ ≤ Ḡ` and a real
  character `ξ` of `Ḡ₀` (`Ḡ₀ = 0` is the rigid case, and every such
  `M` does satisfy all the necessary conditions);
- **(c)** `s = i+1`, and `M = 4s·I_i`;
- **(d)** `s ≥ i+2` cannot occur.

*Proof.* The trace bound gives `s ≤ i+1`. For `s ≤ i` the argument
above runs as far as `λ_min(C) ≥ −s/w = −1`, where the interlacing
step is now silent; then `C + I` is PSD with unit
diagonal and integer entries — the Gram of unit vectors with inner
products in `{0, ±1}` — so `±`-parallelism is an equivalence relation
on `Ḡ`, its class of `0` is a subgroup `Ḡ₀` (by `Ḡ`-invariance), and
`C(c) = ξ(c)` on `Ḡ₀∖0` with `ξ` multiplicative into `{±1}`. Then
`Ĉ(χ) = |Ḡ₀|·[χ|_{Ḡ₀} = ξ] − 1`, which is never `0` when `Ḡ₀ ≠ 0`;
since `Ĉ = 0` off `S`, a nonzero `C` forces `S = Ĝ̄`, i.e. `s = i` —
that is (a) and (b). Conversely each `M` in (b) has
`M̂(χ) = 4i + 4 − 4|Ḡ₀|·[χ|_{Ḡ₀} = ξ] ∈ [4, 4i+4] = [4, N/w]`. For
`s = i+1` all `i` eigenvalues sit at the top of the window, so
`M = 4s·I`. ∎

**Sharpness — PROVEN-BY-CERTIFICATE (cert 16).** Case (b) is
realised *by seeds*. On `G = ℤ₃²` with `K = {(0,b)}` and
`(s,i,w) = (3,3,3)`, the quadruple

```
++-++-++- ,  ++-++-+-+ ,  ++-+-+-++ ,  ++--+++-+     (index (a,b) ↦ 3a+b)
```

has `Σ_q PAF_q = −12` on `K∖0` and `+4` off `K` — that is (H2) with
`M = 16I₃ − 4J₃`, spectrum `{4,16,16}`, **not** a projector multiple
— and a `12×12` `±1` matrix `Q` with `QQᵀ = I₄⊗M` exists. Cert 16
rebuilds all of this from nothing: the seeds by a meet-in-the-middle
over **all** 512 sequences per seed, with no prefilter, which finds
**490 212** ordered solutions on `ℤ₃²` and **none** on `ℤ₉`; the `Q`
by depth-first search, re-multiplied entry by entry. **So `w > s`
cannot be weakened to `w ≥ s`.**

**Exhaustion — PROVEN-BY-CERTIFICATE (cert 16).** The boundary
classification is exhausted at small parameters by **two exact routes
that share no code** — a `Fraction`-LDLᵀ sieve over every
`Ḡ`-invariant `M : Ḡ → 4ℤ` obeying the necessary conditions, and the
explicit construction of the list (a)/(b)/(c) from Ramanujan sums —
which are required to produce the same set at every cell. The default
run covers **22 boundary cells** (`w = s`), over every abelian `Ḡ` of
order `i`, from 45 `(s,i,w)` triples with `s, i ≤ 6`; `--wide` adds
`(7,8,7)` on `ℤ₈`, making 23. The same sweep confirms rigidity at
**26** cells with `w > s` — **20** of them in the range
`s < w ≤ 2s` that Theorem E did not reach — and records **11** cells
below the boundary, where non-projector Grams pass every necessary
condition.

**What this changes here.** The `(1,2,2)`, `(1,4,2)`, `(3,4,6)` and
`(5,6,10)` cells lie in the range `s < w ≤ 2s` and are now covered by
the theorem rather than by measurement; cert 05's exhaustive
`(1,2,2)` result is a corollary of it. Nothing changes at `v = 523`,
where all four live cells have `w > 2s` and Theorem E already
applied. Below the boundary (`w < s`) no rigidity theorem of this
type exists — cert 16 exhibits the non-projector survivors there —
and none is claimed. Nothing is said about (H3) or (H4) at
`(3,3,3)`: there `w = 3 ≤ 2s`, so D3 does not apply and the corner
need not be Hadamard.

### 1.8 Theorem F (border kits) — PROVEN; the `(3,4)` census — PROVEN-BY-CERTIFICATE

Setting: the house form at `i = s+1` (so `S = Ĝ̄∖{1}` and
`M = (4s+4)I − 4J`), `w > 2s`, and a complete instance `(x, E, P, Q)`
on `(G, K, ρ)`. Write `V_S ⊂ ℝ^i` for the span of the nontrivial
characters — the sum-zero vectors — and `Σ̄ = N·I − w·M`.

**Theorem F (structure).**

- **(a)** `QᵀQ = 4i·I_{4s}`: the `4s` columns of `Q` are mutually
  orthogonal and lie blockwise in `V_S` (each `i`-block of each
  column sums to `0`).
- **(b)** `P(I₄ ⊗ P_{S^c}) = 0` and `PᵀP = QQᵀ = I₄⊗M`: the rows of
  `P` are likewise orthogonal and blockwise in `V_S`. In particular
  `Pᵀ` is itself an admissible column table (*dual tables*), and the
  four block-constant rows `h_J ⊗ 1_i`, `h ∈ H(4)`, complete `P` to
  a Hadamard matrix of order `4i` (*forced completion*).
- **(c)** `ĈᵀĈ = ĈĈᵀ = I₄⊗Σ̄`, and `Ĉ` commutes with `I₄ ⊗ P_S`.
- **(d)** `E = −(1/4i)·P·Ĉᵀ·Q`.
- **(e)** Conversely, given any `P, Q` with (a) and (b) and any table
  `σ` whose compressed profile is `Σ̄`, the matrix
  `E := −(1/4i)PĈᵀQ` satisfies `EEᵀ = 4s·I` **automatically**, and
  (H3) then reads `4s + 4iw = N`, an identity. So **(H3) and (H4)
  together reduce to the single condition `E ∈ {±1}^{4s×4s}`.**
- **(f)** (H4) depends on `σ` only through the spectrum
  `(σ̂_q(χ))_{q, χ∈S}` — the **S-part**. Its kernel in `σ` is
  `{σ : σ̂_q(χ) = 0 ∀χ ∈ S}`, of dimension `4(i−s) = 4`: the per-seed
  constants.
- **(g)** (Transport.) A kit `(E,P,Q)` valid for `(G, K, ρ, x)` is
  valid for every `(G′, K′, ρ′, x′)` with `G′/K′ ≅ G/K`,
  `κ′(ρ′) = κ(ρ)` and the same S-part — **at any `w′`**, hence at any
  order `N′ = 4(w′i + s)`. Only the parity `σ′_q(c) ≡ w′ (mod 2)` and
  `Σ_q r′_q² = N′` constrain the row sums.

*Proof.* (a) `QQᵀ = I₄⊗M` has nonzero spectrum `{(4i)^{4s}}` by
Theorem E, and `rank Q = 4s`; `QᵀQ` is `4s×4s` with the same nonzero
spectrum. The block balance is `Q_Iᵀ1 = 0`, from `1ᵀM1 = i·M̂(1) = 0`.
(d) Right-multiply (H4) by `Q` and use (a). (c) `Ĉᵀ` is the standard
array of the seeds `(σ₀∘(−), −σ₁, −σ₂, −σ₃)` at the same `κ(ρ)`, so
Lemma 2 gives `ĈᵀĈ = I₄ ⊗ dev(Σ_q PAF_{σ_q}) = I₄⊗Σ̄`; developed
blocks lie in the commutative group algebra and `R v_χ = χ(κρ)·v_χ̄`
with `S` conjugation-closed, so `Ĉ` commutes with `I₄⊗P_S`.
(e) and (b): `EEᵀ = (1/16i²)·PĈᵀ(I₄⊗M)ĈPᵀ =
(1/4i)·P(I₄⊗P_S)ĈᵀĈPᵀ = (1/4i)·P(I₄⊗P_S)(I₄⊗Σ̄)Pᵀ`, and `Σ̄` acts on
`V_S` as `N − 4iw = 4s`; hence `EEᵀ = (s/i)·P(I₄⊗P_S)Pᵀ`. Forward,
(H3) with D3 gives `EEᵀ = 4sI` and `PPᵀ = 4iI`, so
`P(I₄⊗P_S)Pᵀ = PPᵀ` and `P(I₄⊗P_{S^c}) = 0` (both sides PSD), whence
`PᵀP = 4i·Π_{V_S} = I₄⊗M`. Conversely, with (b),
`EEᵀ = (s/i)·4i·I = 4sI`. Forced completion: the vectors `h_J ⊗ 1_i`
are orthogonal to each other and to every row of `P`, of squared
norm `4i`. (f) `Ĉ` is linear in `σ`, and the `V_{S^c}`-component of
`σ` yields blocks whose rows and columns are `S^c`-vectors,
annihilated by `P`. (g) (a), (b) and (H1) do not mention `w`; (H3) is
the identity of (e); (H4) sees `σ̂|_S` only. ∎

*(At `(s,i) = (1,1)`, `S = {1}` is all of `Ĝ̄`, the kernel in (f) is
`0`-dimensional, and the row sums enter (H4) directly — Theorem D's
(D-e). Theorem F holds there with (f) read as vacuous.)*

**MEASURED (cert 17 [D]).** Clauses (a)–(f) are re-multiplied entry
by entry on the **seven** banked coset-border records with `s ≥ 1`
and `w > 2s` — 668, 716, 1676, 1772 at `(1,1)`; 1916 at `(3,4)`;
1388 at `(5,6)`; 1436 at `(7,8)` — **7 / 7**, together with the
superblock balance of `P` and the (H4)-kernel at `i = s+1`.

**Theorem F (the `(3,4)` census).** Let `(s,i) = (3,4)`,
`Ḡ ∈ {ℤ₄, ℤ₂²}`, `w` even, and assume (H1),(H2) with the house Gram.
Then:

- **(i)** exactly one seed is coset-balanced (`σ_{q*}(c) = r_{q*}/4`,
  with `r_{q*} ≡ 0 (mod 8)` at `w = 130`) and the other three have
  `|σ̂_q(χ)| = 2` at every nontrivial `χ` — the **silent-seed lemma**;
- **(ii)** there are exactly **2048** admissible S-parts per group;
- **(iii)** **every** S-part admits a kit `(E,P,Q)` at **every**
  `κ(ρ)`;
- **(iv)** at `N = 2092` the coset-sum tables form a shell of
  `4192 × 64 = 2048 × 131 = 268 288` tables per group.

*Proof.* (i),(ii): Corollary E2 gives `Σ_q |σ̂_q(χ)|² = 4s = 12` at
each nontrivial `χ`, with `σ̂_q(χ) ∈ 2ℤ` (or `2ℤ[i]`) because each
coset has `w` elements; integrality and parity of `σ = (r + T)/4`
force `T_q(c) (mod 8)` to be constant in `c`; the finite enumeration
returns 2048 per group, all of the silent-seed shape (cert 17 [A],
brute force over the spectra). (iii) is the exhaustive census
(cert 17 [B]). (iv) is Jacobi's `r₄(523) = 4192`, with the mod-8
shape forced by `523 ≡ 3 (mod 8)` — one coordinate `≡ 0 (8)`, three
`≡ 2 (4)` — each shell vector compatible with exactly 64 S-parts and
each S-part with exactly 131 shell vectors. ∎

**Corollary (the `(3,4)` cell is one-layer) — PROVEN.** Let `G` be
abelian of order `4w` with `w > 6`, let `K ≤ G` have index 4, let
`ρ ∈ G`, and let `x₀,…,x₃ : G → {±1}` satisfy
`Σ_q PAF_q = 4n·δ₀ − 12·[K∖0] + 4·[∉K]`. Then a Hadamard matrix of
order `16w + 12` with the bordered-GS layout exists, and its border
is computable from the 16 coset sums. At `N = 2092` (`w = 130`):
**every** house-profile quadruple on `ℤ₈×ℤ₅×ℤ₁₃`, `ℤ₂×ℤ₄×ℤ₆₅` or
`ℤ₂³×ℤ₆₅` with an index-4 subgroup **is** an `H(2092)` — the border
is never the obstruction.

**The census, and what it costs to believe.** Cert 17's engine is a
second implementation: it fixes `Q` from a list of 301 candidate
12-cliques, filters the 648 block-balanced sign representatives of
`P`-rows against `V = (4Ĉ⁰)ᵀQ`, and cliques among the survivors —
where the source laboratory's engine backtracks over `P` first and
looks `Q` up. The default path censuses a fixed deterministic
256-class sample (**256 / 256**), re-verifying each kit exactly at
`w = 130` with a full `σ`; `--full` censuses all
`2 × 4 × 2048 = 16 384` classes, and has been run in this
repository — 2048 / 2048 on each of the eight `(group, κ(ρ))`
sweeps, **16 384 / 16 384**, at the pinned digest, which
`certs/17-border-kits-34/NOTES.md` records. Controls, all through
`verify/verify.py`: four from-scratch `H(44)` at `w = 2`, whose seed
quadruples are found by a full meet-in-the-middle over all 256
sequences on `ℤ₈` and `ℤ₂×ℤ₄`
(and on `ℤ₂³` with `|K| = 2`, which carries none); and `H(124)` at
`w = 7` twice — once with the order-1916 record's kit transported
verbatim, which is Theorem F(g) exercised across `w`, and once with
this engine's own kit.

**What is not claimed.** Nothing about the **seed layer**: no
house-profile quadruple on a group of order 520 is known, and at 2092
the corollary is an implication with an unmet hypothesis. Nothing
about the *number* or *structure* of kits per S-part — existence
only. Nothing at `w ≤ 2s` beyond the `H(44)` controls, which show
that the same kit shape *works* there, not that every kit there has
this shape.

---

## 2. Movement II — the instances

Everything in this section is replayed in the strict sense — nothing
here is audited from a bank. Each cert directory under `certs/`
rebuilds its matrices from the small banked data in `data/`,
re-checks the theorem hypotheses in exact stdlib integer
arithmetic, hands every matrix to `verify/verify.py` (the trust
chain), and compares the canonical SHA-256 against the digest pinned
in its `NOTES.md`. Large matrices are regenerated, never committed.
(Movement III is where the word needs care: the exact 4-profiles
there are **audited** on the default path and **recomputed** only
under `--full`. §3.4 says which is which.)

### 2.1 The twelve public records — PROVEN-BY-CERTIFICATE (full replay)

Seed data encoding twelve Hadamard matrices at orders 668, 716,
892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948 and 1964 was
posted publicly on 2026-08-12 (the sign-stream); expanded matrices
for ten of the twelve appeared in a third-party GitHub repository
created and pushed the same day, 2026-08-12 UTC (see
`PROVENANCE.md`). The parameter records banked in
`data/payload-records.json` — group, seeds, reflection, coset
divisors, and border tables — were **decoded here from those public
artifacts**; they are the posting team's mathematical content, not
ours, and no priority of any kind is claimed on the records
themselves, on the decode, or on existence at these twelve orders.
What this repository adds at these
orders is the verification: cert 01 checks every hypothesis
(H1)–(H4) and every derived law of Theorem C literally on each
record, assembles the matrix, verifies it, and matches the canonical
digest — **12/12 green**. The records saturate the classification of
§1.3: five at `(0,1)`, four at `(1,1)`, and one each at `(3,4)`
(order 1916), `(5,6)` (order 1388), `(7,8)` (order 1436) — corners
of sizes 12, 20, 28 respectively, as D3 forces.

### 2.2 The matrices constructed here at `s ≤ 1` (certs 02–05)

| cert | artifact(s) | what it instantiates |
| --- | --- | --- |
| 02 | `H(668)`, `H(716)`, `H(1676)`, `H(1772)` in the `i = 2` frame | Lemma T (§1.4): the four `s=1` records twisted and re-bordered at `i = 2`; each verifies green and is a different artifact from the decoded one. At all four orders the stronger statement is the theorem of §3.4, §3.6, §3.7 and §3.8: the twisted matrix is a different equivalence class, not merely a different file. |
| 03 | `H(52)` on `G = ℤ₂×ℤ₂×ℤ₃` | Theorem D's gate: a from-scratch `s=1, i=2` instance on a non-cyclic group, `ε = +1` branch. |
| 04 | `H(76)` on a non-cyclic group | a bordered-GS instance on non-cyclic `G` whose seed quadruple is invariant under a **non-scalar multiplier subgroup `Γ ≤ Aut(G)`** (order 2, 12 orbits, re-derived from scratch in the cert) — built as a search gate. Its Theorem-A Gram is the house form; the lab record's name "M" for `Γ` is a letter collision with the Gram and is not used here. |
| 05 | two `H(20)` instances | the hypothesis boundary `w = 2s`: D3's *hypothesis* fails, yet its *conclusion* holds — both corners are Hadamard, and the cert proves exhaustively that **no** non-Hadamard corner exists at `(s,i,w) = (1,2,2)`, so `w > 2s` is sufficient but not necessary there. The pair also exercises an arbitrary subgroup `K` (the diagonal of `ℤ₂×ℤ₂`, not a coordinate kernel) and both array orientations. |

No novelty of existence is claimed for any order in this table —
every order here is long settled. The artifacts instantiate the
theorems; that is their entire role. (The four `i = 2` matrices of
cert 02 can also be produced by Theorem D's transport (D-e); the two
assemblies agree byte for byte.)

### 2.3 The even-`s` branch: four `(2,4)` matrices — PROVEN-BY-CERTIFICATE

Theorem E's **general branch** (§1.2.1) allows `Ḡ = ℤ₄` with
`S = {χ, χ³}`, the two faithful characters, so `s = 2` and

```
M  =  4·(χ + χ³)  =  (8, 0, −8, 0).
```

`s` is **even**, so this cell lies outside the house branch that
Theorem C classifies, and outside every cell instantiated in §2.1 and
§2.2. (H2) reads
`Σ_q PAF_q = 4n·δ₀ − 8·[K∖0] + 0·[odd cosets] + 8·[κ⁻¹(2)]`; seeds
exist in quantity (exhaustively at `n = 12`, by local search at
`n = 20`).

Cert 18 rebuilds four instances from
`data/cell24-records.json` — seeds, `Q′`, `P′`, `ρ` — recomputing the
coset sums, the S-part and the corner `E = −(1/16)PĈᵀQ` from the
seeds rather than reading them, checking (H1)–(H4) and the Σ̄-law
literally, assembling through its own block-explicit assembler, and
handing each matrix to `verify/verify.py`:

| `n` | `w` | `ρ` | order | canonical SHA-256 |
| --- | --- | --- | --- | --- |
| 20 | 5 | 0 | **88** | `942b3f32fcd75e72a64f92d9c294b0d0cedbbd0965fe5a14213c30b8b66ffc8a` |
| 20 | 5 | 1 | **88** | `4cae47d1c5054a86ca48154c5e9cd99845294be4d275dd8560bf6b05fe5f08e7` |
| 12 | 3 | 0 | **56** | `ad67ee2c9d1f4d0343b250824dec301759e097d72d949f1b9c91f22cab026b85` |
| 12 | 3 | 3 | **56** | `fa82808c8dbc0245f3d427312975183fb947438bdcafba56107ee2280d6e4aff` |

**Two structural claims, each a finite exact check.** (i) None of the
four is Hadamard-equivalent to `H(2) ⊗ H′` — the only Kronecker
factorisation an order 88 or 56 admits. If `H = D₁Π₁(H₂⊗H′)Π₂D₂`
then the pointwise product of two rows of `H` is `±` a column
permutation of a pointwise product of two rows of `H₂⊗H′`, and in
`H₂⊗H′` the rows `[h|h]` and `[h|−h]` pair off into `N/2` disjoint
pairs with the *same* pointwise product; so anything equivalent to a
Kronecker product has a sign-normalised pointwise-product class of
size `≥ N/2`. These four have largest class **4**, against `N/2 = 44`
and `28`. Controls in the same run: Sylvester `H(16)` gives 8 (the
test passes a genuine Kronecker product), a scrambled equivalent copy
of it gives 8 again (the statistic is an invariant), and Paley
`H(12)` gives 1. (ii) The common kernel of `S` in `Ḡ = ℤ₄` is `{0}`,
so `S` generates the dual and the instance is not a collapse of a
smaller-index construction.

These are the first coset-border Hadamard matrices with even `s`
built anywhere in this repository, and they are not disguised
`s ≤ 1` objects. **No novelty of existence is claimed at 56 or 88** —
those orders were settled long ago. Nor is any *literature* claim
made for the even-`s` construction type: §4.2's single hedged
sentence bounds exactly the sources §4.1 enumerates, no search was
added for it, and none is implied here.

### 2.4 The `(2,4)` border, and Theorem 3 — PROVEN + PROVEN-BY-CERTIFICATE

**Proposition (the `(2,4)` border).** Let `(s,i) = (2,4)`, `Ḡ = ℤ₄`,
`M = (8,0,−8,0)`, and let `x` satisfy (H2).

- **(a)** (H1) holds **iff** `Q[4I+c+2] = −Q[4I+c]` and the eight
  rows `Q[4I+c]`, `c ∈ {0,1}`, form a Hadamard matrix `Q′ ∈ H(8)`.
- **(b)** If `w > 2s = 4`, then in every complete instance the row
  table is anti-periodic too — `P[r][4J+c+2] = −P[r][4J+c]`, with
  `P′ ∈ H(8)` — and `E = −(1/16)PĈᵀQ ∈ H(8)`.
- **(c)** For **every** admissible S-part
  (`τ_q = (σ_q(0)−σ_q(2), σ_q(1)−σ_q(3))` with even entries and
  `Σ_q |τ_q|² = 4s = 8` — there are 112 of them), **every** `κ(ρ)`
  and **every** `Q′ ∈ H(8)`, an anti-periodic `P` with `E ∈ {±1}`
  exists. **So the border is never the obstruction at `(2,4)`, at
  any `w`.**

*Proof.* (a) `M(2) = −8 = −4s` makes rows `c` and `c+2` antipodal;
what remains of the Gram conditions on the eight rows `c ∈ {0,1}` is
orthogonality. (b) Theorem F(b) with `S^c = {1, χ²}`: the rows of `P`
lie blockwise in `V_S = {p₀+p₂ = 0, p₁+p₃ = 0}`, and `PPᵀ = 16I` is
`P′P′ᵀ = 8I`; `E` is then Theorem F(d),(e). (c) is a census: the
right-orbits of `H(8)` under signed column permutations number
`151 200·2¹⁵/(8!·2⁸) = 480`, and `112 × 4 = 448` classes of
`(S-part, κ(ρ))` give `215 040` classes in all; cert 18 [4] exhibits
a kit for **215 040 / 215 040**. Every kit is `w`-free, since (H3)
reads `8 + 16w = 4(4w+2)`, an identity in `w`. ∎

*(Repair recorded, and it matters: (b) is a **necessity** statement
and needs `w > 4`. The two `H(56)` instances of §2.3 have `w = 3`,
where the anti-periodic kit is what they use and is sufficient — but
is not shown forced.)*

Cert 18's census is a second implementation of the border layer at
this cell. It builds the 480 orbit representatives from the 30 affine
`AG(3,2)` structures on eight labels times 16 plane-sign classes,
cross-checks the count against an independent backtracking count of
normalised `H(8)` (**151 200**), and checks that 3 000
deterministically generated labelled `H(8)` all land in the
representative set. The admissible anti-periodic rows for a given
`Q′` are characterised exactly — `a(p) = 2·v·Q′ᵀ` for some
`v ∈ {±1}⁸`, and then `E`'s row is `−v` — which is what makes an
exhaustive 215 040-class census cost seconds of standard-library
Python. Two hundred of the kits are then multiplied out in exact
integers with a **full `σ`** (an arbitrary `S^c` part added) at
`w = 5` and at `w = 130`.

**Theorem 3.** Assume (H1), (H2), `s ≥ 1`, and `w > s` — so
`M = 4i·P_S` by Theorem E′ (§1.7). If `S` contains **no real
character**, then `4 | i`, `s` is even, and `N ≡ 0 (mod 8)`.
Consequently, at `N ≡ 4 (mod 8)` — in particular at `N = 2092` —
every realisable Gram has a live real character.

*Proof.* Complex conjugation lies in the Galois group and fixes no
element of `S`, so it pairs `S` off: `s` is even. `1 ∉ S` gives
`M̂(1) = 0`, so `|Q_Iᵀ1|² = 1ᵀM1 = i·M̂(1) = 0`: every column of the
block `Q_I` sums to zero over its `i` rows, so `i` is even. A group
`Ḡ` of even order has a nontrivial real character `χ₀`, also `∉ S`,
so `Q_Iᵀv_{χ₀} = 0` likewise; adding, `Q_Iᵀ(1 + v_{χ₀}) = 0` says
every column of `Q_I` sums to zero over `ker χ₀`, a subgroup of order
`i/2` — so `i/2` is even. Finally `N = 4(wi + s)` with `4 | i` and
`s` even is `≡ 0 (mod 8)`. ∎

**Small-case exhaustion (cert 18 [3]).** For every abelian `Ḡ` of
order `i ≤ 16` and every Galois-stable `S` with no real character,
`|S|` is even in every case; and where `4 ∤ i` and `s ≤ 2` — the
budget for an exhaustive depth-first search over `{±1}^{i×4s}` — no
block `Q₀` with `Q₀Q₀ᵀ = M_S` exists, in **9 / 9** cases
(`i = 3, 6, 6, 9, 9, 9, 9, 9, 15`). At `i = 4`, where `4 | i` and the
theorem forbids nothing, the same search finds one — the `(2,4)` cell
of §2.3.

**Remark (the general branch at 2092).** Theorem 3 does *not* say
that at `N ≡ 4 (mod 8)` the general branch collapses onto the house
branch. Cert 12 already exhibits general-branch cells outside the
house identity at other orders, and the source laboratory records
`(H1)`-realisable Grams at `N = 2092` itself — `(3,8)` on `ℤ₂³` and
`ℤ₄×ℤ₂`, `(7,12)` on `ℤ₁₂`, `(11,16)` on `ℤ₁₆` — each with a
Galois-stable `S` containing a real character and `w > 2s`. Those are
not re-derived in this repository, so nothing is claimed for them
here beyond the observation that Theorem 3 does not exclude them;
their seed layers are untested in any case.

---

## 3. Movement III — existence plus separation

**Definition.** `H ~ H'` iff `H' = D_r P_r H P_c D_c` with `P`
permutation matrices and `D` diagonal `±1`. Transpose is **not** in
the group; the invariants below are computed on both sides (rows and
columns) wherever a transpose-extended statement is made, and where
transpose changes a verdict it is called out. The
**transpose-extended** relation is `A ≈ B` iff `A ~ B` or `A ~ Bᵀ`;
refuting it takes two refutations, since `A ~ Bᵀ ⟺ Aᵀ ~ B`.
Transpose-extended statements are now made for all four classes at
668 (§3.4, certs 08 and 15), all three at 716 (§3.6, cert 15) and
all three at 1676 (§3.7, cert 21). Orders 2060 (§3.5) and 1772
(§3.8) are the two that remain **row-side only**, and they say so —
their three-class statements included: at each, the transposed
profiles are a separate leg of the same campaign (at 1772 the `blas`
legs are in and the `bits` legs still running). No
separation statement of any kind is made about the four `(2,4)`
matrices of §2.3: they are existence witnesses for the even-`s`
branch, and beyond cert 18's not-Kronecker test no invariant of
theirs has been computed — not against each other, and not against
the known `H(56)` and `H(88)`.

### 3.1 The invariants

For rows `i ≠ j` let `u_ij[c] = H[i][c]·H[j][c]`, and for a 4-subset
`{i,j,k,l}` let `T4 = Σ_c H[i][c]H[j][c]H[k][c]H[l][c] = ⟨u_ij, u_kl⟩`
(the 4-profile statistic).

- **I1/I2 — row-pair (column-pair) collision profiles**: the multiset
  of block sizes of the partition of row pairs by `u_ij = ±u_kl`.
  Invariant because column negation fixes each `u` exactly, row
  negation flips a global sign, and permutations relabel. A block of
  size `> 1` is a disjoint pair with `|T4| = n`, so I1 is the extreme
  stratum of the 4-profile at `O(n²)` cost.
- **I3 — difference-code dimension.** Over `F₂`, with `r_i` the rows,
  `V = span{r_i + r_j}` and `W = V + ⟨1⟩`: **only `dim W` is
  invariant.**
  > **Trap 1 (documented, machine-demonstrated).** `dim V` is *not*
  > an invariant — row negation moves it (Sylvester `H(4)`: 2 → 3;
  > `H(64)`: 6 → 7 under a random signed permutation). On the 668
  > pair `dim V` *differs* (666 vs 667) and the difference is
  > worthless; `dim W = 667` on both. Cert 06's controls demonstrate
  > the jump so the trap stays closed.
- **I4 — extreme strata, exact.** `|T4|` is invariant; pigeonhole
  blocking enumerates the strata `|T4| ≥ n − 2(B−1)` exhaustively.
- **I5 — the 4-profile**: the histogram of `|T4|` over all `C(n,4)`
  4-subsets. Exact mode is a proven invariant computed exactly;
  sampled mode is COMPUTATIONAL-EVIDENCE only and **can never yield
  INEQUIVALENT**.
  > **Trap 2 (documented, measured).** The *signed* `T4` histogram is
  > **not** an invariant — only `|T4|` is. Measured on a provably
  > equivalent pair at order 1916, the signed histogram manufactures
  > a false separation at `d = 7` (20 bins with `d > 4`; `d` is the
  > raw standardized delta of §3.5); folded, the same data reads
  > exactly zero. The signed statistic is uncalibrated in both
  > directions (it also reads `d = 36.6` on the 668 pair where the
  > folded truth is `d = 1.8`). Every histogram in this repository
  > is folded at measurement.
- **I6 — `p`-rank**, `p` an odd prime dividing `n` (all four factors
  of an equivalence are invertible over `F_p`).

**Lemma (second moment).** For any Hadamard matrix of order `n`,
`Σ_{4-subsets} T4² = n³(n−1)(n−2)/24`.

*Proof.* `Σ_{(i,j,k,l) ∈ [n]⁴} T4² = Σ_{c,c'} (Σ_i H[i][c]H[i][c'])⁴
= n⁵` by column orthogonality. Degenerate tuples contribute
`n·n² + 3n(n−1)·n²`; the remaining `24·Σ_{4-subsets}` gives
`24·Σ T4² = n⁵ − 3n⁴ + 2n³`. ∎

The lemma cuts both ways: the second moment can never separate
anything (it is a function of `n` alone), and it is a sharp
correctness check on any computed profile. Every exact-profile run
in cert 06 asserts it programmatically.

### 3.2 Why the second moment is not an accident — the cheap-invariant wall

Any invariant statistic must be even in the row signs and even in
the column signs, and degree-2 statistics die by orthogonality; the
second-moment lemma shows why one natural aggregate of the
4-profile is forced by `n` alone. No universal minimality claim is
made for the full histogram — cheaper invariants can and do resolve
*other* pairs. The point is pair-specific: on the 668 pair below,
*every* cheap invariant — I1–I4, I6, `dim W`, the dual
weight enumerator — returns identical values (the source
laboratory's measurements; cert 06 recomputes the `F₂`-rank pair
here and audits the banked exact profiles, which `--full`
recomputes), and `2·10⁷` samples of I5 cannot tell the
pair apart (`max d = 1.8`, the raw standardized delta of §3.5); only
the exact 4-profile separates. Matching invariants are never evidence of
equivalence, and this pair is the proof.

### 3.3 Order 1916: the twist that is not a new matrix — PROVEN-BY-CERTIFICATE

The character-twisted sibling of the `(3,4)` record at order 1916
(the non-house-Gram witness of §1.2) satisfies
`H_twisted = S · H_house · S` with
`S = diag(I_{4s}, I_4 ⊗ diag(ψ(g)))` — 952 rows and the same 952
columns negated, no permutation. All `1916²` cells are
machine-checked — and not just for the banked character: cert 09
sweeps the **entire character group** of `Ḡ = ℤ₂×ℤ₂` at this order
(where `ρ = 0`, so `ψ(ρ) = 1` for every `ψ`), and all three
nontrivial twists produce non-house Grams with the *same* Fourier
spectrum as the house form, each verifying green, and each equal to
`S·H·S` cell for cell. This is the `ψ(ρ) = 1` case of §1.4's
proposition doing exactly what it promises: the twist refutes
Gram-forcing as a literal statement while manufacturing nothing
new.

The sweep is also exhaustive, which cert 09 could not say on its
own. By Corollary E1 (§1.2.1) the admissible Grams at
`(s,i) = (3,4)` over `Ḡ = ℤ₂×ℤ₂` are exactly the four
real-character twists of the house form; cert 12 enumerates them
from first principles and finds cert 09's four banked Grams to be
that list, entry for entry. There is nothing at this cell for the
sweep to have missed.

### 3.4 The theorem at order 668 — PROVEN + PROVEN-BY-CERTIFICATE

*(The invariance of the folded 4-profile and the implication
"different profiles ⟹ inequivalent" are paper proofs; the eight
exact profiles — four matrices and their four transposes — the
twelve comparisons the transpose-extended statement needs and the
28 the eight-matrix statement needs are certified machine
computations, certs 06, 08, 13, 15 and 19. The
composite label reflects that split. Where those computations ran
matters, so it is stated here: the `C(668,4)` enumerations were made
in the source laboratory and at the desk, and banked. A default run
of cert 06, 08, 13, 15 or 19 **audits** those banks — the
file digests, the declared matrix identity against the matrix
rebuilt in that same run, the forced congruence, the total, the
second moment, and bin-for-bin agreement of two independent
implementations — and the `--full` flag **recomputes** the profiles
here from the rebuilt rows. Cert 08's `--full` has been run to
completion inside the repository, and cert 06's on its `blas` leg;
each cert's `NOTES.md` records which legs ran and when.)*

> **Theorem.** Order 668 carries at least **three** Hadamard
> equivalence classes. They are exhibited by the decoded `(1,1)`
> record `H`, its Lemma-T `i = 2` rebuild `H'` (cert 02), and the
> paired-Hall-switch matrix `H★` of the public preprint discussed
> below: the three exact 4-profiles are pairwise different, in
> 26 (`H` vs `H'`), 27 (`H` vs `H★`) and 27 (`H'` vs `H★`) of their
> 80 bins — and the theorem holds under the transpose-extended
> relation as well: all six comparisons separate (against the
> transposes: 50, 49 and 50 bins; both transposed profiles computed
> here populate 79 bins).

> **Theorem (the fourth class, 2026-09-02; cert 13).** Let `H″` be
> `H` with its **twelve off-diagonal core blocks negated** and the
> border unchanged. Then `H″` is Hadamard, and it is inequivalent to
> each of `H`, `H'`, `H★`: its exact 4-profile differs from theirs in
> **27, 27 and 26** of the 80 bins (two arithmetics, bin for bin;
> totals and second moment exact; run under the pre-registration
> `REGISTRATION-668-orientation.md`, flushed first). **Order 668
> therefore carries at least four Hadamard equivalence classes.**
> The transposed profile of `H″` was computed on 2026-09-02 under
> `experiments/pr0042/REGISTRATION.md` and separates `(H″)ᵀ` from
> each of `H`, `H'`, `H★` in **50, 50 and 49** of the 80 bins, so the
> fourth class is **transpose-robust** too and the four-class
> statement holds under the transpose-extended relation (cert 15).

*What `H″` is.* `S·H″·S`, with `S = diag(I₄, diag(1,−1,−1,−1)⊗I_n)`, is
the same seeds and border assembled in the **alternate Goethals–Seidel
orientation** (§1.0: the six transposed blocks negated) with the border
strips signed by superblock (`P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`) —
an identity of sign patterns, checked cell by cell in cert 13. And
`D·H_t·D = H″` for `D = diag(I₄, I₄⊗diag(ψ(g)))` and `H_t` the (D-e)
transport instance on the Lemma-T twisted seeds, so the twist and the
orientation switch coincide *up to the border*; that they land in
different classes (`H″ ≁ H'`) is exactly because `H'`'s border is a
Hall-type strip switch of the decoded border (§3.4 below), not its
transport. Two lessons the theorem states rather than suggests: the
GS orientation, a "convention" in §1.0, is **not a gauge for Hadamard
equivalence** at 668; and every bordered GS record found at any order
comes with a second candidate class for free, by negating twelve
blocks. Nothing is claimed here beyond 668, nor that four is the
count; the same switch at 716, 1676, 2060 and 1772 has since been
computed and is §3.6 (cert 14), §3.7 (cert 20), §3.5 (cert 22) and
§3.8 (cert 23).

**Priority, stated first.** An anonymous preprint (hosted at
hadamard-668.vercel.app; its PDF `CreationDate` field reads
2026-08-13 04:14 UTC+05:30, i.e. 2026-08-12 22:44 UTC — a
compilation timestamp, which is not by itself evidence of when the
page became public; it was public when retrieved and verified
firsthand on 2026-08-31) is, as far as our search located, the
first published statement that order 668 carries at least two
Hadamard equivalence classes, and it precedes this repository
either way. Its
`H` is **byte-identical to the decoded record banked here** — an
independent decode of the same public data, border included, which
corroborates ours and confirms that neither party originated the
construction data independently of the public posting; neither claims
priority over it. Its second matrix `H★` is a 1,328-entry paired
Hall switch of `H`, rebuilt here from its published data and
verified (both its SHA-256 digests reproduce). Its separating
statistic `Φ_M` — the correlation profile of a distinguished type-1
quadruple — is sound: its invariance uses the uniqueness of that
quadruple, stated and exactly computed in the preprint as its
Lemma 3, and this repository independently reproduces and certifies
that computation (a distinguished 4-subset is exactly one achieving
`|T4| = 660`, and the exact profile's `660` bin has count **1** on
every matrix here — cert 08 checks it directly). On this order,
however, `Φ_M` is blind where the full 4-profile is not: it is
**bin-for-bin identical on `H` and `H'`**, a pair §3.4 proves
inequivalent. The "at least two
classes" statement is the preprint's; the third class, the Lemma-T
construction that produces it, and the pairwise exact-profile
separations are this repository's.

The separating computation is the exact 4-profile: all
`C(668,4) = 8 222 179 035` 4-subsets, on each matrix, by two
independent implementations (a float32 BLAS Gram of the pair-vector
matrix, and a packed-`uint64` popcount path; opposite bit packings;
both validated bin-for-bin against straight `O(C(n,4))` stdlib
enumeration on five small control matrices), with a third
(canonical-split) implementation agreeing bin-for-bin on the
decoded record (a source-laboratory run, cited in cert 06's notes
rather than banked). All three matrices populate the **same 80
bins** —
the support does not separate them, and neither does the extreme
tail (top bins agree exactly). **The bulk separates them** — for
`H` vs `H'`, 26 of the 80 bin counts differ, e.g.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| decoded | 2 073 064 058 | 1 852 054 148 | 1 491 070 735 | 1 091 442 371 | 732 009 734 | 452 971 620 | 260 220 030 | 139 599 212 |
| rebuilt | 2 073 109 602 | 1 852 009 274 | 1 491 079 303 | 1 091 478 493 | 731 928 178 | 452 961 444 | 260 257 738 | 139 618 428 |

and the `H★` comparisons differ in 27 bins against each. All three
profiles total `C(668,4)` exactly, every pairwise difference vector
sums to zero, and all three hit the second-moment identity
`5 517 193 410 096` to the unit; the first moments — which nothing
forces — differ. The largest `H`-vs-`H'` bin discrepancy is
`1.1·10⁻⁴` of its bin, invisible to any sample of practical size (a
null sampled comparison at `2·10⁷` draws was already in hand and is
worth very little; the exact computation is what settles it).

On transpose: the theorem is **transpose-robust**, and so is the
fourth class. Refuting `A ≈ B` under the transpose-extended relation
needs both `A ≁ B` and `A ≁ Bᵀ`, and all twelve comparisons separate
exactly (certs 08 and 15):

```
H  vs H'  : 26 bins      H  vs (H')ᵀ  : 50 bins
H  vs H★  : 27 bins      H  vs (H★)ᵀ  : 49 bins
H' vs H★  : 27 bins      H' vs (H★)ᵀ  : 50 bins
H  vs H″  : 27 bins      H  vs (H″)ᵀ  : 50 bins
H' vs H″  : 27 bins      H' vs (H″)ᵀ  : 50 bins
H★ vs H″  : 26 bins      H★ vs (H″)ᵀ  : 49 bins
```

(the transposes are genuinely different objects — all four
transposed profiles computed here populate 79 bins where the
originals populate 80, each losing the same bin, `|T4| = 644`).
Order 668 carries at least **four** equivalence classes under either
convention.

> **Theorem (the record and its own transpose, 2026-09-02; cert 19).**
> `H` is **not** equivalent to `Hᵀ`: their exact profiles differ in
> **49** of the 80 bins. Consequently the eight matrices `H, H', H★,
> H″` and their four transposes are **pairwise** inequivalent — all 28
> comparisons separate, the closest pair by 24 bins — so **order 668
> carries at least eight Hadamard equivalence classes under plain
> Hadamard equivalence.** Under the transpose-extended relation the
> count is **unchanged**: exactly **four** classes among those eight,
> since `X ≈ Xᵀ` holds by definition and the twelve refutations above
> already separate the four. The house counts four.

Eight is a lower bound from the artifacts in hand, not a count of the
classes at 668, and four transposes are not four constructions: the
eight matrices come from four objects and one involution. What makes
the comparison a real question rather than a restatement is that the
transpose is **not** in the equivalence group, so `profile(Hᵀ)` is the
same invariant computed on a different matrix; that the bordered GS
family is closed under transposition — so `Hᵀ` is another record of
the same family — is source-laboratory intel, cited in cert 19's notes
and neither proved nor relied on here. Nothing general follows: a
Sylvester matrix is symmetric, and cert 19 measures Paley `H(20)`'s
profile agreeing with its transpose's — agreement being a failure to
separate, never a proof of equivalence. What is settled is one matrix
at one order.

Consequences. The Lemma-T construction at `ψ(ρ) = −1` genuinely
leaves the equivalence class at this order, and the Hall switch
leaves both. The corresponding statement at **716 has since been
made** — §3.6, cert 11: 27 of 87 bins differ, two independent
implementations, so the twist leaves the class at a second order,
and since 2026-09-02 that order's three-class statement holds under
the transpose-extended relation too (cert 15). At **1676 the
statement has since been made** — §3.7, cert 20, 2026-09-02: 68 of
142 bins differ, two arithmetics — and since cert 21 that order's
three-class statement holds under the transpose-extended relation
too; nothing in this section rests on it. At **1772 the statement
has since been made** too — §3.8, cert 23, 2026-09-03: 57 of 89 bins
differ, two arithmetics — row-side, and nothing in this section
rests on that either. The `Φ_M` blindness above is recorded beside the two traps of
§3.1 as the working reason this repository pays for the full
`Θ(n⁴)` invariant: every cheaper statistic tested — including a
sound published one — fails to see at least one true separation at
this order.

### 3.5 The theorem at order 2060 — PROVEN + PROVEN-BY-CERTIFICATE

The formerly-open order 2060 carries three classes among public+banked
artifacts: the publicly posted matrix, the plain GS-array
realisation over the same decoded seed (the `×104` character twist
relates the two seeds), and that realisation in the **alternate
Goethals–Seidel orientation**. Every cheap exact invariant we computed
agrees across them; the exact 4-profile does not.

**Theorem (the pair).** The two matrices are Hadamard-inequivalent;
order 2060 carries at least two Hadamard equivalence classes. *Proof.* `|T4|`
is an equivalence invariant (§3.1, I4). The exact folded 4-profiles
over all `C(2060,4) = 748 155 697 135` row 4-subsets, computed under
the pre-registration `REGISTRATION-2060-exact.md` (flushed before
the first run), **differ in 146 of the 147 union bins and in their
supports — 145 bins for the plain matrix against 133 for the posted
one**; totals and the second-moment identity `Σ T4² = n³(n−1)(n−2)/24`
hold exactly on both, and the two arithmetics (`blas` float32
accumulation of integers below `2²⁴`; `bits` pure-integer popcount)
agree bin for bin on each matrix. ∎ Certificate: cert 07 in EXACT
mode, whose four-file bank (`data/sep2060-exact-{blas,bits}-{plain,
gist}.json`, schema `sep2060-exact-profile/1`, each bound to its
matrix's canonical digest) is pinned byte-for-byte in `run.py` and
was accepted through the fail-closed predicate whose 22 negative
controls the certificate replays on every run. The largest bin
deficit is `≈ 6.2·10⁸` at `|T4| = 4` (0.6 % of that bin) — three orders
of magnitude above anything the sampled statistic below could
resolve, which is why the sampled section is retained as the record
of how the question was first asked. Runs: blas pair 9.5 h wall at 3
threads (2026-09-01), bits pair 22.4 h + a bit-identical resume from
checkpoint after a machine reboot (2026-09-02).

> **Theorem (the third class, 2026-09-03; cert 22).** Let `P` be the
> plain GS realisation and `H″` be `P` with its **twelve off-diagonal
> `515×515` blocks negated** — the orientation switch, in its
> **unbordered** form, `2060` being the degenerate `s = 0` layer where
> the array is a plain `4×4` GS array of circulant blocks and there is
> no border to leave alone (exactly `12·515² = 3 182 700` cells
> change). Then `H″` is Hadamard, and it is inequivalent to each of
> `P` and the posted matrix `G`: its exact 4-profile differs from
> `P`'s in **107 of the 145 bins they share** and from `G`'s in
> **146 of the 147** bins of their union (two arithmetics, bin for
> bin; totals and second moment exact). **Order 2060 therefore
> carries at least three Hadamard equivalence classes.**

*What `H″` is.* `S·H″·S`, with `S = diag(1,−1,−1,−1)⊗I₅₁₅`, is the
same four seeds assembled in the **alternate Goethals–Seidel
orientation** (§1.0: the six transposed blocks negated). There is no
border-signing term here, `s = 0` having emptied it, and because the
seeds are in hand cert 22 checks the identity **twice**: as a
sign-pattern identity cell by cell, as at 668 (§3.4), 716 (§3.6) and
1676 (§3.7), and by **assembling the alternate array directly from the
same raw seeds**, verifying it (canonical digest `40e1d1c8…ca64398d`)
and comparing it to `S·H″·S` cell for cell. The second check exhibits
`H″` as a signed conjugation of an independently built GS array — an
element of the equivalence group on each side — so the theorem is
about *orientation*, not about an arbitrary sign pattern. Hence, at
the founding order too, the GS orientation is **not a gauge for
Hadamard equivalence**: one seed quadruple, three classes, from two
constructions plus the switch. That is **five orders with the same
verdict and no general statement** (668, 716, 1676, 2060, and 1772
since §3.8); nothing is
claimed about `ψ(ρ) = −1` at 2060, the `×104` twist above acting on
the *column index* of each circulant (`data/sep2060-records.json`,
`twist`) rather than on the seed values as Lemma T's `ψ`-twist does,
so no analogue of the 668/716/1676 twist theorem is asserted here.

`H″` and `P` populate the **same** 145 bins — the switch moves counts,
not support — while `G` lacks fourteen of them and populates two of
its own, which is why the comparisons with `G` read *146 of 147*.
`|T4| = 1108` is the one bin where all three agree, at 30 counts
apiece. Two features distinguish `H″` vs `P` from the separations at
the bordered orders: the largest discrepancy is only `3 769 014` at
`|T4| = 12` (`3.7·10⁻⁵` of that bin, against `6.2·10⁸` and 0.6 % for
the pair) — **invisible to any sample of practical size** — and the
**extreme tail does separate**, the top bin `|T4| = 1236` reading 12
against 6, where at 668, 716, 1676 and 1772 the tail always agreed. A
36-bin band, every `|T4|` in `[868, 1180]`, agrees exactly; 106 of
the 108 bins below `868` differ. Certificate: cert 22, whose two new
banks (`data/sep2060-orient-exact-{blas,bits}.json`, schema
`exact-4-profile/1`) are pinned byte-for-byte in its `run.py`
alongside cert 07's four, each bound to the canonical digest of a
matrix that run rebuilds. `dim V` reads 2058 / 2058 / **2059** across
`P`, `G`, `H″` and is worthless; the invariant `dim W` is 2059 on all
three (Trap 1, §3.1). Enumeration: 6 631 s (`blas`) and 17 655 s
(`bits`) on 16 rented threads, under
`experiments/pr0042/REGISTRATION.md`, flushed before the matrix was
built.

*The sampled statistic, kept as history.* The **sampled** 4-profiles
differed first.

*What the statistic is.* Per bin, with counts `p` and `q` at equal
sample size and `p + q ≥ 200`, the reported quantity is the **raw
standardized delta** `d = |p − q| / √(p + q)` — the deviate of a
Poisson common-mean approximation. It is a heuristic z-score, **not
a calibrated sigma**, and this note does not treat it as one. The
draws *are* paired: one seeded 4-subset stream, depending only on
the seed and on `n`, is evaluated on both matrices. The denominator,
however, is the *unpaired* variance, and the per-bin discordance
counts a paired variance needs were never banked, so the correction
is not recoverable from this repository's data. `105` per-bin
comparisons were scanned across the two samplers' row and column
sides, with **no multiplicity control**. Cert 07's `NOTES.md` gives
the full accounting.

At `2·10⁷` paired draws the difference resolves at `max d = 7.4`
(row side, 21 bins with `d > 4`) and `7.1` (column side, 13 bins),
in a coherent monotone pattern across consecutive bins — more
peaked, heavier-tailed — whose deficits and excesses balance as they
must. The independent second sampler, at `3·10⁶` draws, does not by
itself reach `d = 4` (as `√N` scaling predicts); its corroboration
is sign agreement on 29 of the 34 bins the first sampler resolves —
the 29 where the second has enough mass to compare. The calibration
that keeps this honest: the same sampled statistic reads null
(`max d ≤ 2.2`, zero bins with `d > 4`) on the 668 pair that §3.4
**proves** inequivalent. No equivalence exists in the block-affine
family (BOUNDED-NEGATIVE-SEARCH: exhaustive over that family, silent
beyond it). The sampled statistic alone was labelled
**COMPUTATIONAL-EVIDENCE of inequivalence — not a proof**, and this
note never wrote it as one; the exact 4-profile above, the named
upgrade path, has now been run and is the proof. It was never run on
`H″`, and it would not have resolved `H″` against `P` if it had been.
Label for §3.5: **PROVEN + PROVEN-BY-CERTIFICATE** (cert 07 in exact mode for the
pair, cert 22 for the third class). **Row-side only: no
transpose-extended statement is made at 2060** — the transposed
profiles `H_2060-orient-T` and `H_2060-plain-T` are pending legs of
the same campaign. Since 2026-09-03 §3.8 shares that posture at
1772, and those two are the only orders in this note whose
separation statements are row-side. **Not claimed here:** anything
about any of the three matrices versus its own transpose (cert 19
decides that question at 668 only); any general theorem about
orientation; how many classes order 2060 has; and no priority or
novelty claim of any kind — the order was settled by the publicly
posted matrix, which is `G` itself (§4, `PROVENANCE.md`).

### 3.6 The theorem at order 716 — PROVEN + PROVEN-BY-CERTIFICATE

*(The same split as §3.4: the invariance of the folded 4-profile
and the implication "different profiles ⟹ inequivalent" are paper
proofs; the twelve exact profiles are certified machine computations
— certs 11, 14 and 15. Where they ran matters, so it is stated here:
the `C(716,4)` enumerations were made in the source laboratory on
2026-09-01 (the pair), 2026-09-02 (`H″`) and 2026-09-02 (the three
transposes), each under a pre-registration flushed before the
matrices were built, and banked. A default run of cert 11, 14 or 15 **audits** those banks;
`--full` recomputes the profiles here from the rebuilt rows, by a
third arithmetic route. Each cert's `NOTES.md` records which legs
ran.)*

> **Theorem.** Order 716 carries at least **two** Hadamard
> equivalence classes, exhibited by the decoded `(1,1)` record and
> its Lemma-T `i = 2` rebuild (cert 02): the two exact 4-profiles
> populate the **same 87 bins** and **27 of the 87 bin counts
> differ**.

> **Theorem (the third class, 2026-09-02; cert 14).** Let `H″` be
> the decoded record `H` with its **twelve off-diagonal core blocks
> negated** and the border unchanged. Then `H″` is Hadamard, and it
> is inequivalent to each of `H` and `H'`: its exact 4-profile
> differs from theirs in **27** and **25** of the 87 bins (two
> arithmetics, bin for bin; totals and second moment exact; run
> under the pre-registration `REGISTRATION-716-orientation.md`,
> flushed first). **Order 716 therefore carries at least three
> Hadamard equivalence classes.**

> **Theorem (transpose-extended, 2026-09-02; cert 15).** The three
> transposed profiles were computed under
> `experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC, before
> the matrices were built). Every one of the twelve
> transpose-related comparisons separates, so each pair carries
> **both** refutations the relation needs — and both routes, since
> all three transposes are banked: `H` vs `(H')ᵀ` 57 bins and `Hᵀ`
> vs `H'` 56; `H` vs `(H″)ᵀ` 57 and `Hᵀ` vs `H″` 56; `H'` vs
> `(H″)ᵀ` 57 and `(H')ᵀ` vs `H″` 57. **The three classes at 716
> therefore hold under the transpose-extended relation** — the
> statement the two theorems above withheld.

*What `H″` is.* `S·H″·S`, with `S = diag(I₄, diag(1,−1,−1,−1)⊗I_n)`
and `n = 178`, is the same seeds and border assembled in the
**alternate Goethals–Seidel orientation** (§1.0: the six transposed
blocks negated) with the border strips signed by superblock
(`P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`) — an identity of sign
patterns, checked cell by cell in cert 14, as at 668 (§3.4). So the
GS orientation, a "convention" in §1.0, is **not a gauge for
Hadamard equivalence** at 716 either, and the `ψ(ρ) = −1` twist and
the orientation switch land in different classes here as they do at
668. That is **four orders with the same verdict and no general
statement** (668, 716, and — since 2026-09-02 — 1676, §3.7, joined by
1772 in §3.8 on 2026-09-03): nothing here rests on those, and
nothing is said to the effect that three is the count at 716.

The separating computation is the exact 4-profile over all
`C(716,4) = 10 859 143 295` row 4-subsets, on each matrix, by two
independent implementations (a float32 BLAS Gram of the pair-vector
matrix, and a packed-`uint64` popcount path) that agree bin for bin
on each matrix. Both profiles total `C(716,4)` exactly and hit the
second-moment identity `7 807 861 101 040` to the unit; the
difference vector sums to zero, and the first moment — which
nothing forces — does not: `Σ |T4|·Δ = −279 888`.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| decoded | 2 650 505 561 | 2 383 887 265 | 1 944 278 842 | 1 450 375 604 | 996 815 228 | 634 922 458 | 377 261 905 | 210 018 619 |
| rebuilt | 2 650 421 375 | 2 383 982 899 | 1 944 334 154 | 1 450 404 118 | 996 701 290 | 634 954 430 | 377 258 529 | 209 998 559 |
| Δ | −84 186 | +95 634 | +55 312 | +28 514 | **−113 938** | +31 972 | −3 376 | −20 060 |

As at 668, the **bulk** separates and the tail does not: every one
of the 59 bins from `|T4| = 228` to `708` agrees exactly, and the
largest discrepancy — `113 938` at `|T4| = 36` — is `1.1·10⁻⁴` of
its bin, invisible to any sample of practical size. `dim V` differs
across this pair too (714 vs 715) and is worthless; the invariant
`dim W` is `715` on both (Trap 1, §3.1).

Consequence: the Lemma-T construction at `ψ(ρ) = −1` provably
leaves the equivalence class at a **second** order, so §1.4's
question had two proven instances when this section was written
(§3.7 makes it three); the general statement stays unclaimed. A remark, recorded and not headlined: under **plain**
Hadamard equivalence the six objects `H, H', H″` and their
transposes are pairwise inequivalent (all fifteen comparisons
separate, the closest pair by 25 bins), so three constructions
exhibit six classes here; the count stated is the three that survive
**either** convention. **Not claimed here:** nothing about how many
classes order 716 has; the sibling order 1772 was a separate
computation and now carries its own statement in §3.8 (cert 23), on
which nothing here depends. No priority or novelty claim is made at 716: this is a
statement about how many classes are on the table among the
artifacts banked here, not about existence. The sibling order 1676
has since been computed and is §3.7 (cert 20, 2026-09-02); nothing
in this section rests on it.

### 3.7 The theorem at order 1676 — PROVEN + PROVEN-BY-CERTIFICATE

*(The same split as §3.4 and §3.6: the invariance of the folded
4-profile and the implication "different profiles ⟹ inequivalent"
are paper proofs; the ten exact profiles are certified machine
computations — certs 20 and 21. Where they ran matters, so it is
stated here: the five `C(1676,4)` enumerations were made in the
source laboratory on 2026-09-02, on a rented 16-thread machine,
under `experiments/pr0042/REGISTRATION.md` — flushed 10:17 UTC,
before any matrix it governs was built — and banked. The matrices
were built and verified **at the desk**, each through
`verify/verify.py`, before anything was uploaded; the rented machine
enumerated and nothing else. A default run of cert 20 or 21
**audits** those banks. Unlike certs 06, 08, 11, 13, 14, 15 and 19,
neither has an **in-repo `--full` leg**: at this order one
recomputation is of order 6–7
hours — 52× cert 14's 716 leg on the source laboratory's measured
sub-`n⁵` scaling, 70× and ≈ 7.8 h on the `Θ(n⁵)` law §3.4 quotes —
and the `blas` route wants about 9.4 GB, so the flag is offered and
priced but has not been run. Their `NOTES.md` say so, and say which
estimate they are using.)*

> **Theorem (the twist, a third instance).** Order 1676 carries at
> least **two** Hadamard equivalence classes, exhibited by the
> decoded `(1,1)` record `H` and its Lemma-T `i = 2` rebuild `H′`
> (cert 02): the two exact 4-profiles populate the **same 142 bins**
> and **68 of the 142 bin counts differ**. So the Lemma-T
> construction at `ψ(ρ) = −1` provably leaves the equivalence class
> at a **third** order (§1.4).

> **Theorem (the third class).** Let `H″` be `H` with its **twelve
> off-diagonal core blocks negated** and the border unchanged. Then
> `H″` is Hadamard, and it is inequivalent to each of `H` and `H′`:
> its exact 4-profile differs from theirs in **70** and **66** of
> the 142 bins (two arithmetics, bin for bin; totals and second
> moment exact). **Order 1676 therefore carries at least three
> Hadamard equivalence classes.**

> **Theorem (transpose-extended, 2026-09-02; cert 21).** The
> transposed profiles of `H′` and `H″` were computed under the same
> `experiments/pr0042/REGISTRATION.md`. Each of the six
> transpose-related comparisons separates in **139 of the 144** bins
> of the union support, so every pair carries **both** refutations
> the relation needs: `H` vs `(H′)ᵀ` 139 alongside `H` vs `H′` 68;
> `H` vs `(H″)ᵀ` 139 alongside `H` vs `H″` 70; `H′` vs `(H″)ᵀ` 139
> alongside `H′` vs `H″` 66, with `(H′)ᵀ` vs `H″` 139 as the second
> route. **The three classes at 1676 therefore hold under the
> transpose-extended relation** — the statement the two theorems
> above withheld. `Hᵀ` was **not** profiled at this order, so
> nothing is said about `H` versus its own transpose here (cert 19
> decides that question at 668 only), and each pair whose `A` is `H`
> is refuted by its `A` vs `Bᵀ` column alone.

*What `H″` is.* `S·H″·S`, with `S = diag(I₄, diag(1,−1,−1,−1)⊗I_n)`
and `n = 418`, is the same seeds and border assembled in the
**alternate Goethals–Seidel orientation** (§1.0: the six transposed
blocks negated) with the border strips signed by superblock
(`P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`) — an identity of sign
patterns, checked cell by cell in cert 20, as at 668 (§3.4) and 716
(§3.6). So the GS orientation, a "convention" in §1.0, is **not a
gauge for Hadamard equivalence** at 1676 either, and the
`ψ(ρ) = −1` twist and the orientation switch land in different
classes here as they do at 668 and 716. That is **three orders with
the same verdict and no general statement** — four, since §3.8 added
1772 on 2026-09-03; nothing here rests on that, and nothing is said
to the effect that three is the count at 1676.

The separating computation is the exact 4-profile over all
`C(1676,4) = 327 588 749 775` row 4-subsets, on each matrix, by two
independent implementations (a float32 BLAS Gram of the pair-vector
matrix, and a packed-`uint64` popcount path) that agree bin for bin
on each matrix. All three profiles total `C(1676,4)` exactly and hit
the second-moment identity `550 023 273 154 800` to the unit; every
difference vector sums to zero, and the first moments — which
nothing forces — do not: `Σ |T4|·Δ = −20 505 840` (`H` vs `H′`),
`+11 699 520` (`H` vs `H″`), `+32 205 360` (`H′` vs `H″`).

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H` | 55 147 198 687 | 52 329 871 149 | 47 215 439 663 | 40 653 198 021 | 33 547 376 838 | 26 657 993 721 | 20 481 183 809 | 15 276 024 237 |
| `H′` | 55 146 545 415 | 52 329 311 849 | 47 218 249 473 | 40 653 475 393 | 33 547 530 716 | 26 657 747 145 | 20 479 582 049 | 15 275 249 189 |
| `H″` | 55 145 104 333 | 52 329 711 391 | 47 218 534 903 | 40 653 673 809 | 33 547 436 380 | 26 658 420 825 | 20 479 041 633 | 15 275 675 349 |

As at 668 and 716, the **bulk** separates and the tail does not:
every bin above `|T4| = 564` agrees on both comparisons with `H`,
and every bin above `|T4| = 540` on `H′` vs `H″`. The largest
discrepancies — `2 809 810`, `3 095 240` and `1 441 082`, all in the
first three bins — are of order `10⁻⁵` of their bins, invisible to
any sample of practical size. `dim V` differs across the pair here
too (1674 vs 1675) and is worthless; the invariant `dim W` is
`1675` on all five objects §3.7 now profiles (Trap 1, §3.1).

*The transposed profiles* (cert 21) populate 142 bins apiece too,
but **not the same 142**: unlike 668 and 716, where each transpose
drops one bin, at 1676 the two transposes **swap two** — `|T4| =
948` and `964` are populated in every original (counts 30 and 40)
and in neither transpose, and `|T4| = 1180` and `1204` in both
transposes (30 and 30) and in no original. The union support is
therefore **144** bins, which is why every original-versus-transpose
comparison reads *139 of 144*. The bulk separates there too: the top
differing bin is `|T4| = 1204` and the four bins above it agree, the
largest discrepancies (`5.0`–`5.8·10⁶`, mostly at `|T4| = 36`) being
of order `10⁻⁴` of their bins. A remark, recorded and not headlined:
under **plain** Hadamard equivalence the five matrices profiled here
— `H, H′, H″, (H′)ᵀ, (H″)ᵀ` — are pairwise inequivalent (all ten
comparisons separate, the closest pair by 66 bins), so **at least
five** classes are exhibited; the count stated is the three that
survive **either** convention. Since `Hᵀ` was never profiled, that
plain count is a lower bound from five matrices, not the analogue of
§3.4's eight.

**Not claimed here:** anything about `H` versus `Hᵀ` at 1676, whose
profile was never computed; how many classes order 1676 has; and no
priority or novelty claim of any kind — this is a statement about
how many classes are on the table among the artifacts banked here,
not about existence at 1676. The same three constructions exist at
1772; that order was a separate computation and now carries its own
statement in §3.8 (cert 23), on which nothing above depends.

### 3.8 The theorem at order 1772 — PROVEN + PROVEN-BY-CERTIFICATE

*(The same split as §3.4, §3.6 and §3.7: the invariance of the
folded 4-profile and the implication "different profiles ⟹
inequivalent" are paper proofs; the three exact profiles are
certified machine computations — cert 23. Where they ran matters, so
it is stated here: the three `C(1772,4)` enumerations were made in
the source laboratory between 2026-09-02 23:33 UTC and 2026-09-03
11:24 UTC, on a rented 16-thread machine, under
`experiments/pr0042/REGISTRATION.md` — flushed 2026-09-02 10:17 UTC,
before any matrix it governs was built, and naming these three
objects and this decision rule in advance — and banked. The matrices
were built and verified **at the desk**, each through
`verify/verify.py`, before anything was uploaded; the rented machine
enumerated and nothing else. A default run of cert 23 **audits**
those banks. Like certs 20, 21 and 22 it has no **in-repo `--full`
leg**: at this order one recomputation is of order 7–8 hours — 68×
cert 14's 716 leg on the source laboratory's measured sub-`n⁵`
scaling, 93× and ≈ 10.3 h on the `Θ(n⁵)` law §3.4 quotes — and the
`blas` route wants about 11.1 GB, so the flag is offered and priced
but has not been run. Its `NOTES.md` says so, and says which
estimate it is using.)*

> **Theorem (the twist, a fourth instance).** Order 1772 carries at
> least **two** Hadamard equivalence classes, exhibited by the
> decoded `(1,1)` record `H` and its Lemma-T `i = 2` rebuild `H′`
> (cert 02): the two exact 4-profiles populate the **same 89 bins**
> and **57 of the 89 bin counts differ**. So the Lemma-T
> construction at `ψ(ρ) = −1` provably leaves the equivalence class
> at a **fourth** order (§1.4) — and 1772 was the last of cert 02's
> four decoded `(1,1)` orders left unclaimed.

> **Theorem (the third class).** Let `H″` be `H` with its **twelve
> off-diagonal core blocks negated** and the border unchanged. Then
> `H″` is Hadamard, and it is inequivalent to each of `H` and `H′`:
> its exact 4-profile differs from theirs in **58** and **53** of
> the 89 bins (two arithmetics, bin for bin; totals and second
> moment exact). **Order 1772 therefore carries at least three
> Hadamard equivalence classes.**

*What `H″` is.* `S·H″·S`, with `S = diag(I₄, diag(1,−1,−1,−1)⊗I_n)`
and `n = 442`, is the same seeds and border assembled in the
**alternate Goethals–Seidel orientation** (§1.0: the six transposed
blocks negated) with the border strips signed by superblock
(`P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`) — an identity of sign
patterns, checked cell by cell in cert 23 over all `1772²` cells, as
at 668 (§3.4), 716 (§3.6) and 1676 (§3.7). So the GS orientation, a
"convention" in §1.0, is **not a gauge for Hadamard equivalence** at
1772 either, and the `ψ(ρ) = −1` twist and the orientation switch
land in different classes here as they do at the other three
bordered orders. With §3.5's unbordered instance that is **five
orders with the same verdict and no general statement**: nothing is
claimed at any order not computed, and nothing to the effect that
three is the count at 1772.

The separating computation is the exact 4-profile over all
`C(1772,4) = 409 422 905 815` row 4-subsets, on each matrix, by two
independent implementations (a float32 BLAS Gram of the pair-vector
matrix, and a packed-`uint64` popcount path) that agree bin for bin
on each matrix. All three profiles total `C(1772,4)` exactly and hit
the second-moment identity `726 727 740 809 840` to the unit; every
difference vector sums to zero, and the first moments — which
nothing forces — do not: `Σ |T4|·Δ = +55 798 800` (`H` vs `H′`),
`+72 493 792` (`H` vs `H″`), `+16 694 992` (`H′` vs `H″`).

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H` | 61 934 029 130 | 59 717 426 100 | 55 520 798 861 | 49 778 171 281 | 43 034 414 080 | 35 886 330 792 | 28 860 579 439 | 22 388 208 089 |
| `H′` | 61 931 468 288 | 59 716 764 308 | 55 522 617 277 | 49 780 050 525 | 43 034 542 248 | 35 884 341 470 | 28 861 416 579 | 22 388 284 353 |
| `H″` | 61 931 541 520 | 59 716 317 614 | 55 523 091 925 | 49 779 015 109 | 43 035 240 878 | 35 883 263 724 | 28 862 137 199 | 22 388 244 191 |

As at 668, 716 and 1676, the **bulk** separates and the tail does
not: every bin above `|T4| = 476` agrees on both comparisons with
`H`, and every bin above `|T4| = 452` on `H′` vs `H″`. Below those
tops the agreement is nearly empty — three bins, two bins and four
bins respectively, which cert 23 enumerates and asserts. The largest
discrepancies — `2 560 842`, `3 067 068` and `1 077 746` — are of
order `10⁻⁵` of their bins, invisible to any sample of practical
size. `dim V` reads 1770 / 1771 / 1770 across the three and is
worthless; the invariant `dim W` is `1771` on all three (Trap 1,
§3.1).

**Row-side only.** The transposed profiles at 1772 — `(H′)ᵀ` and
`(H″)ᵀ` — are separate legs of the same campaign: their `blas` legs
are in, their `bits` legs were still running when cert 23 was
written, and this repository certifies a profile only once the two
arithmetics agree bin for bin, so they are **not banked** and
**nothing at 1772 is claimed under the transpose-extended
relation**. That is cert 20's caveat at 1676, which cert 21 later
discharged there; a later certificate will do the same here. 1772
and 2060 are now the only two row-side statements in this note.

**Not claimed here:** anything about `H` versus `Hᵀ` at 1772, whose
profile was never computed (cert 19 decides that question at 668 and
at 668 only); how many classes order 1772 has; any general theorem
about `ψ(ρ) = −1` or about orientation, four instances being four
instances; and no priority or novelty claim of any kind — this is a
statement about how many classes are on the table among the
artifacts banked here, not about existence at 1772, which is long
settled.

---

## 4. Prior art and credit

Every credit below was fixed by a firsthand read of the named source
(the full reading trail, with per-source verdicts, is in the lab
record of 2026-08-31). Where a fact is attributed by a source we did
not read, that is said explicitly.

### 4.1 The classical spine

**The four-block array.** The Goethals–Seidel array — four
group-developed blocks with a back-diagonal reflection — is
**Goethals–Seidel 1970** (*A skew Hadamard matrix of order 36*,
J. Austral. Math. Soc. 11, 343–344), the reference every firsthand
source we read gives for it (Wallis–Whiteman 1972 ref [4];
Wallis–Street–Wallis 1972 ref [45]; Spence 1975 (CJM) ref [3];
Georgiou–Koukouvinos–Seberry 2002 ref [40]). Their 1967 paper
(*Orthogonal matrices with zero diagonal*, Canad. J. Math. 19,
1001–1010 — read in full) contains the two-block normal form only
(Thm 2.3, p. 1004, necessary direction, Paley matrices; the paper
numbers two theorems "3.4" — cite it by page), and width-one borders
only.

**The general-abelian setting is Wallis–Whiteman's.** *Some classes
of Hadamard matrices with constant diagonal* (Bull. Austral. Math.
Soc. 7 (1972) 233–249 — read in full) builds the type-1/type-2
incidence calculus over any additive abelian group (Lemmas 2–9,
Corollaries 4/6/10) precisely to carry the array there: Theorem 11
(p. 242) is the unbordered array over general abelian `G`, and our
Lemma 1 is that calculus in developed-matrix form. Nothing about the
abelian setting is claimed here.

**The `s = 1` border is classical, three times over.** Width-4
border strips, constant on the whole group, around the four-block
core: **Wallis–Whiteman 1972, Theorem 12** (p. 243; odd `|G|`, from
`4-{2m+1; m; 2(m−1)}` SDS), restated with the border blocks written
out in **Wallis–Street–Wallis 1972** (LNM 292 — Part 4 read in
full), Thm 4.17 p. 334 and Thm 8.44 pp. 393–394, where the identity
`4 = Σ(2kᵢ−v)²` (p. 393) already pins the feeding SDS block sizes at
this width; and **Spence 1975** (*Hadamard matrices from relative
difference sets*, JCTA 19, 287–300 — read in full), Theorem 2.1
(p. 289), the even-`|G| = 2v` cyclic sibling with blocks from
Elliott–Butson relative difference sets, which Spence himself
introduces as "similar to one of Wallis and Whiteman". (Spence wrote
three related papers; the other two — *Skew-Hadamard matrices of the
Goethals-Seidel type*, CJM 27 (1975), and *Skew Hadamard matrices of
order 2(q+1)*, Discrete Math. 18 (1977) — are unbordered and are not
the paper meant here. His `s` in JCTA Thm 3.2 is a 2-power doubling
parameter, unrelated to our border half-width `s`.) The width-2
sibling on a two-block core is older still: **Blatt–Szekeres 1969**
(Canad. J. Math. 21, 1319–1322 — read in full; received May 1968)
states it over *any* abelian group of order `2m+1`, needed there
precisely because their order-52 difference sets live in `ℤ₅×ℤ₅`
(group order 25) and, by their machine search, in no cyclic group
of order 25 — and
closes "there seems to be no obvious generalization of this
construction"; it reappears as Wallis–Street–Wallis Thm 4.4 (p. 321)
and Wallis–Whiteman Lemma 15. (The WSW book's own reference list
misdates Blatt–Szekeres to CJM 22 (1970); the paper is 21 (1969).) Modern sources treat the `s=1`
array as standard: Momihara–Xiang (arXiv:1801.08776 — read in full)
reprint both bordered arrays, crediting exactly Wallis–Street–Wallis
Thm 4.4 and Wallis–Whiteman 1972.

**The one prior width-parameterised border.** Wallis–Street–Wallis
1972, **Lemma 7.8** (p. 361), borders an `H[4t,4,t]` Hadamard-array
core with strips of width `4t`. Its width parameter is the array
repetition number `t`, not a subgroup index; its strips are constant
on all of `G`; it forces symmetric generating sets (so it cannot
preserve skewness); and it comes with no characterisation and no
classification — the book itself lists both the array parameter
(Question 5, p. 444) and the bordered family's feeding sets
(Question 16, p. 445) as open problems.

**Compression is Đoković–Kotsireas's.** The device behind our
Lemma 3 and the Σ̄ law — summing a sequence over the cosets of a
subgroup and the induced PAF identity — is their *compression*:
arXiv:1302.0571 (Def. 3, Thm 3; cyclic) and, at full abelian
generality, arXiv:1801.07627 (§7, Lemma 3 and Thm 4 — a paper whose
§10 already treats Goethals–Seidel quadruples over abelian groups).
The technique's own genealogy is "subgroup contraction"
(cf. Lumsden–Kotsireas–Bright, arXiv:2408.15611, and the difference-
set literature they cite). What is ours in Lemma 3 is only the
matrix-level bookkeeping that the compressed core is again a
Goethals–Seidel array over `G/K`, reflection included.

**The two-tier profile is a divisible-difference-family condition.**
The profile `4nδ₀ − 4s[K∖0] + 4[∉K]` says exactly that the four base
blocks form a `(G, K, {k_i}, λ, μ)`-**divisible difference family**
in the sense of Momihara–Yamada (arXiv:1212.3021, Def. 1.1) — we use
their language, not new terminology. Their Theorem 4.1 is also the
closest prior *coset-structured* border in the sources enumerated
here: border
rows constant on cosets of a subgroup `N` — around a **two-block**
core, with width tied to `2|N|` by the construction, targeting
symmetric Hadamard matrices of order `n²`, and with no converse
(their Problem 4.3 asks for what a characterisation would provide).
More broadly, coset structure imposed on the **blocks** of a family
is thoroughly classical (cyclotomic classes; multiplier-invariant
blocks); every border strip in the sources enumerated here, that one
exception aside, is constant on the whole group.

**Classification prior art.** Đoković–Kotsireas (arXiv:1802.00556)
*classify* Goethals–Seidel difference families — cyclic, unbordered,
by symmetry type of the blocks, exhaustively for odd `v < 50`, with
an iff on parameter sets (their Prop. 3) — an axis orthogonal to
Theorem C, which classifies border width against subgroup index and
is uniform in `v`. The surveys read in full or at stated coverage —
Seberry–Yamada 1992 (cover to cover; its §8 M-structure borders are
width one *per plug-in block*, Miyamoto's mechanism, not a border
around the assembled array), Seberry–Yamada 2020 (§3.6, §5.10.1–4,
and Appendix A complete from the A.1 heading — every construction
key table A.1–A.15 plus Table A.17 as captured), the *Handbook of
Discrete and Combinatorial Mathematics* (Rosen et al., 2nd ed. 2018;
Chapter 12 read in full with a machine-verified sweep of all 1,615
pages — its only bordered construction is the width-1 Paley border,
p. 902; note its §12.1–12.3 are the Colbourn–Dinitz editors'
self-declared *condensation* of their 2007 Handbook, so this read
does not discharge that item), Georgiou–Koukouvinos–
Seberry 2002 (in full; the word "border" does not occur in it), and
the Cati–Pasechnik database (whose bordered entries are Paley I/II,
Spence's Theorem 2.1, and the Blatt–Szekeres attribution above) —
contain no parameterised border width, no characterisation, and no
such classification. Wallis–Whiteman's Theorem 16 (a wide
conference-matrix Kronecker border) and Spence's Figure 2 (a width-8
border tracking an eight-block doubled core, entries constant on the
whole group) are different devices, named here so they are not
mistaken for the object below.

### 4.2 What is claimed, and the one hedged sentence

Ours are: the border whose **strips are constant on the cosets of a
proper subgroup** `K ≤ G` of index `i ≥ 2` (realised at
`s = 3, 5, 7`), on a single four-block Goethals–Seidel core; the
**iff** (Theorem A); the **classification** (Theorem C: in the
house-Gram branch, under `w > 2s`, and with `n ≥ 3` for the three
cells — `i = s+1`, `s` odd, beyond `(0,1)` and `(1,1)`); Theorem D
and the index-2 collapse; the
twist propositions; and the separation theorem of §3.4. The
`s ≤ 1` layer is classical as credited above, and the parameter
records instantiating `s = 3, 5, 7` are the announcing team's
(§2.1, `PROVENANCE.md`) — the theorems are ours, the instances at
those parameters are not claimed as new objects.

The novelty statement, exactly once and bounded: **in the sources
enumerated in §4.1, each read firsthand at the stated coverage, we
did not locate the coset-constant border form with `i ≥ 2`, the
if-and-only-if characterisation, or the width-against-index
classification; this closes those sources and nothing more.** The
nearest prior objects located — Wallis–Street–Wallis Lemma 7.8 and
Momihara–Yamada Theorem 4.1 — are distinguished above. One named
source could not be obtained and stays outside the enumeration: the
Craigen–Kharaghani Handbook chapter (V.1, 2007).

### 4.3 The public record around the instances

The twelve records derive from the announcement of 2026-08-12 and
its public artifacts (`PROVENANCE.md` carries the pinned chain); the
order-2060 artifact is Schneider's. Two independent Lean
formalizations of a single `H(668)` from the same public data exist:
the comparator of Paul-Lez (registered 2026-08-17), whose
`Challenge.lean` independently exhibits the bordered structure at
668 — circulants, a width-4 border, per-block-constant strips, the
same 4×4 corner — and Ramos–Hulak–de Queiroz (2026-08-28), which
verifies one supplied matrix. Neither states the general
construction theorem or the classification; a formalised Theorem
A/B remains a distinct object, planned separately. The anonymous
vercel-668 preprint (§3.4) independently decoded the same record,
border included, and is — as far as our search located — the first
published statement of an inequivalence at the order; the bounded
priority statement is §3.4. Public s=1 readings of the
same data (independent posts of 2026-08-31) are noted in
`PROVENANCE.md` as corroboration that the construction is broadly
understood.

---

## 5. Labels

| claim | label |
| --- | --- |
| Theorem A (iff); Theorem B; Lemmas 1–3; Theorem C + the classification — within the house-Gram branch, under `w > 2s`, and with `n ≥ 3` for the three surviving cells | **PROVEN** (paper-grade; not machine-checked) |
| Theorem D (D-a, D-a′, D-b … D-e); the index-2 collapse corollary | **PROVEN** (paper-grade) + **MEASURED** ((D-b)/(D-c) machine-validated on the four `i = 2` instances, cert 02; all four clauses (D-a)–(D-d) on the `H(52)` gate, cert 03) + **PROVEN-BY-CERTIFICATE** (the `768²` census of (D-e), and the compressed-block identities (D-a′)/(D-d) rest on; cert 10) |
| Lemma T; the `ψ(ρ)=1` conjugation proposition | **PROVEN** |
| the twelve public records satisfy every hypothesis | **MEASURED** (cert 01, 12/12) |
| the twelve assembled matrices are Hadamard and match their pinned digests | **PROVEN-BY-CERTIFICATE** (cert 01) |
| the twelve matrices constructed here — eight at `s ≤ 1` and in the `i = 2` frame (certs 02–05), and the four even-`s` `(2,4)` matrices of §2.3 (cert 18) | **PROVEN-BY-CERTIFICATE** |
| at `(s,i) = (1,2)` the genuine branch's Gram is house (the alternative is the degenerate `i=1` rewriting) | **PROVEN** (Theorem D-a); that the degenerate branch really is the `i = 1` construction — same borders, same matrices — is (D-a′) |
| **Theorem E** (Gram rigidity): under (H1)+(H2) with `s ≥ 1` and `w > 2s`, every admissible `M` is `4i·P_S` for a conjugation- and Galois-stable `S ⊆ Ĝ̄` of size `s`, so the spectrum is `{0^(i−s), (4i)^s}` at every cell | **PROVEN** (paper-grade; §1.2.1) |
| at `i = s+1` under `w > 2s` the admissible Grams are exactly the real-character twists of the house form (Corollary E1) — this is the former "`M` house form forced for `s ≥ 2` (up to twist)", now with the twist family proved complete | **PROVEN** (paper-grade) **+ MEASURED** (the complete `(3,4)` classification, both groups of order 4, all six Grams `(H1)`-realizable; cert 12) |
| `i ≤ s+1` (Theorem C's D1) fails in the general branch: (H1)–(H4) instances exist at `(s,i) = (1,11)` and beyond, carrying the Gram `4i·P_trivial` as Theorem E predicts | **PROVEN-BY-CERTIFICATE** (eight general-branch witnesses at orders 1676 and 1772; cert 12). That they contain nothing new *in general* is **NOT CLAIMED** — the (D-a′) collapse at general `i` is open |
| **Theorem E′** (§1.7): the mod-4 lemma `PAF_x(t) ≡ n (mod 4)`; Gram rigidity under `w > s`, not merely `w > 2s`; and the classification (a)–(d) of the `w = s` boundary | **PROVEN** (paper-grade) **+ PROVEN-BY-CERTIFICATE** (the classification exhausted by two exact routes sharing no code at 22 boundary cells, rigidity confirmed at 26 cells with `w > s` — 20 of them in the range `s < w ≤ 2s` Theorem E did not reach — and 11 cells recorded below the boundary; cert 16) |
| the `(3,3,3)` escape is realised by seeds, so `w > s` cannot be weakened to `w ≥ s` | **PROVEN-BY-CERTIFICATE** (seeds rebuilt by a full meet-in-the-middle over all 512 sequences per seed — 490 212 ordered solutions on `ℤ₃²`, none on `ℤ₉` — and `Q` by depth-first search; cert 16) |
| **Theorem F** (§1.8) (a)–(g): the border-kit structure, `E = −(1/4i)PĈᵀQ`, (H3)+(H4) `⟺ E ∈ {±1}`, S-part dependence, transport across `w` | **PROVEN** (paper-grade) **+ MEASURED** (all of (a)–(f) re-multiplied on 7 / 7 banked coset-border records with `s ≥ 1`, `w > 2s`; cert 17 [D]) |
| at `(3,4)`: exactly 2048 admissible S-parts per group, all of the silent-seed shape; the `N = 2092` shell `4192 × 64 = 2048 × 131 = 268 288` | **PROVEN-BY-CERTIFICATE** (exhaustive brute force over the spectra, a second implementation; cert 17 [A]) |
| at `(3,4)`: every S-part admits a border kit at every `κ(ρ)` (16 384 classes) | **PROVEN-BY-CERTIFICATE** (a deterministic 256-class sample re-verified exactly at `w = 130` on the default path; the full 16 384-class census on `--full`, run in this repository at the pinned digest — cert 17) |
| the `(3,4)` cell at 2092 is **one-layer**: a house-profile quadruple on a group of order 520 with an index-4 subgroup *is* an `H(2092)` | **PROVEN** (Theorem A + the census). The hypothesis is unmet: no such quadruple is known, and the seed layer is **OPEN** |
| the four even-`s` `(2,4)` matrices — `H(88)` ×2, `H(56)` ×2 — are Hadamard, are not equivalent to `H(2) ⊗ H′`, and are not a collapse | **PROVEN-BY-CERTIFICATE** (cert 18; digests pinned three ways). No novelty of existence at 56 or 88 is claimed |
| the `(2,4)` border proposition (a) and (b) — with (b), the necessity of the anti-periodic row table, holding only for `w > 2s = 4` | **PROVEN** (§2.4). At the `H(56)` instances (`w = 3`) the anti-periodic kit is sufficient, not shown forced |
| the `(2,4)` border census: every `(S-part, κ(ρ), Q′)` class admits a kit, 215 040 / 215 040 | **PROVEN-BY-CERTIFICATE** (exhaustive on the default path, with the 480 `H(8)` right-orbits counted two ways; cert 18 [4]) |
| **Theorem 3** (§2.4): if `S` contains no real character then `4 \| i`, `s` is even and `N ≡ 0 (mod 8)`; so at `N ≡ 4 (mod 8)` every realisable Gram has a live real character | **PROVEN** (paper-grade) **+ PROVEN-BY-CERTIFICATE** (the small-case exhaustion, 9 / 9 cases with `4 ∤ i`, `s ≤ 2` empty by exhaustive DFS; cert 18 [3]) |
| the two-valued Parseval law `Σ_q \|x̂_q(χ∘κ)\|² ∈ {4s, N}` (Corollary E2), generalising (D-d)'s `Σ_q δ_q² = 4` | **PROVEN** |
| all three nontrivial quotient-character twists at order 1916 are diagonally conjugate to the house instance, and the four Grams are the complete admissible list at that cell | **PROVEN-BY-CERTIFICATE** (cert 09) **+ PROVEN** (Corollary E1) **+ MEASURED** (banked set `=` census set; cert 12) |
| **order 668 carries at least four equivalence classes** | **PROVEN** (profile invariance + separation implication) **+ PROVEN-BY-CERTIFICATE** (the exact profiles: `H`, `H'`, `H★` pairwise 26/27/27 differing bins, cert 06/08; the orientation switch `H″` against them 27/27/26, cert 13; two independent implementations per matrix) |
| the **four**-class theorem at 668 under the transpose-extended relation | **PROVEN + PROVEN-BY-CERTIFICATE** (all twelve comparisons separate: the first three pairs in cert 08, the three pairs with `H″` from the new `(H″)ᵀ` profile — 50/50/49 bins — in cert 15) |
| `H ≁ Hᵀ` at 668 (49 of 80 bins), and hence **order 668 carries at least eight equivalence classes** under plain Hadamard equivalence — the four matrices and their four transposes, pairwise separated | **PROVEN** (profile invariance + separation implication) **+ PROVEN-BY-CERTIFICATE** (all 28 comparisons separate, 24–50 bins, two independent implementations per matrix, from the `Hᵀ` profile banked 2026-09-02; cert 19). Under the transpose-extended relation the count is **unchanged at four** — `X ≈ Xᵀ` by definition — and the house counts four. Eight is a lower bound, and four transposes are not four constructions |
| the vercel-668 preprint's `H★` and both its digests reproduce; its `Φ_M` is a valid invariant | **MEASURED** + **PROVEN-BY-CERTIFICATE** (rebuilt firsthand; its Lemma-3 uniqueness reproduced in cert 08 via the `660` bin) |
| `Φ_M` is blind to the `H` vs `H'` separation | **MEASURED** (bin-for-bin identical on a proven-inequivalent pair) |
| **order 716 carries at least three equivalence classes** | **PROVEN** (profile invariance + separation implication) **+ PROVEN-BY-CERTIFICATE** (the exact profiles: `H` vs `H'` 27 of 87 differing bins, cert 11; the orientation switch `H″` against them 27 and 25, cert 14; two independent implementations per matrix) |
| the three-class theorem at 716 under the transpose-extended relation | **PROVEN + PROVEN-BY-CERTIFICATE** (all twelve transpose-related comparisons separate, 26–57 bins each, from the three transposed profiles banked 2026-09-02; cert 15). Under plain equivalence the six objects `H, H', H″` and their transposes are pairwise inequivalent — a remark, not the count |
| **order 1676 carries at least three equivalence classes** | **PROVEN** (profile invariance + separation implication) **+ PROVEN-BY-CERTIFICATE** (the exact profiles: `H` vs `H'` 68 of 142 differing bins — the `ψ(ρ) = −1` twist's **third** proven instance; the orientation switch `H″` against them 70 and 66; two independent implementations per matrix; cert 20). Cert 20 was **row-side only**; cert 21 discharges that caveat, below. Neither has an in-repo `--full` leg: at this order one is ≈ 6–7 h on the measured sub-`n⁵` scaling (≈ 7.8 h on the `Θ(n⁵)` law) and the `blas` route wants ≈ 9.4 GB |
| the three-class theorem at 1676 under the transpose-extended relation | **PROVEN + PROVEN-BY-CERTIFICATE** (all six transpose-related comparisons separate, 139 of the 144 bins of the union support each, from the two transposed profiles banked 2026-09-02; cert 21). Under plain equivalence the five profiled matrices `H, H', H″, (H')ᵀ, (H″)ᵀ` are pairwise inequivalent — at least five classes exhibited, a remark and not the count. `Hᵀ` was **not** profiled at 1676, so nothing is said about `H` versus its own transpose there |
| **order 1772 carries at least three equivalence classes** | **PROVEN** (profile invariance + separation implication) **+ PROVEN-BY-CERTIFICATE** (the exact profiles: `H` vs `H'` 57 of 89 differing bins — the `ψ(ρ) = −1` twist's **fourth** proven instance, and the last of cert 02's four decoded `(1,1)` orders; the orientation switch `H″` against them 58 and 53; two independent implementations per matrix; cert 23, 2026-09-03). **Row-side only**: the transposed 1772 profiles are a pending leg of the same campaign (`blas` in, `bits` still running), so nothing is claimed under the transpose-extended relation, and nothing about `H` versus `Hᵀ` there. Cert 23 has no in-repo `--full` leg: one is ≈ 7–8 h on the measured scaling (≈ 10.3 h on the `Θ(n⁵)` law) and the `blas` route wants ≈ 11.1 GB |
| the order-2060 **pair** is inequivalent — the posted matrix against the plain GS realisation of the same seed — so **order 2060 carries at least two equivalence classes** | **PROVEN** (profile invariance + separation implication) **+ PROVEN-BY-CERTIFICATE** (the exact profiles, 146 of 147 differing bins, supports 145 vs 133, two independent implementations per matrix; cert 07 exact mode). Row-side only. The earlier sampled statistic remains recorded as COMPUTATIONAL-EVIDENCE |
| **order 2060 carries at least three equivalence classes** | **PROVEN** (profile invariance + separation implication) **+ PROVEN-BY-CERTIFICATE** (the orientation switch `H″` of the plain GS realisation — unbordered, `s = 0`, twelve `515`-blocks negated — against `P` in 107 of the 145 bins they share and against the posted `G` in 146 of 147; two independent implementations per matrix; cert 22, 2026-09-03). **Row-side only**, as §3.8's 1772 statement is; those two are the only such statements here: the transposed 2060 profiles are pending legs of the same campaign, so nothing is claimed under the transpose-extended relation, and nothing about any of the three versus its own transpose. Cert 22 has no in-repo `--full` leg: one is ≈ 15 h on the measured scaling (≈ 22 h on the `Θ(n⁵)` law) and the `blas` route wants ≈ 17.5 GB |
| the `s ≥ 2` coset-border novelty statement | **BOUNDED-NEGATIVE-SEARCH** (§4; closes exactly the enumerated sources) |
