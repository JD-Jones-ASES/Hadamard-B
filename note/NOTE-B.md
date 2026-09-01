# NOTE-B — Bordered Goethals–Seidel arrays: theorems, instances, and a separation at order 668

This note is the mathematics of the repository, in three movements:

- **I. The theorems** (§1) — an exact characterisation of when a
  Goethals–Seidel array over a finite abelian group extends to a
  Hadamard matrix through a *coset border* of width `4s`; the forced
  parameter classification; and the complete resolution of the
  `s = 1` border system.
- **II. The instances** (§2) — twelve publicly posted matrices
  re-verified through the theorem's hypotheses, and five matrices
  constructed here, every one replayable from `certs/`.
- **III. Existence plus separation** (§3) — at order 668, two
  Hadamard matrices that are **provably inequivalent**, separated by
  an exactly computed invariant; and the invariant theory that makes
  the separation honest, including two documented traps.

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
classical back-diagonal and `ρ = 0` is group inversion. When `|G|` is
even, `ρ` is **not** removable by relabelling (a shift `g ↦ g+α` sends
`ρ ↦ ρ−2α` and `2` is not invertible).

**The core.** `C = GS(x_0,x_1,x_2,x_3; ρ)` is the `4G × 4G` array

```
       [  A     BR     CR     DR  ]
       [ -BR     A    DᵀR   -CᵀR  ]        A = dev(x_0), B = dev(x_1),
       [ -CR   -DᵀR    A     BᵀR  ]        C = dev(x_2), D = dev(x_3)
       [ -DR    CᵀR  -BᵀR     A   ]
```

We call this orientation *standard*; negating the six transposed
blocks gives the other valid orientation (used, e.g., by SageMath).
The proofs below are orientation-symmetric, but **every theorem in
this note is stated, and every checker in this repository assembles,
the standard orientation only**; the analogues for the other
orientation are routine and are not derived here.

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
- **(H2)** `Σ PAF(t) = 4n·δ_{t,0} − M(κ(t))` for all `t ≠ 0`, and
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
(c) `X R Yᵀ = Y R Xᵀ` (both equal `f(g+h)` with
`f(t) = Σ_u x(t−u)y(u)`);
(d) `X R Y = Yᵀ R Xᵀ` (both equal `Σ_u x(ρ−g−h+u)y(u)`).
Each is a one-line change of summation variable; (c) and (d) are
where `R` being the *reflection* `k+h = ρ` (rather than an arbitrary
involution) is used. ∎

**Lemma 2 (core Gram).** `C Cᵀ = I_4 ⊗ Σ`, where
`Σ[g,h] = Σ PAF(g−h)`.

*Proof.* Diagonal blocks: `R Rᵀ = I` collapses every reversed block,
and `XᵀX = XXᵀ` (Lemma 1b) makes all four diagonal blocks equal
`Σ_q X_q X_qᵀ = Σ`. Off-diagonal blocks: each is a sum of four
products which pair up as `−A R Bᵀ + B R Aᵀ` (zero by 1c),
`CD − DC` (zero by 1a), `B Cᵀ − Cᵀ B` (zero by 1b),
`−A R D + Dᵀ R Aᵀ` (zero by 1d) — one such pair for each of the six
off-diagonal block positions. Nothing but Lemma 1 is used, so **no
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
equivalent to (H1)+(H2). This also shows the **two-tier shape of the
profile is forced by the ansatz**: constant on `K∖{0}` because the
diagonal of a `±1` Gram is constant, and constant on each other coset
because `Q Qᵀ` cannot see anything finer than a coset. This is the
sentence the whole construction turns on.

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
matrices `M`. What is known: under `w > 2s`, any admissible `M` is
PSD, `Ḡ`-invariant, has diagonal `4s` and rank exactly `s`, and its
nonzero eigenvalues lie in `(4i−2s+2, 4i+2)` and sum to `4si`; the
house value (all equal to `4i`) is the equality case of that
averaging bound. At `(s,i) = (1,2)` the house form **is** forced
(Theorem D below). For `s ≥ 2` forcedness is **open**, and it is
open *up to a character twist at best*: a character twist `x_q ↦ ψx_q`
with `ψ² = 1`, `ψ|_K = 1` produces valid instances whose Gram is
non-house (an explicit `N = 1916` witness is banked; §3.3 explains
why, when `ψ(ρ) = 1`, that witness is a diagonal conjugation of the
house instance and not a new matrix). Any exhaustive treatment of
admissible profiles at `i ≥ 4` must therefore quotient by the twist.
Label of "`M` house is forced, `s ≥ 2`": **CONJECTURE**.

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

So the twist **never manufactures a new matrix when `ψ(ρ) = 1`**,
however non-house its Gram. When `ψ(ρ) = −1` the argument breaks:
the conjugation produces the array with all twelve off-diagonal
blocks negated, and that sign pattern is not of rank one
(`s_I s_J = −1` for all `I ≠ J` is impossible since
`(s₁s₂)(s₁s₃)(s₂s₃) = (s₁s₂s₃)² = +1`). At order 668 — where
`G = ℤ₁₆₆`, `ρ = −1`, and the unique `ψ` has `ψ(ρ) = −1` — the
twisted matrix is in fact **provably a new equivalence class**
(§3.4). Whether a twist with `ψ(ρ) = −1` always leaves the
equivalence class is not claimed in general; 668 is the one proven
case, and the sibling orders 716, 1676, 1772 remain unclaimed.

### 1.5 Theorem D (the `s = 1, i = 2` border system, resolved) — PROVEN

Setting: `G` abelian of even order `n`; `K ≤ G` of index `2`;
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
construction written with `i = 2` bookkeeping). So the house-`M`
question of §1.2 is **closed at `(s,i) = (1,2)`**; it stays open for
`s ≥ 2`.

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
row-sum vector. Every admissible `d` (all 24 with `Σd² = 4`) admits
a border: exhaustively, over all `768²` ordered pairs of 4×4
Hadamard `(p,U)`, every `±2e_j` admits exactly `½` of the pairs and
every `(±1,±1,±1,±1)` exactly `¾` (mechanism: an orthogonal basis of
`{±1}⁴` is one of the two weight-parity cosets of `F₂³`, and `E` is
`±1` iff the row-classes of `p·Λ(d)ᵀ/2` and the column-classes of
`U` are the different coset). And the `i = 1` universal border
transports verbatim: with `P[r][2J+c] = (−1)^c·P₁[r][J]`,
`Q[2I] = Q₁[I]`, `Q[2I+1] = −Q₁[I]`, `E = E₁`, an `i = 1` border is
a valid `i = 2` border whenever `d = (2,0,0,0)`.

**One sentence:** *at `s = 1` the `i = 1` and `i = 2` border systems
are the same 4×4 system `E Uᵀ + p Λ(·)ᵀ = 0`; only the argument
changes, from the row-sum vector `r` to the twisted vector `d`.*

Theorem D is stated for the standard orientation only, and its
sufficiency direction inherits from Theorem A. Every clause is
machine-validated on the four `i = 2` instances of §2, and the
theorem is *gated*: a brand-new `s = 1, i = 2` instance on the
non-cyclic group `G = ℤ₂×ℤ₂×ℤ₃` (the `ε = +1` branch, which the
four instances do not exercise) was built from scratch and verifies
green — `H(52)`, cert 04. Label: PROVEN (paper-grade) + MEASURED
(4/4 clauses) + PROVEN-BY-CERTIFICATE (the gate).

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
`i = 2` cell is a relabelling of `i = 1` whenever the index-2
subgroup is unique, and its border system is the `i = 1` system in
different coordinates. The underlying Parseval mechanism is that the
trivial character's constraint is tight at `i = 1` while at `i = 2`
the binding constraint is the character trivial on `K` but not on
`G`, which contributes `Σ_q δ_q² = 4n − 4(s+1)w + 4s`.)

---

## 2. Movement II — the instances

Everything in this section is replayable: each cert directory under
`certs/` rebuilds its matrices from the small banked data in
`data/`, re-checks the theorem hypotheses in exact stdlib integer
arithmetic, hands every matrix to `verify/verify.py` (the trust
chain), and compares the canonical SHA-256 against the digest pinned
in its `NOTES.md`. Large matrices are regenerated, never committed.

### 2.1 The twelve public records — PROVEN-BY-CERTIFICATE (replay)

Twelve Hadamard matrices at orders 668, 716, 892, 1132, 1244, 1388,
1436, 1676, 1772, 1916, 1948 and 1964 were posted publicly on
2026-08-12 (the sign-stream; expanded matrices on GitHub from
2026-08-13; see `PROVENANCE.md`). The parameter records banked in
`data/payload-records.json` — group, seeds, reflection, coset
divisors, and border tables — were **decoded here from those public
artifacts**; they are the posting team's mathematical content, not
ours, and no priority of any kind is claimed over them or over
anything derived from them. What this repository adds at these
orders is the verification: cert 01 checks every hypothesis
(H1)–(H4) and every derived law of Theorem C literally on each
record, assembles the matrix, verifies it, and matches the canonical
digest — **12/12 green**. The records saturate the classification of
§1.3: five at `(0,1)`, four at `(1,1)`, and one each at `(3,4)`
(order 1916), `(5,6)` (order 1388), `(7,8)` (order 1436) — corners
of sizes 12, 20, 28 respectively, as D3 forces.

### 2.2 The matrices constructed here

| cert | artifact(s) | what it instantiates |
| --- | --- | --- |
| 02 | `H(668)`, `H(716)`, `H(1676)`, `H(1772)` in the `i = 2` frame | Lemma T (§1.4): the four `s=1` records twisted and re-bordered at `i = 2`; each verifies green and is a different artifact from the decoded one. At 668 the stronger statement is §3.4's theorem. |
| 03 | `H(52)` on `G = ℤ₂×ℤ₂×ℤ₃` | Theorem D's gate: a from-scratch `s=1, i=2` instance on a non-cyclic group, `ε = +1` branch. |
| 04 | `H(76)` on a non-cyclic group | a bordered-GS instance on non-cyclic `G` with non-scalar `M`, built as a search gate; the first such the lab located anywhere. |
| 05 | two `H(20)` instances | the hypothesis boundary: bordered instances outside `w > 2s`, where D3's forcing fails and the corner need not be Hadamard — kept as the documented edge of Theorem C. |

No novelty of existence is claimed for any order in this table —
every order here is long settled. The artifacts instantiate the
theorems; that is their entire role. (The four `i = 2` matrices of
cert 02 can also be produced by Theorem D's transport (D-e); the two
assemblies agree byte for byte.)

---

## 3. Movement III — existence plus separation

**Definition.** `H ~ H'` iff `H' = D_r P_r H P_c D_c` with `P`
permutation matrices and `D` diagonal `±1`. Transpose is **not** in
the group; every invariant below is computed on both sides (rows and
columns), and where transpose changes a verdict it is called out.

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
  > a `7σ` false separation (20 bins over `4σ`); folded, the same
  > data reads exactly zero. The signed statistic is uncalibrated in
  > both directions (it also reads `36.6σ` on the 668 pair where the
  > folded truth is `1.8σ`). Every histogram in this repository is
  > folded at measurement.
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
the column signs; degree-2 statistics die by orthogonality, and
every aggregate degree-4 statistic collapses to a function of `n`
(the lemma above is the clean form). The first statistic with
resolution is the 4-profile histogram itself, at `Θ(n⁴)`. On the 668
pair below, *every* cheap invariant — I1–I4, I6, `dim W`, the dual
weight enumerator — returns identical values, and `2·10⁷` samples of
I5 cannot tell the pair apart (`max |z| = 1.8`); only the exact
4-profile separates. Matching invariants are never evidence of
equivalence, and this pair is the proof.

### 3.3 Order 1916: the twist that is not a new matrix — PROVEN-BY-CERTIFICATE

The character-twisted sibling of the `(3,4)` record at order 1916
(the non-house-Gram witness of §1.2) satisfies
`H_twisted = S · H_house · S` with
`S = diag(I_{4s}, I_4 ⊗ diag(ψ(g)))` — 952 rows and the same 952
columns negated, no permutation. All `1916²` cells are
machine-checked. This is the `ψ(ρ) = 1` case of §1.4's proposition
doing exactly what it promises: the twist refutes Gram-forcing as a
literal statement while manufacturing nothing new.

### 3.4 The theorem at order 668 — PROVEN

> **Theorem.** The two Hadamard matrices of order 668 banked here —
> the decoded `(1,1)` record and its Lemma-T `i = 2` rebuild — are
> Hadamard-**inequivalent**.

The separating computation is the exact 4-profile: all
`C(668,4) = 8 222 179 035` 4-subsets, on each matrix, by two
independent implementations (a float32 BLAS Gram of the pair-vector
matrix, and a packed-`uint64` popcount path; opposite bit packings;
both validated bin-for-bin against straight `O(C(n,4))` stdlib
enumeration on five small control matrices). Both matrices populate
the **same 80 bins** — the support does not separate them, and
neither does the extreme tail (top bins agree exactly). **The bulk
separates them: 26 of the 80 bin counts differ**, e.g.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| decoded | 2 073 064 058 | 1 852 054 148 | 1 491 070 735 | 1 091 442 371 | 732 009 734 | 452 971 620 | 260 220 030 | 139 599 212 |
| rebuilt | 2 073 109 602 | 1 852 009 274 | 1 491 079 303 | 1 091 478 493 | 731 928 178 | 452 961 444 | 260 257 738 | 139 618 428 |

Both profiles total `C(668,4)` exactly, the differences sum to zero,
and both hit the second-moment identity `5 517 193 410 096` to the
unit; the first moment — which nothing forces — differs. The largest
bin discrepancy is `1.1·10⁻⁴` of its bin, invisible to any sample of
practical size (a null sampled comparison at `2·10⁷` draws was
already in hand and is worth very little; the exact computation is
what settles it). The separation also survives adding transpose to
the group.

Consequences. Order 668 carries at least two Hadamard equivalence
classes among the artifacts banked here; the Lemma-T construction at
`ψ(ρ) = −1` genuinely leaves the equivalence class at this order.
The corresponding statement at 716, 1676 and 1772 is **not claimed**
(the same computation costs ≈1.4×, ≈98× and ≈130× the 668 run and
has not been made).

### 3.5 The pair at order 2060 — COMPUTATIONAL-EVIDENCE

The formerly-open order 2060 carries two candidate classes among
public+banked artifacts: the publicly posted matrix and the plain
GS-array realisation over the same decoded seed (the `×104`
character twist relates the two seeds). Every exact invariant we
computed agrees between them; the **sampled** 4-profiles differ
systematically and reproducibly (two seeds, both sides, `4–7σ` in a
coherent monotone pattern across consecutive bins; the shape
difference — more peaked, heavier-tailed — balances as it must), and
no equivalence exists in the block-affine family
(BOUNDED-NEGATIVE-SEARCH: exhaustive over that family, silent beyond
it). Label: **COMPUTATIONAL-EVIDENCE of inequivalence — not a
proof**, and this note does not write it as one. The exact 4-profile
at 2060 is priced (≈11–22 core-hours, memory-aware enumeration) and
is the named upgrade path.

---

## 4. Prior art and credit

*[SECTION FINALIZED AFTER SOURCE ADJUDICATION — placeholder pending
the firsthand-read lanes of 2026-08-31; see intel in the lab repo.]*

---

## 5. Labels

| claim | label |
| --- | --- |
| Theorem A (iff); Theorem B; Lemmas 1–3; Theorem C + classification | **PROVEN** (paper-grade; not machine-checked) |
| Theorem D (D-a … D-e); the index-2 collapse corollary | **PROVEN** (paper-grade) + **MEASURED** (every clause machine-validated 4/4) |
| Lemma T; the `ψ(ρ)=1` conjugation proposition | **PROVEN** |
| the twelve public records satisfy every hypothesis | **MEASURED** (cert 01, 12/12) |
| the twelve assembled matrices are Hadamard and match their pinned digests | **PROVEN-BY-CERTIFICATE** (cert 01) |
| the five matrices constructed here (certs 02–05) | **PROVEN-BY-CERTIFICATE** |
| `M` house form forced at `(s,i) = (1,2)` | **PROVEN** (Theorem D-a) |
| `M` house form forced for `s ≥ 2` (up to twist) | **CONJECTURE** |
| the 1916 pair is equivalent (explicit witness) | **PROVEN-BY-CERTIFICATE** |
| **the 668 pair is inequivalent** | **PROVEN** (exact 4-profile, two independent implementations; a third implementation exists in the lab record) |
| the 2060 pair is inequivalent | **COMPUTATIONAL-EVIDENCE** (sampled profiles; block-affine family exhausted) |
| the `s ≥ 2` coset-border novelty statement | **BOUNDED-NEGATIVE-SEARCH** (§4; closes exactly the enumerated sources) |
