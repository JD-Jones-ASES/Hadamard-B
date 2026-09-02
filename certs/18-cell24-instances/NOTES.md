# cert 18 — the even-`s` branch: four `(2,4)` matrices, the border census, Theorem 3

**Label: PROVEN-BY-CERTIFICATE** (the four matrices, the two structural
claims, the border census, the Theorem-3 exhaustion) **+ PROVEN** (the border
proposition (a),(b) and Theorem 3 are paper-grade proofs,
`note/NOTE-B.md` §2.4). Default run:
`python certs/18-cell24-instances/run.py` from the repository root. Standard
library only, exact integers only. **5.8 s, 97 checks, exit 0.** Nothing here
is audited from a bank: every claim is recomputed in the run. There is no
`--full` — the census is exhaustive on the default path.

---

## The cell

Theorem E's **general branch** (`note/NOTE-B.md` §1.2.1) admits `Ḡ = ℤ₄`,
`S = {χ, χ³}` (the two faithful characters), `s = 2`, and the Gram

```
M = 4·(χ + χ³) = (8, 0, −8, 0).
```

`s` is **even**, so the cell lies outside the house branch that Theorem C
classifies. (H2) reads `Σ_q PAF_q = 4n·δ₀ − 8·[K∖0] + 0·[odd cosets] +
8·[κ⁻¹(2)]`.

## The four matrices — `note/NOTE-B.md` §2.3

| record | `n` | `w` | `ρ` | order | canonical SHA-256 |
| --- | --- | --- | --- | --- | --- |
| `H88_cell24_n20_seed0_rho0` | 20 | 5 | 0 | **88** | `942b3f32fcd75e72a64f92d9c294b0d0cedbbd0965fe5a14213c30b8b66ffc8a` |
| `H88_cell24_n20_seed1_rho1` | 20 | 5 | 1 | **88** | `4cae47d1c5054a86ca48154c5e9cd99845294be4d275dd8560bf6b05fe5f08e7` |
| `H56_cell24_n12_seed0_rho0` | 12 | 3 | 0 | **56** | `ad67ee2c9d1f4d0343b250824dec301759e097d72d949f1b9c91f22cab026b85` |
| `H56_cell24_n12_seed2_rho3` | 12 | 3 | 3 | **56** | `fa82808c8dbc0245f3d427312975183fb947438bdcafba56107ee2280d6e4aff` |

Each digest is pinned **three ways** and all three must agree: computed in
this run from the assembled rows, reported by `verify/verify.py` in its
verdict line, and carried by `data/cell24-records.json`.

These are the first coset-border Hadamard matrices with **even `s`** known to
this laboratory. **The orders were never open, and no novelty of existence is
claimed at 56 or 88** — what is new here is the construction type.

## The data file — `data/cell24-records.json`

SHA-256 `9727b392940d416d3f25dca5d51d2db71cd499bc73c3b8dc4efd22801180f179`,
6 198 bytes; pinned in `run.py`. Schema `cell24-records/1`, one entry per
matrix:

| field | meaning |
| --- | --- |
| `name` | the record's name (also the temporary file name during verification) |
| `order`, `group`, `s`, `coset_divisors`, `w` | `N = 4(n+s)`, `G = ℤ_n`, `s = 2`, `[4]` (so `K = 4ℤ_n`, `i = 4`), `w = n/4` |
| `r_shift`, `rho_bar` | the Goethals–Seidel reflection `ρ ∈ G`, and `κ(ρ) ∈ Ḡ = ℤ₄` |
| `gs_variant` | `"standard"` — the orientation of `note/NOTE-B.md` §1.0 |
| `seeds` | four `±1` sequences on `G`, tape order `g = 0 … n−1` |
| `coset_sums` | `σ_q(c) = Σ_{g ≡ c (4)} x_q(g)` |
| `tau` | the S-part `τ_q = (σ_q(0) − σ_q(2), σ_q(1) − σ_q(3))` |
| `col_table_8` | `Q′`, the `8×8` Hadamard half of the column table |
| `row_table_8` | `P′`, the `8×8` Hadamard half of the row table |
| `corner` | `E` as produced upstream |
| `canonical_sha256` | the digest `verify/verify.py` reports |

**Everything except `seeds`, `row_table_8`, `col_table_8` and `r_shift` is
redundant, and is there to be checked rather than believed.** `coset_sums`
and `tau` are recomputed from `seeds`; `corner` is **not used** — `E` is
recomputed here as `−(1/16)·P·Ĉᵀ·Q` from this certificate's own
Goethals–Seidel array of the coset sums over `ℤ₄` (the compression lemma) and
required to equal the banked value; `canonical_sha256` is compared against
the digest of the matrix this run assembles.

## What `run.py` checks (exit 0 iff every check passes)

**[0]** the data file's SHA-256; the cell declaration; four records with the
expected names.

**[1] the four instances**, each: the record's declared parameters are
self-consistent; the coset sums and the S-part recomputed from the seeds
equal the banked ones; (H2) lag by lag; Parseval at the live characters
(`Σ_q |σ̂_q(χ)|² = N − w·M̂(χ) = 8`); the row-sum law (`Σ_q r_q² = N`);
`Q′, P′ ∈ H(8)`; the anti-periodic doubling to `Q` (16×8) and `P` (8×16);
(H1) `QQᵀ = I₄⊗M`; `PPᵀ = 16I`; **`E` recomputed** and matched; `E ∈ H(8)`;
(H3) `EEᵀ + wPPᵀ = NI`; (H4) `EQᵀ + PĈᵀ = 0`; the Σ̄-law
`ĈĈᵀ = ĈᵀĈ = I₄ ⊗ dev(Σ̄)` (`Σ̄ = [48,0,40,0]` at `n = 20`, `[32,0,24,0]` at
`n = 12`); assembly by this certificate's own block-explicit assembler;
`verify/verify.py` exit 0 with `HADAMARD order=N`; the digest pinned three
ways. The matrix file is deleted afterwards; nothing is committed.

**[2] the two structural claims**, with controls run first. (i) *Not a
Kronecker product.* If `H = D₁Π₁(H₂ ⊗ H′)Π₂D₂` then the pointwise product of
two rows of `H` is `±` a column permutation of a pointwise product of two
rows of `H₂ ⊗ H′` (row signs give the global sign; column signs cancel in
products). In `H₂ ⊗ H′` the rows `[h|h]` and `[h|−h]` pair off into `N/2`
disjoint pairs whose pointwise product is the same vector. So anything
equivalent to a Kronecker product has a sign-normalised pointwise-product
class of size `≥ N/2`. All four instances have largest class **4** (against
`N/2 = 44` and `28`). Controls: Sylvester `H(16)` gives 8 `≥` 8 (the test
*passes* a genuine Kronecker product); a scrambled Hadamard-equivalent copy
gives the same 8 (the statistic really is an invariant); Paley `H(12)`, which
admits no Kronecker factorisation, gives 1. (ii) *Not a collapse.* The common
kernel of `S = {χ, χ³}` in `Ḡ = ℤ₄` is `{0}`, so `S` generates the dual and
the instance does not descend to a smaller-index construction.

**[3] Theorem 3, small cases.** For every abelian `Ḡ` of order `i ≤ 16` and
every Galois-stable `S` with no real character: `|S|` is even in every case
(the theorem's first clause). Where `4 ∤ i` and `s ≤ 2` (the DFS budget:
`4s ≤ 8` columns), an **exhaustive** depth-first search for one block
`Q₀ ∈ {±1}^{i×4s}` with `Q₀Q₀ᵀ = M_S` finds nothing — **9 cases, all empty**:
`i = 3` (`ℤ₃`), `i = 6` (`ℤ₆`, ×2), `i = 9` (`ℤ₉`, and `ℤ₃²` ×4), `i = 15`
(`ℤ₁₅`). At `i = 4`, where `4 | i` and the theorem forbids nothing, the same
search **does** find `Q₀` in 5 nodes — the `(2,4)` cell, whose instances this
certificate builds.

**[4] the `(2,4)` border census, exhaustive.**

* The **480** right-orbit representatives of `H(8)` are built
  deterministically from the **30** affine `AG(3,2)` structures on eight
  labels (one per permutation class, enumerated over the `7!` permutations
  fixing label 0) times `2⁴ = 16` plane-sign classes. They are checked to be
  Hadamard, pairwise distinct as right-orbits (480 canonical forms), and to
  catch 3 000 deterministically generated labelled `H(8)`. The count **480**
  is re-derived independently: a backtracking count of **normalised** `H(8)`
  (first row and column all `+`) gives **151 200**, so the labelled total is
  `151 200 · 2¹⁵ = 4 954 521 600` and the orbit count is
  `T/(8!·2⁸) = 480`. Two routes to the same number.
* The **112** S-parts are enumerated from Parseval (`Σ_q |τ_q|² = 4s = 8`)
  and parity, both signs kept; with `κ(ρ)` that is **448** classes per `Q′`,
  and **112 × 4 × 480 = 215 040** classes in all.
* For each class the admissible anti-periodic rows are found and an
  **8-clique in their orthogonality graph** is searched for. **215 040 /
  215 040** classes admit a kit — the border is never the obstruction at
  `(2,4)`, at any `w`. The histogram of admissible rows per class is
  `{8: 114 688, 16: 86 016, 32: 14 336}`.
* **200** kits, chosen by a fixed deterministic stride, are re-verified in
  exact integers with a **full `σ`** — an arbitrary `S^c` part (a constant
  plus `b·(−1)^c` per seed) added to the S-part — at **both** `w = 5` and
  `w = 130`: (H1), `E ∈ H(8)`, `PPᵀ = 16I`, (H4), and (H3) at each width.
  And (H3) is checked to be **`w`-free**: it reads `8 + 16w = 4(4w+2)`, an
  identity, which is why one kit serves every order.

### How the admissible rows are characterised

`p` (anti-periodic, `p[0] = +1`) is admissible for `Q′` iff
`w := a(p)·Q′ ∈ {±16}⁸`, where `a(p)` is the 8-vector
`a[2I+j] = u[4I+j] − u[4I+j+2]` and `u = Ĉ⁰p`. Multiplying on the right by
`Q′ᵀ` and using `Q′Q′ᵀ = 8I`, that condition is **exactly**
`a(p) = 2·v·Q′ᵀ` for some `v ∈ {±1}⁸`, and then `E`'s row is `−v`. So the
admissible `a` for a given `Q′` form a set of 256 vectors, computed once per
`Q′` and looked up 128 times per class. This is what makes the exhaustive
census cost 3.6 s of pure Python instead of hours, and it is an identity, not
an approximation: the sampled re-verification multiplies the resulting kits
out in full.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **5.8 s** (97 checks, exit 0; measured here 2026-09-02, Python 3.14, one worker) |
| clauses [0]–[3] (the four matrices, structure, Theorem 3) | 0.5 s |
| clause [4]: the normalised-`H(8)` count | 0.6 s |
| clause [4]: the 215 040-class census | 3.6 s |
| clause [4]: 200 exact re-verifications at two widths | 0.6 s |

## What is NOT claimed

* **Nothing about `H(2092)`.** The `(2,4)` cell does not land at 2092:
  `4(4w+2) = 2092` has no integer solution. These are existence witnesses for
  the even-`s` branch, at orders that were never open.
* **No novelty of the orders.** Only the construction type is new to this
  laboratory (a coset border with even `s`); whether it is new to the
  literature is not searched here.
* **Nothing about Hadamard equivalence** between the four, or to other known
  `H(56)`/`H(88)`. Only "not equivalent to a Kronecker product with `H(2)`"
  is proved, and by one invariant.
* **The anti-periodic row table is not shown FORCED at `w ≤ 4`.** Border
  proposition (b) is a necessity statement and uses `w > 2s = 4`. The two
  `H(56)` instances have `w = 3`, where the anti-periodic kit is *sufficient*
  — it is what they use — but not shown necessary.
* The census settles the **border** layer only. It says nothing about which
  `(2,4)` seed quadruples exist at which `n`.
* Theorem 3's small-case run is an **exhaustion of 9 cases**, not a proof;
  the proof is `note/NOTE-B.md` §2.4. Cases with `s > 2` are outside the DFS
  budget (`4s > 8` columns) and are not run.
* The `S`-part census counts `τ` and `−τ` separately (112 S-parts, 448
  classes with `κ(ρ)`). A census that identifies them would report half as
  many classes for the same statement.

## How to re-run

```
python verify/verify.py --selftest
python certs/18-cell24-instances/run.py
```
