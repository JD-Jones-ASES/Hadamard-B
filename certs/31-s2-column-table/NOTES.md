# cert 31 — the general-branch border layer at `(s,i) = (7,12)`

**Label: PROVEN-BY-CERTIFICATE.** Default run:
`python certs/31-s2-column-table/run.py` from the repository root. Standard
library only, exact integers and `Fraction`s only, **no floats anywhere**, no
network, nothing imported or opened outside this repository. **1.1 s, 50
checks, exit 0.** There is no `--full`: every enumeration here is exhaustive
on the default path.

The mathematics is `note/NOTE-B.md` §1.12 (Theorem F in the general branch,
the row-factorisation lemma, the dual form) and §2.6.

---

## The setting

`Ḡ = ℤ₁₂`, `i = 12`, `s = 7`, so `4i = 48` and `4s = 28`. The Galois-stable
size-7 supports of `Ẑ₁₂` containing the trivial character are exactly four
(cert 30 [B]):

```
S1 = {0,1,3,5,7,9,11}    S2 = {0,1,4,5,7,8,11}
S3 = {0,1,2,5,7,10,11}   S4 = {0,2,3,4,8,9,10}
```

and `M_S = 4i·P_S` is the only admissible Gram wherever `w > s`
(Theorem E′, §1.7). By Theorem F(a) in the general branch (§1.12), an
**admissible column table** is a `Q ∈ {±1}^{48×28}` with

* **(H1)** `QQᵀ = I₄ ⊗ M_S`,
* **(a)** `QᵀQ = 48·I₂₈`,
* **(blk)** every 12-block of every column lies in `V_S`,

and (a) + (blk) already imply (H1). So the whole vocabulary of the column
layer is the **admissible block alphabet** `B_S = {±1}¹² ∩ V_S`, and it is
small.

## What this certificate establishes

**[A] the alphabet, decided twice.** Over all `2¹² = 4096` sign vectors, by
two membership tests that share no code —

* *cyclotomic divisibility*: `Φ_d(x) | v(x)` for every character order `d`
  outside `S` (legitimate exactly because `S` is Galois-stable, so it is a
  union of full order-classes);
* *the projector fixed point*: `circ(M_S/4)·v = 12·v`, i.e. `P_S v = v`.

The two sets agree at every support, and

| support | `\|B_S\|` | span rank |
| --- | --- | --- |
| `S1` | 66 | 7 |
| `S2` | 32 | 7 |
| `S3` | **2** | **1** |
| `S4` | 24 | 7 |

**[B] the `S3` kill — `(7,12)/S3` admits no column table at any order,
PROVEN-BY-CERTIFICATE.** At `S3` the only
`±1` vectors of length 12 with Fourier support inside `S3` are the **two
constants**. Three independent obstructions follow, and any one of them
suffices:

1. *counting.* Every column is then `(ε₀𝟙, ε₁𝟙, ε₂𝟙, ε₃𝟙)` with
   `ε ∈ {±1}⁴`, so at most `2⁴ = 16` distinct columns exist at all —
   against the `4s = 28` pairwise orthogonal ones required.
2. *rank.* Those 16 columns span a space of dimension **4**, and
   `rank Q = 28`.
3. *Gram.* Inside one superblock any two rows are equal or antipodal, so
   their inner product is `±28`, while (H1) demands `M_{S3}(c)` for `c ≠ 0`,
   whose values are `{−8, −4, 8}` — none of them `±28`.

So **(H1) has no solution at `(7,12)/S3`**: at any `σ`, any `w`, any
`κ(ρ)`, and therefore at **every** order `N = 4(12w + 7)`, not only at 2092.
The `χ₆`-twin support `6 + S3 = {1,4,5,6,7,8,11}` dies with it by the
Lemma-T twist (§1.4) — **checked here, not asserted**: `v ↦ χ₆·v` is a
bijection of `{±1}¹²` carrying `V_{S3}` onto `V_{6+S3}`, so that alphabet is
the two *alternating* vectors, `M_{6+S3}(c) = (−1)^c·M_{S3}(c)` takes no
value `±28` off zero, and the same three obstructions apply. This is the
first general-branch border kill in this repository.

**[C] the artifact at `S2`.** `data/q-7_12-S2.json` holds a `4 × 28 × 12`
array of `±1` blocks. Every block is required to lie in the alphabet computed
in [A]; the assembled `Q` is checked for `QᵀQ = 48 I₂₈` and
`QQᵀ = I₄ ⊗ M_{S2}` entry by entry, with
`M_{S2} = (28, 0, 8, 12, −8, 0, −4, 0, −8, 12, 8, 0)` recomputed here by
Ramanujan sums. So **(H1) has a solution at `(7,12)/S2`.**

**[D] the forced per-superblock multiset.** Writing `R_I` for the `12 × 28`
superblock-`I` slice, (H1) says `R_I R_Iᵀ = circ(M)`, i.e.
`Σ_j m_j b_j b_jᵀ = circ(M)` over the **16 negation classes** of `B_{S2}`
(1 constant, 3 of type A — functions of `c mod 3` — and 12 of type B). The 16
outer products are linearly independent, so together with `Σ_j m_j = 28` the
system has a **unique** rational solution, found here by exact `Fraction`
elimination; it is integral and nonnegative:

> every superblock of every admissible `Q` carries the constant class once,
> each type-A class once, and each type-B class twice — `1 + 3 + 24 = 28`.

The artifact realises exactly that multiset in all four superblocks. Because
the multiset is forced, column permutations and negations put superblock 0
into one canonical arrangement, and a search over superblocks 1–3 is
**exhaustive** — which is what makes a negative at this cell a theorem rather
than a search failure.

**[E] the row-factorisation lemma, concrete.** With `σ` the `(7,12)` table
`B2` of cert 30 — whose support *is* `S2` — and
`Ĉ = GS(σ; κ(ρ))` the `48 × 48` Goethals–Seidel array over `ℤ₁₂` in standard
orientation, this run checks at **every** `κ(ρ) ∈ ℤ₁₂`:

* `ĈᵀĈ = ĈĈᵀ = I₄ ⊗ circ(F)`, `F = 2092·δ₀ − w·M` (Theorem F(c));
* `Ĉ` commutes with `I₄ ⊗ P_S`;
* `W := ĈᵀQ` satisfies `W Wᵀ = 4i·4s·(I₄ ⊗ P_S) = 28·(I₄ ⊗ M)`.

Hence every admissible row `p` has `‖Ĉp‖² = 16is = 1344` **exactly**, so the
corresponding row of `E` has squared norm exactly `4s = 28` whatever `p` is,
and a kit exists at a fixed `(Q, κ(ρ), σ)` **iff** the flat set contains `4s`
mutually orthogonal admissible rows — a condition on each row alone.

**[F] the dual form, on this repository's own kits.**

> `p` is flat **⟺** `Ĉp ∈ Q·{±1}^{4s}`; equivalently a kit is a Hadamard `E`
> of order `4s` with `P = −(1/4s)·E·Qᵀ·Ĉ ∈ {±1}`.

Checked

* at `(2,4)` on **all four** banked records of cert 18, *exhaustively over
  all 256 admissible rows*: `Ĉp` lies in `col(Q)` with `‖Ĉp‖² = 16is` for
  every one; and on the four kits themselves, `Ĉp = −Qe` with `e` the corner
  row, together with `P = −(1/4s)EQᵀĈ`;
* at `(3,4)` on kits found by this certificate's own engine — a fixed
  deterministic list of 40 candidate `Q`s, three S-parts at two `κ(ρ)` each —
  where `E ∈ H(12)`, `E ∈ {±1}`, (H4) holds, and both dual identities hold on
  every kit found. **5 of the 6 points yield a kit**; the sixth is a miss of
  the candidate list, not a statement about that point — this clause is a
  check on the identities, not a census, and cert 17 is where the `(3,4)` kit
  census lives.

## What is **not** claimed

* **No kit at `(7,12)`.** (H1) is the **column** layer only. Whether the `S2`
  table extends to a border kit `(E, P, Q)` is a different question and is
  not answered here; the source laboratory reports **no** border kit anywhere
  in the general branch at 2092, and none is exhibited in this repository.
* Nothing about **seeds** at any `(7,12)` support; nothing about `H(2092)`.
* Nothing about `S1` or `S4` beyond their alphabet sizes and ranks. Their
  column tables exist (the source laboratory constructed them); this
  certificate does not carry them.
* The `S3` kill is a statement about the **column table**, not about the
  σ-layer: `S3`'s σ-table is perfectly valid and is banked in cert 30 as
  `B3_7_12_43`. A valid compressed shell does not imply a valid border.
* Nothing about the *number* of admissible `Q` at `S2`, and nothing about
  which of them are equivalent.

## Controls that can fail

Four negatives and one positive, all through the same acceptance path:

| control | must |
| --- | --- |
| one entry of one block flipped | be rejected — the block leaves the alphabet **and** `QᵀQ = 48I` fails |
| one block replaced by another **admissible** block | be rejected by (H1) while the alphabet test still passes — the two tests are independent |
| two columns permuted inside one superblock only | be rejected by (H1) |
| a support that is not Galois-stable | be rejected outright by the Ramanujan path, not silently evaluated |
| a whole column of `Q` negated | **pass** — a genuine symmetry, so the tests above are not vacuously strict |

## Pinned digests

| file | sha256 |
| --- | --- |
| `data/q-7_12-S2.json` | `e23e1ab73feffa2d948aa0e91777b2799f6a0c51b0b13e495dc25775b13caffb` |
| `data/general-branch-sigma-tables.json` (read by [E]) | `d2c37aace3bf016dc245a26c3fe7a2bd4aa524aafc874ba717c7e562d98c4a97` |
| `data/cell24-records.json` (read by [F]) | `9727b392940d416d3f25dca5d51d2db71cd499bc73c3b8dc4efd22801180f179` |

No matrix is generated and no digest of a matrix is pinned: the `(3,4)` kits
of [F] depend on the candidate-`Q` list, and pinning them would bind the
certificate to a search order rather than to a statement.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`** | **1.07 s** (50 checks, exit 0; measured here 2026-09-05, Python 3.14, one worker) |
| [A], the four alphabets over `2¹²` vectors, twice each | 0.2 s |
| [E], the twelve `κ(ρ)` at `48 × 48` | 0.4 s |
| [F], the `(2,4)` records and the `(3,4)` kit search | 0.4 s |

## Provenance

The `S2` column table was found in the source laboratory on 2026-09-04, after
an earlier depth-first search there had failed to find one and had recorded
`(7,12)/S2` as OPEN under a **BOUNDED-NEGATIVE-SEARCH** label — a search
failure, correctly labelled, which the forced-multiset lemma then made
unnecessary. Nothing about that search is trusted here: a witness certifies
itself, and this run re-derives the alphabet, re-checks every inner product,
and re-solves the multiset system from scratch. See `PROVENANCE.md`.

## How to re-run

```
python verify/verify.py --selftest
python certs/30-general-branch-sigma/run.py
python certs/31-s2-column-table/run.py
```
