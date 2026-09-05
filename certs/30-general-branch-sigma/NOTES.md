# cert 30 — the general-branch σ-layer at `N = 2092`

**Label: PROVEN-BY-CERTIFICATE**, for exactly this: *these six coset-sum
tables exist and carry these profiles, and the support enumerations behind
them are complete.* Default run:
`python certs/30-general-branch-sigma/run.py` from the repository root.
Standard library only, exact integers only, **no floats anywhere**, no
network, nothing imported or opened outside this repository. **0.1 s, 59
checks, exit 0.** There is no `--full`: everything here is small and
everything here is run.

The mathematics is `note/NOTE-B.md` §2.6, and the theorems it leans on are
§1.7 (Theorem E′) and §2.4 (Theorem 3).

---

## The statement

`N = 2092 = 4·523`. A **general-branch cell** is a quotient `Ḡ` of order `i`,
a fibre size `w = |K|` with `N = 4(wi + s)`, and a Galois-stable spectral
support `S ⊆ Ĝ̄` with `|S| = s`. Since `w > 2s` at every cell here, Theorem E
— a fortiori Theorem E′ — makes

> `M_S(c) = 4 Σ_{χ∈S} χ(c) = 4i·P_S`,  `F(c) = N·δ₀(c) − w·M_S(c)`

the **only** admissible Gram: a conclusion, not an ansatz. The **σ-layer** of
the cell is the set of integral coset-sum tables `σ = (σ₀, σ₁, σ₂, σ₃)` on
`Ḡ` with

1. `Σ_q PAF_{σ_q}(c) = F(c)` for every `c ∈ Ḡ`;
2. `σ_q(c) ≡ w (mod 2)` and `|σ_q(c)| ≤ w` — each entry is a sum of `w`
   signs;
3. `Σ_{q,c} σ_q(c)² = F(0)` and `Σ_q (Σ_c σ_q(c))² = 4s`.

**The σ-layer of every advertised general branch at 2092 is non-empty**, by
explicit artifact:

| id | cell `(s,i,w)` | `Ḡ` | support `S` | row sums | budget on `S` |
| --- | --- | --- | --- | --- | --- |
| `A_11_16_32` | (11, 16, 32) | `ℤ₁₆` | orders 1, 4, 16 | (6, 2, 2, 0) | 44 |
| `B1_7_12_43` | (7, 12, 43) | `ℤ₁₂` | orders 1, 4, 12 | (4, 2, 2, 2) | 28 |
| `B2_7_12_43` | (7, 12, 43) | `ℤ₁₂` | orders 1, 3, 12 | (4, 2, 2, 2) | 28 |
| `B3_7_12_43` | (7, 12, 43) | `ℤ₁₂` | orders 1, 6, 12 | (4, 2, 2, 2) | 28 |
| `B4_7_12_43` | (7, 12, 43) | `ℤ₁₂` | orders 1, 3, 4, 6 | (4, 2, 2, 2) | 28 |
| `C_3_8_65` | (3, 8, 65) | `𝔽₂³` | `{0,1,2}` binary | (0, −2, −2, −2) | 12 |

and the Fourier budget is `N = 2092` off `S` in every one — the two-valued
Parseval law of Corollary E2, here at general cells.

**The support enumerations are complete.** Exactly **one** Galois-stable
size-11 subset of `Ẑ₁₆` contains the trivial character — table A's — and
exactly **four** Galois-stable size-7 subsets of `Ẑ₁₂` do — B1…B4's. So the
`(7,12)` cell is populated on *every one of its rational branches*, not on a
lucky one. (`Ẑ₁₂` has Galois orbits of sizes 1, 1, 2, 2, 2, 4; exactly eight
unions total 7, four containing `χ₀` and four their `χ₆`-twists.)

**The `(3,8)` cell is not an independent door.** All **56** size-3 subsets of
`𝔽̂₂³` and all **12** Galois-stable size-3 subsets of `Ẑ₄ × Ẑ₂` are
real-character twists of `W ∖ {1}` for a subgroup `W` of order 4, so every
such cell coarsens onto a house `(3,4)` profile. Exhibited concretely: table
`C` twisted by the real character 3 and compressed along the line
`J = {0,4}` is `[[−8,−8,−8,−8],[−4,−2,−2,−2],[4,6,6,6],[4,6,6,6]]`, with
aggregate PAF `[532, 520, 520, 520]` — the house `(3,4)` compressed profile
`[4w+12, 4w, 4w, 4w]` at `w = 130`, which is the width of `N = 2092` in that
cell (§1.8).

## What is **not** claimed

* **Nothing about seeds.** A compressed σ-table is a *necessary* shadow of a
  seed quadruple on a group of order `wi`, never a sufficient one. **No seed
  quadruple is known at any of these cells**, and none is exhibited here.
* **Nothing about borders.** That is cert 31, and its answer at `(7,12)` is
  partial: an admissible column table at `S2`, a kill at `S3`, and no border
  kit anywhere.
* **Nothing about `H(2092)`.** Neither existence nor non-existence.
* **Nothing about the number of tables** in any σ-layer. One per support is
  exhibited; the census is not attempted and would be astronomically large.
* Nothing about cells this file does not list. The three cells are the ones
  the source laboratory advertised; that they are the only general-branch
  cells at 2092 is **not** claimed here.

## What `run.py` checks (exit 0 iff every check passes)

**[0]** `data/general-branch-sigma-tables.json`, SHA-256 pinned in `run.py`.

**[A] the six tables.** For each: `M_S` and `F` are **re-derived from
`(S, i, w)` alone** and only then compared with the listed ones; the
aggregate PAF is computed twice and must equal `F`; parity, box, row sums,
total norm `Σ σ² = F(0)`, the row-sum shell `Σ_q r_q² = 4s`; and the Fourier
budget, computed by a third route, must be `4s` on `S` and `N` off it.

**[B] support uniqueness.** Complete enumeration of Galois-stable supports of
the given size containing `χ₀`, on `Ẑ₁₆` and `Ẑ₁₂`, matched against the
banked ones.

**[C] the `(3,8)` collapse and the transport.** The 56 + 12 size-3 supports;
the twisted, compressed table; its aggregate PAF against the house profile.

**[D] the theorems the cells sit under.** `N = 4(wi + s)`; `w > 2s`;
`M/(4i)` idempotent with trace `4is` (so the spectrum is `{0^{i−s}, (4i)^s}`,
Theorem E); and `S` contains a real character, as Theorem 3 (§2.4) requires
at `N ≡ 4 (mod 8)`.

**[E] controls that can fail.** Four perturbed tables through the identical
acceptance path — one entry moved by `+2`, one support element swapped
(breaking Galois stability), a row sum falsified, `w` changed to 30 — each
must be **rejected**, and the run reports which condition caught it. The
unperturbed table is then run again through the same path and must pass.

## Two implementations, and why a bad support is rejected rather than evaluated

| quantity | implementation A | implementation B |
| --- | --- | --- |
| `PAF_σ` | direct double sum `Σ_j σ_j σ_{j+t}` | row 0 of `G Gᵀ` for the group matrix `G[a][b] = σ(b−a)` (or `σ(a⊕b)`) — a matrix is built and multiplied, with no autocorrelation formula anywhere on that path |
| `M_S` on `ℤ_i` | root-of-unity sums reduced modulo `Φ_i(x)` | **Ramanujan sums** `c_d(c) = μ(d/g)·φ(d)/φ(d/g)`, `g = gcd(c,d)`, over the order-classes of `S` |
| `M_S` on `𝔽₂³` | direct character sum | a fast Walsh–Hadamard butterfly on `1_S` |
| the aggregate | `Σ_q PAF` against `F` | the Fourier budget, `4s` on `S` and `N` off it |

A support that is **not** Galois-stable makes implementation A's `M_S`
irrational — the reduction mod `Φ_i` leaves a remainder of degree `> 0` — and
makes implementation B's order-class decomposition fail outright. So it is
**rejected**, not silently evaluated at some rounded value. That is what the
`CTRL-support-swapped` control exercises.

## Pinned digests

| file | sha256 |
| --- | --- |
| `data/general-branch-sigma-tables.json` | `d2c37aace3bf016dc245a26c3fe7a2bd4aa524aafc874ba717c7e562d98c4a97` |

Nothing else is pinned: this certificate generates no matrix, and every
number in it is recomputed on every run.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`** | **0.10 s** (59 checks, exit 0; measured here 2026-09-05, Python 3.14, one worker) |

## Provenance

The six tables were produced by the **Sol station**, an external reviewer, in
a consult reply of 2026-09-03, and were banked verbatim in the source
laboratory before any checking code existed. Everything asserted about them
here is re-derived from the support `S` alone; the station's own `M` and `F`
vectors are read only to be compared against what this run computed, and its
conclusions are not inputs to any decision here. See `PROVENANCE.md`.

## How to re-run

```
python verify/verify.py --selftest
python certs/30-general-branch-sigma/run.py
python certs/30-general-branch-sigma/run.py --verbose
```
