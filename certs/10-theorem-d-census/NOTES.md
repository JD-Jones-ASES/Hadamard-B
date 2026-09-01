# cert 10 — the Theorem-D census: `768²` border pairs, and the compressed-block identities

## The claim

> **The census.** Let `p, U` range over the 4×4 Hadamard matrices and
> let `d ∈ ℤ⁴` satisfy `Σ_q d_q² = 4`. Over all `768²` ordered pairs
> `(p, U)`, the corner `E = −¼·p·Λ(d)ᵀ·U` of NOTE-B.md §1.5 **(D-d)**
> has all entries `±1`
>
> - for exactly `768²/2 = 294 912` pairs, for each of the eight even
>   arguments `d = ±2e_j`;
> - for exactly `3·768²/4 = 442 368` pairs, for each of the sixteen odd
>   arguments `d = (±1,±1,±1,±1)`.
>
> Every one of the 24 arguments therefore admits a border, and the
> **(D-e)** mechanism predicts the outcome of every single one of the
> `24 × 768² = 14 155 776` triples correctly.
>
> **The identities.** Over `Ḡ = ℤ₂`, for both `ε = +1` and `ε = −1`,
> the compressed Goethals–Seidel array `Ĉ` satisfies
>
> ```
> Ĉ[(I,1),(J,c′)] − Ĉ[(I,0),(J,c′)] = −(−1)^{c′}·Λ(d)[I][J]
> Ĉ[(I,0),(J,0)] − Ĉ[(I,0),(J,1)]  =              Λ(d)[I][J]
> Ĉ[(I,c),(J,0)] + Ĉ[(I,c),(J,1)]  =              Λ(r)[I][J]
> ```
>
> identically in the eight coset sums — the row-difference and row-sum
> tables that **(D-a′)** turns into the collapse of the degenerate
> branch, and the column-difference table that **(D-d)** turns into the
> closed form.

## Why this artifact exists

NOTE-B.md §1.5 **(D-e)** reported the census as a source-lab
measurement, with no certificate in this repository replaying it. The
census is small enough to enumerate completely, so it is enumerated
here, from nothing.

**This cert banks no data and reads no data.** There is no `data/`
input, no matrix, and no pinned digest — every number in the output is
produced by enumeration inside the run. `tools/bordered_gs.py` is
imported (unmodified) only to build `Ĉ` for the identity check, so the
identities are checked against the same array builder the other certs
use rather than against a private re-implementation.

## What is checked, and why each check is complete

1. **`Λ(y)Λ(y)ᵀ = (Σ_q y_q²)·I₄`.** Each entry of the difference is a
   polynomial in `y₀…y₃` of degree ≤ 2 overall and ≤ 2 in each
   variable, so vanishing on the `5⁴` grid `{−2,…,2}⁴` (5 > 2 values per
   variable) forces it to vanish on `ℤ⁴`. This is a **proof**, not a
   sample.
2. **The 24 admissible arguments.** `Σ_q y_q² = 4` forces `|y_q| ≤ 2`,
   so the box `{−2,…,2}⁴` is exhaustive; the cert checks that the
   solutions are exactly the eight `±2e_j` and the sixteen
   `(±1,±1,±1,±1)`.
3. **The compressed-block identities.** Both sides are linear in the
   eight coset sums `σ_q(c)`, so vanishing on `{−1,0,1}⁸` (3 > 1 values
   per variable) forces them to vanish on `ℤ⁸`. Both reflection
   branches `κ(ρ) ∈ {0,1}` are swept, i.e. `ε = ±1`. Again a proof.
4. **768 Hadamard matrices.** All `2¹⁶` sign matrices of order 4 are
   enumerated; exactly 768 are Hadamard. The cert also checks the two
   facts the mechanism uses: the four rows of each lie in a single
   weight-parity class, so do the four columns, and both classes split
   `384/384`.
5. **The census**, over all `24 × 768² = 14 155 776` triples, by
   **three evaluators** which must agree pair by pair, not merely in
   total:

   | evaluator | what it does |
   | --- | --- |
   | `naive` | tests `⟨W_r, u⟩ ∈ {±4}` for every row `r` of `W = p·Λ(d)ᵀ` and every column `u` of `U`, in integers, pair by pair — the entrywise `±1` test on `E = −¼·p·Λ(d)ᵀ·U`, without forming `E` or dividing by 4 |
   | `table` | precomputes `⟨W_r, v⟩` against all sixteen sign vectors `v`, then tests the four columns of `U` by bitmask |
   | `mechanism` | the structural prediction of **(D-e)**: rows of `W/2` are either `±2e_k` spikes (accepting every `U`) or sign vectors of one weight-parity class (accepting exactly the `U` whose columns carry the other class) |

   These are **not three independent arithmetic routes**. `W` is
   computed once and all three read it; `naive` and `table` share the
   acceptance test as well and differ only in loop order and
   memoisation, and `mechanism` is structurally different but reads the
   same `W`. What their agreement rules out is a bookkeeping error, not
   an error in `W`. The arithmetic itself is cross-checked by the
   interpolation identities of 1 and 3 and by the exact aggregate
   counts, which a corrupted `Λ` moves — that is what the negative
   controls below exercise.

   Measured: **0** naive-vs-table disagreements, **0** mechanism
   mispredictions, **0** pairs on which the mechanism fails to apply.
6. **The corner is Hadamard.** On a deterministic stride through the
   accepted triples (1 152 of them), `E` is built explicitly and
   `E Eᵀ = 4I₄` is checked — the "given that, `EEᵀ = 4I` is automatic"
   of (D-d), exhibited rather than asserted.

## Why the mechanism needs two cases

The sign-vector case — the rows of `W = p·Λ(d)ᵀ/2` are sign vectors —
gives `½`, and it is the whole story for the even arguments, but it
cannot give `¾`. The second case is real and is measured here: for an
**odd** `d`, `Λ(d)` is
itself a 4×4 Hadamard matrix, and for exactly `384` of the `768` choices
of `p` (namely `p = S·Λ(d)` with `S` a signed permutation) every row of
`W` is a spike `±2e_k`, which admits **every** `U`. Hence
`½·½ + ½·1 = ¾`. Measured over all 24 arguments: `6 144` of the
`24 × 768` pairs `(d, p)` are spike (`= 16 × 384`, odd arguments only)
and `12 288` are sign.

## Negative controls

Both are single-cell corruptions of the `Λ` table, run through the same
census code as the real table. Each must move the census; the cert fails
if either leaves it unchanged.

| | control | what it must move | result |
| --- | --- | --- | --- |
| **C1** | sign flip at `Λ[1][3]` (`−y₂ ↦ +y₂`) | at least `d = (1,1,1,1)` | fires: the census moves on 16 of the 24 arguments (all sixteen odd ones, `442 368 → 0`); the eight even arguments are untouched, because flipping one sign still leaves every `Λ(±2e_j)` twice a signed permutation matrix, which is all the even census sees |
| **C2** | variable substitution at `Λ[0][0]` (`y₀ ↦ y₁`) | at least `d = (2,0,0,0)` | fires: the census moves on 12 of the 24 arguments, including four even ones (`294 912 → 0`) |

Neither control is silent, and neither is trivially total: each leaves
part of the census at its true value, so the cert is checking a shape,
not just an exception.

`python certs/10-theorem-d-census/run.py --negative-control` (or
`--negative-control=index`) installs the corruption into the main census
instead of the control slot. That run **must fail**, with a non-zero
exit code — a pass would mean the census does not see the `Λ` table. It
does fail: 17 problems for C1, 13 for C2.

## Honesty labels

| part of the claim | label |
| --- | --- |
| `Λ(y)Λ(y)ᵀ = (Σ_q y_q²)I₄` | **PROVEN** by a complete finite check (interpolation grid; see 1 above) |
| the compressed-block identities over `Ḡ = ℤ₂`, both `ε` branches | **PROVEN** by a complete finite check (interpolation grid; see 3 above) |
| exactly 768 of the `2¹⁶` order-4 sign matrices are Hadamard; rows and columns each split `384/384` by weight parity | **PROVEN** by a complete finite check (full enumeration) |
| the 24 census counts, `768²/2` and `3·768²/4` | **MEASURED** (exact integer arithmetic) and, being a complete finite enumeration replayed here, **PROVEN-BY-CERTIFICATE** — all `14 155 776` triples, no banked input. The label rests on the completeness of the enumeration and on the interpolation-proved identities its arithmetic uses; the three evaluators' agreement rules out bookkeeping error only |
| the **(D-e)** mechanism as a per-pair predicate, spike case included | **PROVEN** by a complete finite check (0 mispredictions on `14 155 776` pairs) |
| Theorem D itself, and the (D-a′) collapse | **PROVEN** (paper-grade; NOTE-B.md §1.5). This cert replays (D-e)'s census and the identities (D-a′) and (D-d) rest on; it is not the theorem's proof. |
| what the census says about *realizable* border data | nothing beyond the even half — seed-derived `d` is always even (§1.5 (D-e)), so the sixteen odd arguments are abstract arguments of the same 4×4 equation and are censused only for completeness |

## How to re-run

From the repository root, on bare Python 3.9 or newer, stdlib only, no
network:

```
python verify/verify.py --selftest
python certs/10-theorem-d-census/run.py
```

Runtime ≈ 20–30 s, machine-dependent — measured here at 21.4 s, of which
the three-evaluator census over all `14 155 776` triples is 18.7 s; the
two negative-control censuses run table-only and cost ≈ 3 s together.
Nothing is written anywhere, inside the repository or outside it.

Exit code 0 iff every identity held, the enumeration returned 768, all
24 census counts hit their exact expected values, the three evaluators
agreed on every pair, and both negative controls fired.
