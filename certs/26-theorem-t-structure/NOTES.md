# cert 26 — Theorem T: the family is closed under transposition

**Label: PROVEN** (the identities; `note/NOTE-B.md` §1.9) **+
COMPUTATIONAL-EVIDENCE** (the small-order transposition census). Default run:
`python certs/26-theorem-t-structure/run.py` from the repository root.
Standard library only, **≈ 16.9 s**, exit 0, **46 checks**. `--full` is offered,
priced, and **has been run in this repository** — 2026-09-05, **45.4 s**, 64
checks, `ALL CHECKS PASS`.

This certificate adds no new matrix to the repository and no new separation. It
adds the *structural* fact that the transposed matrices the separation certs
measure are members of the same family, and it checks the sign bookkeeping that
makes that true **entrywise**, on arbitrary data and on every banked
coset-border record.

---

## The theorem

> **Theorem T (`note/NOTE-B.md` §1.9).** Write
> `BGS(x₀,x₁,x₂,x₃; ρ; E,P,Q)` for the house-orientation bordered array of
> Theorem A, `S′ = diag(−1,−1,+1,+1) ⊗ I_i` acting on the table index
> `iI + c`, `S₀₁ = diag(−1,−1,+1,+1) ⊗ I_n` acting on the core, and
> `D̃ = diag(I_{4s}, S₀₁)`. Then
>
> `Hᵀ = D̃ · BGS(x₀∘(−), −x₁, x₂, x₃; ρ; Eᵀ, (S′Q)ᵀ, (PS′)ᵀ) · D̃`,
>
> a house-orientation member of the family. And, **exactly and with no
> conjugation**, `(H″)ᵀ = BGS(x₀∘(−), x₁, x₂, x₃; ρ; Eᵀ, Qᵀ, Pᵀ)`, where `H″`
> is the orientation switch of `H`.

**The theorem is branch-free.** Its proof uses nothing beyond `H Hᵀ = N I`: the
right-hand side *is* `D̃ Hᵀ D̃`, `D̃` is a `±1` diagonal, so it is Hadamard
exactly when `H` is; it is a house instance on the data displayed; and the
*only if* half of Theorem A hands back (H1)–(H4) for it in one line. In
particular the route through `EEᵀ = 4s·I` — which holds only in the house
branch under `w > 2s` — is **not** used. Theorem T holds at **every** `(s,i)`
and **every** `w`.

## Remark R, and both sides of it

The `(H″)ᵀ` identity is an identity of **signed arrays**, whether or not either
side is Hadamard. The orientation switch *with the border unchanged* is
Hadamard **only when `P` annihilates the off-diagonal part of `Ĉ`** — automatic
at `(1,1)`, where the row sums are `(±2,0,0,0)` and `Ĉ = ±2I`, and **false in
general**. Clause [4] evaluates the condition on every banked record and then
puts `H″` through `verify/verify.py`:

| record | cell | `P·offdiag(Ĉ)ᵀ = 0` | `verify.py` on `H″` |
| --- | --- | --- | --- |
| 668, 716 | `(1,1)` | **yes** | `HADAMARD` |
| 1676, 1772 (`--full`) | `(1,1)` | **yes** | `HADAMARD` |
| 1916 | `(3,4)` | **no** | `FAIL — non-orthogonal row pairs` |
| 1388 (`--full`) | `(5,6)` | **no** | `FAIL` |
| 1436 (`--full`) | `(7,8)` | **no** | `FAIL` |

A remark that only ever confirmed itself would be worth nothing; this one is
exercised in **both** directions, and the certificate asserts the agreement
record by record rather than asserting the good half. Wherever `H″` or `(H″)ᵀ`
is called a Hadamard matrix in `note/NOTE-B.md`, this condition is what is
being used.

## The evidence chain

**[0]** SHA-256 pins on the four banked data files read. **[1]** the block
facts (F1),(F2) entrywise on 30 deterministically generated cores over ten
groups, cyclic and not, with arbitrary seeds — the identities are seed-free:
every off-diagonal block symmetric, block `(J,I) = −` block `(I,J)`,
`Cᵀ = C^{sw}(x₀∘(−), x₁, x₂, x₃)` (only the **type-1** seed reverses), and
`C^{sw}(y) = S_{0j}·C(…,−y_j,…)·S_{0j}` for `j = 1, 2, 3`. **[2](a)** Theorem T
and the switch identity entrywise on **12 arbitrary bordered shapes** — random
corners and random tables, satisfying none of (H1)–(H4), neither side Hadamard.
**[2](b),[3]** the same two identities entrywise on the banked coset-border
records: 668, 716 at `(1,1)` and 1916 at `(3,4)` by default; 1388 at `(5,6)`,
1436 at `(7,8)`, 1676 and 1772 at `(1,1)` under `--full`. **[4]** remark R,
both sides, through the trust chain. **[5]** the **T-images** of 668, 716 and
1916 assembled here from the right-hand side of the theorem, handed to
`verify/verify.py`, and pinned. **[6]** the small-order census.

## Pinned digests

The three T-images built in clause [5], canonical SHA-256 (the digest
`verify/verify.py` reports):

| order | canonical SHA-256 of `BGS(the T-image data)` |
| --- | --- |
| 668 | `6396100a41b75a2ddbf308396f8ec15c1cbab8ae56decdb007eb5e04f2bca2ba` |
| 716 | `1bff41d81d80ab63e60660eb36cd31d3d2e909a900bd239e6db589342c833b40` |
| 1916 | `2cfb31f52f6f4cd612716a19ab87b06a0f161f0d3c250e6a88b2bc6b1b031b3d` |

These are **not** the digests certs 15, 19, 21, 24 and 25 carry for their
transposed objects. Those certificates bank `Mᵀ` itself; this one banks
`D̃ Mᵀ D̃`, the house instance. `D̃` is a `±1` diagonal, so it changes neither
the equivalence class nor the `|T4|` profile — only the canonical digest.
Recorded so that the two families of pins are never confused for each other.
The source laboratory recorded the same three digests to eight hex digits
(`6396100a…`, `1bff41d8…`, `2cfb31f5…`) when it first checked these
identities; all three prefixes reproduce here, from this repository's own
records and this certificate's own assembler.

## The small-order census — COMPUTATIONAL-EVIDENCE

Eleven instances, all built here: the five Goethals–Seidel arrays `GS(28)`,
`GS(36)`, `GS(44)`, `GS(52)`, `GS(60)` on **Williamson** seeds — all four seeds
symmetric — found by an exhaustive meet-in-the-middle over the symmetric `±1`
sequences on `ℤ_t` (deterministic: the first quadruple in lexicographic order);
the order-52 gate instance at `(1,2)` (`data/h52-gate.json`, cert 03); the
order-76 non-scalar record at `(1,1)` (`data/h76-nonscalar.json`, cert 04); and
the four `(2,4)` matrices of cert 18 (`data/cell24-records.json`), the two
`H(88)` under `--full`. Every instance is put through `verify/verify.py` before
it is profiled.

Each is decided by the exact `|T4|` 4-profile of `H` against that of `Hᵀ` over
all `C(N,4)` row 4-subsets. That profile is a Hadamard-equivalence invariant
(§3.1, invariant **I5**) and the transpose of a Hadamard matrix is Hadamard, so
a **different** profile is a **proof** of `H ≁ Hᵀ`.

| instance | cell | differing bins / union | verdict |
| --- | --- | --- | --- |
| `GS(28)` Williamson(7) | `s = 0` | 0 / 4 | undecided |
| `GS(36)` Williamson(9) | `s = 0` | 0 / 4 | undecided |
| `GS(44)` Williamson(11) | `s = 0` | **4 / 5** | `H ≁ Hᵀ` |
| `GS(52)` Williamson(13) | `s = 0` | 0 / 6 | undecided |
| `GS(60)` Williamson(15) | `s = 0` | **5 / 7** | `H ≁ Hᵀ` |
| `H(52)` gate | `(1,2)` | **4 / 6** | `H ≁ Hᵀ` |
| `H(76)` non-scalar | `(1,1)` | **4 / 9** | `H ≁ Hᵀ` |
| `H(56)` seed0 ρ̄0 | `(2,4)` | 0 / 7 | undecided |
| `H(56)` seed2 ρ̄3 | `(2,4)` | **8 / 8** | `H ≁ Hᵀ` |
| `H(88)` seed0 ρ̄0 (`--full`) | `(2,4)` | **8 / 10** | `H ≁ Hᵀ` |
| `H(88)` seed1 ρ̄1 (`--full`) | `(2,4)` | **9 / 11** | `H ≁ Hᵀ` |

**Totals: 11 instances, 7 proved inequivalent to their own transpose, 4
undecided by the profile** (9 and 5 on the default path).

**"Undecided" is not "equivalent."** This certificate runs **no isomorphism
search**. An agreeing invariant is a failure to separate, never a proof of
equivalence, and the run's own output says so instance by instance. The four
undecided rows are open questions here, not fixed points.

**What the census is evidence for.** That **all four seeds symmetric does not
suffice** for `H ~ Hᵀ`: two of the five Williamson-seeded arrays are proved
inequivalent to their transposes, which is enough to refute the sufficiency
claim outright. Beyond that it is evidence, at eleven instances, that
`H ~ Hᵀ` is not the normal state of affairs in this family. It is **not** a
theorem, and §1.9 states the conjecture as a conjecture.

## What is cited and not replayed

The source laboratory ran a wider census: **44 instances** across the cells
`s = 0`, `(1,1)`, `(1,2)`, `(3,4)`, `(2,4)`, each decided exactly — by the
4-profile, or, where the profiles agreed, by an explicit signed permutation
found by a finder-side individualisation-refinement search on the Hadamard
graph. Totals **16 equivalent, 28 inequivalent, 0 undecided**; at orders `≥ 44`
exactly **one** equivalent (`H56_cell24_n12_seed0_rho0`) against **24**
inequivalent, the three Williamson-seeded GS arrays at 44, 52 and 60 among the
inequivalent.

That search is **finder-side**: a `False` from it rests on the completeness of
the refinement, so it is not in this repository's trust chain and is **not**
replayed here. It is recorded as a **source-laboratory measurement**. Nothing in
this certificate or in `note/NOTE-B.md` §1.9 rests on it.

Two points of contact are worth recording, because they are corroboration
rather than repetition. The one instance the laboratory found equivalent to its
transpose at orders `≥ 44` is `H56_cell24_n12_seed0_rho0` — and that is exactly
the `H(56)` this certificate's profile leaves **undecided**, as it must if the
two are equivalent. And where this certificate proves `H ≁ Hᵀ` at `GS(44)` and
`GS(60)` on its own Williamson quadruples, the laboratory reached the same
verdict at those orders on quadruples of its own. `GS(52)` is the one place the
two part company: the profile does not decide this certificate's Williamson(13)
quadruple, where the laboratory's search decided a different one.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **≈ 16.9 s** (exit 0, 46 checks; measured 2026-09-05 on the desk) |
| **`run.py --full`** | **45.4 s** (exit 0, 64 checks; **run here**, 2026-09-05) |
| the entrywise identity at order 1916 | ≈ 5 s (two matrices of `1916²` cells built and compared) |
| `verify/verify.py` on an order-1916 matrix | ≈ 1.5 s |
| the `C(88,4) = 2 331 890` profile, per matrix | ≈ 4 s |

`--full` is cheap enough to run, and was run; there is no hours-long leg
anywhere in this certificate. What `--full` buys is the four remaining banked
records — including the `(5,6)` and `(7,8)` cells, where remark R's condition
fails and `H″` is correctly rejected — and the two `H(88)`.

## What is NOT claimed

* **No structural criterion for `H ~ Hᵀ`.** Theorem T reduces the question to
  whether the class of `H` is a fixed point of an involution on the parameter
  set; no sufficient condition is known, and *all four seeds symmetric does not
  suffice* — this certificate refutes that at `GS(44)` and `GS(60)`.
* **An agreeing profile is not equivalence.** No isomorphism search is run
  here; four census rows are undecided and are reported as undecided.
* **Nothing about the separations of §3.** The class counts at 668, 716, 1676,
  1772 and 2060 are certs 08, 15, 19, 20, 21, 24 and 25's, and are **cited, not
  re-banked**. Theorem T supplies the structural reason the transposes are
  still members of the family; it is not needed for the separations themselves,
  and none of them depends on it.
* **Nothing about `H(2092)`.** Theorem T is a closure statement about the
  family, not an existence statement at any order.
* **The general orientation switch is not claimed Hadamard.** Remark R is the
  whole of what is asserted, and it is asserted with its failures.

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/26-theorem-t-structure/run.py
python certs/26-theorem-t-structure/run.py --full
```

Standard library only, no network, no numpy anywhere — including under
`--full`. Exit code 0 iff every check passed.

## Provenance

Theorem T and its sign bookkeeping are the source laboratory's
(`Hadamard-2060`, skeptic pass II and III; decisions **D-062** and **D-068**),
with the branch-free argument and the restoration of remark R being the repairs
the desk audit imposed at adoption (D-068 items 8(a) and 10). This certificate
is written here: it shares no code with the laboratory's
`skeptic-pass/transpose_structure.py`, reads only this repository's own `data/`,
and re-derives every identity from the block array. The 44-instance census is
the laboratory's and is cited above, not replayed. Credit is to stations, as
everywhere in this repository.
