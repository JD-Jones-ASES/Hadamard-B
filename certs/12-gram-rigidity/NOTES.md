# cert 12 — Gram rigidity: the complete `(3,4)` classification, and the general-branch witnesses

**Replay:** `python certs/12-gram-rigidity/run.py` from the repository root.
Standard library only, exact integers only, no floating point, no network,
about **2 s**, exit 0. Nothing is written anywhere.

---

## The claim

`note/NOTE-B.md` §1.2.1 proves **Theorem E**: under (H1)+(H2) with `s ≥ 1`
and `w > 2s`, every admissible Gram is `M = 4i·P_S` for a conjugation- and
Galois-stable `S ⊆ Ĝ̄` of size `s`, so its spectrum is `{0^(i−s), (4i)^s}`;
and **Corollary E1**: at `i = s+1` the admissible Grams are *exactly* the
`|Ĝ̄[2]|` real-quotient-character twists of the house form
`M_house = (4s+4)I_i − 4J_i`. This certificate measures both edges of that
statement.

> **A — the complete classification at `(s,i) = (3,4)`.** Over both abelian
> groups of order 4, the Grams satisfying the conditions provable *without*
> Theorem E are exactly **2** for `Ḡ = ℤ₄` and **4** for `Ḡ = ℤ₂×ℤ₂`; every
> one has spectrum `{16,16,16,0}`; the lists are precisely the
> real-character twist orbits of the house form; and every member is
> **(H1)-realizable** by an explicit `±1` table `Q ∈ {±1}^{16×12}`. The four
> Grams cert 09 sweeps at order 1916 are the whole `ℤ₂×ℤ₂` list.
>
> **B — the general branch.** Theorem C's **D1 (`i ≤ s+1`) is not a
> consequence of (H1)–(H4)**. The decoded `(s,i) = (1,1)` records at orders
> 1676 (`G = ℤ₄₁₈`) and 1772 (`G = ℤ₄₄₂`), re-read on the subgroups of index
> `i ∈ {11, 19, 22, 38}` and `{13, 17, 26, 34}` — every index `i > 1` with
> `w = n/i > 2s = 2` — satisfy all four hypotheses with Gram `M = 4J_i`,
> and `i > s+1 = 2` in every case. Each Gram is `4i·P_{trivial}`, exactly
> Theorem E's prediction. **8/8.**

## Why this artifact exists

Two gaps, one at each end of Theorem E.

At the small end, cert 09 sweeps the character group of `Ḡ` at order 1916 and
banks four Grams, but its own `NOTES.md` said, correctly at the time, that it
"does not show there are no other non-house Grams". Theorem E says there are
none. That is a paper proof; this cert is the measurement that the proof and
the bank agree at the one cell where a bank exists — and, since `(3,4)` is
the smallest open cell of §1.3's surviving family, the measurement is the
whole classification there.

At the large end, Theorem E holds at *every* cell `(s,i)`, not only at
`i = s+1`, which raises the question of whether the extra cells are empty.
They are not, and part B exhibits them. The consequence is a correction of
scope, not of content: `note/NOTE-B.md` §1.3's classification was already
stated inside the house branch, and part B is *why* that qualifier is
necessary rather than cautious.

## Part A — what is checked, and why the enumeration is complete

The candidate set is `M : Ḡ → ℤ` with the six conditions, each with its
provenance in the note. `w > 2s` is assumed throughout (at order 1916,
`w = 119`).

| | condition at `(3,4)` | where it comes from |
| --- | --- | --- |
| **A1** | `M(0) = 4s = 12` | Theorem A: the diagonal of a Gram of `±1` rows of length `4s`. Forced by the ansatz. |
| **A2** | `M(−c) = M(c)` | `M` symmetric. |
| **A3** | `M(c) ∈ ℤ`, `\|M(c)\| ≤ 12`, `M(c)` even | `M(c) = ⟨u,v⟩ = 12 − 2d` for two `±1` rows of length 12 (Theorem E's step (L2)). **Only mod 2**: `M(c) ≡ 0 (mod 4)` is *not* forced and is nowhere assumed. |
| **A4** | `M` positive semidefinite | it is a Gram matrix. |
| **A5** | `rank M ≤ s = 3` | the four superblock row-spans are mutually orthogonal subspaces of `ℝ^{4s}`, each of dimension `rank M`. |
| **A6** | every eigenvalue `< 4i + 2 = 18` | the Parseval window of Theorem E's step (L1): `λ ≤ N/w = 4i + 4s/w`, and `w > 2s`. |

**The box is exhaustive.** A1 fixes `M(0) = 12`; A3 bounds the off-origin
entries by 12 in absolute value; so `{−12,…,12}` on each free entry covers
every Gram of `±1` rows of length 12 and nothing is enumerated away. A2 is
imposed by construction on `ℤ₄` (`M(1) = M(3)`) and is automatic on
`ℤ₂×ℤ₂`; the run re-checks it on every candidate anyway. This makes A1–A6 an
**over-approximation** of the Grams that really occur, which is the right
direction for a positive verdict: if even the over-approximation is inside
the twist orbit, so is the truth.

**The eigenvalues are verified, not trusted.** They are computed from the
character formulas — `λ_t = Σ_c (−1)^{popcount(t∧c)}M(c)` on `ℤ₂×ℤ₂`, and
`(M₀+2a+b, M₀−b, M₀−2a+b, M₀−b)` with `a = M(1) = M(3)`, `b = M(2)` on `ℤ₄`
— and then, **per matrix**, the monic integer characteristic polynomial
built from the power sums `tr(M^k)`, `k = 1..4`, by Newton's identities is
compared against the polynomial with exactly those roots. A monic integer
polynomial determines its root multiset, so no transcribed formula is
load-bearing. Every division in Newton's identities is exact over `ℤ`.

**What each hypothesis does** (measured in the same pass):

| filter | `ℤ₄` | `ℤ₂×ℤ₂` |
| --- | --- | --- |
| A1–A2 + integer entries + PSD + `rank ≤ 3`, **no window** | 48 | 1 154 |
| the same, entries even (A3) | 24 | 290 |
| all of A1–A6 but with A3's parity **dropped** | **2** | **4** |
| all of A1–A6 | **2** | **4** |

Two readings, both asserted by the run. The Parseval window A6 is doing
essentially all of the cutting (`1 154 → 4`). And the entry-parity clause of
A3 is *not needed* for the classification: dropping it changes nothing. That
is not luck — it is step 4 of the note's argument, where "even, hence
`≤ 4i`" already follows from the window plus the common parity of the
eigenvalues.

**The witnesses.** `data/gram34-witnesses.json` banks six `±1` tables
`Q ∈ {±1}^{16×12}`, one per admissible Gram, rows indexed `k = iI + c`. The
run re-checks all `16 × 16 = 256` inner products of `QQᵀ` against `I₄ ⊗ M`
in exact integers, and checks that the six cover the six admissible Grams
exactly. So the classification is not merely of abstract Grams: **(H1) prunes
nothing further** — every admissible `M` really does carry a `±1` column
table.

**The 1916 cross-check.** The four Grams in `data/n1916-twist.json`'s
`instances` block — cert 09's character-group sweep at order 1916, where
`Ḡ = ℤ₂×ℤ₂` and `w = 119 > 2s` — are compared as a set against the `ℤ₂×ℤ₂`
admissible list. They are equal. Cert 09's sweep is therefore an exhaustive
classification of the admissible Grams at that cell, not merely a sweep of
one construction.

## Part B — what is checked, and what it does and does not mean

Take a decoded `(s,i) = (1,1)` record on `G` with border data `(E, P, Q)`,
`Q` a `4×4` matrix, and a subgroup `K ≤ G` of index `idx`. Re-read the same
instance at index `i = idx`:

```
Q'[(I,c)] := Q[I]  for every c ∈ Ḡ        (row repetition)
P'[r][iJ + c] := P[r][J]                  (column repetition)
M := 4 on all of Ḡ,  i.e. M = 4·J_i
```

The run checks, in exact integer arithmetic and from the record alone:

1. **`Q` is `4×4` Hadamard** (`QQᵀ = 4I₄`) — the `i = 1` Gram.
2. **(H1)** `Q'Q'ᵀ = I₄ ⊗ M` over all `(4i)²` entries. Row repetition sends
   the `i = 1` Gram `4I₄` to `I₄ ⊗ 4J_i`, so the Gram of the re-reading is
   `4J_i`, and *not* the house form `(4s+4)I_i − 4J_i`.
3. **(H2)** `Σ PAF(t) = −M(κt) = −4` for every `t ≠ 0`, and `4n` at `t = 0`.
   The aggregate PAF is recomputed over all of `G` from the definition, once
   per order, and read against the fibres of `κ`.
4. **(H3)** `EEᵀ + w·P'P'ᵀ = N·I_{4s}`. Column repetition gives
   `P'P'ᵀ = i·PPᵀ` and `w·i = n`, so this is the `i = 1` condition
   unchanged; it is checked as stated, not argued away.
5. **(H4)** `EQ'ᵀ + P'Ĉ'ᵀ = 0`, with `Ĉ' = GS(σ'_0,…,σ'_3; κ(ρ))` the
   compressed Goethals–Seidel array over `Ḡ = ℤ_i`, built by
   `tools/bordered_gs.py`'s `gs_array` in the **standard orientation** from
   the record's own coset sums. This is the one real computation: the
   coupling is a `4s × 4i` system that has to hold class by class, and there
   is no reason from the `i = 1` data alone that it should. It does, for all
   eight indices.
6. **Theorem E's prediction.** The `i×i` Gram matrix is `4J_i` exactly; its
   rank is `1 = s`; and its power sums `tr(M^k) = (4i)^k` for `k = 1..4`,
   which for a PSD matrix of rank 1 pins the spectrum to `{4i, 0^(i−1)}`.
   `4J_i = 4i·P_S` for `S = {trivial}`, and the run enumerates the
   Galois-stable singletons of `ℤ̂_i` — the `a ∈ ℤ_i` fixed by every
   `σ_k : χ_a ↦ χ_{ka}`, `gcd(k,i) = 1` — confirming that for `i` prime
   (11, 19, 13, 17) the trivial character is the *only* one available, so
   Theorem E leaves no other Gram possible at those indices.
7. **`i > s+1` and `w > 2s`** in every case: `11, 19, 22, 38 > 2` at 1676
   with `w = 38, 22, 19, 11`; `13, 17, 26, 34 > 2` at 1772 with
   `w = 34, 26, 17, 13`.
8. **The re-reading is inert.** `P̃[r,(J,h)] = P'[r][iJ + κ(h)] = P[r][J]`
   and `Q̃[(I,g),c] = Q'[iI + κ(g)][c] = Q[I][c]` — no dependence on `h` or
   `g` — so the border strips, and hence the assembled matrix, are exactly
   the `i = 1` ones. Checked on a spread of group elements.

**What this proves.** `i ≤ s+1` is not derivable from (H1)–(H4): here are
eight instances that satisfy all four with `i` as large as 38 and `s = 1`.
Theorem C's D1 uses the house form, and the classification corollary of
`note/NOTE-B.md` §1.3 is correspondingly stated inside the house branch.
Part B is the certificate that the qualifier is load-bearing.

**What this does not prove.** Nothing new exists. Item 8 says the general
branch here carries no border and no matrix the classical `i = 1`
construction does not already carry, which is the `s = 1` shape of §1.5's
(D-a′) collapse read at general `i`. The *general* collapse — that every
`i > s+1` instance is an `i = 1` instance in `i`-fold coordinates — is **not
proved**, here or in the note; it would need the (D-a′) border argument run
at general `i`. What is proved is the negative: D1 is a house-branch fact.

## Negative controls

Both must fire, and the cert fails if either is silent.

| | control | what it must move | result |
| --- | --- | --- | --- |
| **C1** | the Parseval window corrupted from `< 4i+2 = 18` to `< 4i+6 = 22` | the Part A census | fires: `ℤ₄` `2 → 4`, `ℤ₂×ℤ₂` `4 → 40`. Not a trivially total failure — the true four survive inside the forty, and 36 of the survivors carry spectra other than `{16,16,16,0}` (e.g. `{20,20,8,0}`), which is exactly the shape the window exists to exclude. Parity also stops being redundant (`112` vs `40`), so the control moves the *structure* of the census, not just its size. |
| **C2** | the Part B Gram perturbed to `M(1) := −4`, the house value, everything else unchanged | **(H1)** at 1676, index 11 | fires: (H1) fails, and so does (H2) — the perturbed profile no longer matches the record's aggregate PAF on that fibre. |

`python certs/12-gram-rigidity/run.py --negative-control[=window|gram]`
installs the corruption into the main path instead of the control slot. That
run **must fail**, with a non-zero exit code — a pass would mean the cert does
not see what it claims to check. It does fail: 11 problems for `window`,
8 for `gram` (one per general-branch index).

## Provenance

The mathematics replayed here — Theorem E, its corollaries, the `(3,4)`
enumeration, and the general-branch witnesses — was **derived and verified at
the source laboratory, 2026-09-01**, and adjudicated there by an independent
skeptic pass before any of it reached this repository. What crosses the
boundary is small and self-certifying:

* the six `Q` tables of `data/gram34-witnesses.json`, banked from the source
  lane's exhaustive `(H1)`-realizability search (`16×12`, two verdict-exact
  symmetry reductions). **Nothing about that search is trusted here**: a
  witness certifies itself, and this cert re-checks every inner product. The
  search's other outputs — the solution counts, the orbit sizes — are *not*
  banked and are not used;
* the eight `(order, index)` pairs of part B, which are re-derived from
  `data/payload-records.json` and `tools/bordered_gs.py` and from nothing
  else. The claim that they satisfy (H1)–(H4) is not banked; it is
  recomputed.

Everything else in this cert — the census, the eigenvalues, the control
table, the 1916 comparison, both negative controls — is enumerated inside
the run.

`data/gram34-witnesses.json` is SHA-256 pinned in `run.py`.
`data/payload-records.json` and `data/n1916-twist.json` are shared with certs
01, 06, 08 and 09 and are bound there by the canonical digests of the
matrices they produce; they are not file-pinned here, on purpose.

| object | sha256 |
| --- | --- |
| `data/gram34-witnesses.json` | `abaf4728e8ba5cd737024b9ab319640c8c634e497884a444b683fb5ee4b93307` |

## Honesty labels

| part of the claim | label |
| --- | --- |
| Theorem E and Corollary E1 themselves | **PROVEN** (paper-grade; `note/NOTE-B.md` §1.2.1). This cert measures their content at `(3,4)` and exhibits the general-branch cells they permit; it is not their proof. |
| the complete admissible list at `(s,i) = (3,4)`, both groups of order 4, under `w > 2s`: 2 and 4, every spectrum `{16,16,16,0}`, both lists equal to the real-character twist orbit | **MEASURED** (exact integer arithmetic) and, being a complete finite enumeration replayed here over an exhaustive candidate box, **PROVEN-BY-CERTIFICATE**. The label rests on the completeness of the box and on the per-matrix verification of the eigenvalues against `tr(M^k)`. |
| `(H1)`-realizability of all six admissible Grams | **PROVEN-BY-CERTIFICATE** (six banked witnesses, every inner product re-checked). The *number* of solutions `Q` is **NOT CLAIMED** — the source search did not compute it and it is not banked. |
| cert 09's four Grams are the complete admissible list at order 1916 | **PROVEN** (Corollary E1) **+ MEASURED** (banked set `=` census set, here) |
| eight `(H1)`–`(H4)` instances with `i > s+1`, at orders 1676 and 1772, with Gram `4i·P_trivial` | **PROVEN-BY-CERTIFICATE** (every hypothesis recomputed from the records; nothing audited) |
| hence: `i ≤ s+1` (Theorem C's D1) is a house-branch fact, not a consequence of (H1)–(H4) | **PROVEN-BY-CERTIFICATE** (the eight witnesses are the counterexample) |
| that the general branch at `i > s+1` contains nothing new *in general* | **NOT CLAIMED.** Item 8 shows it for these eight re-readings only; the (D-a′)-style collapse at general `i` is open. |
| anything about Hadamard equivalence, about order 668 or 716, or about the twist at `ψ(ρ) = −1` | **NOT CLAIMED.** Theorem E is a statement about `M`, hence about seeds and column tables, and this cert is a statement about Theorem E. |

## How to re-run

From the repository root, on bare Python 3.9 or newer, stdlib only, no
network:

```
python verify/verify.py --selftest
python certs/12-gram-rigidity/run.py
python certs/12-gram-rigidity/run.py --negative-control=window   # must FAIL
python certs/12-gram-rigidity/run.py --negative-control=gram     # must FAIL
```

Runtime ≈ 2 s, machine-dependent — measured here at 2.1 s, of which a full
census pass (16 250 candidates over the two groups, a characteristic
polynomial each) is ≈ 0.8 s and there are two of them, the main one and C1's;
the two aggregate-PAF recomputations over `ℤ₄₁₈` and `ℤ₄₄₂` are ≈ 0.3 s
together, and the nine general-branch re-reads (eight plus C2's) ≈ 0.1 s. No
numpy anywhere; no flag adds any. Nothing is written, inside the repository
or outside it.

Exit code 0 iff the two censuses returned exactly their twist orbits with
the house spectrum, the control table matched, all six witnesses verified and
covered the lists, the 1916 sets agreed, all eight general-branch instances
passed every hypothesis and Theorem E's prediction, and both negative
controls fired.
