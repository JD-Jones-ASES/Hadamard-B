# cert 16 — Theorem E′ (Gram rigidity under `w > s`) and the `w = s` boundary

**Label: PROVEN** (the theorems are paper-grade proofs, `note/NOTE-B.md` §1.7)
**+ PROVEN-BY-CERTIFICATE** (the `(3,3,3)` escape, and the exhaustion of the
boundary classification at the cells run). Default run:
`python certs/16-theorem-eprime-boundary/run.py` from the repository root.
Standard library only, exact integers only, **no data file at all** — nothing
here is audited from a bank, and nothing here reads anything. **2.6 s, 59
checks, exit 0.**

---

## The theorems

> **Lemma (mod 4).** For every `x : G → {±1}` and every `t ≠ 0`,
> `PAF_x(t) ≡ n (mod 4)`. Hence under (H2), `M(c) ∈ 4ℤ` for every `c ∈ Ḡ`.

> **Theorem E′.** Assume (H1), (H2), `s ≥ 1` and **`w > s`**. Then
> `S = {χ : M̂(χ) ≠ 0}` has `|S| = s`, is closed under conjugation and under
> the Galois action, and `M = 4·Σ_{χ∈S} χ = 4i·P_S`, with spectrum
> `{0^(i−s), (4i)^s}`.

> **Theorem E′ (boundary).** Let `w = s` and let `M` satisfy the necessary
> conditions (`M(0) = 4s`, `M ∈ 4ℤ`, `|M(c)| ≤ 4s`, PSD, `M̂ ≤ N/w`,
> `rank ≤ s`). Then exactly one of: (a) `s < i` and `M = 4i·P_S`; (b) `s = i`
> and `M = 4i·I − 4C` with `C(c) = ξ(c)·[c ∈ Ḡ₀∖0]` for a subgroup `Ḡ₀ ≤ Ḡ`
> and a real character `ξ` of `Ḡ₀`; (c) `s = i+1` and `M = 4s·I`;
> (d) `s ≥ i+2` cannot occur.

The proofs are in `note/NOTE-B.md` §1.7. Theorem E (§1.2.1) needed `w > 2s`;
the mod-4 lemma halves the lock. **This certificate proves neither theorem.**
It carries the three finite things around them: a control on the lemma, the
sharpness witness rebuilt from nothing, and the boundary classification
exhausted at small parameters by two routes that share no code.

## What `run.py` checks (exit 0 iff every check passes)

**[A] the mod-4 lemma, as a control.** `PAF_x(t) ≡ n (mod 4)` on **3 990**
random `(group, sequence, lag)` triples over eleven abelian groups
(`ℤ₇, ℤ₈, ℤ₁₂, ℤ₂²×ℤ₃, ℤ₃², ℤ₂×ℤ₆, ℤ₄², ℤ₂⁴, ℤ₅², ℤ₉, ℤ₂×ℤ₉`), zero
violations; and the proof's own statement — along each cycle of `u ↦ u+t` the
number of sign changes of `x` is even — instantiated on `ℤ₁₂`.

**[B] the `(3,3,3)` escape, rebuilt from nothing.** Every `±1` quadruple on
`G = ℤ₃²`, `K = {(0,b)}`, whose aggregate PAF is `−12` on `K∖0` and `+4` off
`K`, found by a meet-in-the-middle over **all 512** sequences per seed — no
coset-sum prefilter, so nothing is assumed about the shape of a solution.
**490 212** ordered `(pair,pair)` solutions; the same enumeration on `ℤ₉`
with `K = ⟨3⟩` returns **0**. One witness is re-checked lag by lag, and so is
the witness quoted in `note/NOTE-B.md` §1.7
(`++-++-++-, ++-++-+-+, ++-+-+-++, ++--+++-+`, row sums `(3,3,3,3)`). A
`12×12` `±1` matrix `Q` with `QQᵀ = I₄ ⊗ (16I₃ − 4J₃)` is found by
depth-first search and re-multiplied entry by entry; `16I − 4J` is PSD of
rank `3 = s = i` by exact `Fraction` LDLᵀ — spectrum `{4,16,16}`, not a
projector multiple. So `w > s` cannot be weakened to `w ≥ s`.

**[C] the boundary classification, two routes, no shared code.**

* *route A (sieve).* Every `Ḡ`-invariant symmetric `M : Ḡ → 4ℤ` with
  `M(0) = 4s` and `|M(c)| ≤ 4s` is enumerated, and kept iff `M` is PSD,
  `rank M ≤ s`, and `N·I − w·M` is PSD (the Parseval window). Every decision
  is an exact `Fraction` LDLᵀ. **No floats anywhere.**
* *route B (construction).* The classification's list (a)/(b)/(c) built
  explicitly — Galois-orbit character sums evaluated as **Ramanujan sums**,
  so every entry is an exact integer; the `s = i` family enumerated over all
  subgroups `Ḡ₀ ≤ Ḡ` and all real characters `ξ` of each.

At every `w = s` cell the two **sets** must coincide, element for element; a
disagreement prints the symmetric difference and fails the run. At every
`w > s` cell every survivor must be a projector form (Theorem E′ proper).
The `w < s` cells are recorded and nothing is asserted at them.

**[D] global checks.** The sieve is not vacuous below the boundary
(non-projector survivors exist there); on the boundary the escapes occur
**exactly** at `s = i` and never with `s < i`; the cell tally is consistent.

## The cells run, and the counts

Default: **45 `(s,i,w)` triples**, over **every** abelian `Ḡ` of order `i`,
giving **59 `(s,i,w,Ḡ)` cells** — **22** on the boundary (`w = s`, two routes
agree), **26** with `w > s` (rigid: **20** in the range `s < w ≤ 2s` that
Theorem E did not reach, 6 with `w > 2s` where it already applied), **11**
below it. `--wide` adds `(7,8,w)` on `ℤ₈` for `w ∈ {6,7,8,15}`: **63** cells,
**23 / 28 / 12**, 62 checks.

Selected survivor counts, all from the run:

| cell | `Ḡ` | candidates | survivors | escapes |
| --- | --- | --- | --- | --- |
| `(2,2,2)` | `ℤ₂` | 5 | 3 | 2 |
| `(3,3,3)` | `ℤ₃` | 7 | 2 | 1 |
| `(4,4,4)` | `ℤ₄` | 81 | 5 | 4 |
| `(4,4,4)` | `ℤ₂²` | 729 | 11 | 10 |
| `(5,5,5)` | `ℤ₅` | 121 | 2 | 1 |
| `(6,6,6)` | `ℤ₆` | 2 197 | 6 | 5 |
| `(3,4,3)` | `ℤ₄` / `ℤ₂²` | 49 / 343 | 2 / 4 | 0 / 0 |
| `(5,6,5)` | `ℤ₆` | 1 331 | 2 | 0 |
| `(5,3,5)` | `ℤ₃` | 11 | 0 | 0 |
| `(2,4,2)` | `ℤ₄` / `ℤ₂²` | 25 / 125 | 2 / 6 | 0 / 0 |
| `(7,8,7)` (`--wide`) | `ℤ₈` | 50 625 | 2 | 0 |

Every boundary cell with `s < i` is rigid; every escape sits at `s = i` (or
at `s = i+1`, where `M = 4s·I` is counted as non-projector because no
projector form of size `s > i` exists).

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **2.6 s** (59 checks, exit 0; measured here 2026-09-02, Python 3.14, one worker) |
| **`run.py --wide`** | **46.8 s** (62 checks, exit 0; the `(7,8,·)` sieve is 15⁴ = 50 625 candidates per cell, ≈ 11 s each) |

## What is NOT claimed

* **Neither theorem is machine-checked.** Theorem E′ and the boundary
  classification are paper proofs (`note/NOTE-B.md` §1.7). This certificate
  exhausts the classification *at the cells it runs* and rebuilds the
  sharpness witness; it does not verify the proofs.
* Nothing about **(H3)/(H4) at `(3,3,3)`**. The escape is a statement about
  seeds and column table only. Whether a bordered Hadamard matrix of order
  `4(9+3) = 48` exists with this `M` is not asked here — `w = 3 ≤ 2s`, so
  Theorem C's D3 does not apply and the corner need not be Hadamard.
* Nothing about `w = s` **beyond the necessary-condition level**. The
  classification lists what `(H1)+(H2)`'s consequences allow and exhibits
  one realised escape; it does not say every listed form is realised by
  seeds.
* Nothing new about `v = 523`: all four live cells there have `w > 2s`,
  where Theorem E already applied. The `333`-cell degenerate tail
  (`w ≤ 9 < s`) is *below* the boundary — this certificate shows the sieve
  is non-vacuous there, i.e. that **no** rigidity theorem of this type
  exists in that regime; it neither closes nor opens any of those cells.
* The `(7,8,·)` cells on `ℤ₂×ℤ₄` and `ℤ₂³` are **not swept** (the sieve has
  `15⁶` and `15⁷` candidates there). The classification at those cells rests
  on the proof, not on this enumeration.
* The `490 212` count and the `0` on `ℤ₉` are **MEASURED** by one exhaustive
  implementation here. They agree with the source laboratory's count, which
  used a coset-sum prefilter; this run uses none.

## How to re-run

```
python verify/verify.py --selftest
python certs/16-theorem-eprime-boundary/run.py
python certs/16-theorem-eprime-boundary/run.py --wide
```
