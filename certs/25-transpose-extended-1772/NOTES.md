# cert 25 — the three classes at order 1772 hold with the transpose in the group

**Label: PROVEN.** Default run:
`python certs/25-transpose-extended-1772/run.py` from the repository root.
Standard library only, **≈ 7.9 s**, exit 0, **202 checks**. That run **audits a
banked exact computation**; `--full` is offered and priced below and **has not
been run in this repository at this order** — cert 23's position at 1772,
unchanged and for the same reasons. The trust boundary is the one certs 15, 19,
20, 21, 22, 23 and 24 draw, and is set out below.

This certificate adds no new matrix and no new construction. It adds **two
exact 4-profiles** — the transposes of the Lemma-T rebuild and of the
orientation switch at order 1772 — and with them the statement cert 23
explicitly withheld: the **row-side only** caveat at this order is discharged.
With cert 24 at 2060, published the same day, **no separation statement in
`note/NOTE-B.md` is row-side any longer**; every class count in the note holds
under the transpose-extended relation.

---

## The relation, and why each pair needs two refutations

`A ~ B` is Hadamard equivalence: `B = D_r P_r A P_c D_c`. The
**transpose-extended** relation is

> `A ≈ B`  iff  `A ~ B`  or  `A ~ Bᵀ`.

Since `~` is symmetric and transposition is an involution, `A ~ Bᵀ ⟺ Aᵀ ~ B`,
so refuting `A ≈ B` takes **two** refutations: `profile(A) ≠ profile(B)`, and
`profile(A) ≠ profile(Bᵀ)` — or, reaching the same statement from the other
side, `profile(Aᵀ) ≠ profile(B)`. The multiset `{|T4|}` over all `C(n,4)` row
4-subsets is a Hadamard-equivalence invariant (`note/NOTE-B.md` §3.1, invariant
**I5**), and the transpose of a Hadamard matrix is Hadamard, so `profile(Bᵀ)`
is that same invariant computed on that other matrix. Nothing else is needed,
and nothing weaker will do: transpose is **not** in the equivalence group, so a
row-side separation says nothing about the transpose-extended relation on its
own — which is exactly why cert 23 stopped where it did.

## The theorem

> **Theorem (1772, transpose-extended).** Let `H` be the decoded
> `(s,i) = (1,1)` bordered Goethals–Seidel record at order 1772 (certs 01/23),
> `H′` its Lemma-T `i = 2` rebuild (certs 02/23), and `H″` the orientation
> switch of `H` — the twelve off-diagonal core blocks negated, the `4×4`
> border, the four row strips and the four column strips unchanged (cert 23).
> Then `H`, `H′` and `H″` are pairwise inequivalent **under the
> transpose-extended relation**. **Order 1772 therefore carries at least three
> Hadamard equivalence classes with the transpose in the group.**

*Proof.* Every leg is an exact `|T4|` 4-profile comparison over all
`C(1772,4) = 409 422 905 815` row 4-subsets, in two arithmetics that agree bin
for bin. `H ≈ H′` is refuted by `H ≁ H′` (**57 of 89**, cert 23) together with
`H ≁ (H′)ᵀ` (**91 of 92**, new here); `H ≈ H″` by `H ≁ H″` (**58**, cert 23) and
`H ≁ (H″)ᵀ` (**91**, new); `H′ ≈ H″` by `H′ ≁ H″` (**53**, cert 23) and
`H′ ≁ (H″)ᵀ` (**91**, new) — and, by the other route, `(H′)ᵀ ≁ H″` (**91**,
new). An invariant that differs is a separation. ∎

**Remark, recorded and not headlined.** Under **plain** Hadamard equivalence
the five matrices profiled at this order — `H, H′, H″, (H′)ᵀ, (H″)ᵀ` — are
**pairwise inequivalent**: all ten profile comparisons separate, the least
separated pair by **52** bins. **At least five classes are therefore exhibited
at 1772** by three constructions and transposition. **The house counts three.**
The transpose-extended relation is the one under which a matrix and its
transpose are the same object, and three is the count that survives *either*
convention — so it is the number stated in the theorem, the note and the
README. The remark is recorded because it is true and because it is what makes
the transpose-extended statement non-trivial: `H′` and `H″` are each
inequivalent to their own transposes here (91 bins apiece), so the transpose
legs were genuinely open until they were computed. Certs 15, 19, 21 and 24
record the same kind of remark at 716 (six matrices), 668 (eight), 1676 (five)
and 2060 (five).

## What is **not** profiled, and what therefore is not said

**`Hᵀ` — the decoded record's own transpose at order 1772 — was never
enumerated.** It is not in the campaign's object list
(`experiments/pr0042/REGISTRATION.md` §2 names `H_1772-twisted-T` and
`H_1772-orient-T` and no third 1772 transpose), it is not banked here, and
`run.py` does not build or pin it. Three consequences, all of them stated in
the run's own output:

* **Nothing is claimed about `H` vs `Hᵀ` at 1772.** Cert 19 decided that
  question at **668 only** (`H ≁ Hᵀ`, 49 of 80 bins). Whether it goes the same
  way here is **open**, and no analogue is asserted.
* The plain-equivalence count above is "**at least five exhibited**", not
  eight. Five matrices are profiled; a sixth (`Hᵀ`) exists and is unprofiled,
  so nothing is said about how it sits against the other five.
* In the verdict table the third column, `Aᵀ vs B`, reads **n/a** for the two
  pairs whose `A` is `H`. Those pairs are refuted by their middle column
  (`A` vs `Bᵀ`) alone, which is all the relation requires; the third column is
  the redundant crossing, available only where `Aᵀ` is banked. `run.py` asserts
  that the `n/a` entries are exactly the ones with `A = H`.

This is the same posture as cert 21 at 1676, object for object.

## The trust boundary — what a default run does and does not establish

The two new `C(1772,4)` enumerations **were not run inside this repository.**
They ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py` (the engine unchanged since the 2060
registration, and the same one certs 11, 13, 14, 15, 19, 20, 21, 22, 23 and 24
rest on), numpy, **16 threads on a rented `c2d-highcpu-16` (`prof42-1`,
`us-east1-b`), 2026-09-03** — under the pre-registration
`experiments/pr0042/REGISTRATION.md`, **flushed 10:17 UTC on 2026-09-02, before
any matrix it governs was built**. That registration fixes the objects (§2),
the decision rules (§4: "differs in any bin ⟹ inequivalent"; "equal in every
bin ⟹ MEASURED, nothing proved"; and, in as many words, `H ≁ H′ᵀ`,
`H ≁ H″ᵀ`, `H′ ≁ H″ᵀ` ⟹ the three-class statement holds under the
transpose-extended relation) and the kill criteria (§5: `blas ≠ bits` in any
bin is a hard stop; a builder digest mismatch ⟹ nothing uploaded) in advance.

The division of labour matters and is stated: **the matrices were built and
verified at the desk** (`experiments/pr0042/build_matrices.py`, every digest in
its `manifest.json`, each through this repository's `verify/verify.py`) and
only then uploaded; the rented machine **enumerated and did nothing else** — it
never assembled or verified a matrix.

**A default `run.py` establishes:**

* the ten bank files are byte-for-byte the ones pinned in `run.py` (SHA-256
  compared in code);
* `H` rebuilds from the banked record through the full master-theorem
  hypothesis re-check (H0–H4, D1/D3/D5, `Σ̄`, the compression-lemma
  cross-check), passes `verify/verify.py`, and carries cert 01's canonical
  digest — computed in-process;
* `H″` is formed from those rows by negating exactly `12·n² = 2 344 368` cells,
  satisfies the alternate-orientation identity `S·H″·S = H_alt` cell by cell
  over all `1772²` cells, and carries cert 23's digest;
* `H′` rebuilds with its seeds **re-derived here** as the `ψ`-twist of the
  decoded seeds (`ψ(g) = (−1)^g` on `ℤ₄₄₂`; `ρ = 441` is odd, so
  `ψ(ρ) = −1`), compared character for character against the bank, so cert 02's
  shared record is bound to `payload-records.json` by computation rather than
  by a file pin;
* `(H′)ᵀ` and `(H″)ᵀ` are **transposed in-process** from those rows, pass
  `verify/verify.py`, and carry the two new pins — so each new bank's declared
  matrix is bound to a matrix **this run built**, not merely to a name;
* transposition is checked as an involution and, on both objects, as a genuine
  move: `(Mᵀ)ᵀ = M` cell for cell, and `Mᵀ ≠ M` in **all 1 772** rows;
* the five matrices carry five distinct canonical digests, established before
  any profile is opened;
* every banked profile satisfies the forced identities — 89 bins on each
  original and 88 on each transpose, all `≡ 4 (mod 8)` (`1772 ≡ 4 mod 8`), every
  key canonical in `[0, 1772]` and every count a positive integer, total
  `C(1772,4) = 409 422 905 815`, second moment
  `n³(n−1)(n−2)/24 = 726 727 740 809 840` — recomputed **and** compared against
  the `second_moment`, `second_moment_want`, `total`, `n` and `C_n_4` fields the
  bank declares, alongside its `schema`, its `folded` field (the **signed** `T4`
  histogram is not an invariant — §3.1), its `impl`, its `engine` and its
  `producer_filename`;
* each bank's declared matrix name and canonical digest against the in-process
  digest of the matrix rebuilt in **this** run;
* `blas == bits` bin for bin on each of the **five** matrices, so every leg of
  the theorem, the transpose legs included, rests on two independent
  implementations;
* all **ten** pair comparisons, with their counts, union support sizes and
  first moments, in **both** arithmetics; that on every one of the six
  original-versus-transpose legs the agreeing bins are **exactly** `[1764]`;
  and the transpose-extended verdicts derived in code from those counts, with
  the second route shown wherever the transpose of `A` is banked;
* controls C0–C6b.

**A default run does not establish that the banked histograms were computed
from the matrices `run.py` rebuilt.** They are *producer-banked*: the digest
each carries is the one the engine recorded against the file it enumerated, and
it equals the digest this certificate pins and the digest computed in-process
here. A self-declared digest is metadata, not a computation. `--full` is what
would close that gap — and at this order it has not been run (see *Runtimes*).

## The evidence chain

**[0]** ten file pins. **[1]** rebuild `H` (hypotheses H0–H4, D1/D3/D5, `Σ̄`,
compression lemma), verify, pin; check the layout `n = |G| = 442`, `s = 1`,
`N = 4(n+s)`; form `H″`, count the negated cells, check `S·H″·S = H_alt`,
verify, pin; re-derive the `ψ`-twist seeds and rebuild `H′`, verify, pin;
transpose `H′` and `H″`, check the involution and that the transpose is a
genuine move, verify, pin; five distinct digests. `Hᵀ` deliberately absent.
**[1b]** control C6 — the dim-`V` trap on the real objects. **[2]** ten banks
audited in exact integers; header, schema, folding and matrix identity per
bank; `blas == bits` ×5. **[3]** the ten pair comparisons, each asserted to the
exact bin count and union size, deltas summing to zero, a non-zero first
moment; then the single agreeing bin on each mixed leg; then, for the three
legs the theorem rests on, the tail that does *not* separate and the first
eight divergent bins printed. **[4]** the transpose-extended verdicts, derived
from those counts, then the plain-equivalence remark and the explicit
`Hᵀ`-is-absent assertion. **[5]** controls C0–C6b. **[6]** `--full`:
`certs/06-668-separation/full_recompute.py` imported (not copied), smoke-tested
on the forced profile of Sylvester `H(128)`, then run on the verified
transposed rows — offered, priced, **not run here**.

**Controls.**

| | control | result |
| --- | --- | --- |
| **C0** | every control matrix is what it claims | the GS condition `Σ_q PAF_q(g) = 4v·[g=0]` re-verified for both GS seed quadruples; all five controls checked Hadamard by brute force; the orientation switch applied to GS `H(28)` is still Hadamard, moves exactly `12·7²` cells, and satisfies the `S·H″·S = H_alt` sign identity cell by cell |
| **C1** | five small Hadamard matrices profiled two ways, one of them the route `--full` takes | straight `O(C(n,4))` enumeration `==` the pair-vector/popcount route on Sylvester `H(8)`, `H(16)`, Paley I `H(20)`, GS `H(28)`, `H(36)`; the two Sylvester profiles match their **forced** values; every bin `≡ n (mod 8)` |
| **C2** | the transposed-profile route, on matrices small enough for straight enumeration | Sylvester `H(8)` and `H(16)` are **symmetric**, so `profile(Mᵀ) = profile(M)` is *forced* — and holds; Paley I `H(20)` is not symmetric and is **MEASURED**, never asserted (its profiles happen to agree) |
| **C3** | the transposes are genuinely different objects | at 1772 they **drop four** of the originals' bins and **add three**, where at 1676 the swap is two for two (cert 21) and at 668/716 each transpose simply drops one (cert 15). Also: `\|T4\| = 1764` is the isolated top bin of all five supports, at **1** count apiece, and `(H′)ᵀ` vs `(H″)ᵀ` stops separating at `\|T4\| = 428` with the 34 bins above it agreeing |
| **C4** | the comparator in the null direction | every banked profile against itself: 0 differing bins, 10 times |
| **C5** | *negative* control: a **total-preserving** corruption of the new `(H″)ᵀ` bank (one count moved from `\|T4\| = 4` to `\|T4\| = 1764`), which only the second-moment identity can catch | the assert fires; the corrupted profile still totals `C(1772,4)` |
| **C6** | the dim-`V` trap on the real objects | `dim W` = **1771** on all five — an invariant that separates **none** of the ten pairs clause [4] proves inequivalent; `dim V` reads 1770 / 1771 / 1770 / 1771 / 1770 on `H`, `H′`, `H″`, `(H′)ᵀ`, `(H″)ᵀ` and is worthless (Trap 1, §3.1) |
| **C6b** | the same trap demonstrated | under a seeded signed row negation of Sylvester `H(16)`, `dim V` moves 4 → 5, `dim W` does not move, and the `|T4|` profile does not move either |

## The separations

All ten pairs, as `run.py` prints them; the counts are identical in both
arithmetics.

| pair | differing | union | `Σ \|T4\|·Δ` | source |
| --- | --- | --- | --- | --- |
| `H` vs `H′` | 57 | 89 | +55 798 800 | cert 23 |
| `H` vs `H″` | 58 | 89 | +72 493 792 | cert 23 |
| `H′` vs `H″` | 53 | 89 | +16 694 992 | cert 23 |
| `H` vs `(H′)ᵀ` | **91** | 92 | +33 942 976 | new |
| `H` vs `(H″)ᵀ` | **91** | 92 | +56 189 584 | new |
| `H′` vs `(H′)ᵀ` | **91** | 92 | −21 855 824 | new |
| `H′` vs `(H″)ᵀ` | **91** | 92 | +390 784 | new |
| `H″` vs `(H′)ᵀ` | **91** | 92 | −38 550 816 | new |
| `H″` vs `(H″)ᵀ` | **91** | 92 | −16 304 208 | new |
| `(H′)ᵀ` vs `(H″)ᵀ` | **52** | 88 | +22 246 608 | new |

The verdicts, derived from those counts:

| pair | `A` vs `B` | `A` vs `Bᵀ` | `Aᵀ` vs `B` | verdict |
| --- | --- | --- | --- | --- |
| `H` vs `H′` | 57 | **91** | n/a | **SEPARATED** |
| `H` vs `H″` | 58 | **91** | n/a | **SEPARATED** |
| `H′` vs `H″` | 53 | **91** | 91 | **SEPARATED** |

The two `n/a` entries are `Hᵀ`, which was never profiled; each of those pairs
is refuted by its middle column alone, which is all the relation asks. The
third row shows both routes, and both hold.

### The support: four dropped, three added — and an isolated top bin

At 668 and 716 every transpose populated **one bin fewer** than its original,
and it was the same bin at each order (cert 15, control C3). At 1676 the count
did not change and the transposes **swapped two** bins with the originals
(cert 21, control C3). **At 1772 it is neither:** the three originals share one
support of **89** bins, the two transposes share another of **88**, and the two
differ by **four dropped and three added**.

| bin | `H` | `H′` | `H″` | `(H′)ᵀ` | `(H″)ᵀ` |
| --- | --- | --- | --- | --- | --- |
| `\|T4\| = 620` | — | — | — | 24 | 24 |
| `\|T4\| = 636` | 24 | 24 | 24 | — | — |
| `\|T4\| = 668` | 6 | 6 | 6 | — | — |
| `\|T4\| = 708` | 18 | 18 | 18 | — | — |
| `\|T4\| = 724` | — | — | — | 16 | 16 |
| `\|T4\| = 772` | 12 | 12 | 12 | — | — |
| `\|T4\| = 916` | — | — | — | 4 | 4 |

So the union support of the five profiles is **92** bins, which is why every
original-versus-transpose comparison is reported as *91 of 92*, the
same-side originals as *n of 89* and the transpose pair as *52 of 88*. Those
seven sparse bins carry counts between 4 and 24 out of `4.1·10¹¹` 4-subsets;
the separations do not rest on them — the bulk does — but they are a clean
structural fact and `run.py` asserts them.

**The isolated top bin.** All five profiles have their maximum at
`|T4| = 1764 = n − 8`, at **one** 4-subset each, with a wide gap below it: the
next populated bin down is `|T4| = 772` in an original and `916` in a
transpose. That single bin is **the one bin every original-versus-transpose
comparison agrees on** — 91 of the 92 union bins differ, and the 92nd is 1764 —
so "the bulk separates and the extreme tail does not" is at this order the
sharpest it can be: exactly one bin of agreement, and it is the tail.
`run.py` asserts both halves (the agreeing set is exactly `[1764]`, and the top
differing bin is `|T4| = 916`).

### The three new legs the theorem rests on

`run.py` prints the first eight divergent bins of each; the full 91-bin lists
are the banked JSONs.

**`H` vs `(H′)ᵀ`** — 91 of 92 bins differ; largest `|Δ| = 5 778 020` at
`|T4| = 20`, i.e. `1.04·10⁻⁴` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H` | 61 934 029 130 | 59 717 426 100 | 55 520 798 861 | 49 778 171 281 | 43 034 414 080 | 35 886 330 792 |
| `(H′)ᵀ` | 61 931 496 909 | 59 716 876 356 | 55 526 576 881 | 49 775 031 462 | 43 037 400 888 | 35 883 182 292 |
| Δ | −2 532 221 | −549 744 | **+5 778 020** | −3 139 819 | +2 986 808 | −3 148 500 |

**`H` vs `(H″)ᵀ`** — 91 of 92; largest `|Δ| = 6 177 634` at `|T4| = 20`, i.e.
`1.11·10⁻⁴` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H` | 61 934 029 130 | 59 717 426 100 | 55 520 798 861 | 49 778 171 281 | 43 034 414 080 | 35 886 330 792 |
| `(H″)ᵀ` | 61 931 975 703 | 59 714 957 052 | 55 526 976 495 | 49 775 901 272 | 43 037 560 744 | 35 881 849 506 |
| Δ | −2 053 427 | −2 469 048 | **+6 177 634** | −2 270 009 | +3 146 664 | −4 481 286 |

**`H′` vs `(H″)ᵀ`** — 91 of 92; largest `|Δ| = 4 359 218` at `|T4| = 20`, i.e.
`7.85·10⁻⁵` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H′` | 61 931 468 288 | 59 716 764 308 | 55 522 617 277 | 49 780 050 525 | 43 034 542 248 | 35 884 341 470 |
| `(H″)ᵀ` | 61 931 975 703 | 59 714 957 052 | 55 526 976 495 | 49 775 901 272 | 43 037 560 744 | 35 881 849 506 |
| Δ | +507 415 | −1 807 256 | **+4 359 218** | −4 149 253 | +3 018 496 | −2 491 964 |

The fourth leg, `(H′)ᵀ` vs `H″`, is the other route for the same pair: 91 bins,
largest `|Δ| = 4 523 043` at `|T4| = 52`. The other three new legs are
`H′` vs `(H′)ᵀ` (largest `|Δ| = 5 019 063` at `|T4| = 28`), `H″` vs `(H″)ᵀ`
(`3 884 570` at `|T4| = 20`) and `(H′)ᵀ` vs `(H″)ᵀ` (`1 919 304` at
`|T4| = 12`, `3.21·10⁻⁵` of that bin). Every difference vector sums to zero, as
it must; the first moment, which nothing forces, does not — and all ten are
non-zero. As everywhere else in this repository, the largest discrepancies are
of order `10⁻⁴`–`10⁻⁵` of their bins: invisible to any sample of practical
size, and nothing cheaper than the exact profile would have found them.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports) — all
five rebuilt in clause [1] of every run:

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `H` decoded `(1,1)` record | `1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2` | yes (cert 01's pin) |
| `H′` Lemma-T `i = 2` rebuild | `82484769a28ac93201f208ca3256bfd491f8edc8bb5e0309764c4f609a113378` | yes (cert 02's pin; seeds re-derived) |
| `H″` orientation switch | `7f1fae050def5b9b7bdc491c05b24551465cbea8d3d9482a9cd23c98ba607e53` | yes (cert 23's pin; formed from `H`) |
| `(H′)ᵀ` | `471f705168cc87b4a1256625ce345a0308f356c5a4e9b1807154163078773238` | yes (new) |
| `(H″)ᵀ` | `0dffc98fbb6e290a6592bf4e253cf7ab973add8be8b642af83434386a2568864` | yes (new) |

The two transposed digests are the desk's own, recorded in
`experiments/pr0042/manifest.json` when the matrices were built and verified
there, and re-derived in-process by every run of this certificate. **`Hᵀ` is
not in this table**, on purpose: no profile of it exists.

**Banked files** (SHA-256 of the file bytes, compared in `run.py`). The two
`-T-` pairs are this certificate's own, banked 2026-09-03 under cert 25; the
six others are cert 23's, reused verbatim and re-pinned here:

| file | SHA-256 | banked for |
| --- | --- | --- |
| `data/sep1772-twisted-T-exact-blas.json` | `20507070c8cf28702bc9093b2f5324e0736744959cfed3c7f44ea9900df3d101` | cert 25 |
| `data/sep1772-twisted-T-exact-bits.json` | `673d558086ac938e75ba28c2a2b240edf3737482a07714cdaa33aa419b1300dc` | cert 25 |
| `data/sep1772-orient-T-exact-blas.json` | `960c9a1893accfd9e29eb44febd19226340f4db612188090daaa5a50d222e0b7` | cert 25 |
| `data/sep1772-orient-T-exact-bits.json` | `e2901d9976067293c12794b6cf7cf004af9f044acd2b92d7f1abbad25374e94a` | cert 25 |
| `data/sep1772-decoded-exact-blas.json` | `5985d5f9e1a7ceb54d12dc65e7d5179412eaffb78d94ac3ecf8366db5edbc0d4` | cert 23 |
| `data/sep1772-decoded-exact-bits.json` | `f4b2522d2b8d8ff06ed0195f72051ef28e1f47d7aa691d0afd425a1c84cf98e6` | cert 23 |
| `data/sep1772-twisted-exact-blas.json` | `3b010bf5406916977f060a033990427633bebff0891895b72a0e943afde0f76e` | cert 23 |
| `data/sep1772-twisted-exact-bits.json` | `481ba9f78be1c1d0d2072cfad67598d30c74423ec12e12702350dbba6fb72b35` | cert 23 |
| `data/sep1772-orient-exact-blas.json` | `f58bb4d7db1106950e6506899ffc1329787740cc35b4e0c48574b1d6f06bbea7` | cert 23 |
| `data/sep1772-orient-exact-bits.json` | `57ea8694738d3a6bffe5fa56bf7bab5c6e0d354425b7678f57f6b0d691472f5d` | cert 23 |

`data/payload-records.json` and `data/twisted-i2-records.json` are not
file-pinned here, for cert 23's reason: they are shared with certs 01, 02 and
23, and the binding pin on each is the canonical digest of the matrix it
produces, checked in clause [1] — reinforced, for the twisted record, by the
`ψ`-twist re-derivation that binds it to `payload-records.json` outright.

## Provenance of the two new profiles

Produced by the source laboratory's unchanged engine
(`Hadamard-2060`, `experiments/inequiv/exact_profile_big.py`) on the rented
`c2d-highcpu-16` `prof42-1` (`us-east1-b`), 16 threads, 2026-09-03, under
`experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC 2026-09-02;
Amendment 1 ~11:05 UTC re-priced the campaign). The `orient-T` `blas` leg
landed 03:20 UTC and the `twisted-T` `blas` leg by 04:24 UTC, which completed
the five 1772 `blas` legs; the `orient-T` `bits` leg landed 13:48 UTC and the
`twisted-T` `bits` leg 16:07 UTC, the last of the instance's worklist, after
which it self-deleted. The seconds and peak resident sets below are the values
those runs recorded in the JSONs themselves.

Banked into `data/` by the laboratory's `experiments/pr0042/bank.py` under
`--cert 25 --date 2026-09-03`, which refuses unless both implementations audit,
agree bin for bin, and declare the digest `manifest.json` pins for that matrix.
It adds the seven header fields each bank's own `banked_note` names — `schema`,
`description`, `matrix`, `matrix_canonical_sha256`, `producer_filename`,
`arithmetic`, `banked_note` — and **nothing numeric**: every numeric field is
the producer's own output, unaltered. `run.py` checks `schema`, `matrix`,
`producer_filename`, `matrix_canonical_sha256`, `matrix_sha256`, `impl`,
`folded`, `engine`, `n`, `C_n_4`, `total`, `second_moment` and
`second_moment_want` against values it recomputes or rebuilds.

Credit, as everywhere in this repository, is to stations: the theorems and the
certificate are this repository's; the engine, the pre-registration, the matrix
building and verification, and the rented enumeration are the source
laboratory's; `payload-records.json` encodes the posting team's mathematical
content (`PROVENANCE.md`), and no priority of any kind is claimed on it.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 7.9 s** (exit 0, 202 checks; measured 2026-09-03 on the desk, three runs, 7.9 / 7.9 / 7.8 s) |
| rebuild + full hypothesis re-check + `verify.py`, per matrix | 1.5 s (0.67 s of it the hypothesis re-check) |
| the orientation switch (0.27 s), the cell count (0.09 s) and the `S·H″·S` identity (0.72 s) at `N = 1772` | ≈ 1.1 s |
| forming a transpose at `N = 1772` | 0.02 s (the `verify.py` call that follows dominates) |
| `dim V` / `dim W`, per 1772 matrix | ≈ 0.15 s (five of them) |
| producing `sep1772-orient-T-exact-{blas,bits}.json` | 3 401.2 s (peak 788.0 MB) / 8 511.5 s (peak 109.6 MB) |
| producing `sep1772-twisted-T-exact-{blas,bits}.json` | 3 405.4 s (peak 787.9 MB) / 8 501.5 s (peak 109.5 MB) |
| `run.py --full` at 1772 | **not run in this repository** — see below |

**Why no `--full` leg was run here.** The flag is offered and wired exactly as
in certs 11, 13, 14, 15, 19, 20, 21, 22, 23 and 24 —
`certs/06-668-separation/full_recompute.py` imported by a `sys.path` insert
rather than copied, BLAS threads capped at three before numpy loads,
smoke-tested against the forced profile of Sylvester `H(128)` first — and
`--matrix` and `--impl` select the leg. The price is cert 23's, unchanged,
because the order and the module are the same: one 1772 leg is about **68×**
the 716 leg the same module took in this repository (cert 14, 400.3 s), i.e.
roughly **7–8 hours** for a single `blas` matrix at three threads on this desk.
That 68× is the source laboratory's **measured** sub-`n⁵` scaling — its
desk-measured 716→2060 ratio of 137
(`experiments/pr0042/REGISTRATION.md`, Amendment 1) implies the exponent
`4.66`, and `(1772/716)^4.66 = 68`. On the `Θ(n⁵)` law quoted elsewhere here
the same leg is `(1772/716)⁵ = 93×`, i.e. ≈ **10.3 h**. Both are estimates and
both say *hours*; the smaller, measured one is quoted so that the decision not
to run is not defended with an inflated price. The `blas` route also
materialises a `C(n,2) × n` pair matrix — `1 569 106 × 1 772`, which is 2.78 GB
as `int8` and **11.1 GB** as the `float32` copy `_profile_blas` makes — past
this desk's memory; the `bits` route (a `1 569 106 × 28` `uint64` packing) is
the tractable one at this order. So this certificate's default verdict is
*banked exact computation audited*, and it says so in its own output. Certs 06,
08, 11, 13, 14, 15 and 19 each have an in-repo `--full` `blas` leg on the
record; certs 20, 21, 22, 23, 24 and 25 do not, and the word *replayed* is not
used of them.

## What is NOT claimed

* **The default run recomputes nothing.** Only `--full` would bind a bank to a
  matrix by computation, and **no `--full` leg has been run in this repository
  at 1772**, for cert 23 or for this certificate.
* **Nothing about `H` vs `Hᵀ` at 1772.** `Hᵀ` was never enumerated. Cert 19
  decided the analogous question at 668 and at 668 only; nothing here extends
  it, and the plain-equivalence count is "at least five exhibited", not eight.
* **No general theorem.** That the `ψ(ρ) = −1` twist and the orientation switch
  land in different classes — and stay apart with the transpose in the group —
  is now known at four bordered orders (668, 716, 1676, 1772), and the
  orientation switch is a class of its own at the unbordered 2060 as well
  (cert 24). Four instances are four instances.
* **No claim that three is the count** at 1772. It is a lower bound exhibited by
  the matrices in hand, under either convention.
* **No novelty or priority claim of any kind** at 1772. This certificate counts
  classes among the artifacts banked here; it says nothing about who first
  exhibited a Hadamard matrix of this order, which is long settled. The decoded
  record is a public artifact (`PROVENANCE.md`).
* **Matching invariants prove nothing.** `dim W` agrees on all five objects the
  4-profile proves pairwise inequivalent, and Paley `H(20)`'s profile agrees
  with its transpose's. Read every "agrees" in this repository as "did not
  separate", never as "the same".

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/25-transpose-extended-1772/run.py
```

The default path is standard library only, no network, and is the whole
certificate. The `--full` paths below import numpy (finder-side only, never in
the trust chain) and cap BLAS threads at three; at this order each is hours,
and the `blas` variants want memory this desk does not have. **None of them has
been run here.**

```
python certs/25-transpose-extended-1772/run.py --full --impl bits --matrix or-T
python certs/25-transpose-extended-1772/run.py --full
```

Exit code 0 iff every check passed.
