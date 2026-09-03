# cert 24 — the three classes at order 2060 hold with the transpose in the group

**Label: PROVEN.** Default run:
`python certs/24-transpose-extended-2060/run.py` from the repository root.
Standard library only, **≈ 10 s**, exit 0, **219 checks**. That run **audits a
banked exact computation**; `--full` is offered and priced below and **has not
been run in this repository at this order** — cert 22's position at 2060,
unchanged and for the same reasons. The trust boundary is the one certs 15, 19,
20, 21, 22 and 23 draw, and is set out below.

This certificate adds no new matrix and no new construction. It adds **two
exact 4-profiles** — the transposes of the plain Goethals–Seidel array and of
the unbordered orientation switch at order 2060 — and with them the statement
cert 22 explicitly withheld: the **row-side only** caveat at this order is
discharged. With cert 25 at 1772, published the same day, **no separation
statement in `note/NOTE-B.md` is row-side any longer**; every class count in
the note holds under the transpose-extended relation.

---

## The relation, and why each pair needs two refutations

`A ~ B` is Hadamard equivalence: `B = D_r P_r A P_c D_c`. The
**transpose-extended** relation is

> `A ≈ B`  iff  `A ~ B`  or  `A ~ Bᵀ`.

Refuting `A ≈ B` means refuting **both** disjuncts. The first is
`profile(A) ≠ profile(B)`. For the second, `~` is symmetric and transposition
is an involution, so `A ~ Bᵀ ⟺ Aᵀ ~ B`; **either** of
`profile(A) ≠ profile(Bᵀ)` and `profile(Aᵀ) ≠ profile(B)` therefore refutes it,
and one of the two suffices. The multiset `{|T4|}` over all `C(n,4)` row
4-subsets is a Hadamard-equivalence invariant (`note/NOTE-B.md` §3.1, invariant
**I5**), and the transpose of a Hadamard matrix is Hadamard, so each of those is
that same invariant computed on that other matrix. Nothing weaker will do:
transpose is **not** in the equivalence group, so a row-side separation says
nothing about the transpose-extended relation on its own — which is exactly why
cert 22 stopped where it did.

**At this order the two routes are not interchangeable**, because `Gᵀ` is not
banked (below). The pair `{P, G}` is carried only by the `Aᵀ vs B` route, the
pair `{G, H″}` only by the `A vs Bᵀ` route, and `{P, H″}` by both. Every pair
has at least one, which is all the relation asks; `run.py` asserts exactly which
entries are `n/a` and that each `n/a` pair is carried by its other route. This
is the one structural difference from cert 21 at 1676, where the middle column
was available for every pair.

## The theorem

> **Theorem (2060, transpose-extended).** Let `P` be `2060-plain`, the plain
> Goethals–Seidel array over the raw `ℤ₅₁₅` seed; `G` be `2060-gist`, the
> `×104`-twisted array that is byte-for-byte the publicly posted `H(2060)`
> (both cert 07); and `H″` be `P` with its **twelve off-diagonal `515`-blocks
> negated** — the orientation switch in its **unbordered** form, 2060 being the
> degenerate `s = 0` layer where the array is a plain `4×4` GS array of
> circulant blocks and there is no border to leave alone (cert 22). Then `P`,
> `G` and `H″` are pairwise inequivalent **under the transpose-extended
> relation**. **Order 2060 therefore carries at least three Hadamard
> equivalence classes with the transpose in the group.**

*Proof.* Every leg is an exact `|T4|` 4-profile comparison over all
`C(2060,4) = 748 155 697 135` row 4-subsets, in two arithmetics that agree bin
for bin.

* `P ≈ G` is refuted by `P ≁ G` (**146 of 147**, cert 07) together with
  `Pᵀ ≁ G` (**134 of 134**, new here);
* `P ≈ H″` by `P ≁ H″` (**107 of 145**, cert 22) and `P ≁ (H″)ᵀ` (**145 of
  146**, new) — with `Pᵀ ≁ H″` (**145 of 146**, new) as the redundant second
  route;
* `G ≈ H″` by `G ≁ H″` (**146 of 147**, cert 22) and `G ≁ (H″)ᵀ` (**134 of
  134**, new).

An invariant that differs is a separation. ∎

**Remark, recorded and not headlined.** Under **plain** Hadamard equivalence
the five matrices profiled at this order — `P, G, H″, Pᵀ, (H″)ᵀ` — are
**pairwise inequivalent**: all ten profile comparisons separate, the least
separated pair by **92** bins. **At least five classes are therefore exhibited
at 2060** by two constructions, the orientation switch and transposition.
**The house counts three.** The transpose-extended relation is the one under
which a matrix and its transpose are the same object, and three is the count
that survives *either* convention — so it is the number stated in the theorem,
the note and the README. The remark is recorded because it is true and because
it is what makes the transpose-extended statement non-trivial: `P` and `H″` are
each inequivalent to their own transposes here (145 bins apiece), so the
transpose legs were genuinely open until they were computed. Certs 15, 19 and 21
record the same kind of remark at 716 (six matrices), 668 (eight) and 1676
(five).

## What is **not** profiled, and what therefore is not said

**`Gᵀ` — the posted matrix's own transpose — was never enumerated.** It is not
in the campaign's object list (`experiments/pr0042/REGISTRATION.md` §2 names
`H_2060-orient-T` and `H_2060-plain-T` and no third 2060 transpose), it is not
banked here, and `run.py` does not build or pin it. Three consequences, all of
them stated in the run's own output:

* **Nothing is claimed about `G` vs `Gᵀ` at 2060.** Cert 19 decided the
  analogous question at **668 only** (`H ≁ Hᵀ`, 49 of 80 bins). Whether it goes
  the same way here is **open**, and no analogue is asserted.
* The plain-equivalence count above is "**at least five exhibited**", not
  eight. Five matrices are profiled; a sixth (`Gᵀ`) exists and is unprofiled, so
  nothing is said about how it sits against the other five.
* In the verdict table the `A vs Bᵀ` column reads **n/a** for the pair whose
  `B` is `G`, and the `Aᵀ vs B` column reads **n/a** for the pair whose `A` is
  `G`. Each of those pairs is refuted by its other route, which is all the
  relation requires; `run.py` asserts that the `n/a` entries are exactly those
  two and that both pairs are carried.

None of this weakens the theorem: each of the three pairs has both of its
required refutations in hand.

## The trust boundary — what a default run does and does not establish

The two new `C(2060,4)` enumerations **were not run inside this repository.**
They ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py` (the engine unchanged since the 2060
registration, and the same one certs 11, 13, 14, 15, 19, 20, 21, 22 and 23 rest
on), numpy, **16 threads on a rented `c2d-highcpu-16` (`prof42-2`,
`us-east1-b`), 2026-09-03** — under the pre-registration
`experiments/pr0042/REGISTRATION.md`, **flushed 10:17 UTC on 2026-09-02, before
any matrix it governs was built** (Amendment 1, ~11:05 UTC, added the second
instance for the 2060 legs and re-priced the campaign; both predate the legs).
That registration fixes the objects (§2), the decision rules (§4: "differs in
any bin ⟹ inequivalent"; "equal in every bin ⟹ MEASURED, nothing proved"; and,
in as many words, that the two transposes give the transpose-extended form for
all three pairs) and the kill criteria (§5: `blas ≠ bits` in any bin is a hard
stop; a builder digest mismatch ⟹ nothing uploaded) in advance.

The division of labour matters and is stated: **the matrices were built and
verified at the desk** (`experiments/pr0042/build_matrices.py`, every digest in
its `manifest.json`, each through this repository's `verify/verify.py`, the
plain source re-verified against cert 07's pin first) and only then uploaded;
the rented machine **enumerated and did nothing else** — it never assembled or
verified a matrix.

**A default `run.py` establishes:**

* the ten bank files are byte-for-byte the ones pinned in `run.py` (SHA-256
  compared in code);
* `P` and `G` rebuild from `data/sep2060-records.json` — the four normalised
  seeds, the normalising shifts, the `×104` twist and the CRT relabelling —
  after the classical Goethals–Seidel condition `Σ_q PAF_q(t) = 4v·[t=0]` is
  **re-verified on the raw seeds** (the `s = 0` layer of the master theorem,
  checked and not assumed), pass `verify/verify.py`, and carry cert 07's two
  canonical digests — computed in-process;
* `H″` is formed from `P`'s rows by negating exactly `12·515² = 3 182 700`
  cells, satisfies the alternate-orientation identity **twice** — as a sign
  pattern `S·H″·S = H_alt` cell by cell, and against the alternate array
  **assembled from the same raw seeds** and separately verified — and carries
  cert 22's digest;
* `Pᵀ` and `(H″)ᵀ` are **transposed in-process** from those rows, pass
  `verify/verify.py`, and carry the two new pins — so each new bank's declared
  matrix is bound to a matrix **this run built**, not merely to a name;
* transposition is checked as an involution and, on both objects, as a genuine
  move: `(Mᵀ)ᵀ = M` cell for cell, and `Mᵀ ≠ M` in **all 2 060** rows;
* the five matrices carry five distinct canonical digests (the alternate array's
  is a sixth), established before any profile is opened;
* every banked profile satisfies the forced identities — its declared bin count,
  all bins `≡ 4 (mod 8)` (`2060 ≡ 4 mod 8`), every key canonical in `[0, 2060]`
  and every count a positive integer, total `C(2060,4) = 748 155 697 135`,
  second moment `n³(n−1)(n−2)/24 = 1 543 448 476 598 000` — recomputed **and**
  compared against the header fields the bank declares, alongside its `schema`,
  its `folded` field (the **signed** `T4` histogram is not an invariant — §3.1)
  and its `impl`. Two schemas are in play — cert 07's
  `sep2060-exact-profile/1` and the campaign's `exact-4-profile/1` — and each is
  checked on its own terms;
* each bank's declared matrix name and canonical digest against the in-process
  digest of the matrix rebuilt in **this** run;
* `blas == bits` bin for bin on each of the **five** matrices, so every leg of
  the theorem, the transpose legs included, rests on two independent
  implementations;
* all **ten** pair comparisons, with their counts, union support sizes and
  first moments, in **both** arithmetics; the support structure; and the
  transpose-extended verdicts derived in code from those counts, with the route
  each pair takes named and the `n/a` entries asserted;
* controls C0–C7b.

**A default run does not establish that the banked histograms were computed
from the matrices `run.py` rebuilt.** They are *producer-banked*: the digest
each carries is the one the engine recorded against the file it enumerated, and
it equals the digest this certificate pins and the digest computed in-process
here. A self-declared digest is metadata, not a computation. `--full` is what
would close that gap — and at this order it has not been run (see *Runtimes*).

## The evidence chain

**[0]** ten file pins. **[1]** rebuild `P` and `G` from the seed record (the GS
condition re-verified on the raw seeds; the declared row sums checked against
the seeds' own; the assembly's seeds checked to be the ones just verified),
verify, pin; form `H″` by the **unbordered** switch, count the negated cells,
check `S·H″·S = H_alt` as a sign pattern **and** against the independently
assembled alternate array, verify both, pin; transpose `P` and `H″`, check the
involution and that the transpose is a genuine move, verify, pin; five distinct
digests plus the alternate array's sixth. `Gᵀ` deliberately absent.
**[1b]** control C7 — the dim-`V` trap on the real objects. **[2]** ten banks
audited in exact integers; headers, schema, folding and matrix identity per
bank, per schema; `blas == bits` ×5. **[3]** the ten pair comparisons, each
asserted to the exact bin count and union size, deltas summing to zero, a
non-zero first moment; then the support structure, the single agreeing bin, the
tail behaviour, and the first eight divergent bins of the three legs the
theorem rests on. **[4]** the transpose-extended verdicts, derived from those
counts, with the two-route logic and the `n/a` assertions, then the
plain-equivalence remark and the explicit `Gᵀ`-is-absent assertion.
**[5]** controls C0–C7b. **[6]** `--full`:
`certs/06-668-separation/full_recompute.py` imported (not copied), smoke-tested
on the forced profile of Sylvester `H(128)`, then run on the verified
transposed rows — offered, priced, **not run here**.

**Controls.**

| | control | result |
| --- | --- | --- |
| **C0** | every control matrix is what it claims | the GS condition `Σ_q PAF_q(t) = 4v·[t=0]` re-verified for both GS seed quadruples; all five controls checked Hadamard by brute force |
| **C1** | five small Hadamard matrices profiled two ways, one of them the route `--full` takes | straight `O(C(n,4))` enumeration `==` the pair-vector/popcount route on Sylvester `H(8)`, `H(16)`, Paley I `H(20)`, GS `H(28)`, `H(36)`; the two Sylvester profiles match their **forced** values; every bin `≡ n (mod 8)` |
| **C2** | the **unbordered** orientation switch, exercised end to end where the Hadamard property can be checked by brute force — `s = 0` there too, so these are the exact structural analogue of the 2060 switch | on GS `H(28)` and `H(36)`: the switch is still Hadamard, moves exactly `12v²` cells, satisfies identity (a) as a sign pattern and identity (b) against the alternate array assembled from the same seeds, the alternate array is Hadamard, and the two `\|T4\|` profiles are equal — so (b) exhibits `H″` as the other orientation rather than a coincidence of counts |
| **C3** | the transposed-profile route, on matrices small enough for straight enumeration | Sylvester `H(8)` and `H(16)` are **symmetric**, so `profile(Mᵀ) = profile(M)` is *forced* — and holds; Paley I `H(20)` is not symmetric and is **MEASURED**, never asserted (its profiles happen to agree) |
| **C4** | the comparator in the null direction | every banked profile against itself: 0 differing bins, 10 times |
| **C5** | *negative* control: a **total-preserving** corruption of the new `Pᵀ` bank (one count moved from `\|T4\| = 4` to `\|T4\| = 1204`), which only the second-moment identity can catch | the assert fires; the corrupted profile still totals `C(2060,4)` |
| **C7** | the dim-`V` trap on the real objects | `dim W` = **2059** on all five — an invariant that separates **none** of the ten pairs clause [4] proves inequivalent; `dim V` reads 2058 / 2058 / 2059 / 2059 / 2058 on `P`, `G`, `H″`, `Pᵀ`, `(H″)ᵀ` and is worthless (Trap 1, §3.1) |
| **C7b** | the same trap demonstrated | under a seeded signed row negation of Sylvester `H(16)`, `dim V` moves 4 → 5, `dim W` does not move, and the `|T4|` profile does not move either |

The support control lives in clause [3] rather than in the control table,
because at this order it is part of the statement: see *The support* below.

## The separations

All ten pairs, as `run.py` prints them; the counts are identical in both
arithmetics.

| pair | differing | union | `Σ \|T4\|·Δ` | source |
| --- | --- | --- | --- | --- |
| `P` vs `G` | 146 | 147 | −54 709 381 744 | cert 07 |
| `P` vs `H″` | 107 | 145 | +242 005 088 | cert 22 |
| `G` vs `H″` | 146 | 147 | +54 951 386 832 | cert 22 |
| `P` vs `Pᵀ` | **145** | 146 | +8 049 832 320 | new |
| `P` vs `(H″)ᵀ` | **145** | 146 | +7 782 382 064 | new |
| `G` vs `Pᵀ` | **134** | **134** | +62 759 214 064 | new |
| `G` vs `(H″)ᵀ` | **134** | **134** | +62 491 763 808 | new |
| `H″` vs `Pᵀ` | **145** | 146 | +7 807 827 232 | new |
| `H″` vs `(H″)ᵀ` | **145** | 146 | +7 540 376 976 | new |
| `Pᵀ` vs `(H″)ᵀ` | **92** | 123 | −267 450 256 | new |

The verdicts, derived from those counts:

| pair | `A` vs `B` | `A` vs `Bᵀ` | `Aᵀ` vs `B` | verdict |
| --- | --- | --- | --- | --- |
| `P` vs `G` | 146 | n/a | **134** | **SEPARATED** |
| `P` vs `H″` | 107 | **145** | 145 | **SEPARATED** |
| `G` vs `H″` | 146 | **134** | n/a | **SEPARATED** |

The bold entry in each row is the refutation of the second disjunct that the
theorem uses. The two `n/a` entries are `Gᵀ`, which was never profiled; each of
those pairs is carried by its other route, and the middle row shows both routes,
both holding.

### The support: twenty-three dropped, one added

At 668 and 716 every transpose populated **one bin fewer** than its original,
and it was the same bin at each order (cert 15, control C3). At 1676 the count
did not change and the transposes **swapped two** bins with the originals
(cert 21). **At 2060 the transposes are the sparsest object here:** `P` and `H″`
populate the same **145** bins, `Pᵀ` and `(H″)ᵀ` the same **123**, and the two
supports differ by **twenty-three dropped and one added**.

| | bins | relation to the originals' support |
| --- | --- | --- |
| `P`, `H″` | 145 | — |
| `Pᵀ`, `(H″)ᵀ` | 123 | 23 of the originals' bins unpopulated; one new bin, `\|T4\| = 1204`, at 6 counts in both transposes |
| `G` | 133 | 14 of the originals' unpopulated; two of its own, `\|T4\| = 892, 908` (cert 22) |

The 23 bins the transposes drop are `|T4| = 844, 860, 876, 924, 940, 956, 972,
988, 1004, 1020, 1036, 1052, 1068, 1076, 1084, 1100, 1116, 1124, 1148, 1156,
1164, 1180` and `1236` — twenty-three of the forty bins the originals populate
at `|T4| ≥ 844`, including the top four and the maximum, carrying counts
between 6 and 72 out of `7.5·10¹¹` 4-subsets. Above `|T4| = 1140` the
transposes populate nothing but `1204`. So an original-versus-transpose
union support is **146** bins (which is why those comparisons read *145 of 146*),
`G`'s union with a transpose is only **134**, and the union of all five supports
is **148**.

Two consequences the run asserts:

* **Exactly one bin agrees** in each of the four original-versus-transpose
  comparisons that involve `P` or `H″`: `|T4| = 900`, at **300** counts on each
  side. `G` reads **1380** in that same bin — which is precisely why `G` against
  either transpose separates in **every** bin of the union, the strongest form
  the comparison can take: **134 of 134**.
* **The extreme tail does separate.** The originals' top bin, `|T4| = 1236`, is
  populated in no transpose, so it differs on every original-versus-transpose
  leg. That matches cert 22's row-side finding at this order (the top bin reads
  12 against 6 for `H″` vs `P`) and is unlike 668, 716, 1676 and 1772, where the
  tail always agreed. On `Pᵀ` vs `(H″)ᵀ`, by contrast, the familiar pattern
  returns: the top differing bin is `|T4| = 788` and the **24** bins above it
  agree exactly.

### The three new legs the theorem rests on

`run.py` prints the first eight divergent bins of each; the full lists are the
banked JSONs.

**`G` vs `Pᵀ`** — **134 of 134** bins differ; largest `|Δ| = 755 660 673` at
`|T4| = 4`, i.e. `7.10·10⁻³` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `G` | 106 441 862 685 | 102 965 772 295 | 96 360 629 035 | 87 289 459 295 | 76 582 817 395 | 65 122 330 205 |
| `Pᵀ` | 105 686 202 012 | 102 338 966 122 | 95 976 526 112 | 87 178 728 849 | 76 719 856 458 | 65 431 766 493 |
| Δ | **−755 660 673** | −626 806 173 | −384 102 923 | −110 730 446 | +137 039 063 | +309 436 288 |

**`P` vs `(H″)ᵀ`** — 145 of 146; largest `|Δ| = 135 208 122` at `|T4| = 4`,
i.e. `1.28·10⁻³` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `P` | 105 825 054 681 | 102 453 207 083 | 96 023 541 023 | 87 167 542 521 | 76 658 478 649 | 65 350 483 547 |
| `(H″)ᵀ` | 105 689 846 559 | 102 342 122 541 | 95 978 082 148 | 87 179 171 099 | 76 718 630 462 | 65 430 536 650 |
| Δ | **−135 208 122** | −111 084 542 | −45 458 875 | +11 628 578 | +60 151 813 | +80 053 103 |

**`G` vs `(H″)ᵀ`** — **134 of 134**; largest `|Δ| = 752 016 126` at `|T4| = 4`,
i.e. `7.07·10⁻³` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `G` | 106 441 862 685 | 102 965 772 295 | 96 360 629 035 | 87 289 459 295 | 76 582 817 395 | 65 122 330 205 |
| `(H″)ᵀ` | 105 689 846 559 | 102 342 122 541 | 95 978 082 148 | 87 179 171 099 | 76 718 630 462 | 65 430 536 650 |
| Δ | **−752 016 126** | −623 649 754 | −382 546 887 | −110 288 196 | +135 813 067 | +308 206 445 |

The redundant fourth leg, `Pᵀ` vs `H″`, is the other route for the pair
`{P, H″}`: 145 of 146, largest `|Δ| = 135 792 463` at `|T4| = 4`. Every
difference vector sums to zero, as it must; the first moment, which nothing
forces, does not — and all ten are non-zero.

The `G`-versus-transpose legs are the largest separations anywhere in this
repository as a fraction of the bins they touch (`134/134`) and in absolute
delta (`7.6·10⁸`, `7.1·10⁻³` of the first bin). That is not a stronger *kind*
of statement than a one-bin separation — an invariant that differs is an
invariant that differs — but it is worth recording that these two legs are not
close calls. The `P`- and `H″`-versus-transpose legs sit at `1.3·10⁻³` of their
first bins, and `Pᵀ` vs `(H″)ᵀ` at `3.45·10⁻⁵`, which is the familiar scale
everywhere else here: invisible to any sample of practical size.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports) — all
five rebuilt in clause [1] of every run:

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `P` plain GS array | `510f89b7b423c85da0c7cada52cb0f62e0415d736d2040d6701f66c4e524cf6a` | yes (cert 07's pin) |
| `G` the posted `H(2060)` | `c7a145d86210740dd3f8ea21ca896a54d6916007a042638f17c8c47f097200f7` | yes (cert 07's pin; the public artifact's digest) |
| `H″` unbordered orientation switch | `4e1891b095b8aafa21176e494038f199b495c96a840bdb003e231c160870b801` | yes (cert 22's pin; formed from `P`) |
| `Pᵀ` | `5a980e3ad69f02fdece4f0aca40f46b9afafe15703e87b82941ff661db1a2960` | yes (new) |
| `(H″)ᵀ` | `8558904d9d61c7547b835c25791da97f7d3e0bc1cd852de082b8475c01b34337` | yes (new) |
| the alternate GS array (clause [1](b)) | `40e1d1c8cd40e94016c453f12e520a8518e7d29b773d3adaae3f484eca64398d` | yes (cert 22's pin; assembled from the seeds) |

The two transposed digests are the desk's own, recorded in
`experiments/pr0042/manifest.json` when the matrices were built and verified
there, and re-derived in-process by every run of this certificate. **`Gᵀ` is
not in this table**, on purpose: no profile of it exists.

**Banked files** (SHA-256 of the file bytes, compared in `run.py`). The two
`-T-` pairs are this certificate's own, banked 2026-09-03 under cert 24; the two
`orient` banks are cert 22's and the four `plain`/`gist` banks cert 07's, all
reused verbatim and re-pinned here at exactly the values those certificates
carry:

| file | SHA-256 | banked for |
| --- | --- | --- |
| `data/sep2060-plain-T-exact-blas.json` | `1f8aa2469f22d3fee20ce3cf3618dbc275c81d4f0a6493271a0532eb77e9ae57` | cert 24 |
| `data/sep2060-plain-T-exact-bits.json` | `8a5389e0b74c0d27c5981ac1ca5605b834517f751e83a81b14d8759aa9e557d6` | cert 24 |
| `data/sep2060-orient-T-exact-blas.json` | `b21f217eaa70dfe01da599cac3c260f4cef9436d398833962bc46224298a771f` | cert 24 |
| `data/sep2060-orient-T-exact-bits.json` | `468386dab0c5b6c7c317cdc4f3694113f1a2d8cfac449f2a19ba618209c56694` | cert 24 |
| `data/sep2060-orient-exact-blas.json` | `38135aef205b4428760dc0439b29196b6776f215397fb34880156c66ee283f00` | cert 22 |
| `data/sep2060-orient-exact-bits.json` | `18dcf4c2e5d603324182eeef45c13e89ef80d5c0c0a5add3a7f90c333e4e5e87` | cert 22 |
| `data/sep2060-exact-blas-plain.json` | `5428aeac7b570fff55975c2b737fae9e8d0b717ec511735b68893e609a0037d8` | cert 07 |
| `data/sep2060-exact-bits-plain.json` | `e6c3af94712d0ba5cf3a3047796ccd474970036fec211b41a5579b7ff892ca49` | cert 07 |
| `data/sep2060-exact-blas-gist.json` | `a20b9a63cd3d93046c251b5c19aabeeac412b8f7933bbafa82d0210320e3aef0` | cert 07 |
| `data/sep2060-exact-bits-gist.json` | `9d8cc4b55c297c7e948df3e7639613a0580fc3e54af9eb12399bc010337f8a93` | cert 07 |

`data/sep2060-records.json` is not file-pinned here, for cert 22's reason: it is
shared with certs 07 and 22, and the binding pin on it is the canonical digest
of each matrix it produces, checked in clause [1].

## Provenance of the two new profiles

Produced by the source laboratory's unchanged engine
(`Hadamard-2060`, `experiments/inequiv/exact_profile_big.py`) on the rented
`c2d-highcpu-16` `prof42-2` (`us-east1-b`, created 18:14 UTC 2026-09-02),
16 threads, under `experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC
2026-09-02; Amendment 1 ~11:05 UTC added this instance and re-priced the
campaign). The `orient-T` `blas` leg landed ~02:50 UTC 2026-09-03 and its `bits`
leg 07:55 UTC; the `plain-T` `blas` leg 09:35 UTC and its `bits` leg 14:32 UTC,
at which point the instance reached the end of its worklist and self-deleted.
The seconds and peak resident sets below are the values those runs recorded in
the JSONs themselves.

Banked into `data/` by the laboratory's `experiments/pr0042/bank.py` under
`--cert 24 --date 2026-09-03`, which refuses unless both implementations audit,
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
laboratory's; `G` is a public artifact and the seed record encodes its
mathematical content (`PROVENANCE.md`), and no priority of any kind is claimed
on it — order 2060 was settled by that posted matrix.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 10 s** (exit 0, 219 checks; measured 2026-09-03 on the desk, three runs, 10.3 / 9.7 / 9.7 s) |
| assemble `P` and `G` from the seed record (0.24 s), GS condition re-verified on the raw seeds (0.08 s) | 0.3 s |
| the unbordered switch (0.18 s), the cell count (0.13 s), identity (a) as a sign pattern (0.46 s), identity (b) against the assembled alternate array (0.26 s) | ≈ 1.0 s |
| forming a transpose at `N = 2060` | 0.03 s (the `verify.py` call that follows dominates) |
| `verify/verify.py` per matrix at `N = 2060`, six of them | ≈ 1.3 s apiece (`manifest.json` records the same) |
| `dim V` / `dim W`, per 2060 matrix | ≈ 0.2 s (five of them) |
| producing `sep2060-orient-T-exact-{blas,bits}.json` | 6 645.1 s (peak 860.4 MB) / 17 685.3 s (peak 131.9 MB) |
| producing `sep2060-plain-T-exact-{blas,bits}.json` | 6 650.4 s (peak 860.2 MB) / 17 678.4 s (peak 131.7 MB) |
| `run.py --full` at 2060 | **not run in this repository** — see below |

**Why no `--full` leg was run here.** The flag is offered and wired exactly as
in certs 11, 13, 14, 15, 19, 20, 21, 22 and 23 —
`certs/06-668-separation/full_recompute.py` imported by a `sys.path` insert
rather than copied, BLAS threads capped at three before numpy loads,
smoke-tested against the forced profile of Sylvester `H(128)` first — and
`--matrix` and `--impl` select the leg. The price is cert 22's, unchanged,
because the order and the module are the same: one 2060 leg is about **137×**
the 716 leg the same module took in this repository (cert 14, 400.3 s), i.e. of
order **15 hours** for a single `blas` matrix at three threads on this desk.
That 137× is the source laboratory's **measured** desk ratio
(`experiments/pr0042/REGISTRATION.md`, Amendment 1), which implies the exponent
`4.66`; on the `Θ(n⁵)` law quoted elsewhere here the same leg is
`(2060/716)⁵ = 197×`, i.e. ≈ **22 h**. Both are estimates and both say *hours*;
the smaller, measured one is quoted so that the decision not to run is not
defended with an inflated price. The `blas` route also materialises a
`C(n,2) × n` pair matrix — `2 120 770 × 2 060`, which is 4.4 GB as `int8` and
**17.5 GB** as the `float32` copy `_profile_blas` makes — far past this desk's
memory; the `bits` route (a `2 120 770 × 33` `uint64` packing) is the tractable
one at this order. So this certificate's default verdict is *banked exact
computation audited*, and it says so in its own output. Certs 06, 08, 11, 13,
14, 15 and 19 each have an in-repo `--full` `blas` leg on the record; certs 20,
21, 22, 23, 24 and 25 do not, and the word *replayed* is not used of them.

## What is NOT claimed

* **The default run recomputes nothing.** Only `--full` would bind a bank to a
  matrix by computation, and **no `--full` leg has been run in this repository
  at 2060**, for cert 22 or for this certificate.
* **Nothing about `G` vs `Gᵀ` at 2060.** `Gᵀ` was never enumerated. Cert 19
  decided the analogous question at 668 and at 668 only; nothing here extends
  it, and the plain-equivalence count is "at least five exhibited", not eight.
* **No general theorem about orientation.** That the GS orientation switch is a
  class of its own — and stays one with the transpose in the group — is now
  known at five orders (668, 716, 1676, 1772, 2060), one of them unbordered.
  Five orders are five orders.
* **No claim that three is the count** at 2060. It is a lower bound exhibited by
  the matrices in hand, under either convention.
* **Nothing about `ψ(ρ) = −1` at 2060.** The `×104` twist that relates `P` and
  `G` acts on the *column index* of each circulant, not on the seed values as
  Lemma T's `ψ`-twist does, so no analogue of the 668/716/1676/1772 twist
  theorem is asserted at this order.
* **No novelty or priority claim of any kind** at 2060. This certificate counts
  classes among the artifacts banked here; it says nothing about existence at
  2060 or about who first exhibited a Hadamard matrix of this order. `G` **is**
  the publicly posted matrix (`PROVENANCE.md`).
* **Matching invariants prove nothing.** `dim W` agrees on all five objects the
  4-profile proves pairwise inequivalent, and Paley `H(20)`'s profile agrees
  with its transpose's. Read every "agrees" in this repository as "did not
  separate", never as "the same".

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/24-transpose-extended-2060/run.py
```

The default path is standard library only, no network, and is the whole
certificate. The `--full` paths below import numpy (finder-side only, never in
the trust chain) and cap BLAS threads at three; at this order each is of order
15 hours, and the `blas` variants want memory this desk does not have. **None of
them has been run here.**

```
python certs/24-transpose-extended-2060/run.py --full --impl bits --matrix or-T
python certs/24-transpose-extended-2060/run.py --full
```

Exit code 0 iff every check passed.
