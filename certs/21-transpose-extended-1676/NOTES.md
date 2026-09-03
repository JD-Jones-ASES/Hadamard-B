# cert 21 — the three classes at order 1676 hold with the transpose in the group

**Label: PROVEN.** Default run:
`python certs/21-transpose-extended-1676/run.py` from the repository root.
Standard library only, **≈ 8.7 s**, exit 0, **173 checks**. That run **audits a
banked exact computation**; `--full` is offered and priced below and **has not
been run in this repository at this order** — cert 20's position at 1676,
unchanged and for the same reasons. The trust boundary is the one certs 15, 19
and 20 draw, and is set out below.

This certificate adds no new matrix and no new construction. It adds **two
exact 4-profiles** — the transposes of the Lemma-T rebuild and of the
orientation switch at order 1676 — and with them the statement cert 20
explicitly withheld: the **row-side only** caveat at this order is discharged.
This was the **first** such discharge in the repository. Certs 24 (2060) and 25
(1772) did the same on 2026-09-03, and with them **no separation statement in
`note/NOTE-B.md` is row-side any longer**; nothing below depends on either.

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
own — which is exactly why cert 20 stopped where it did.

## The theorem

> **Theorem (1676, transpose-extended).** Let `H` be the decoded
> `(s,i) = (1,1)` bordered Goethals–Seidel record at order 1676 (certs 01/20),
> `H′` its Lemma-T `i = 2` rebuild (certs 02/20), and `H″` the orientation
> switch of `H` — the twelve off-diagonal core blocks negated, the `4×4`
> border, the four row strips and the four column strips unchanged (cert 20).
> Then `H`, `H′` and `H″` are pairwise inequivalent **under the
> transpose-extended relation**. **Order 1676 therefore carries at least three
> Hadamard equivalence classes with the transpose in the group.**

*Proof.* Every leg is an exact `|T4|` 4-profile comparison over all
`C(1676,4) = 327 588 749 775` row 4-subsets, in two arithmetics that agree bin
for bin. `H ≈ H′` is refuted by `H ≁ H′` (**68** bins, cert 20) together with
`H ≁ (H′)ᵀ` (**139**, new here); `H ≈ H″` by `H ≁ H″` (**70**, cert 20) and
`H ≁ (H″)ᵀ` (**139**, new); `H′ ≈ H″` by `H′ ≁ H″` (**66**, cert 20) and
`H′ ≁ (H″)ᵀ` (**139**, new) — and, by the other route, `(H′)ᵀ ≁ H″` (**139**,
new). An invariant that differs is a separation. ∎

**Remark, recorded and not headlined.** Under **plain** Hadamard equivalence
the five matrices profiled at this order — `H, H′, H″, (H′)ᵀ, (H″)ᵀ` — are
**pairwise inequivalent**: all ten profile comparisons separate, the least
separated pair by **66** bins. **At least five classes are therefore exhibited
at 1676** by three constructions and transposition.
**The house counts three.** The transpose-extended relation is the one under
which a matrix and its transpose are the same object, and three is the count
that survives *either* convention — so it is the number stated in the theorem,
the note and the README. The remark is recorded because it is true and because
it is what makes the transpose-extended statement non-trivial: `H′` and `H″`
are each inequivalent to their own transposes here (139 bins apiece), so the
transpose legs were genuinely open until they were computed. Certs 15 and 19
record the same kind of remark at 716 (six matrices) and 668 (eight).

## What is **not** profiled, and what therefore is not said

**`Hᵀ` — the decoded record's own transpose at order 1676 — was never
enumerated.** It is not in the campaign's object list
(`experiments/pr0042/REGISTRATION.md` §2 names `H_1676-twisted-T` and
`H_1676-orient-T` and no third 1676 transpose), it is not banked here, and
`run.py` does not build or pin it. Three consequences, all of them stated in
the run's own output:

* **Nothing is claimed about `H` vs `Hᵀ` at 1676.** Cert 19 decided that
  question at **668 only** (`H ≁ Hᵀ`, 49 of 80 bins). Whether it goes the same
  way here is **open**, and no analogue is asserted.
* The plain-equivalence count above is "**at least five exhibited**", not
  eight. Five matrices are profiled; a sixth (`Hᵀ`) exists and is unprofiled,
  so nothing is said about how it sits against the other five.
* In the verdict table the third column, `Aᵀ vs B`, reads **n/a** for the two
  pairs whose `A` is `H`. Those pairs are refuted by their middle column
  (`A` vs `Bᵀ`) alone, which is all the relation requires; the third column is
  the redundant crossing, available only where `Aᵀ` is banked. `run.py`
  asserts that the n/a entries are exactly the ones with `A = H`.

None of this weakens the theorem: each of the three pairs has both of its
required refutations in hand.

## The trust boundary — what a default run does and does not establish

The two new `C(1676,4)` enumerations **were not run inside this repository.**
They ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py` (the engine unchanged since the 2060
registration, and the same one certs 11, 13, 14, 15, 19 and 20 rest on), numpy,
**16 threads on a rented `c2d-highcpu-16` (`prof42-1`, `us-east1-b`),
2026-09-02** — under the pre-registration
`experiments/pr0042/REGISTRATION.md`, **flushed 10:17 UTC, before any matrix it
governs was built**. That registration fixes the objects (§2), the decision
rules (§4: "differs in any bin ⟹ inequivalent"; "equal in every bin ⟹
MEASURED, nothing proved"; and, in as many words, `H ≁ H′ᵀ`, `H ≁ H″ᵀ`,
`H′ ≁ H″ᵀ` ⟹ the three-class statement holds under the transpose-extended
relation) and the kill criteria (§5: `blas ≠ bits` in any bin is a hard stop;
a builder digest mismatch ⟹ nothing uploaded) in advance.

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
* `H″` is formed from those rows by negating exactly `12·n² = 2 096 688` cells,
  satisfies the alternate-orientation identity `S·H″·S = H_alt` cell by cell,
  and carries cert 20's digest;
* `H′` rebuilds with its seeds **re-derived here** as the `ψ`-twist of the
  decoded seeds (`ψ(g) = (−1)^g` on `ℤ₄₁₈`; `ρ = 417` is odd, so
  `ψ(ρ) = −1`), compared character for character against the bank, so cert 02's
  shared record is bound to `payload-records.json` by computation rather than
  by a file pin;
* `(H′)ᵀ` and `(H″)ᵀ` are **transposed in-process** from those rows, pass
  `verify/verify.py`, and carry the two new pins — so each new bank's declared
  matrix is bound to a matrix **this run built**, not merely to a name;
* transposition is checked as an involution and, on both objects, as a genuine
  move: `(Mᵀ)ᵀ = M` cell for cell, and `Mᵀ ≠ M` in **all 1 676** rows;
* the five matrices carry five distinct canonical digests, established before
  any profile is opened;
* every banked profile satisfies the forced identities — 142 bins, all
  `≡ 4 (mod 8)` (`1676 ≡ 4 mod 8`), total `C(1676,4) = 327 588 749 775`, second
  moment `n³(n−1)(n−2)/24 = 550 023 273 154 800` — recomputed **and** compared
  against the `second_moment`, `total`, `n` and `C_n_4` fields the bank
  declares, alongside its `schema`, its `folded` field (the **signed** `T4`
  histogram is not an invariant — §3.1) and its `impl`;
* each bank's declared matrix name and canonical digest against the in-process
  digest of the matrix rebuilt in **this** run;
* `blas == bits` bin for bin on each of the **five** matrices, so every leg of
  the theorem, the transpose legs included, rests on two independent
  implementations;
* all **ten** pair comparisons, with their counts, union support sizes and
  first moments, in **both** arithmetics; and the transpose-extended verdicts
  derived in code from those counts, with the second route shown wherever the
  transpose of `A` is banked;
* controls C0–C6b.

**A default run does not establish that the banked histograms were computed
from the matrices `run.py` rebuilt.** They are *producer-banked*: the digest
each carries is the one the engine recorded against the file it enumerated, and
it equals the digest this certificate pins and the digest computed in-process
here. A self-declared digest is metadata, not a computation. `--full` is what
would close that gap — and at this order it has not been run (see *Runtimes*).

## The evidence chain

**[0]** ten file pins. **[1]** rebuild `H` (hypotheses H0–H4, D1/D3/D5, `Σ̄`,
compression lemma), verify, pin; check the layout `n = |G| = 418`, `s = 1`,
`N = 4(n+s)`; form `H″`, count the negated cells, check `S·H″·S = H_alt`,
verify, pin; re-derive the `ψ`-twist seeds and rebuild `H′`, verify, pin;
transpose `H′` and `H″`, check the involution and that the transpose is a
genuine move, verify, pin; five distinct digests. `Hᵀ` deliberately absent.
**[1b]** control C6 — the dim-`V` trap on the real objects. **[2]** ten banks
audited in exact integers; header, schema, folding and matrix identity per
bank; `blas == bits` ×5. **[3]** the ten pair comparisons, each asserted to the
exact bin count and union size, deltas summing to zero, a non-zero first
moment; for the three legs the theorem rests on, the tail that does *not*
separate and the first eight divergent bins printed. **[4]** the
transpose-extended verdicts, derived from those counts, then the plain-
equivalence remark and the explicit `Hᵀ`-is-absent assertion. **[5]** controls
C0–C6b. **[6]** `--full`: `certs/06-668-separation/full_recompute.py` imported
(not copied), smoke-tested on the forced profile of Sylvester `H(128)`, then
run on the verified transposed rows — offered, priced, **not run here**.

**Controls.**

| | control | result |
| --- | --- | --- |
| **C0** | every control matrix is what it claims | the GS condition `Σ_q PAF_q(g) = 4v·[g=0]` re-verified for both GS seeds; all five controls checked Hadamard by brute force; the orientation switch applied to GS `H(28)` is still Hadamard and moves exactly `12·7²` cells |
| **C1** | five small Hadamard matrices profiled two ways, one of them the route `--full` takes | straight `O(C(n,4))` enumeration `==` the pair-vector/popcount route on Sylvester `H(8)`, `H(16)`, Paley I `H(20)`, GS `H(28)`, `H(36)`; the two Sylvester profiles match their **forced** values; every bin `≡ n (mod 8)` |
| **C2** | the transposed-profile route, on matrices small enough for straight enumeration | Sylvester `H(8)` and `H(16)` are **symmetric**, so `profile(Mᵀ) = profile(M)` is *forced* — and holds; Paley I `H(20)` is not symmetric and is **MEASURED**, never asserted (its profiles happen to agree) |
| **C3** | the transposes are genuinely different objects | at 1676 they do **not** drop a bin as at 668 and 716 (cert 15, C3) — they **swap two**: the three originals share one 142-bin support, the two transposes share another, `\|T4\| = 948, 964` are populated in every original and no transpose (counts 30, 40) and `\|T4\| = 1180, 1204` in every transpose and no original (counts 30, 30). Union support: **144** bins |
| **C4** | the comparator in the null direction | every banked profile against itself: 0 differing bins, 10 times |
| **C5** | *negative* control: a **total-preserving** corruption of the new `(H″)ᵀ` bank (one count moved from `\|T4\| = 4` to `\|T4\| = 1668`), which only the second-moment identity can catch | the assert fires; the corrupted profile still totals `C(1676,4)` |
| **C6** | the dim-`V` trap on the real objects | `dim W` = **1675** on all five — an invariant that separates **none** of the ten pairs clause [4] proves inequivalent; `dim V` reads 1674 / 1675 / 1674 / 1675 / 1674 on `H`, `H′`, `H″`, `(H′)ᵀ`, `(H″)ᵀ` and is worthless (Trap 1, §3.1) |
| **C6b** | the same trap demonstrated | under a seeded signed row negation of Sylvester `H(16)`, `dim V` moves 4 → 5, `dim W` does not move, and the `|T4|` profile does not move either |

## The separations

All ten pairs, as `run.py` prints them; the counts are identical in both
arithmetics.

| pair | differing | union | `Σ \|T4\|·Δ` | source |
| --- | --- | --- | --- | --- |
| `H` vs `H′` | 68 | 142 | −20 505 840 | cert 20 |
| `H` vs `H″` | 70 | 142 | +11 699 520 | cert 20 |
| `H′` vs `H″` | 66 | 142 | +32 205 360 | cert 20 |
| `H` vs `(H′)ᵀ` | **139** | 144 | −62 530 104 | new |
| `H` vs `(H″)ᵀ` | **139** | 144 | −65 482 040 | new |
| `H′` vs `(H′)ᵀ` | **139** | 144 | −42 024 264 | new |
| `H′` vs `(H″)ᵀ` | **139** | 144 | −44 976 200 | new |
| `H″` vs `(H′)ᵀ` | **139** | 144 | −74 229 624 | new |
| `H″` vs `(H″)ᵀ` | **139** | 144 | −77 181 560 | new |
| `(H′)ᵀ` vs `(H″)ᵀ` | **70** | 142 | −2 951 936 | new |

The verdicts, derived from those counts:

| pair | `A` vs `B` | `A` vs `Bᵀ` | `Aᵀ` vs `B` | verdict |
| --- | --- | --- | --- | --- |
| `H` vs `H′` | 68 | 139 | n/a | **SEPARATED** |
| `H` vs `H″` | 70 | 139 | n/a | **SEPARATED** |
| `H′` vs `H″` | 66 | 139 | 139 | **SEPARATED** |

The two `n/a` entries are `Hᵀ`, which was never profiled; each of those pairs
is refuted by its middle column alone, which is all the relation asks. The
third row shows both routes, and both hold.

### The support: a swap, not a drop

At 668 and 716 every transpose populated **one bin fewer** than its original
and it was the same bin at each order (cert 15, control C3). **At 1676 the
count does not change** — every one of the five profiles populates 142 bins —
but the *support* does. The three originals share one 142-bin support and the
two transposes share another, and the two differ in two bins each way:

| bin | `H` | `H′` | `H″` | `(H′)ᵀ` | `(H″)ᵀ` |
| --- | --- | --- | --- | --- | --- |
| `\|T4\| = 948` | 30 | 30 | 30 | — | — |
| `\|T4\| = 964` | 40 | 40 | 40 | — | — |
| `\|T4\| = 1180` | — | — | — | 30 | 30 |
| `\|T4\| = 1204` | — | — | — | 30 | 30 |

So the union support of the five profiles is **144** bins, which is why every
original-versus-transpose comparison is reported as *139 of 144* and every
same-side comparison as *n of 142*. Four of those five sparse bins carry counts
of 30 or 40 out of `3.3·10¹¹` 4-subsets; the separations do not rest on them —
the bulk does — but they are a clean structural fact and `run.py` asserts them.

### The three new legs the theorem rests on

`run.py` prints the first eight divergent bins of each; the full 139-bin lists
are the banked JSONs.

**`H` vs `(H′)ᵀ`** — 139 of 144 bins differ; largest `|Δ| = 5 195 974` at
`|T4| = 36`, i.e. `1.55·10⁻⁴` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H` | 55 147 198 687 | 52 329 871 149 | 47 215 439 663 | 40 653 198 021 | 33 547 376 838 | 26 657 993 721 |
| `(H′)ᵀ` | 55 150 948 176 | 52 326 378 631 | 47 215 018 097 | 40 653 295 278 | 33 552 572 812 | 26 656 410 724 |
| Δ | +3 749 489 | −3 492 518 | −421 566 | +97 257 | **+5 195 974** | −1 582 997 |

**`H` vs `(H″)ᵀ`** — 139 of 144; largest `|Δ| = 5 778 194` at `|T4| = 36`,
i.e. `1.72·10⁻⁴` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H` | 55 147 198 687 | 52 329 871 149 | 47 215 439 663 | 40 653 198 021 | 33 547 376 838 | 26 657 993 721 |
| `(H″)ᵀ` | 55 150 769 956 | 52 326 344 165 | 47 216 280 093 | 40 652 388 966 | 33 553 155 032 | 26 655 531 362 |
| Δ | +3 571 269 | −3 526 984 | +840 430 | −809 055 | **+5 778 194** | −2 462 359 |

**`H′` vs `(H″)ᵀ`** — 139 of 144; largest `|Δ| = 5 624 316` at `|T4| = 36`,
i.e. `1.68·10⁻⁴` of that bin. (The fourth leg, `(H′)ᵀ` vs `H″`, is the other
route for the same pair: 139 bins, largest `|Δ| = 5 843 843` at `|T4| = 4`.)

Every difference vector sums to zero, as it must; the first moment, which
nothing forces, does not — and all ten are non-zero. The **bulk** separates and
the extreme tail does not: on every original-versus-transpose comparison the
top differing bin is `|T4| = 1204` and the four bins above it — `1228`, `1252`,
`1268` and `1668` — agree exactly (as does `1116`, the fifth of the five
agreeing bins); on `(H′)ᵀ` vs `(H″)ᵀ` every bin above `|T4| = 564` agrees, 71
of them. As everywhere else in this repository, the largest discrepancies are
of order `10⁻⁴`–`10⁻⁵` of their bins — invisible to any sample of practical
size, and nothing cheaper than the exact profile would have found them.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports) — all
five rebuilt in clause [1] of every run:

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `H` decoded `(1,1)` record | `8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99` | yes (cert 01's pin) |
| `H′` Lemma-T `i = 2` rebuild | `6a4938371ddbe4ad8bd35f21d7e61dad683b15f8f2ec1c88e88ce579c4907405` | yes (cert 02's pin; seeds re-derived) |
| `H″` orientation switch | `16d1617cc62532b26c010f3b174c741f0b9388089516759834030d9056a84346` | yes (cert 20's pin; formed from `H`) |
| `(H′)ᵀ` | `46f29432da47cb6e106a19d5c8e453e4bf797edb48924bfb0f9664a4af261ce3` | yes (new) |
| `(H″)ᵀ` | `ad870e1c216a76462c4913a667bb747847f64e2d197353d179f26a15ea85303d` | yes (new) |

The two transposed digests are the desk's own, recorded in
`experiments/pr0042/manifest.json` when the matrices were built and verified
there, and re-derived in-process by every run of this certificate. **`Hᵀ` is
not in this table**, on purpose: no profile of it exists.

**Banked files** (SHA-256 of the file bytes, compared in `run.py`). The two
`-T-` pairs are this certificate's own, banked 2026-09-02 under cert 21; the
six others are cert 20's, reused verbatim and re-pinned here:

| file | SHA-256 | banked for |
| --- | --- | --- |
| `data/sep1676-twisted-T-exact-blas.json` | `d951a78a0af94cb3979eec3878583d750ecb2c417946ab4c6159a3a7bfaaa040` | cert 21 |
| `data/sep1676-twisted-T-exact-bits.json` | `2a249d15e2d8c0aaa71d0fbeb503da43175aed7b9466c32dfaa77db1401bad6c` | cert 21 |
| `data/sep1676-orient-T-exact-blas.json` | `4d7506d93e35125d236432d92457d60321c1e1c4f014febc45e5ded0b06b820a` | cert 21 |
| `data/sep1676-orient-T-exact-bits.json` | `46750b71d6400a51ddfde19ec574cc6a8d92c345aabc521f5cd866ca3eb84e41` | cert 21 |
| `data/sep1676-decoded-exact-blas.json` | `57b9a43caf5246de779ad3205a45642c98f7a211be47e4ed12d718fe098781c9` | cert 20 |
| `data/sep1676-decoded-exact-bits.json` | `469e0b0382d479a6d917316246807cadfa1f113bb2bfcd1429ec1712622e7b94` | cert 20 |
| `data/sep1676-twisted-exact-blas.json` | `328ed05c9614a223d95bd35583c83433a8700249c0b50872fbfc8d846e9b5a49` | cert 20 |
| `data/sep1676-twisted-exact-bits.json` | `311ef88606d0543967e2b0cf46aad4f3fb3f1353cc59b018c5369e447c0c2bb1` | cert 20 |
| `data/sep1676-orient-exact-blas.json` | `a83b239695a3bd820de222e829e65a10a5dd66a432858af57cc950eb4ff40be2` | cert 20 |
| `data/sep1676-orient-exact-bits.json` | `af198e51aecd165e8a2a22ee5ece8dfa73d8ddedf314fa94684ee367db14e9d5` | cert 20 |

`data/payload-records.json` and `data/twisted-i2-records.json` are not
file-pinned here, for cert 20's reason: they are shared with certs 01, 02 and
20, and the binding pin on each is the canonical digest of the matrix it
produces, checked in clause [1] — reinforced, for the twisted record, by the
`ψ`-twist re-derivation that binds it to `payload-records.json` outright.

## Provenance of the two new profiles

Produced by the source laboratory's unchanged engine
(`Hadamard-2060`, `experiments/inequiv/exact_profile_big.py`) on the rented
`c2d-highcpu-16` `prof42-1` (`us-east1-b`), 16 threads, 2026-09-02, under
`experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC; Amendment 1 ~11:05 UTC
re-priced the campaign). The two `blas` legs landed 14:18Z; the `orient-T`
`bits` leg 21:40Z; the `twisted-T` `bits` leg was the last of the 1676
worklist. The seconds and peak resident sets below are the values those runs
recorded in the JSONs themselves.

Banked into `data/` by the laboratory's `experiments/pr0042/bank.py` under
`--cert 21`, which refuses unless both implementations audit, agree bin for
bin, and declare the digest `manifest.json` pins for that matrix. It adds the
seven header fields each bank's own `banked_note` names — `schema`,
`description`, `matrix`, `matrix_canonical_sha256`, `producer_filename`,
`arithmetic`, `banked_note` — and **nothing numeric**: every numeric field is
the producer's own output, unaltered. `run.py` checks `schema`, `matrix`,
`matrix_canonical_sha256`, `impl`, `folded`, `n`, `C_n_4`, `total` and
`second_moment` against values it recomputes or rebuilds.

Credit, as everywhere in this repository, is to stations: the theorems and the
certificate are this repository's; the engine, the pre-registration, the
matrix building and verification, and the rented enumeration are the source
laboratory's; `payload-records.json` encodes the posting team's mathematical
content (`PROVENANCE.md`), and no priority of any kind is claimed on it.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 8.7 s** (exit 0, 173 checks; measured 2026-09-02 on the desk, three runs, 8.7 / 8.7 / 8.7 s) |
| rebuild + full hypothesis re-check + `verify.py`, per matrix | 1.4 s |
| the orientation switch, the cell count and the `S·H″·S` identity at `N = 1676` | 0.9 s |
| forming a transpose at `N = 1676` (the `verify.py` call dominates) | < 0.1 s |
| `dim V` / `dim W` on all five 1676 matrices | ≈ 1.9 s |
| producing `sep1676-twisted-T-exact-{blas,bits}.json` | 2 650.8 s (peak 762.2 MB) / 6 661.6 s (peak 104.5 MB) |
| producing `sep1676-orient-T-exact-{blas,bits}.json` | 2 665.8 s (peak 762.3 MB) / 6 660.6 s (peak 104.5 MB) |
| `run.py --full` at 1676 | **not run in this repository** — see below |

**Why no `--full` leg was run here.** The flag is offered and wired exactly as
in certs 11, 13, 14, 15, 19 and 20 — `certs/06-668-separation/full_recompute.py`
imported by a `sys.path` insert rather than copied, BLAS threads capped at
three before numpy loads, smoke-tested against the forced profile of Sylvester
`H(128)` first — and `--matrix` and `--impl` select the leg. The price is
cert 20's, unchanged, because the order and the module are the same: one 1676
leg is about **52×** the 716 leg the same module took in this repository
(cert 14, 400.3 s), i.e. roughly **6–7 hours** for a single `blas` matrix at
three threads on this desk. That 52× is the source laboratory's **measured**
sub-`n⁵` scaling — its desk-measured 716→2060 ratio of 137
(`experiments/pr0042/REGISTRATION.md`, Amendment 1) implies the exponent
`4.66`, and `(1676/716)^4.66 = 52`. On the `Θ(n⁵)` law quoted elsewhere here
the same leg is `(1676/716)⁵ = 70×`, i.e. ≈ **7.8 h**. Both are estimates and
both say *hours*; the smaller, measured one is quoted so that the decision not
to run is not defended with an inflated price. The `blas` route also
materialises a `C(n,2) × n` pair matrix — `1 403 650 × 1 676`, which is 2.35 GB
as `int8` and **9.4 GB** as the `float32` copy `_profile_blas` makes — past
this desk's memory; the `bits` route (a `1 403 650 × 27` `uint64` packing) is
the tractable one at this order. So this certificate's default verdict is
*banked exact computation audited*, and it says so in its own output. Certs 06,
08, 11, 13, 14, 15 and 19 each have an in-repo `--full` `blas` leg on the
record; certs 20 and 21 do not, and the word *replayed* is not used of them.

## What is NOT claimed

* **The default run recomputes nothing.** Only `--full` would bind a bank to a
  matrix by computation, and **no `--full` leg has been run in this repository
  at 1676**, for cert 20 or for this certificate.
* **Nothing about `H` vs `Hᵀ` at 1676.** `Hᵀ` was never enumerated. Cert 19
  decided the analogous question at 668 and at 668 only; nothing here extends
  it, and the plain-equivalence count is "at least five exhibited", not eight.
* **Nothing at 1772 or 2060.** The same registration governs the orientation
  switches and transposed profiles at those orders; they are separate
  computations and carry their own certificates. Order 2060 was, as of this
  certificate, the **only** order in `note/NOTE-B.md` whose separation
  statement was row-side only. *(All of that has since been settled, and
  nothing here rests on any of it: **certs 22** at 2060 and **23** at 1772
  made the row-side three-class statements on 2026-09-03, and **certs 24**
  and **25** discharged both caveats the same day — 2060 by `Pᵀ` vs `G`
  **134 of 134**, `P` vs `(H″)ᵀ` 145 of 146 and `G` vs `(H″)ᵀ` **134 of
  134**; 1772 by **91 of 92** on every original-versus-transpose leg. So
  **1676 was the first order where a row-side caveat was discharged — this
  certificate — and since 2026-09-03 the note has none left.**)*
* **No general theorem.** That the `ψ(ρ) = −1` twist and the orientation switch
  land in different classes — and stay apart with the transpose in the group —
  is now known at three orders (668, 716, 1676). Three orders are three orders.
* **No claim that three is the count** at 1676. It is a lower bound exhibited
  by the matrices in hand, under either convention.
* **No novelty or priority claim of any kind** at 1676. This certificate counts
  classes among the artifacts banked here; it says nothing about who first
  exhibited a Hadamard matrix of this order. The decoded record is a public
  artifact (`PROVENANCE.md`).
* **Matching invariants prove nothing.** `dim W` agrees on all five objects the
  4-profile proves pairwise inequivalent, and Paley `H(20)`'s profile agrees
  with its transpose's. Read every "agrees" in this repository as "did not
  separate", never as "the same".

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/21-transpose-extended-1676/run.py
```

The default path is standard library only, no network, and is the whole
certificate. The `--full` paths below import numpy (finder-side only, never in
the trust chain) and cap BLAS threads at three; at this order each is hours,
and the `blas` variants want memory this desk does not have. **None of them has
been run here.**

```
python certs/21-transpose-extended-1676/run.py --full --impl bits --matrix or-T
python certs/21-transpose-extended-1676/run.py --full
```

Exit code 0 iff every check passed.
