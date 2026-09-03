# cert 22 — order 2060 carries at least three equivalence classes

**Label: PROVEN + PROVEN-BY-CERTIFICATE. ROW-SIDE ONLY.** Default run:
`python certs/22-2060-three-classes/run.py` from the repository root.
Standard library only, **≈ 7.5 s**, exit 0, **137 checks**. That run **audits a
banked exact computation**; `--full` is offered and priced below and **has not
been run in this repository at this order** — certs 20 and 21's position at
1676, for the same reasons and at a larger price. The trust boundary is the one
certs 15, 19, 20 and 21 draw, and is set out below.

This certificate adds **one matrix and one exact 4-profile** at the founding
order: the orientation switch of the plain Goethals–Seidel realisation, in
both arithmetics. With cert 07's pair it upgrades order 2060 from **two**
classes to **three**.

---

## The theorem

> **Theorem (2060, row-side).** Let `P` be `2060-plain`, the plain
> Goethals–Seidel array over the raw `ℤ₅₁₅` seed, and `G` be `2060-gist`, the
> `×104`-twisted array that is byte-for-byte the publicly posted `H(2060)`
> (both cert 07, which proves `P ≁ G`). Let `H″` be `P` with its **twelve
> off-diagonal `515×515` blocks negated** — the orientation switch, in its
> **unbordered** form. Then `H″` is a Hadamard matrix, and it is
> Hadamard-inequivalent to each of `P` and `G`. **Order 2060 therefore carries
> at least three Hadamard equivalence classes.**

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all
`C(2060,4) = 748 155 697 135` row 4-subsets, with
`T4 = Σ_c H[i][c]H[j][c]H[k][c]H[l][c]`, is a Hadamard-equivalence invariant
(`note/NOTE-B.md` §3.1, invariant **I5**: each row negation contributes one
sign to `T4` so `|T4|` is fixed, each column negation contributes `d_c⁴ = 1`,
and permutations relabel). `H″` and `P` populate the **same 145 bins** and
**107** of those bin counts differ; `H″` and `G` have union support **147** and
**146** of those bins differ. (Cert 07's own leg, `P` against `G`, is 146 of
147 and is re-derived here.) Every profile totals `C(2060,4)` exactly and hits
the second moment `n³(n−1)(n−2)/24 = 1 543 448 476 598 000` to the unit, in two
arithmetics that agree bin for bin on each matrix. An invariant that differs is
a separation. ∎

### The unbordered orientation switch

At 668, 716 and 1676 the array is **bordered**: `N = 4(n+s)` with `s = 1`, and
certs 13, 14 and 20 negate the twelve off-diagonal **core** blocks and leave
the `4×4` corner and the eight border strips alone. **At 2060 there is no
border.** The order sits at the degenerate `s = 0` layer of the master theorem
— `N = 4·515`, a plain `4×4` Goethals–Seidel array of circulant `515`-blocks
over `ℤ₅₁₅` — so the switch negates twelve blocks and there is nothing to leave
alone. That is what `experiments/pr0042/REGISTRATION.md` §2 registered
(`H_2060-orient`: "the plain GS realisation (cert 07 pin) with its twelve
off-diagonal `515×515` blocks negated (no border, `s = 0`)"), and it is what
`run.py` builds.

Exactly `12·515² = 3 182 700` cells change.

### What `H″` is — checked two ways

`S·H″·S`, with `S = diag(1,−1,−1,−1) ⊗ I₅₁₅`, is the same four seeds assembled
in the **alternate Goethals–Seidel orientation** — the six transposed blocks
negated (`note/NOTE-B.md` §1.0, where the standard orientation is called a
*convention of this repository*). Certs 13, 14 and 20 check that as a
sign-pattern identity, cell by cell, carrying a border-signing term
(`P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`) that is empty here. Because at 2060
the **seeds themselves are in hand**, this certificate checks the identity
twice:

* **(a)** as the sign-pattern identity of certs 13/14/20 with the border term
  gone: `S[r]·S[c]·H″[r][c] = alt(r,c)·P[r][c]` for every one of the
  `2060² = 4 243 600` cells, with `alt(r,c) = −1` exactly on the six transposed
  blocks;
* **(b)** by **assembling the alternate array directly from the same raw
  seeds** —
  `[[A, BR, CR, DR], [−BR, A, −DᵀR, CᵀR], [−CR, DᵀR, A, −BᵀR], [−DR, −CᵀR, BᵀR, A]]`,
  §1.0's standard array with its six transposed blocks negated — putting it
  through `verify/verify.py` (HADAMARD, canonical digest
  `40e1d1c8…ca64398d`), and comparing it to `S·H″·S` cell for cell.

**(b) is what (a) can only assert.** It exhibits `H″` as a **signed
conjugation of an independently constructed Goethals–Seidel array**: not
"twelve blocks flipped" but the other orientation of the founding seed, and `S`
is a diagonal `±1` on each side, i.e. an element of the equivalence group,
which is why `H″` and the alternate array have the same profile and the theorem
is a statement about orientation rather than about an arbitrary sign pattern.
The same pair of checks is exercised on the GS controls `H(28)` and `H(36)`,
which are unbordered too — the **exact structural analogue** of the 2060 switch,
unlike the bordered controls certs 13/14/20 could offer.

**Remark, as at 668, 716 and 1676.** The GS orientation, a "convention" in
§1.0, is **not a gauge for Hadamard equivalence at the founding order
either**. One seed quadruple; three classes; **two constructions** — the plain
array and the `×104` twist that *is* the public artifact — **plus the
orientation switch**. Four orders now carry that verdict (668, 716, 1676,
2060). Four orders are four orders: **no general theorem** about orientation is
claimed, here or anywhere in this repository.

## Row-side only — the caveat this certificate keeps

**Nothing here is claimed under the transpose-extended relation.** The
transposed profiles at 2060 — `H_2060-orient-T` and `H_2060-plain-T`, both in
the campaign's object list (`REGISTRATION.md` §2) — are **separate legs of the
same campaign, still running when this certificate was written**, and neither
is banked here. Transpose is **not** in the equivalence group
(`note/NOTE-B.md` §3), so a row-side separation says nothing about `A ≈ B` on
its own; refuting `A ≈ B` needs `profile(A) ≠ profile(B)` **and**
`profile(A) ≠ profile(Bᵀ)`, and the second refutation does not exist here for
any of the three pairs.

This is **exactly cert 20's caveat at 1676**, which cert 21 later discharged
there once the two transposed profiles landed. Order 2060 was the only
row-side statement in `note/NOTE-B.md` before this certificate and remains the
only one after it — the count changes from two to three; the caveat does not
move. `run.py` asserts in clause [4] that no transposed 2060 profile is among
the banks it reads, so the caveat cannot be quietly dropped by adding a file.

Nor is anything claimed about `P` versus `Pᵀ`, or `G` versus `Gᵀ`, or `H″`
versus `(H″)ᵀ`: cert 19 decided the analogous question at **668 only**.

## The trust boundary — what a default run does and does not establish

The `C(2060,4)` enumeration of `H″` **was not run inside this repository.** It
ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py` (the engine unchanged since the 2060
registration, and the same one certs 07, 11, 13, 14, 15, 19, 20 and 21 rest
on), numpy, **16 threads on a rented `c2d-highcpu-16` (`prof42-2`,
`us-east1-b`), 2026-09-02/03** — under the pre-registration
`experiments/pr0042/REGISTRATION.md`, **flushed 10:17 UTC on 2026-09-02,
before any matrix it governs was built** (Amendment 1, ~11:05 UTC, added the
second instance for the 2060 legs and re-priced the campaign; both predate the
leg). That registration fixes the object (§2), the decision rules (§4, in as
many words: "`H″ ≁ plain` and `H″ ≁ gist` ⟹ order 2060 carries at least THREE
classes (row-side)"; "differs in any bin ⟹ inequivalent"; "equal in every bin
⟹ MEASURED, nothing proved") and the kill criteria (§5: `blas ≠ bits` in any
bin is a hard stop; a builder digest mismatch ⟹ nothing uploaded) in advance.
The four comparison profiles are cert 07's, produced under
`REGISTRATION-2060-exact.md` and reused here verbatim.

The division of labour matters and is stated: **the matrix was built and
verified at the desk** (`experiments/pr0042/build_matrices.py`, digest in its
`manifest.json`, through this repository's `verify/verify.py`, with the plain
source re-verified against cert 07's pin first) and only then uploaded; the
rented machine **enumerated and did nothing else**.

**A default `run.py` establishes:**

* the six bank files are byte-for-byte the ones pinned in `run.py` (SHA-256
  compared in code);
* the seed record is the order-2060 pair, and its four raw seeds satisfy the
  classical Goethals–Seidel condition `Σ_q PAF_q(t) = 4v·[t = 0]` — the `s = 0`
  layer of the master theorem, **re-verified here on all 515 shifts**, not
  assumed — with the record's declared row sums re-derived from the seeds;
* `P` and `G` rebuild from that record and pass `verify/verify.py` with cert
  07's canonical digests, computed in-process;
* `H″` is formed from `P`'s rows by negating exactly `12·515² = 3 182 700`
  cells, satisfies the alternate-orientation identity **both ways** (the cell-
  by-cell sign pattern, and equality with the independently assembled alternate
  array), and passes `verify/verify.py` with the desk's digest
  `4e1891b0…0870b801`;
* the alternate array itself passes `verify/verify.py` and is a fourth,
  distinct object; the four matrices carry four distinct canonical digests;
* every banked profile satisfies the forced identities — 145 / 145 / 133 bins,
  all `≡ 4 (mod 8)` (`2060 ≡ 4 mod 8`), every key canonical in `[0, 2060]` and
  every count a positive integer, total `C(2060,4) = 748 155 697 135`, second
  moment `1 543 448 476 598 000` — recomputed **and** compared against the
  `second_moment`, `total`, `n` (and, where the schema has it,
  `second_moment_want` and `C_n_4`) fields the bank declares, alongside its
  `schema`, its `folded` field (the **signed** `T4` histogram is not an
  invariant — §3.1), its `impl`, and its own name for the matrix;
* each bank's declared matrix digest against the **in-process** digest of the
  matrix rebuilt in **this** run, and the producer's own `matrix_sha256`
  against that declared digest;
* `blas == bits` bin for bin on each of the **three** matrices;
* all three pair comparisons, with their counts, union support sizes and first
  moments, in **both** arithmetics; the support structure; and the theorem
  derived in code from those counts, with the row-side caveat asserted;
* controls C0–C7.

**A default run does not establish that the banked histograms were computed
from the matrices `run.py` rebuilt.** They are *producer-banked*: the digest
each carries is the one the engine recorded against the file it enumerated, and
it equals the digest this certificate pins and the digest computed in-process
here. A self-declared digest is metadata, not a computation. `--full` is what
would close that gap — and at this order it has not been run (see *Runtimes*).

### Two schemas, each checked on its own terms

The four cert-07 banks predate the campaign and declare
`schema: "sep2060-exact-profile/1"` with `tag`, `producer` and
`producer_matrix_sha256`; the two new ones declare
`schema: "exact-4-profile/1"` with `matrix`, `matrix_sha256`, `C_n_4`,
`second_moment_want`, `arithmetic` and `banked_note`. `run.py` checks each set
against its own schema and refuses either if a field is absent, mistyped or
disagrees with what it recomputes; it does **not** accept a missing header as a
pass.

## The evidence chain

**[0]** six file pins. **[1]** the seed record, the GS condition on the raw
seeds, the declared row sums; rebuild `P` and `G`, verify, pin against cert
07's digests; form `H″`, count the negated cells, check the alternate-
orientation identity as a sign pattern **and** against the independently
assembled alternate array, verify both, pin; four distinct digests. **[1b]**
control C7 — the dim-`V` trap on the real 2060 objects. **[2]** six banks
audited in exact integers; identities, headers, schema, folding, arithmetic,
matrix name and matrix identity per bank; `blas == bits` ×3. **[3]** the three
pair comparisons, each asserted to the exact bin count and union size, deltas
summing to zero, a non-zero first moment; the support structure; where the
invariant does and does not separate; the first eight divergent bins of the two
new legs. **[4]** the theorem derived from those counts, and the row-side
caveat asserted. **[5]** controls C0–C5. **[6]** `--full`:
`certs/06-668-separation/full_recompute.py` imported (not copied), smoke-tested
on the forced profile of Sylvester `H(128)`, then run on the verified rows —
offered, priced, **not run here**.

**Controls.**

| | control | result |
| --- | --- | --- |
| **C0** | every control matrix is what it claims | the GS condition `Σ_q PAF_q(t) = 4v·[t=0]` re-verified for both GS seed sets **and for the 2060 seeds themselves**; all five small controls checked Hadamard by brute force |
| **C1** | five small Hadamard matrices profiled two ways, one of them the route `--full` takes | straight `O(C(n,4))` enumeration `==` the pair-vector/popcount route on Sylvester `H(8)`, `H(16)`, Paley I `H(20)`, GS `H(28)`, `H(36)`; the two Sylvester profiles match their **forced** values; every bin `≡ n (mod 8)` |
| **C2** | the **unbordered** orientation switch, end to end, on GS `H(28)` and `H(36)` — `s = 0` there too, so the exact structural analogue of the 2060 switch | the switch moves exactly `12v²` cells; the result is still Hadamard by brute force; identity **(a)** holds cell by cell; identity **(b)** holds — `S·H″·S` *is* the alternate GS array over the same seeds; that array is Hadamard too; and `profile(H″) = profile(alt)`, as a conjugation inside the group must give |
| **C3** | the comparator in the null direction | every banked profile against itself: 0 differing bins, 6 times |
| **C4** | *negative* control: a **total-preserving** corruption of the new `H″` bank (one count moved from `\|T4\| = 4` to `\|T4\| = 1236`), which only the second-moment identity can catch | the assert fires; the corrupted profile still totals `C(2060,4)` |
| **C5** | the dim-`V` trap demonstrated | under a seeded signed row negation of Sylvester `H(16)`, `dim V` moves 4 → 5, `dim W` does not move, and the `|T4|` profile does not move either |
| **C7** | the dim-`V` trap on the real objects | `dim W` = **2059** on all three — an invariant that separates **none** of the three pairs clause [4] proves inequivalent; `dim V` reads 2058 / 2058 / **2059** on `P`, `G`, `H″` and is worthless (Trap 1, §3.1). Cert 07 measured the same `dim W` on the pair |

## The separations

All three pairs, as `run.py` prints them; the counts are identical in both
arithmetics.

| pair | differing | union | `Σ \|T4\|·Δ` | largest `\|Δ\|` | source |
| --- | --- | --- | --- | --- | --- |
| `H″` vs `P` | **107** | 145 | −242 005 088 | 3 769 014 at `\|T4\| = 12` (`3.7·10⁻⁵` of that bin) | new |
| `H″` vs `G` | **146** | 147 | −54 951 386 832 | 619 868 210 at `\|T4\| = 4` (`5.9·10⁻³`) | new |
| `P` vs `G` | 146 | 147 | −54 709 381 744 | 616 808 004 at `\|T4\| = 4` (`5.8·10⁻³`) | cert 07 |

Every difference vector sums to zero, as it must; the first moment, which
nothing forces, does not — and all three are non-zero.

### The support, and where the invariant separates

`H″` and `P` populate **the same 145 bins**: the orientation switch moves
counts, not support. `G` populates **133** — it lacks fourteen of theirs
(`|T4| = 940, 972, 988, 1004, 1020, 1036, 1052, 1068, 1084, 1100, 1116, 1148,
1164, 1180`) and populates two of its own (`892`, `908`), which is why every
comparison with `G` reads *146 of 147* and `H″` vs `P` reads *107 of 145*.
`|T4| = 1108` is the **one** bin where all three agree, at 30 counts apiece.

`H″` vs `P` has a shape worth stating, because it is **not** the shape of the
other orders:

* the **bulk separates** — 106 of the 108 bins below `|T4| = 868` differ, the
  only exceptions being `|T4| = 644` and `852`;
* a **36-bin band agrees exactly**: every `|T4|` in `[868, 1180]`;
* and then the **extreme tail separates**, which it does not at 668, 716 or
  1676: the top bin `|T4| = 1236` differs, **12** for `H″` against **6** for
  `P`.

So the two profiles are one small-`|T4|` bulk apart plus a single six-count
difference in the highest bin. `run.py` asserts all three facts.

The largest discrepancy between `H″` and `P` is `3 769 014` at `|T4| = 12` —
`3.7·10⁻⁵` of that bin, and the whole separation lives at that scale.
Against `G` the discrepancies are two orders larger (`6.2·10⁸` at `|T4| = 4`,
0.6 % of the bin), which is the difference cert 07's *sampled* statistic could
see at `d = 7.4` and this one cannot: **`H″` vs `P` is invisible to any sample
of practical size**, and nothing cheaper than the exact profile would have
found it.

### The two new legs

`run.py` prints the first eight divergent bins of each; the full lists are the
banked JSONs.

**`H″` vs `P`** — 107 of 145 bins differ.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H″` | 105 821 994 475 | 102 449 438 069 | 96 022 641 731 | 87 168 327 023 | 76 659 086 757 | 65 352 165 541 |
| `P` | 105 825 054 681 | 102 453 207 083 | 96 023 541 023 | 87 167 542 521 | 76 658 478 649 | 65 350 483 547 |
| Δ | +3 060 206 | **+3 769 014** | +899 292 | −784 502 | −608 108 | −1 681 994 |

**`H″` vs `G`** — 146 of 147.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H″` | 105 821 994 475 | 102 449 438 069 | 96 022 641 731 | 87 168 327 023 | 76 659 086 757 | 65 352 165 541 |
| `G` | 106 441 862 685 | 102 965 772 295 | 96 360 629 035 | 87 289 459 295 | 76 582 817 395 | 65 122 330 205 |
| Δ | **+619 868 210** | +516 334 226 | +337 987 304 | +121 132 272 | −76 269 362 | −229 835 336 |

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports) — all
four rebuilt in clause [1] of every run:

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `P` = `2060-plain` | `510f89b7b423c85da0c7cada52cb0f62e0415d736d2040d6701f66c4e524cf6a` | yes (cert 07's pin) |
| `G` = `2060-gist`, the posted `H(2060)` | `c7a145d86210740dd3f8ea21ca896a54d6916007a042638f17c8c47f097200f7` | yes (cert 07's pin) |
| `H″` orientation switch | `4e1891b095b8aafa21176e494038f199b495c96a840bdb003e231c160870b801` | yes (formed from `P`; the desk's `manifest.json` digest) |
| the alternate-orientation array `S·H″·S` | `40e1d1c8cd40e94016c453f12e520a8518e7d29b773d3adaae3f484eca64398d` | yes (assembled from the seeds; no profile of it is banked) |

`P`'s and `G`'s digests were measured before this repository existed, and `G`'s
is additionally the digest of the public artifact itself.

**Banked files** (SHA-256 of the file bytes, compared in `run.py`). The two
`orient` banks are this certificate's own; the four others are cert 07's,
reused verbatim and re-pinned here:

| file | SHA-256 | banked for |
| --- | --- | --- |
| `data/sep2060-orient-exact-blas.json` | `38135aef205b4428760dc0439b29196b6776f215397fb34880156c66ee283f00` | cert 22 |
| `data/sep2060-orient-exact-bits.json` | `18dcf4c2e5d603324182eeef45c13e89ef80d5c0c0a5add3a7f90c333e4e5e87` | cert 22 |
| `data/sep2060-exact-blas-plain.json` | `5428aeac7b570fff55975c2b737fae9e8d0b717ec511735b68893e609a0037d8` | cert 07 |
| `data/sep2060-exact-bits-plain.json` | `e6c3af94712d0ba5cf3a3047796ccd474970036fec211b41a5579b7ff892ca49` | cert 07 |
| `data/sep2060-exact-blas-gist.json` | `a20b9a63cd3d93046c251b5c19aabeeac412b8f7933bbafa82d0210320e3aef0` | cert 07 |
| `data/sep2060-exact-bits-gist.json` | `9d8cc4b55c297c7e948df3e7639613a0580fc3e54af9eb12399bc010337f8a93` | cert 07 |

`data/sep2060-records.json` is **not** file-pinned here, for cert 20 and 21's
reason: it is shared with cert 07, which pins it
(`1c9742fe…a13202bc`), and the binding pin on it is the canonical digest of
each matrix it produces, checked in clause [1] — reinforced by the GS-condition
re-verification, which binds the seeds to the arrays by computation rather than
by a file hash.

Note that the four cert-07 banks are pinned here at exactly the digests cert
07's own `EXACT_FILE_PINS` carries, so the two certificates cannot drift apart
on the shared half of the evidence.

## Provenance of the new profile

Produced by the source laboratory's unchanged engine
(`Hadamard-2060`, `experiments/inequiv/exact_profile_big.py`) on the rented
`c2d-highcpu-16` `prof42-2` (`us-east1-b`, created 2026-09-02 18:14 UTC when
route B's cap freed the quota), 16 threads, under
`experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC 2026-09-02;
Amendment 1 ~11:05 UTC added this instance and re-priced the campaign). The
`blas` leg landed 20:10 UTC on 2026-09-02 and the `bits` leg 01:00 UTC on
2026-09-03; **`blas == bits` bin for bin**, so kill criterion 5 did not fire.
The seconds and peak resident sets below are the values those runs recorded in
the JSONs themselves: **6 631.0 s** at 860.7 MB (`blas`) and **17 654.7 s** at
131.9 MB (`bits`).

Banked into `data/` by the laboratory's `experiments/pr0042/bank.py` under
`--cert 22 --date 2026-09-03` — **2026-09-03 UTC**, the date the banks' own
`banked_note` carries. `bank.py` refuses unless
both implementations audit, agree bin for bin, and declare the digest
`manifest.json` pins for that matrix. It adds the seven header fields each
bank's `banked_note` names — `schema`, `description`, `matrix`,
`matrix_canonical_sha256`, `producer_filename`, `arithmetic`, `banked_note` —
and **nothing numeric**: every numeric field is the producer's own output,
unaltered. `run.py` checks `schema`, `matrix`, `matrix_canonical_sha256`,
`matrix_sha256`, `impl`, `folded`, `n`, `C_n_4`, `total`, `second_moment` and
`second_moment_want` against values it recomputes or rebuilds.

Credit, as everywhere in this repository, is to stations: the theorem and the
certificate are this repository's; the engine, the pre-registration, the matrix
building and verification, and the rented enumeration are the source
laboratory's; **the order-2060 artifact — `G`, byte for byte — is the public
record's**, and `PROVENANCE.md` carries its dated chain. Cert 07's disclosure
language binds here: no priority of any kind is claimed on it, on the decode,
or on existence at order 2060.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 7.5 s** (exit 0, 137 checks; measured 2026-09-03 on the desk, three runs, 7.5 / 7.5 / 7.5 s) |
| the GS condition on the four length-515 seeds, all 515 shifts | 0.3 s |
| assembling both `2060 × 2060` arrays from the seed record | 0.3 s |
| `verify/verify.py` on each of the four matrices (2 120 770 row pairs, exact integers) | ≈ 1.0 s each |
| the orientation switch, the cell count, and both forms of the `S·H″·S` identity at `N = 2060` | ≈ 1.3 s |
| `dim V` / `dim W` on all three 2060 matrices | ≈ 0.7 s |
| producing `sep2060-orient-exact-{blas,bits}.json` (16 rented threads) | 6 631.0 s (peak 860.7 MB) / 17 654.7 s (peak 131.9 MB) |
| producing cert 07's four banks (3 desk threads) | 15 691.0 s / 30 646.0 s (`P`), 13 800.3 s / 31 805.4 s (`G`) |
| `run.py --full` at 2060 | **not run in this repository** — see below |

**Why no `--full` leg was run here.** The flag is offered and wired exactly as
in certs 11, 13, 14, 15, 19, 20 and 21 — `certs/06-668-separation/full_recompute.py`
imported by a `sys.path` insert rather than copied, BLAS threads capped at
three before numpy loads, smoke-tested against the forced profile of Sylvester
`H(128)` first — and `--matrix` and `--impl` select the leg. The price is the
largest in the repository:

* one 2060 leg is about **137×** the 716 leg the same module took here
  (cert 14, 400.3 s), i.e. **≈ 15 h** for a single `blas` matrix at three
  threads on this desk. That 137 is the source laboratory's **measured**
  desk-side 716→2060 ratio (`experiments/pr0042/REGISTRATION.md`,
  Amendment 1), which implies the exponent 4.66;
* on the `Θ(n⁵)` law quoted elsewhere here, `(2060/716)⁵ = 197×`, i.e.
  **≈ 22 h**. Both are estimates and both say *many hours*; the smaller,
  measured one is quoted so that the decision not to run is not defended with
  an inflated price;
* the `blas` route materialises a `C(n,2) × n` pair matrix —
  `2 120 770 × 2 060`, **4.4 GB as `int8`** and **17.5 GB** as the `float32`
  copy `_profile_blas` makes — far past this desk (cert 07's `NOTES.md` gives
  the same 4.4 GB figure). The `bits` route packs the same pairs into 33
  `uint64` words per row, ≈ 560 MB, and is the tractable one at this order —
  and it is also the slower of the two.

The banked legs' own `seconds` fields price the campaign's engine for
comparison: **6 631.0 s `blas` / 17 654.7 s `bits` on 16 rented threads** for
`H″`, against cert 07's **15 691 s / 30 646 s on three desk threads** for `P` —
so even the memory-aware canonical-split engine is 4–9 desk-hours per leg here.
Certs 06, 08, 11, 13, 14, 15 and 19 each have an in-repo `--full` `blas` leg on
the record; certs 20, 21 and 22 do not, and the word *replayed* is not used of
them.

## What is NOT claimed

* **Nothing under the transpose-extended relation.** `H_2060-orient-T` and
  `H_2060-plain-T` are pending legs of the same campaign; neither is banked
  here. Order 2060 remains the only row-side separation statement in
  `note/NOTE-B.md`, exactly as cert 20 left 1676 before cert 21.
* **Nothing about any matrix versus its own transpose at 2060.** Cert 19
  decided that question at 668 and at 668 only.
* **The default run recomputes nothing.** Only `--full` would bind a bank to a
  matrix by computation, and **no `--full` leg has been run in this repository
  at 2060** — nor, for the four cert-07 banks, has one ever been.
* **No general theorem.** That the orientation switch lands in a class of its
  own is now known at four orders (668, 716, 1676, 2060). Four orders are four
  orders; nothing is claimed about orientation in general, and nothing about
  `ψ(ρ) = −1` at 2060 — the `×104` twist relating `P` and `G` acts on the
  **column index** of each circulant
  (`twist(x,k)[i][j] = circ(x)[i][(k·j) mod v]`, the record's own
  convention), not on the seed values as Lemma T's `ψ`-twist does, so no
  analogue of the 668/716/1676 twist theorem is asserted here.
* **No claim that three is the count** at 2060. It is a lower bound exhibited
  by the three matrices in hand.
* **No novelty or priority claim of any kind.** Order 2060 was settled by the
  **publicly posted matrix**, which is `G` itself; this certificate counts
  classes among the artifacts banked here and says nothing about existence at
  2060 or about who first exhibited a Hadamard matrix of this order. Cert 07's
  disclosure language binds; `PROVENANCE.md` carries the dated chain.
* **Matching invariants prove nothing.** `dim W` agrees on all three objects
  the 4-profile proves pairwise inequivalent, and cert 07 records every cheap
  exact invariant agreeing on the pair. Read every "agrees" in this repository
  as "did not separate", never as "the same".
* **The sampled statistic of cert 07 is not evidence here.** It was never run
  on `H″`, and it would not have resolved `H″` against `P` if it had been: the
  largest bin discrepancy is `3.7·10⁻⁵` of its bin.

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/22-2060-three-classes/run.py
```

The default path is standard library only, no network, and is the whole
certificate. The `--full` paths below import numpy (finder-side only, never in
the trust chain) and cap BLAS threads at three; at this order each is many
hours, and the `blas` variants want memory this desk does not have. **None of
them has been run here.**

```
python certs/22-2060-three-classes/run.py --full --impl bits --matrix orient
python certs/22-2060-three-classes/run.py --full
```

Exit code 0 iff every check passed.
