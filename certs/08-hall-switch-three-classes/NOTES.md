# cert 08 — order 668 carries at least three Hadamard equivalence classes

**Label: PROVEN.** Default run:
`python certs/08-hall-switch-three-classes/run.py` from the repository root.
Standard library only, about 1.4 s, 111 checks, exit 0. That run **audits a
banked exact computation**; the word *replay* belongs to `--full`, which
recomputes all five profiles here. The trust boundary is set out in full
below.

---

## The theorem

> **Theorem.** At order 668, let
>
> * `H` be the decoded `(s, i) = (1, 1)` bordered Goethals–Seidel record
>   (`data/payload-records.json`),
> * `H'` be its Lemma-T `i = 2` rebuild (`data/sep668-twisted-record.json`,
>   cert 06),
> * `H★` be the paired Hall switch of `H` — the 1 328 entries named by
>   `data/sep668-hall-switch.json`, negated.
>
> Then `H`, `H'` and `H★` are **pairwise Hadamard-inequivalent**: for no pair
> is there a `D_r P_r · P_c D_c` carrying one to the other. Moreover, adding
> the **transpose** to the group merges no pair either.

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all
`C(668,4) = 8 222 179 035` row 4-subsets, with
`T4 = Σ_c H[i][c]H[j][c]H[k][c]H[l][c]`, is a Hadamard-equivalence invariant
(`note/NOTE-B.md` §3.1, invariant **I5**). The three profiles populate the
**same 80 bins** and differ in **26** (`H` vs `H'`), **27** (`H` vs `H★`) and
**27** (`H'` vs `H★`) of the bin counts. For the transpose-extended relation,
`A ≈ B` requires `A ~ B` **or** `A ~ Bᵀ`, so each pair needs a second
refutation; the exact profiles of `(H')ᵀ` and `(H★)ᵀ` supply them — **50**,
**49** and **50** differing bins. An invariant that differs is a separation. ∎

## The trust boundary — what a default run does and does not establish

The five `C(668,4)` enumerations **were not run inside this repository.** They
ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile.py` (the upper-triangle route) and
`experiments/inequiv/exact_profile_big.py` (the canonical-split route), numpy,
three BLAS threads, `$0`, no VM — and their output was banked into nine JSON
files holding ten profiles: this certificate's five
(`sep668-hall-exact-{blas,bits}`, `sep668-hall-T-exact-{blas,bits}`, and
`sep668-twisted-T-exact`, which carries both implementations) and cert 06's
four. Which producer ran for which bank is recorded in each bank's own
`engine`, `enumeration` and `arithmetic` fields and repeated in the Runtimes
table below.

**A default `run.py` establishes:**

* the eleven bank files are byte-for-byte the ones pinned in `run.py`;
* all five matrices — `H`, `H'`, `H★`, `(H★)ᵀ`, `(H')ᵀ` — rebuild from banked
  parameters, pass `verify/verify.py`, and carry the pinned canonical
  digests, each computed **in-process**;
* the preprint's `X₁…X₄`, its eq. (6) blocks and both of its published
  digests reproduce from `data/payload-records.json` alone;
* every bank that declares a matrix digest declares the digest of the matrix
  **rebuilt in this run** — the comparison is against the in-process digest,
  not against a static string;
* every profile satisfies the forced identities (bins `≡ 4 (mod 8)`, total
  `C(668,4)`, second moment `n³(n−1)(n−2)/24` to the unit), and two
  independent implementations agree bin for bin on all five matrices.

**A default run does not establish that a banked histogram was computed from
the matrix `run.py` rebuilt.** A `matrix_canonical_sha256` field is a
declaration, not a proof of computation. All ten profiles now carry one —
the four cert-06 banks had theirs backfilled 2026-09-01 from digests this
repository already pinned, recorded in each bank's `banked_note` — so every
bank is bound to a matrix rebuilt in this run; a bank declaring none would
still be named in a `[NOTE]` line rather than passed over in silence.

**`--full` is what closes the gap.** `python
certs/08-hall-switch-three-classes/run.py --full` recomputes the exact
4-profile of **all five** matrices here, from the rows clause [1] verified,
and compares each fresh profile against **both** banked implementations of
that matrix bin for bin. Every leg of the theorem — the transpose-extended
legs included — then rests on a computation made inside this repository. The
machinery is cert 06's `full_recompute.py`, imported by `sys.path` insert
rather than copied, so the two certificates cannot drift apart. Default
arithmetic is the float32-BLAS route, which is **exact** here: `T4` is a sum
of 668 signed units, so every partial sum is an integer of absolute value
`≤ 668`, far inside float32's exactly representable integer range `2²⁴` — no
rounding can occur at these sizes. numpy is imported only under this flag, is
finder-side only, and is never in the trust chain; BLAS threads are capped at
three.

**It has been run.** On 2026-08-31, `--full` completed on this desktop in
**1 932.5 s** (32 min) — 123 checks, no failures, exit 0 — and every one of
the five fresh profiles matched **both** of its banked implementations bin for
bin: `decoded` 80 bins/355 s, `twisted` 80/355, `H★` 80/382, `(H★)ᵀ` 79/480,
`(H')ᵀ` 79/359. That was this repository's first in-repo regeneration of all
five banks. `--impl bits` costs roughly five times as much and `--impl both`
six; neither has been run at 668 inside this repository.

One thing `--full` does **not** reproduce: the banks for `H★`, `(H★)ᵀ` and
`(H')ᵀ` were produced by two *different enumerations* (upper-triangle
with `/3`, and the canonical-split bijection). `full_recompute.py` runs the
upper-triangle enumeration on both of its arithmetic paths, so a `--full` run
re-establishes the numbers from the rebuilt matrices; the enumeration
diversity remains a property of the banks, checked on the default path by
`blas == bits`.

Read the default verdict as *banked exact computation audited*; read `--full`
as the replay.

## Priority — stated first, and conceded cleanly

An **anonymous preprint** hosted at `hadamard-668.vercel.app`
(*Two H-Inequivalent Hadamard Matrices of Order 668*; no author on the page,
empty `/Author` in the PDF; "August 2026" on the page; PDF `CreationDate`
`D:20260813041406+05'30'` = 2026-08-12 22:44 UTC; retrieved and verified
firsthand 2026-08-31; PDF SHA-256 `ca4850fc…`) is, **as far as this
laboratory's search located, the first publication of the statement that
order 668 carries at least two Hadamard equivalence classes.** That statement
is theirs. No date archaeology changes it and none is attempted; PDF
timestamps establish compilation, not first public availability, and are
recorded here as evidence, not as fact.

Three things that verification established, and that `run.py` re-derives:

* **Their `H` is this repository's decoded record, byte for byte.** Their four
  subsets `X₁…X₄ ⊂ Z₁₆₆` are exactly the negative supports of the four decoded
  seeds, element for element; their `K`, `T`, `S` are the record's `corner`,
  `row_table` and `col_table`ᵀ. So an outside party independently recovered
  the same `s = 1` border from Alpöge's public artifacts. That **corroborates
  this repository's decode**, and it means neither party originated the
  construction data independently of the public posting; it is the announcing
  team's content, and neither party claims priority over it.
* **Both SHA-256 digests they publish reproduce**, from an independent build,
  once the ambiguous sentence "row-major bytes, encoding +1 by 1 and −1 by 0"
  is pinned to *one byte per entry, `0x01`/`0x00`, 446 224 bytes*. Five
  readings were tried upstream; exactly one matches, on **both** matrices.
  `run.py` reproduces both.
* **Their second matrix is not this repository's second matrix.** Theirs is a
  Hall switch of `H`; this repository's `H'` is a Lemma-T `i = 2` frame
  rebuild — a different construction, and now proven a **third** class.

**What is this repository's:** the third class, the construction that produces
it, and the three pairwise exact-profile separations, standard and
transpose-extended. What is the preprint's: the "≥ 2 classes" statement, and
the matrix `H★`.

**Their statistic is blind where the full profile is not, and measurably
so.** Their
separating statistic `Φ_M` is a canonical *slice* of the `|T4|` 4-profile (the
4-subsets with exactly one row in the distinguished border quadruple). It is a
genuine invariant — its invariance uses the preprint's own **Lemma 3**
(uniqueness of the distinguished quadruple), stated and exactly computed
there; this certificate independently reproduces that computation from the
exact profile's `660` bin (control C3 below) — but on this order it is
**bin-for-bin identical on `H` and `H'`**, a pair this certificate proves
inequivalent. A
published cheap invariant returning "no separation" on a true separation is
the working reason this repository pays for the full `Θ(n⁴)` computation. That
measurement is upstream (`intel/hadamard-668-vercel-intake.md` §4.3), is
**MEASURED** rather than replayed here, and is recorded in `note/NOTE-B.md`
§3.4 beside the `dim V` and signed-`T4` traps.

## The evidence chain

`run.py` runs all of the following before it prints a verdict.

**[0] The bank is pinned.** Eleven banked files, each SHA-256 compared in code.
Five are shared read-only with cert 06 and carry cert 06's values unchanged.

**[1] Five matrices are rebuilt, not assumed.** Both records go through
`tools/bordered_gs.py`'s `check_record`, which re-checks **every hypothesis of
the master theorem** — H0 shape, H1 the two-tier PAF profile, H2 the
corner/row-table budget, H3 the column-table Gram, H4 the coupling (that
module's numbering; see its `LABEL MAPPING` block, which is a permutation of
`note/NOTE-B.md` Theorem A's), the
derived D1/D3/D5 and the Σ̄ law, and the compression-lemma cross-check — and
only then assembles. `H★` is then built by negating the mask, and `(H★)ᵀ` and
`(H')ᵀ` by transposition. **All five** go to `verify/verify.py`, the trust
chain, and each canonical digest is compared against the pin. The matrices are
deleted afterwards; nothing multi-megabyte is committed.

In the same clause, from `data/payload-records.json` alone: `X₁…X₄` are
re-derived as the negative supports of the seeds and compared against a pinned
digest of the preprint's own transcription; `K K^T = T T^T = S S^T = 4I` and
`K S^T = −2T` (the preprint's eq. 6) are re-derived from the record's border;
and both published digests are reproduced.

**[1b] The switch is checked, not trusted.** The mask is rebuilt from the
banked rectangles *and* independently from the closed forms in the record's
own `(s, |G|)`; the two must agree. Then:

| | control | result |
| --- | --- | --- |
| **C1** | the set of positions where the rebuilt `H` and `H★` differ is compared with the mask | equal, 1 328 positions |
| **C2** | the switch is applied twice | the canonical digest of `H` comes back — it is an involution |
| **C3** | the border quadruple `P = {1,2,3,4}` is checked type-1 (minority 4, `\|T4\| = 660`) in `H`, `H'` and `H★` | holds in all three — the hypothesis the preprint's own Lemma 3 needs |

**[2] Eight banked profiles are audited — not recomputed — in exact
integers**, plus the two in the twisted-transpose bank. Per profile: every
populated bin is `≡ 4 (mod 8)`; the counts total `C(668,4)`; and the second
moment equals `n³(n−1)(n−2)/24 = 5 517 193 410 096` — the closed form proved
in `note/NOTE-B.md` §3.1, a 13-digit number hit to the unit and tuned for by
nothing. Where a banked `second_moment` field is present it is checked too.
Where a banked `matrix_canonical_sha256` or `matrix_sha256` field is present
it is compared against the digest of the matrix **rebuilt in this run**,
computed in-process in clause [1], rather than against the static pin at the
top of the script. Comparing against the in-process digest rather than the
static pin leaves no place for a stale constant to hide. All ten profiles
declare such a digest — the four cert-06 profiles had theirs backfilled
2026-09-01, and the check label says so where that is the case — and
`run.py` prints a `[NOTE]` for any bank that declares none. A field that is
present but empty or malformed is a hard failure, not a missing binding.
A declared digest says the bank *names* this matrix. It does not say
the histogram was computed from it; only `--full` does that.

**[3] Two independent implementations agree bin for bin**, on **all five**
matrices — `H`, `H'`, `H★`, `(H★)ᵀ` and `(H')ᵀ`. So every leg of the theorem,
the transpose-extended legs included, rests on two implementations; no
comparison in this certificate rides on a single arithmetic route. `blas` is a
float32 Gram of the pair-vector matrix (exact at these sizes: every entry and
partial sum is an integer below `2²⁴`); `bits` packs rows into `uint64` words
and uses `|T4| = |n − 2·popcount(u_P ⊕ u_Q)|`. On `H★` and `(H★)ᵀ` the two
also use **different enumerations** — the upper-triangle-with-`/3` route
(`exact_profile.py`) against the canonical-split bijection
(`exact_profile_big.py`) — so there nothing but the answer is shared; on
`(H')ᵀ` both run the canonical split and differ in arithmetic alone.

**[4] Controls.**

| | control | result |
| --- | --- | --- |
| **C4** | *negative* control: a banked profile is corrupted in a **total-preserving** way (one count moved from `\|T4\| = 4` to `\|T4\| = 660`), so only the second-moment identity can catch it, and `audit()` is required to raise | the assert fires |
| **C5** | the comparator is exercised in the null direction: every banked profile against itself | 0 differing bins, five times |
| **C6** | the `\|T4\|` profile is invariant: a deterministically seeded (`20260831`) signed row **and** column permutation of Sylvester `H(16)` and Paley I `H(20)`, checked Hadamard and checked different from the original | the profile does not move |
| **C7** | the transposed-profile route, on matrices small enough for straight `O(C(n,4))` enumeration: Sylvester `H(8)` and `H(16)` are symmetric, so `profile(Mᵀ) = profile(M)` is **forced** — and holds; Paley I `H(20)` is not symmetric and is **measured**, not asserted (its profiles happen to agree) | forced cases hold; the non-forced case is reported, never asserted |

C7 exists because clause 4 rests on transposed profiles. At 668 the route is
not vacuous: `(H★)ᵀ` populates **79** bins where `H★` populates 80 — bin 644
is present in `H★` and empty in `(H★)ᵀ` — and `(H')ᵀ` likewise populates 79.

**[6] `--full` — the replay, off by default.** All five profiles recomputed
here from the rebuilt rows and compared against both banked implementations
bin for bin — run 2026-08-31, all ten comparisons exact; see *The trust
boundary* above. `run.py` first re-derives the canonical digest of the rows it
is about to enumerate, so the recomputation is demonstrably on the matrix
`verify.py` accepted, and smoke-tests the ported numpy path against the
**forced** profile of Sylvester `H(128)` before spending the half hour.

## The separation, under the standard relation

All three share the same 80 bins, `≡ 4 (mod 8)` from 4 to 660 (`620`, `636`,
`652` empty in all three). The *support* does not separate them, and neither
does the extreme tail — `604: 19, 612: 18, 628: 2, 644: 1, 660: 1` in all
three. The **bulk** separates them.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 |
| --- | --- | --- | --- | --- | --- |
| `H` decoded | 2 073 064 058 | 1 852 054 148 | 1 491 070 735 | 1 091 442 371 | 732 009 734 |
| `H'` rebuild | 2 073 109 602 | 1 852 009 274 | 1 491 079 303 | 1 091 478 493 | 731 928 178 |
| `H★` Hall | 2 073 010 592 | 1 852 087 772 | 1 491 075 635 | 1 091 423 423 | 732 034 330 |

| pair | differing bins | same support | `Σ \|T4\|·Δ` |
| --- | --- | --- | --- |
| `H` vs `H'` | **26** of 80 | yes | −306 848 |
| `H` vs `H★` | **27** of 80 | yes | +712 400 |
| `H'` vs `H★` | **27** of 80 | yes | +1 019 248 |

`run.py` prints all 27 of the `H`-vs-`H★` bins. Every difference vector sums
to zero, as it must; the **first** moment, which nothing forces, does not. The
largest `H`-vs-`H★` discrepancy is 53 466 at `|T4| = 4`, i.e. `2.6·10⁻⁵` of
its bin — invisible to any sample of practical size. The `2·10⁷`-draw sampled
comparison of the `H`/`H'` pair banked in
`data/sep668-sampled-histograms.json` reads `max |z| = 1.8` and sees nothing
at all. **A null sampled comparison is worth very little.**

## The separation, under the transpose-extended relation

`A ≈ B` requires `A ~ B` **or** `A ~ Bᵀ`, so each pair needs two refutations.
Both are exact-profile comparisons; the second column is the one the
`(H')ᵀ` run was bought for.

| pair | vs `B` | vs `Bᵀ` | verdict |
| --- | --- | --- | --- |
| `H` vs `H'` | 26 | **50** | **SEPARATED** |
| `H` vs `H★` | 27 | **49** | **SEPARATED** |
| `H'` vs `H★` | 27 | **50** | **SEPARATED** |

So the three-class theorem is **transpose-robust**: all three pairs stay
apart. `run.py` also checks the redundant crossings that must hold if the
above do — `(H')ᵀ` vs `H★` (49) and `H'` vs `(H★)ᵀ` (50).

## Pinned digests

**Matrices** (canonical SHA-256 of the `+/-` serialisation — the digest
`verify/verify.py` reports; all five re-derived by `run.py` from banked
parameters alone):

| matrix | canonical SHA-256 |
| --- | --- |
| `H`, decoded `(1,1)` record | `bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0` |
| `H'`, Lemma-T `i = 2` rebuild | `600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3` |
| `H★`, paired Hall switch | `7f6af1d9e9c80fae52c6f178e298144d209c01f6c103db119e6c12425aa1a722` |
| `(H★)ᵀ` | `565a9ca5a9db739f74215474364202a4d35fd691542b18e8ace8bcbb3c190c65` |
| `(H')ᵀ` | `32afde351e1f44aa4236cb3c406fbbcd59c5ab2cc3e02c8380e08680db2d19d7` |

**The preprint's own published digests** (one byte per entry, `0x01`/`0x00`,
row-major, 446 224 bytes), reproduced by `run.py`:

| matrix | published SHA-256 |
| --- | --- |
| their `H` | `19f435b31eb6561dd97356b761e9b0f174824e42c215c33fbafc20f9b1b20744` |
| their `H★` | `08428eb6779986c50ee2997584daa5c87e00520cc1e6f62d7793dc519f85ea45` |

**Banked files** (SHA-256 of the file bytes, compared in `run.py`):

| file | SHA-256 |
| --- | --- |
| `data/sep668-hall-switch.json` | `13efd2402b8394c62c901af4f7cfbec7b2e474832dd3055c6b9e9e220b351c85` |
| `data/sep668-hall-exact-blas.json` | `35e716ecb43bb6190d5dd6f4160e0bc2bed4f61a3aacf07a36ff9d190810c154` |
| `data/sep668-hall-exact-bits.json` | `a6f703b499d98995f6446a1aed671284c47e99cfe869f3ce8dc8b5fd9394accb` |
| `data/sep668-hall-T-exact-blas.json` | `151fb5d6e70cf56d6a1c2aa124a597a837bca0ecf5d64958b43a34c05383e0db` |
| `data/sep668-hall-T-exact-bits.json` | `48fdb26f8b1ee5135ed278ec866e204c1ab47df168c043fabff8699c0f4fd8bb` |
| `data/sep668-twisted-T-exact.json` | `38355274ec61d33fcd96e24255e4a7b02874150cd914fdfb928d28cee751fc4a` |
| `data/sep668-twisted-record.json` | `fe8154179ba2ebfe097c82e468368cdc8a070548555bb10140949af0560611fb` |
| `data/sep668-exact-blas-decoded.json` | `370fffe6c2f5dc53c09d3b74f8c09dd2bc2a39a1ac2b27fb5167ab4d3559387b` |
| `data/sep668-exact-bits-decoded.json` | `7bace61441f17b5e95fff433bdc5939da212e2b8735e8738d7ed3078fae456b7` |
| `data/sep668-exact-blas-twisted.json` | `8526b3cfa7938a9af334e23f722b1c215ffd1e318c0c713ecc3da1b91f5b3afe` |
| `data/sep668-exact-bits-twisted.json` | `f40bbb8c3906d6fc7374e3e04c2b68eaf29393e50b5662f69eee2426ed3f1e9a` |

`data/payload-records.json` is not file-pinned here, on purpose: it is shared
with certs 01 and 06, and the binding pin on it is the canonical digest of the
matrix it produces, which is checked above.

The four `*-exact-*.json` files banked by this certificate carry a
`banked_note` field naming the header fields added at banking time —
`schema`, `description`, `matrix`, `matrix_canonical_sha256`,
`producer_filename`, `engine`, `enumeration`, `arithmetic`. **Every numeric
field is the producer's own output, unaltered.** The added
`matrix_canonical_sha256` binds each profile to the matrix `run.py` rebuilds,
so a profile can no longer drift away from the object it describes;
`producer_filename` keeps the producing run's own filename on the record.
`data/sep668-twisted-T-exact.json` differs in shape: it holds **both**
implementations in one file, under `implementations.bits` and
`implementations.blas`, each with the producer's `matrix_sha256`,
`peak_rss_mb`, `threads` and block geometry.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 1.4 s** (111 checks, exit 0) |
| **`run.py --full`, all five recomputed (`blas`)** | **1 932.5 s** (123 checks, exit 0; measured 2026-08-31) |
| `--full`, per matrix | 355 s, 355 s, 382 s, 480 s, 359 s |
| rebuild + hypothesis re-check + `verify.py`, per 668 matrix | ≈ 0.15 s |
| producing `data/sep668-hall-exact-blas.json` | 365.0 s, upper-triangle route |
| producing `data/sep668-hall-exact-bits.json` | 168.8 s, canonical split, 70.8 MB peak |
| producing `data/sep668-hall-T-exact-blas.json` | 327.5 s, upper-triangle route |
| producing `data/sep668-hall-T-exact-bits.json` | 204.4 s, canonical split, 70.7 MB peak |
| producing `data/sep668-twisted-T-exact.json`, `bits` | 206.8 s, 70.9 MB peak |
| producing `data/sep668-twisted-T-exact.json`, `blas` | 106.1 s, 409.6 MB peak |

All profile runs: single desktop process, `OMP_NUM_THREADS = 3`, no VM, `$0`.
The profile producers are numpy and therefore **finder-side only** — they are
never in the trust chain, which is `verify/verify.py` and stdlib integers. The
seconds on the `producing …` rows are the values those upstream runs recorded
in the JSONs themselves; the `--full` rows are wall-clock from the in-repo run
of 2026-08-31, on a desktop that was also carrying another 3-thread job. That
run's check count is what it printed on the day; the 2026-09-01 digest
backfill added four checks to every path, so a fresh run reports four more.

## What is NOT claimed

* **The default run does not claim to have recomputed anything.** The five
  `C(668,4)` enumerations ran upstream, in `Hadamard-2060`'s
  `experiments/inequiv/exact_profile.py` and `exact_profile_big.py`; the
  default path audits their banked output and checks the declared matrix
  digests against the matrices it rebuilt. Only `--full` recomputes.
* **Nothing about orders 716, 1676 or 1772.** The same constructions exist
  there and the corresponding exact computation costs roughly `1.4×`, `98×`
  and `130×` the 668 run. *This* certificate says nothing about them. The
  `716` computation has since been made and carries its own certificates —
  certs 11, 14 and 15, NOTE-B.md §3.6 (three classes, and the statement holds
  under the transpose-extended relation); `1676` and `1772` have not been
  made, so nothing is said there.
* **No claim of priority on "≥ 2 classes at 668."** That is the preprint's;
  see above.
* **No claim of novelty of existence at 668.** Order 668 was settled by the
  publicly posted matrix.
* **No claim that these are the only three classes** at order 668.
* **The preprint's 16-mask family is not verified here.** Its claim that all
  16 switch masks are Hadamard and that a pinned colored-graph
  canonicalization splits them into exactly two H-classes by mask parity is
  **REPORTED-FROM-SOURCE**: this lab has no `nauty` and did not check it.
* **`Φ_M`'s blindness on the `H`/`H'` pair is MEASURED upstream, not replayed
  here.** `run.py` computes no `Φ`.
* **Matching invariants prove nothing.** Every cheap invariant tested on these
  matrices — `dim W`, the collision profiles, the exact extreme strata,
  `rank₁₆₇`, the dual weight enumerator, `2·10⁷` sampled 4-subsets, and a
  *published* invariant, `Φ_M` — returned identical values on at least one
  pair that is in fact inequivalent. Read every "agrees" in this repository as
  "did not separate", never as "the same".
