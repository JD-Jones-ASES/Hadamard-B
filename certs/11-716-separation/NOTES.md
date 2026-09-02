# cert 11 — order 716 carries at least two Hadamard equivalence classes

**Label: PROVEN.** Default run: `python certs/11-716-separation/run.py` from
the repository root. Standard library only, about one second, exit 0. That
run **audits a banked exact computation**; the word *replay* belongs to
`--full`, which recomputes both profiles here. The trust boundary is set out
in full below.

---

## The theorem

> **Theorem.** Let `H₁` be the decoded `(s, i) = (1, 1)` bordered
> Goethals–Seidel record at order 716 (`data/payload-records.json`) and let
> `H₂` be its Lemma-T `i = 2` rebuild (`data/twisted-i2-records.json`,
> order 716 — cert 02). Then `H₁` and `H₂` are **not Hadamard-equivalent**:
> there is no `H₂ = D_r P_r H₁ P_c D_c` with `P` permutation matrices and
> `D` diagonal `±1`.

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all
`C(716,4) = 10 859 143 295` row 4-subsets, with
`T4 = Σ_c H[i][c]H[j][c]H[k][c]H[l][c]`, is a Hadamard-equivalence
invariant. The two profiles populate the **same 87 bins**, and **27 of the
87 bin counts differ**. An invariant that differs is a separation. ∎

**Consequence.** Order 716 carries at least two Hadamard equivalence
classes, and the Lemma-T construction at `ψ(ρ) = −1` provably leaves the
equivalence class at a **second** order — 668 (cert 06) was the first, and
`note/NOTE-B.md` §1.4 previously left 716, 1676 and 1772 unclaimed.

The computation is finite and exact; the two profiles, in two
implementations each, are banked in `data/`.

## The trust boundary — what a default run does and does not establish

The `C(716,4)` enumeration **was not run inside this repository.** It ran in
the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py`, numpy, three threads, on
**2026-09-01**, under the pre-registration
`experiments/inequiv/REGISTRATION-716-exact.md`, which was **flushed before
the matrices were built** and fixed the decision rule, the two audit
identities and five kill criteria in advance. Its output was banked into
`data/sep716-exact-{blas,bits}-{decoded,twisted}.json`; the seconds each run
recorded in its own JSON are in the Runtimes table below.

**A default `run.py` establishes:**

* the four bank files are byte-for-byte the ones pinned in `run.py` (SHA-256
  compared in code);
* both matrices rebuild from the banked records through the full
  master-theorem hypothesis re-check, pass `verify/verify.py`, and carry the
  pinned canonical digests — that digest is computed **in-process** and is
  what the bank-identity check is compared against;
* the twisted record — cert 02's shared bank, not a file of this
  certificate's own — has seeds that are re-derived here as the `ψ`-twist of
  the decoded seeds and compared character for character, so the two records
  are bound to each other by computation rather than by a file pin;
* each of the four banks declares `matrix_canonical_sha256`, and each
  declaration is compared against that in-process digest — so no bank can
  drift onto a different object;
* every banked profile satisfies the forced identities: bins `≡ 4 (mod 8)`
  (`716 ≡ 4 mod 8`), total `C(716,4)`, second moment
  `n³(n−1)(n−2)/24 = 7 807 861 101 040` to the unit;
* the two independent banked implementations agree bin for bin on each
  matrix;
* the separation itself — 27 of 87 bins differ — follows from those banked
  numbers.

**A default run does not establish that a banked histogram was computed from
the matrix `run.py` rebuilt.** All four banks are *producer-banked*: the
`matrix_canonical_sha256` each carries is the digest the producing engine
itself recorded against the matrix file it enumerated, and it equals the
canonical digest this repository already pins (certs 01, 02). No numeric
field was altered at banking, and each bank's own `banked_note` says exactly
which fields were added. It remains a *declaration*: a self-declared digest
is metadata, not a proof of computation.

**`--full` is what closes the gap.** It recomputes both profiles here, from
the rows clause [1] verified — the identity of those rows is re-derived
immediately before the enumeration — and compares each fresh profile against
**both** banked implementations bin for bin. At this order it is also a
genuinely **third arithmetic route**: the banks came from the
canonical-split engine, while `full_recompute.py` enumerates the upper
triangle of `U Uᵀ` and divides out the threefold overcount. It has been run:
on 2026-09-01, `--full --impl blas` matched both banks bin for bin (see
*Runtimes*). Read the default verdict as *banked exact computation audited*;
read `--full` as the replay.

## Why the |T4| profile is an invariant

`note/NOTE-B.md` §3.1, invariant **I5**. In one line: a row negation
contributes exactly one sign to `T4` so `|T4|` is fixed; a column negation
contributes `d_c⁴ = 1`; row and column permutations relabel the 4-subsets
and the coordinates. Hence the *multiset* of `|T4|` values is carried to the
corresponding multiset by any element of the equivalence group.

Two warnings from §3.1 that this certificate honours mechanically:

* **the signed `T4` histogram is not an invariant** — only `|T4|` is. No
  signed histogram is banked or compared anywhere in this repository; every
  banked 716 profile declares its folding in its own `folded` field.
* **`dim V` is not an invariant** — only `dim W = V + ⟨1⟩` is. On this pair
  `dim V` **differs** (714 vs 715) and the difference is worthless;
  `run.py` measures it (control C4) precisely so it is on the record as a
  trap and never as evidence. `dim W = 715` on both.

## The evidence chain

`run.py` runs all of the following before it prints a verdict.

**[0] The bank is pinned.** Every banked profile file's SHA-256 is compared
in code.

**[1] Both matrices are rebuilt, not assumed.** Each record goes through
`tools/bordered_gs.py`'s `check_record`, which re-checks **every hypothesis
of the master theorem** — H0 shape, H1 the two-tier PAF profile, H2 the
corner/row-table budget, H3 the column-table Gram, H4 the coupling (that
module's numbering; see its `LABEL MAPPING` block), the derived D1/D3/D5 and
the Σ̄ law, and the compression-lemma cross-check of the assembled core
against the `G/K` Goethals–Seidel array — and only then assembles. Each
assembled matrix is handed to `verify/verify.py`, the trust chain, and its
canonical digest compared against the pin. Before that, the twisted record's
four seeds are re-derived as `x'_q = ψ·x_q` with `ψ(g) = (−1)^g` on `ℤ₁₇₈`
and compared character for character against the bank, and its group and
reflection shift are compared against the decoded record's. The matrices are
deleted afterwards; nothing multi-megabyte is committed.

**[2] Each banked profile is audited — not recomputed — in exact
integers.** Per profile: every populated bin is `≡ 4 (mod 8)`; the counts
total `C(716,4)`; and the second moment equals
`n³(n−1)(n−2)/24 = 7 807 861 101 040` — the closed form proved in
`note/NOTE-B.md` §3.1, a 13-digit number hit to the unit by both matrices
and neither tuned for. All four 716 banks declare `second_moment`
themselves, so `run.py` recomputes it from the profile *and* checks the
declared value. All four declare `matrix_canonical_sha256`, and `run.py`
compares each declaration against the digest of the matrix rebuilt in the
same run; a bank that declared none would be named out loud in a `[NOTE]`
line.

**[3] Two independent implementations agree bin for bin.** `blas` is a
float32 Gram of the pair-vector matrix (exact at these sizes: every entry
and every partial sum is an integer below `2²⁴`); `bits` packs rows into
`uint64` words and uses `|T4| = |n − 2·popcount(u_P ⊕ u_Q)|`. They share no
arithmetic. Both hit `C(716,4)` and the second moment, and agree on all 87
bins of both matrices — kill criterion 4 of the pre-registration, clear.

**[4] Controls.**

| | control | result |
| --- | --- | --- |
| **C0** | every control matrix is checked to *be* Hadamard before anything is asked of it | pass, 5 matrices |
| **C1** | the full `\|T4\|` profile of Sylvester `H(8)`, Sylvester `H(16)`, Paley `H(20)`, and two plain Goethals–Seidel arrays `H(28)` and `H(36)` built by this repository's own assembler at the degenerate `s = 0` layer (three- and four-bin profiles, so the agreement is not vacuous), each computed **twice** — by straight `O(C(n,4))` enumeration over the `±1` entries, and by the pair-vector / Gram-triangle route `--full` takes | agree bin for bin on all five; both routes hit `C(n,4)` and the second moment; the two Sylvester cases additionally match their **forced** profile `{0: C(n,4) − n(n−1)(n−2)/24, n: n(n−1)(n−2)/24}`, so C1 is a positive control with a predicted answer |
| **C2** | *negative* control: a banked 716 profile is corrupted in a **total-preserving** way (one count moved from `\|T4\| = 4` to `\|T4\| = 708`), so only the second-moment identity can catch it, and `audit()` is required to raise | the assert fires |
| **C3** | the `dim V` trap on Sylvester `H(16)` under a deterministically seeded (`20260901`) random signed row negation | `dim V` moves `4 → 5`; `dim W` stays `5`; the `\|T4\|` profile is unchanged |
| **C4** | the same trap on the real objects | `dim V` 714 (decoded) vs 715 (twisted) — **differs and is worthless**; `dim W = 715` on both |

The C1 route-2 code is the same bookkeeping `--full` depends on (drop the
diagonal, drop the `n·C(n−1,2)` index-sharing pairs, divide by three), run
where straight enumeration can check it.

**[5] `--full` — the replay.** `python certs/11-716-separation/run.py --full`
recomputes both profiles **here**, from the rows clause [1] rebuilt and
`verify.py` accepted, with numpy, and compares each fresh profile against
**both** banked implementations bin for bin. The recompute module is
`certs/06-668-separation/full_recompute.py`, imported by a `sys.path` insert
rather than copied — exactly as cert 08 imports it — so the three
certificates cannot drift apart. `--impl blas` restricts it to the float32
route (about a quarter of an hour rather than an hour and a half); that
route is exact — each dot product is a sum of 716 signed units, so every
partial sum is an integer of absolute value `≤ 716`, far inside float32's
exactly representable integer range `2²⁴`. numpy is imported **only** under
this flag, is finder-side only, and is never in the trust chain; BLAS
threads are capped at three, set in `run.py` before numpy loads. Before
spending the time, `--full` smoke-tests each ported path against the forced
profile of Sylvester `H(128)`, which needs two `uint64` words per row.

## The separation

Same 87 bins, `≡ 4 (mod 8)`. The *support* does not separate them and
neither does the tail: **every one of the 59 bins from `|T4| = 228` to
`708` agrees exactly**. The **bulk** separates them — all 27 differing bins
lie at `|T4| ≤ 220`. The first eight:

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| decoded | 2 650 505 561 | 2 383 887 265 | 1 944 278 842 | 1 450 375 604 | 996 815 228 | 634 922 458 | 377 261 905 | 210 018 619 |
| rebuild | 2 650 421 375 | 2 383 982 899 | 1 944 334 154 | 1 450 404 118 | 996 701 290 | 634 954 430 | 377 258 529 | 209 998 559 |
| Δ | −84 186 | +95 634 | +55 312 | +28 514 | **−113 938** | +31 972 | −3 376 | −20 060 |

`run.py` prints all 27. The differences sum to zero, as they must. The
**first** moment, which nothing forces, does not: `Σ |T4|·Δ = −279 888`.

Nothing cheaper could have found this. The largest discrepancy, 113 938 at
`|T4| = 36`, is `1.1·10⁻⁴` of its bin — the same order of invisibility as
at 668, where a `2·10⁷`-draw sampled comparison of the analogous pair read
`max d = 1.8` and saw nothing at all. **A null sampled comparison is worth
very little**; no sampled comparison of the 716 pair is banked in this
repository, and none would have added anything.

## Why the rebuild leaves the class

`G = ℤ₁₇₈` and `ρ = 177 = −1`, so the only character with `ψ² = 1` has
`ψ(ρ) = −1`. `note/NOTE-B.md` §1.4's proposition — a character twist with
`ψ(ρ) = +1` is a diagonal conjugation `S H S`, and therefore manufactures
nothing — **does not apply**. The `i = 2` border is a genuinely different
border, and the exact profile confirms the construction really does produce
a new equivalence class here, as it does at 668. Whether it always does at
`ψ(ρ) = −1` is still **not claimed**: there are now three proven instances
(668, 716, 1676 — cert 20, 2026-09-02) and no general argument.

## Pinned digests

**Matrices** (canonical SHA-256 of the `+/-` serialisation, the digest
`verify/verify.py` reports; both re-derived by `run.py` from the banked
parameters alone):

| matrix | canonical SHA-256 |
| --- | --- |
| decoded `(1,1)` record | `3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6` |
| Lemma-T `i = 2` rebuild | `6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7` |

**Banked files** (SHA-256 of the file bytes, compared in `run.py`):

| file | SHA-256 |
| --- | --- |
| `data/sep716-exact-blas-decoded.json` | `80ee1e151ec1f759d7213d500623603716b9afa6fc382a385ce6970efac35a6b` |
| `data/sep716-exact-bits-decoded.json` | `a0d5b3a65b83c39c905ec2a1d3b25ca1c58e0106b76aaa6eb3b2feee3748aeed` |
| `data/sep716-exact-blas-twisted.json` | `e2076eb890557775edaccda3c9dcbab7e585f5181e0da1bb5914913bf0749b46` |
| `data/sep716-exact-bits-twisted.json` | `c385773a7a2bf4506b94752406be787bac3d85b8a9ee18bb23668080c6afe7bc` |

`data/payload-records.json` and `data/twisted-i2-records.json` are not
file-pinned here on purpose: they are shared with certs 01 and 02, and the
binding pin on each is the canonical digest of the matrix it produces, which
is checked in clause [1] — reinforced, for the twisted record, by the
`ψ`-twist re-derivation that binds it to `payload-records.json` outright.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 1.1 s** (73 checks, exit 0) |
| rebuild + full hypothesis re-check + `verify.py`, per matrix | 0.3–0.4 s |
| `dim V` / `dim W` on both 716 matrices | < 0.2 s |
| controls C1–C3 | 0.2 s |
| producing the banked profiles, `blas` | 124.8 s (decoded), 183.9 s (twisted) |
| producing the banked profiles, `bits` | 228.0 s (decoded), 243.2 s (twisted) |
| **`run.py --full --impl blas`, both recomputed** | **868.5 s** (80 checks, exit 0; measured 2026-09-01) |
| `run.py --full` (both matrices, both paths) | ≈ 1.5 h, **not yet run in this repository** |

The four banked profiles were produced upstream on a shared desktop with
three threads (a second 3-thread job of another lane ran concurrently); the
`producing …` seconds are the values those runs recorded in the JSONs
themselves, and their peak resident sets — 422 MB (`blas`), 73 MB (`bits`)
— are recorded there too.

**Which legs ran where.** The source-laboratory run covered **all four**
legs: `blas` and `bits`, on both matrices, agreeing bin for bin on each.
Inside this repository, **`--full --impl blas` has been run and is green**:
on 2026-09-01 it completed in 868.5 s — 80 checks, no failures, exit 0 —
and each fresh profile matched **both** banked implementations bin for bin
(`decoded` 87 bins in 456 s, `twisted` 87 bins in 411 s). That is the first
in-repo regeneration of these four banks, and — because `full_recompute.py`
takes the upper-triangle route while the banks came from the
canonical-split engine — it is also a third arithmetic route agreeing with
the two the certificate rests on. The `bits` path at order 716 has **not**
been run inside this repository; upstream it cost 228.0 s and 243.2 s. The
ported `bits` code is exercised here against both stdlib routes on the five
C1 controls and, under `--full --impl bits` or `--impl both`, against the
forced Sylvester profile at `n = 128`.

## What is NOT claimed

* **The default run does not claim to have recomputed anything.** The
  `C(716,4)` enumeration ran upstream, in `Hadamard-2060`'s
  `experiments/inequiv/exact_profile_big.py`; the default path audits its
  banked output. Only `--full` recomputes, and only `--full` binds a bank to
  a matrix by computation.
* **Nothing about the transpose-extended relation at 716.** Unlike order
  668 — where cert 08 computes the transposed profiles and all six
  comparisons separate — **the transposed profiles at 716 had not been
  computed when this certificate was written**. Refuting `A ≈ B` under the
  transpose-extended relation needs both `A ≁ B` and `A ≁ Bᵀ`, and only the
  first is in hand here; the statement of this certificate is row-side only.
  The three transposes *have* since been computed and banked — cert 15,
  2026-09-02 — and the transpose-extended statement at 716 now holds; this
  certificate does not rest on it.
* **Nothing about orders 1676 or 1772.** The same construction gives
  Lemma-T `i = 2` rebuilds there, and the corresponding exact computation
  costs roughly `98×` and `130×` the 668 run (the profile is `Θ(n⁵)` in the
  pair-vector Gram). The `1676` computation has since been made — cert 20,
  2026-09-02, NOTE-B.md §3.7 — and nothing in this certificate depends on
  it; `1772` has not been made, so nothing is said there.
* **No claim of novelty of existence at 716**, and no priority claim of any
  kind. Order 716 is long settled by the publicly posted matrix; this
  certificate is about *how many classes* are on the table among the
  artifacts banked here.
* **No claim that these are the only two classes** at order 716.
* **No general statement about `ψ(ρ) = −1`.** Three orders now have it
  proven (668, 716, 1676); that is three instances, not a theorem.
* **Matching invariants prove nothing.** `dim W` agrees on this pair, as
  every cheap invariant did at 668. Read every "agrees" in this repository
  as "did not separate", never as "the same". Beyond `dim V`/`dim W`, no
  cheap-invariant sweep of the 716 pair was made here — cheap agreement at
  716 is neither claimed nor needed.

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/11-716-separation/run.py
python certs/11-716-separation/run.py --full --impl blas
python certs/11-716-separation/run.py --full
```

The default path is standard library only, no network. The last two import
numpy (finder-side only, never in the trust chain) and cap BLAS threads at
three. Exit code 0 iff every check passed.
