# cert 06 — order 668 carries at least two Hadamard equivalence classes

**Label: PROVEN.** Replay: `python certs/06-668-separation/run.py` from the
repository root. Standard library only, about one second, exit 0.

---

## The theorem

> **Theorem.** Let `H₁` be the decoded `(s, i) = (1, 1)` bordered
> Goethals–Seidel record at order 668 (`data/payload-records.json`) and let
> `H₂` be its Lemma-T `i = 2` rebuild (`data/sep668-twisted-record.json`).
> Then `H₁` and `H₂` are **not Hadamard-equivalent**: there is no
> `H₂ = D_r P_r H₁ P_c D_c` with `P` permutation matrices and `D` diagonal
> `±1`.

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all
`C(668,4) = 8 222 179 035` row 4-subsets, with
`T4 = Σ_c H[i][c]H[j][c]H[k][c]H[l][c]`, is a Hadamard-equivalence
invariant. The two profiles populate the **same 80 bins**, and **26 of the
80 bin counts differ**. An invariant that differs is a separation. ∎

The computation is finite, exact, and replayed by `run.py`; the two
profiles are banked in `data/`.

## Why the |T4| profile is an invariant

`note/NOTE-B.md` §3.1, invariant **I5**. In one line: a row negation
contributes exactly one sign to `T4` so `|T4|` is fixed; a column negation
contributes `d_c⁴ = 1`; row and column permutations relabel the 4-subsets
and the coordinates. Hence the *multiset* of `|T4|` values is carried to
the corresponding multiset by any element of the equivalence group.

Two warnings from §3.1 that this certificate honours mechanically:

* **the signed `T4` histogram is not an invariant** — only `|T4|` is. No
  signed histogram is banked or compared anywhere in this repository, and
  cert 07's control C2 exhibits a signed profile moving under a legal row
  negation while the folded one stays fixed.
* **`dim V` is not an invariant** — only `dim W = V + ⟨1⟩` is. On this very
  pair `dim V` **differs** (666 vs 667) and the difference is worthless;
  `run.py` measures it (control C4) precisely so it is on the record as a
  trap and never as evidence. `dim W = 667` on both.

## The evidence chain

`run.py` runs all of the following before it prints a verdict.

**[0] The bank is pinned.** Every banked file's SHA-256 is compared in code.

**[1] Both matrices are rebuilt, not assumed.** Each record goes through
`tools/bordered_gs.py`'s `check_record`, which re-checks **every hypothesis
of the master theorem** — H0 shape, H1 the two-tier PAF profile, H2 the
corner/row-table budget, H3 the column-table Gram, H4 the coupling (that
module's numbering; see its `LABEL MAPPING` block, which is a permutation of
`note/NOTE-B.md` Theorem A's), the
derived D1/D3/D5 and the Σ̄ law, and the compression-lemma cross-check of
the assembled core against the `G/K` Goethals–Seidel array — and only then
assembles. Each assembled matrix is handed to `verify/verify.py`, the trust
chain, and its canonical digest compared against the pin. The matrices are
deleted afterwards; nothing multi-megabyte is committed.

**[2] Each banked profile is audited in exact integers.** Per profile:
every populated bin is `≡ 4 (mod 8)`; the counts total `C(668,4)`; and the
second moment equals `n³(n−1)(n−2)/24 = 5 517 193 410 096` — the closed
form proved in `note/NOTE-B.md` §3.1, a 13-digit number hit to the unit by
both matrices and neither tuned for. Three of the four banked JSONs predate
the addition of a `second_moment` field upstream, so `run.py` **recomputes
the second moment from the profile itself** rather than trusting a banked
number; where the field is present it is checked too.

**[3] Two independent implementations agree bin for bin.** `blas` is a
float32 Gram of the pair-vector matrix (exact at these sizes: every entry
and every partial sum is an integer below `2²⁴`); `bits` packs rows into
`uint64` words and uses `|T4| = |n − 2·popcount(u_P ⊕ u_Q)|`. They share no
arithmetic. Both hit `C(668,4)` and the second moment, and agree on all 80
bins of both matrices.

**[4] Controls.**

| | control | result |
| --- | --- | --- |
| **C0** | every control matrix is checked to *be* Hadamard before anything is asked of it | pass, 5 matrices |
| **C1** | the full `\|T4\|` profile of Sylvester `H(8)`, Sylvester `H(16)`, Paley `H(20)`, and two plain Goethals–Seidel arrays `H(28)` and `H(36)` built by this repository's own assembler at the degenerate `s = 0` layer (three- and four-bin profiles, so the agreement is not vacuous), each computed **twice** — by straight `O(C(n,4))` enumeration over the `±1` entries, and by the pair-vector / Gram-triangle route the banked profiles use | agree bin for bin on all five; both routes hit `C(n,4)` and the second moment; the two Sylvester cases additionally match their **forced** profile `{0: C(n,4) − n(n−1)(n−2)/24, n: n(n−1)(n−2)/24}`, so C1 is a positive control with a predicted answer |
| **C2** | *negative* control: a banked profile is corrupted in a **total-preserving** way (one count moved from `\|T4\| = 4` to `\|T4\| = 660`), so only the second-moment identity can catch it, and `audit()` is required to raise | the assert fires |
| **C3** | the `dim V` trap on Sylvester `H(16)` under a deterministically seeded (`20260831`) random signed row negation | `dim V` moves `4 → 5`; `dim W` stays `5`; the `\|T4\|` profile is unchanged |
| **C4** | the same trap on the real objects | `dim V` 666 (decoded) vs 667 (twisted) — **differs and is worthless**; `dim W = 667` on both |

The C1 route-2 code is the same bookkeeping the banked numpy profiles
depend on (drop the diagonal, drop the `n·C(n−1,2)` index-sharing pairs,
divide by three), run where straight enumeration can check it. That is what
makes the banked profiles auditable rather than merely reproducible.

**[5] `--full`.** `python certs/06-668-separation/run.py --full` recomputes
both profiles from the rebuilt matrices with numpy
(`full_recompute.py`, both arithmetic paths) and compares them to the bank
bin for bin. numpy is imported **only** under this flag, is finder-side
only, and is never in the trust chain; BLAS threads are capped at three,
set in `run.py` before numpy loads. Before spending the hour, `--full`
smoke-tests both ported paths against the forced profile of Sylvester
`H(128)`, which needs two `uint64` words per row and so exercises the
multi-word packing that the small controls cannot.

## The separation

Same 80 bins, `≡ 4 (mod 8)` from 4 to 660 (`620`, `636`, `652` empty in
both). The *support* does not separate them and neither does the extreme
tail — `604: 19, 612: 18, 628: 2, 644: 1, 660: 1` in both. The **bulk**
separates them.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| decoded | 2 073 064 058 | 1 852 054 148 | 1 491 070 735 | 1 091 442 371 | 732 009 734 | 452 971 620 | 260 220 030 | 139 599 212 |
| rebuild | 2 073 109 602 | 1 852 009 274 | 1 491 079 303 | 1 091 478 493 | 731 928 178 | 452 961 444 | 260 257 738 | 139 618 428 |
| Δ | +45 544 | −44 874 | +8 568 | +36 122 | **−81 556** | −10 176 | +37 708 | +19 216 |

`run.py` prints all 26. The differences sum to zero, as they must. The
**first** moment, which nothing forces, does not: `Σ |T4|·Δ = −306 848`.

Nothing cheaper could have found this. The largest discrepancy, 81 556 at
`|T4| = 36`, is `1.1·10⁻⁴` of its bin; a `2·10⁷`-draw sampled comparison of
this pair reads `max |z| = 1.8` with zero bins over `4σ`, i.e. it sees
nothing at all (that null reading is banked in
`data/sep668-sampled-histograms.json` and used as cert 07's calibration
control). **A null sampled comparison is worth very little**, and this pair
is the proof of it.

## Why the rebuild leaves the class

`G = Z₁₆₆` and `ρ = 165`, so the only character with `ψ² = 1` has
`ψ(ρ) = −1`. `note/NOTE-B.md` §1.4's proposition — a character twist with
`ψ(ρ) = +1` is a diagonal conjugation `S H S`, and therefore manufactures
nothing — **does not apply**. The `i = 2` border is a genuinely different
border, and the exact profile confirms the construction really does produce
a new equivalence class here.

## Pinned digests

**Matrices** (canonical SHA-256 of the `+/-` serialisation, the digest
`verify/verify.py` reports; both re-derived by `run.py` from the banked
parameters alone):

| matrix | canonical SHA-256 |
| --- | --- |
| decoded `(1,1)` record | `bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0` |
| Lemma-T `i = 2` rebuild | `600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3` |

**Banked files** (SHA-256 of the file bytes, compared in `run.py`):

| file | SHA-256 |
| --- | --- |
| `data/sep668-twisted-record.json` | `fe8154179ba2ebfe097c82e468368cdc8a070548555bb10140949af0560611fb` |
| `data/sep668-exact-blas-decoded.json` | `22df5ce9fcd6eb307f56981c507bb46b2a18b79861d903349dc13458a6dffcbf` |
| `data/sep668-exact-bits-decoded.json` | `0bafbf8219d33b9c74786700106aeba3086bbf577ee02bcda43768f35978fdd8` |
| `data/sep668-exact-blas-twisted.json` | `c4d8db3ba40cf8c5a244607032dab6b66d878b8fe6b98784351f7b8ae70e5a17` |
| `data/sep668-exact-bits-twisted.json` | `91d154d05ccea87a6fa98a02b4fcbf275dc6b4025650116941647216a69faf5a` |

`data/payload-records.json` is not file-pinned here on purpose: it is
shared with cert 01, and the binding pin on it is the canonical digest of
the matrix it produces, which is checked above.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 1.0 s** (66 checks, exit 0) |
| rebuild + full hypothesis re-check + `verify.py`, per matrix | 0.3 s |
| `dim V` / `dim W` on both 668 matrices | < 0.1 s |
| controls C1–C3 | 0.2 s |
| producing the banked profiles, `blas` | 276.4 s (decoded), 274.6 s (twisted) |
| producing the banked profiles, `bits` | 1 479.4 s (decoded), 1 465.0 s (twisted) |
| `run.py --full` (recompute both, both paths) | ≈ 1 h |

The four banked profiles were produced upstream on a shared desktop with
`OMP_NUM_THREADS = 3`; the seconds above are the values the runs recorded
in the JSONs themselves. `--full` was **not** executed at order 668 while
this certificate was written, to avoid contending with a concurrent local
2060 run; the ported paths were validated against both stdlib routes on
the five C1 controls and against the forced Sylvester profile at
`n = 8, 16, 128`.

**Third-engine corroboration** (observed, not banked): a separate
memory-aware engine upstream (`exact_profile_big.py`, canonical-split
enumeration, 32-way, 3 threads, 410 MB peak) recomputed the decoded
profile in 84.2 s and returned it **bin for bin identical** to the banked
`blas` and `bits` profiles, on the matrix whose digest it independently
records as `bdeb5059…`. That is a third arithmetic route agreeing with the
two the certificate relies on. It is cited rather than banked because it is
an in-flight artifact of another lane's run.

## What is NOT claimed

* **Nothing about orders 716, 1676 or 1772.** The same construction gives
  Lemma-T `i = 2` rebuilds at those orders, and the corresponding exact
  computation costs roughly `1.4×`, `98×` and `130×` the 668 run
  (the profile is `Θ(n⁵)` in the pair-vector Gram). It has not been made,
  so nothing is said.
* **No claim of novelty of existence at 668.** Order 668 was settled by
  the publicly posted matrix; this certificate is about *how many classes*
  are on the table among the artifacts banked here, not about existence.
* **No claim that these are the only two classes** at order 668.
* **Matching invariants prove nothing.** Every cheap invariant on this pair
  — `dim W`, the collision profiles, the exact extreme strata, `rank₁₆₇`,
  the dual weight enumerator, and `2·10⁷` sampled 4-subsets — returned
  *identical* values on a pair that is in fact inequivalent. Read every
  "agrees" in this repository as "did not separate", never as "the same".
  Of that list only `dim V`/`dim W` and the sampled comparison are replayed
  by `run.py`; the I1/I2 collision profiles, the I4 extreme strata,
  `rank₁₆₇`, the dual weight enumerator and the block-affine exhaustion are
  **MEASURED upstream, not replayed here**.
