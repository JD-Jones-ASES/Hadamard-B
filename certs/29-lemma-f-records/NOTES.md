# cert 29 — Lemma F, and no odd prime acts on the 668 / 716 records

**Label: PROVEN** (Lemma F (i)–(iv), the admissible-prime lists, the 16-block
shape at 523, and the implication from the class sizes; `note/NOTE-B.md` §1.11)
**+ PROVEN-BY-CERTIFICATE** (the per-row `|T4|` class partition at the record
orders, computed twice — see *The trust boundary* below). Default run:
`python certs/29-lemma-f-records/run.py` from the repository root. Standard
library only, **≈ 1.3 s**, exit 0, **40 checks**. That run **audits a banked
exact computation**; `--full` is offered and priced below and **has not been run
in this repository**.

---

## Lemma F

> **Lemma F (`note/NOTE-B.md` §1.11).** Let `H` be Hadamard of **order `N`**
> (the note fixes `n = |G|` and `N = 4v`; this certificate uses `N` throughout)
> and let `(P, Q)` be an automorphism — signed permutation matrices with
> `P H Qᵀ = H` — of odd prime order `p`. Then
>
> * **(i)** the signs are removable by diagonal conjugation: on a `p`-cycle the
>   product of the signs is `+1` because `P^p = I`, so a diagonal makes that
>   cycle unsigned; on a fixed point the sign `e` satisfies `e^p = e = 1`;
> * **(ii)** `#`fixed rows `= #`fixed columns, since `P = H Q Hᵀ / N` gives
>   `tr P = tr Q`;
> * **(iii)** if `(P,Q) ≠ (I,I)` then `f := #`fixed rows `≤ N/p`; and if
>   `f ≤ m := (N−f)/p` then `f ≤ N/(p+1)`;
> * **(iv)** the `m` orbits give `p×p` circulant blocks, the fixed rows are
>   constant on column orbits, and `AAᵀ + p·BBᵀ = N·I_f`.

*(iii)'s argument, in one line: if `f > m` then `rank(BBᵀ) ≤ m < f`, so
`NI − AAᵀ` is singular and `AAᵀ` has eigenvalue `N` with multiplicity `≥ f−m`;
`tr(AAᵀ) = f²` then gives `p f² − N(p+1) f + N² ≥ 0`, whose roots are `N/p` and
`N`, and `f = N` is the identity.*

**At `N = 2092`, `p = 523`.** `523 | 2092`, so `f ≡ 0 (mod 523)`, and (iii)
gives `f ≤ 4`: hence `f = 0`, `m = 4` — a **16-block circulant array**.

**The repair, and it matters.** That is *all* it is. "The all-type-1 Williamson
array" is a **4-seed subfamily**; Lemma F forces the block shape and no more.
The Goethals–Seidel array is not of that shape — its off-diagonal blocks are
back-circulants — and its translation automorphisms would need `2a = 0`, which
no odd `|G|` has.

**The admissible primes.** An odd `p` can act only if some `f` with
`f ≡ N (mod p)` and `0 ≤ f ≤ N/p` exists (`f = 0` only when `p | N`). At 2092
that leaves **22** odd primes and excludes **293**, among them 131 and 349. A
constant-strip multi-block border with four fixed rows needs `4 | b`, since
`|(AAᵀ)_{rr′}| ≤ 4 < p` forces `(BBᵀ)_{rr′} = 0` (every odd prime but 3) and
makes `B` a `4×b` partial Hadamard matrix.

## The records

For row `i` let `π(i)` be the multiset of `|T4(i,j,k,l)|` over 3-subsets of the
other rows. A signed row/column permutation preserves `π`, so the `π`-classes
are an invariant partition of the rows. If an automorphism of odd prime order
`p` exists, every class is a union of `p`-cycles and fixed rows, so it
contributes at least `|C| mod p` fixed rows, while (iii) caps the total at
`N/p`. Hence

```
Σ_C (|C| mod p)  >  N/p     ⟹     no automorphism of order p.
```

`|T3|` is **not** invariant under column signs and is not used anywhere here.

| matrix | canonical SHA-256 | classes | sizes |
| --- | --- | --- | --- |
| `H(668)` | `bdeb5059…a387b0` | 336 | 4 singletons + 332 `τ`-pairs |
| `H′(668)` | `600849b0…2008a3` | 336 | 4 + 332 |
| `H″(668)` | `af1c285c…2953c7` | 336 | 4 + 332 |
| `H(716)` | `3adcb1bb…664ee6` | 360 | 4 + 356 |

The four singletons are the **border rows**; `τ` — the shift by `n/2` on every
superblock, `n = |G|`, so 83 at 668 and 89 at 716 — is an automorphism, with
border sign `+1` on the `(1,1)` records and `−1` on `H′`, and it pairs the rest.
Every class has size 1 or 2, so the forced fixed-row count is `N` itself, which
exceeds `N/p` for every odd `p`.

> **Theorem.** No automorphism of odd prime order acts on `H(668)`, `H′(668)`,
> `H″(668)` or `H(716)`. `Aut` is a **2-group** at each; every automorphism
> fixes the four border rows up to sign and acts within `τ`-pairs.

## The trust boundary — where the class partition was computed

The **implication** — those class sizes ⟹ no odd prime acts — is PROVEN, and
this certificate applies it to the pinned sizes on every run. Its **input**, the
per-row `|T4|` partition at the record orders, is a heavy computation, and it
has been made **twice**, both times in the **source laboratory**
(`Hadamard-2060`):

1. by a **numpy pair-histogram finder**
   (`skeptic-pass/c9b_lemma_f_records.py`, 2026-09-02), which produced the class
   counts banked in `pins.json`;
2. by **this script's route B under `--full`** — standard library, exact
   integers, no numpy — run on all four records as detached single-core
   processes on 2026-09-03/04, `ALL CHECKS PASS` in each case
   (`certs/0028-lemma-f-records/replay-2026-09-03/`):

| record | `--full` wall, one core |
| --- | --- |
| `H(668)` | 9 044.0 s |
| `H′(668)` | 9 062.2 s |
| `H″(668)` | 9 049.1 s |
| `H(716)` | 11 273.6 s |
| **total** | **38 428.9 s ≈ 10.7 core-hours** (≈ 3.1 h wall on four workers) |

That second run is the independent implementation **D-008** asks for, and it is
why the records claim is **PROVEN-BY-CERTIFICATE** rather than
COMPUTATIONAL-EVIDENCE. Before it, the claim rested on one implementation and
the source laboratory's own certificate said so: decision **D-068** item 5
recorded the label as COMPUTATIONAL-EVIDENCE *until `--full` has run*, and the
desk's audit had caught the draft calling the numpy finder "the second
implementation already run" — inverted.

**It has not been run in this repository.** `--full` here is the in-repo replay:
offered, wired, priced. Say **"banked exact computation audited"** of a default
run; the word *replayed* belongs to `--full`, and at these orders it belongs to
the source laboratory's run and not to this one. This is the posture certs 20,
21, 22, 23, 24 and 25 take with their own heavy legs.

## The evidence chain

**[A] Lemma F on small matrices, and the arithmetic at 2092.** (i)–(iv)
exercised on Sylvester `H(8)` (automorphisms of orders 7 and 3) and Paley
`H(12)` (orders 11, 5, and a signed Möbius map of order 3), each time checking
the fixed-row bound, the circulant orbit blocks and `AAᵀ + pBBᵀ = N I_f`. Then
22 admissible / 293 excluded odd primes at 2092, with 131 and 349 among the
excluded; `523 ⟹ f = 0`; and `4 | b` for a `4×b` partial Hadamard matrix,
`b = 2..12`.

**[B] Two `|T4|` implementations on controls (D-008).** Route A (direct
`C(N−1,3)` per row) and route B (the pair-histogram) agree on Sylvester `H(8)`,
Paley `H(12)`, and the `GS(28)` array on Williamson(7) seeds. The two
Aut-transitive controls return **one** class, as they must; `GS(28)` returns
more than one, which is the case the implication is about. Note what this does
and does not establish: the two implementations agree **on the controls**, at
orders where both are cheap. At the record orders the second implementation is
`--full`, above.

**[C] The four records, reassembled here, and the implication.** Both banked
data files are **SHA-256 file-pinned**. Each of `H(668)`, `H′(668)`, `H″(668)`,
`H(716)` is rebuilt in this run and its canonical digest matched against the
pin. `τ` is exhibited as an automorphism of each, with the border sign the
structure predicts. The twisted record's assembled digest is additionally
matched against the `pinned_sha256` the data file itself carries. Then the
implication is applied to the pinned class sizes, and the pinned partition is
checked for internal consistency — the sizes total `N`, the class count is the
pinned one, and the number of singletons is the pinned four.

**[D] `--full [ --matrix NAME ]`.** Recomputes route B on the named record (or
all four) **in this repository** and requires the pinned class sizes, the four
border singletons, the `τ`-pairing, and an empty surviving-prime list.

## What the port rewired

The source laboratory's certificate read three lab-only inputs. This one reads
only this repository's `data/`:

| lab input | here |
| --- | --- |
| `experiments/pr0023/alpoge_full_decode.json` | `data/payload-records.json` (same record schema; file-pinned) |
| `experiments/bordered_gs_theorem/twist_report.json` (the twisted border) | `data/twisted-i2-records.json`, which banks the twisted instance **whole** — seeds and border together (cert 02's record; file-pinned) |
| the committed pass-II output files, audited for their text | **dropped.** A text audit of a finder's own log is not evidence this repository can stand behind. In its place: the `ψ`-twist of the decoded seeds is **re-derived here** and required to equal the twisted record's seeds, the assembled digest is matched against the data file's own `pinned_sha256`, and the pinned partition is checked for internal consistency. `pins.json` keeps the digests, class counts and sizes, and carries a `source` field naming both computations and their walls. |

Nothing else changed. The mathematics, the two `|T4|` routes, the Lemma-F
checks and the `--full` path are the laboratory's code unaltered.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **≈ 1.3 s** (exit 0, 40 checks; measured 2026-09-05 on the desk) |
| `--full`, per 668 record | ≈ 2.5 h on one core (**not run here**; the laboratory's walls are tabled above) |
| `--full`, `H(716)` | ≈ 3.1 h on one core (**not run here**) |
| `--full`, all four | ≈ **10.7 core-hours** |

## What is NOT claimed

* **The 2-part of `Aut` is not computed.** `Aut = ⟨−I, τ⟩` is *consistent* with
  the row classes at all four records and is **not proved**. Every automorphism
  fixes the border rows up to sign and acts within `τ`-pairs; how many there are
  is open.
* **`H′(716)` and `H″(716)` were never run.** The statement covers the four
  records named and no others.
* **The default run recomputes nothing** at the record orders. Only `--full`
  would, and it has not been run here.
* **Nothing about `Aut` of a hypothetical `H(2092)`.** Lemma F constrains what
  an odd-prime automorphism of one would have to look like; it does not say one
  exists, and §1.11 asserts nothing about existence at 2092.
* **`|T3|` is not used.** It is not invariant under column signs, and no claim
  rests on it.

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/29-lemma-f-records/run.py
```

The default path is standard library only, no network, no numpy, and is the
whole certificate. The `--full` paths are hours each and are also standard
library only:

```
python certs/29-lemma-f-records/run.py --full --matrix "H(668)"
python certs/29-lemma-f-records/run.py --full
```

Exit code 0 iff every check passed.

## Provenance

Lemma F, the records computation and this certificate's code are the source
laboratory's (`Hadamard-2060`, `certs/0028-lemma-f-records`, filed by a Cursor
cloud lane under skeptic pass III and adopted as decision **D-068** after a desk
replay and an independent audit; the `n`/`N` alignment with `note/NOTE-B.md`
§1.0, the 16-block repair, and the COMPUTATIONAL-EVIDENCE label pending
`--full` are that audit's repairs, D-068 items 8(d), 8(e) and 10). The `--full`
run that upgraded the label was made at the source laboratory's desk on
2026-09-03/04 and is recorded in that certificate's `NOTES.md` replay record and
in the dated desk note appended to `NOTE-B-draft-2.md`. Credit is to stations,
as everywhere in this repository.
