# cert 23 — order 1772 carries at least three Hadamard equivalence classes

**Label: PROVEN + PROVEN-BY-CERTIFICATE.** Default run:
`python certs/23-1772-three-classes/run.py` from the repository root.
Standard library only, about six seconds, exit 0, 132 checks. That run
**audits a banked exact computation**; `--full` is offered and priced below
and **has not been run in this repository at this order**. The trust
boundary is the one certs 11, 13, 14, 20 and 22 draw, and is set out below.
**Row-side only**: the transposed profiles at 1772 are a separate leg of the
same campaign — their `blas` legs are in, their `bits` legs were still
running when this certificate was written, and this repository certifies a
profile only once `blas` and `bits` agree bin for bin — so they are not
banked here and **nothing at all is claimed under the transpose-extended
relation at this order**. That is exactly cert 20's caveat at 1676, which
cert 21 later discharged there; a later certificate will do the same here.

---

## The theorems

> **Theorem 1 (the twist, a fourth instance).** Let `H` be the decoded
> `(s, i) = (1, 1)` bordered Goethals–Seidel record at order 1772
> (`data/payload-records.json`; cert 01) and let `H'` be its Lemma-T
> `i = 2` rebuild (`data/twisted-i2-records.json`, order 1772; cert 02).
> Then `H` and `H'` are **not Hadamard-equivalent**: there is no
> `H' = D_r P_r H P_c D_c` with `P` permutation matrices and `D` diagonal
> `±1`.

> **Theorem 2 (three classes).** Let `H″` be `H` with its **twelve
> off-diagonal core blocks negated** and the `4×4` corner, the four row
> strips and the four column strips unchanged. Then `H″` is a Hadamard
> matrix, and it is **not Hadamard-equivalent** to either of `H` and `H'`.
> With Theorem 1, **order 1772 carries at least three Hadamard equivalence
> classes.**

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all
`C(1772,4) = 409 422 905 815` row 4-subsets is a Hadamard-equivalence
invariant (`note/NOTE-B.md` §3.1, invariant I5). All three profiles populate
the same 89 bins; `H` differs from `H'` in **57** of them, `H″` from `H` in
**58**, and `H″` from `H'` in **53**. An invariant that differs is a
separation. ∎

**Consequence.** The Lemma-T construction at `ψ(ρ) = −1` provably leaves the
equivalence class at a **fourth** order — 668 (cert 06), 716 (cert 11) and
1676 (cert 20) were the first three. Four instances; `note/NOTE-B.md`
§1.4's question is **still not answered in general**, and no general theorem
is claimed here. 1772 was the last decoded `(1,1)` order left unclaimed:
cert 02's list of four orders is now complete.

## What `H″` is, and what the theorem says about orientation

`note/NOTE-B.md` §1.0 fixes a *standard* Goethals–Seidel orientation and
calls the other one (the six transposed blocks negated) a convention. Put
`S = diag(I₄, diag(1,−1,−1,−1) ⊗ I_n)`, `n = 442`. Then `S·H″·S` is exactly
the same seeds and border assembled in the alternate orientation with the
border strips signed by superblock — `P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`
— an identity of sign patterns that `run.py` checks cell by cell over all
`1772² = 3 139 984` cells (clause [1]). So the theorem says: **at 1772, as at
668, 716, 1676 and 2060, the GS orientation is not a gauge for Hadamard
equivalence** — the two orientations of one record are two classes — and
every bordered GS record found at any order carries a second candidate class
for free, by negating twelve blocks.

**Corollary (twist versus orientation, at four orders).** The `ψ(ρ) = −1`
Lemma-T twist and the orientation switch land in **different** classes at
1772 (`H″ ≁ H'`, 53 bins) exactly as they do at 1676 (cert 20, 66 bins), at
716 (cert 14, 25 bins) and at 668 (cert 13, 27 bins). **No structural
account of why is claimed here**: at 1772 the statement rests on the profiles
alone. Cert 13 records such an account at 668 — the two moves agree on the
core and differ on the border — from a transport identity checked in the
source laboratory; nothing corresponding was computed or checked at 1772, so
nothing corresponding is asserted. That is now four orders with the same
verdict, and **no general statement is made**: nothing here says the twist
and the orientation switch differ at every order, or at any order not
computed.

## The trust boundary — what a default run does and does not establish

The three `C(1772,4)` enumerations **were not run inside this repository.**
They ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py`, the engine unchanged since the
2060 registration, numpy, **16 threads on a rented `c2d-highcpu-16`
(`prof42-1`, `us-east1-b`), 2026-09-02 23:33Z – 2026-09-03 11:24Z** — under
the pre-registration `experiments/pr0042/REGISTRATION.md`, **flushed
2026-09-02 10:17 UTC, before any matrix it governs was built**, which fixed
the decision rules and the kill criteria in advance (`blas ≠ bits` in any
bin: hard stop, no claim; a builder digest mismatch: nothing uploaded). §2
of that registration names `H_1772-decoded`, `H_1772-twisted` and
`H_1772-orient` in as many words, and §4 fixes the rule "`H ≁ H'` ⟹ the
Lemma-T twist leaves the class at that order; `H″ ≁ H` and `H″ ≁ H'` ⟹ the
order carries at least THREE classes (row-side)". Its output is banked in
`data/sep1772-{decoded,twisted,orient}-exact-{blas,bits}.json`.

The matrices themselves were **built and verified at the desk**
(`build_matrices.py`, digests in `manifest.json`, each through this
repository's `verify/verify.py`), and the two pinned ones reproduced their
cert-01 and cert-02 pins before anything was uploaded; **the rented machine
enumerated and nothing else.** It never assembled or verified a matrix.

**A default `run.py` establishes:**

* the six bank files are byte-for-byte the ones pinned in `run.py` (SHA-256
  compared in code);
* `H` rebuilds from the banked record through the full master-theorem
  hypothesis re-check (H0–H4, D1/D3/D5, Σ̄, the compression-lemma
  cross-check), passes `verify/verify.py` over all `1 569 106` row pairs,
  and carries cert 01's canonical digest `1852e951…77236ba2` — computed
  in-process;
* `H″` is formed from those rows by negating exactly `12·n² = 2 344 368`
  cells, passes `verify/verify.py`, and carries the canonical digest
  `7f1fae05…8ba607e53`; the alternate-orientation identity holds cell by
  cell;
* `H'` rebuilds too, with its seeds **re-derived here** as the `ψ`-twist of
  the decoded seeds (`ψ(g) = (−1)^g` on `ℤ₄₄₂`; `ρ = 441` is odd, so
  `ψ(ρ) = −1`), compared character for character against the bank, so cert
  02's shared record is bound to `payload-records.json` by computation
  rather than by a file pin;
* the three matrices carry three **distinct** canonical digests;
* each of the six banks declares `matrix_canonical_sha256`, and for all
  three matrices it is compared against the in-process digest of the matrix
  rebuilt in this same run — and against the producer's own `matrix_sha256`
  field;
* every banked profile satisfies the forced identities: 89 bins, all
  `≡ 4 (mod 8)` (`1772 ≡ 4 mod 8`), every key canonical in `[0, 1772]` and
  every count a positive integer, total `C(1772,4) = 409 422 905 815`,
  second moment `n³(n−1)(n−2)/24 = 726 727 740 809 840` — recomputed here
  and also compared against the `second_moment`, `second_moment_want`,
  `total`, `n` and `C_n_4` fields the bank declares, along with its
  `schema`, `folded`, `impl`, `matrix`, `producer_filename` and `engine`;
* `blas == bits` bin for bin on each of the three matrices;
* the three separations — 57, 58 and 53 bins — **in both arithmetics**, with
  their difference vectors summing to zero, their first moments non-zero,
  the tail that does not separate, and the agreeing bins below the top
  divergent one enumerated;
* the theorems themselves, re-derived in clause [4] from the counts rather
  than restated, together with the assertion that no transposed 1772 profile
  is banked here.

**A default run does not establish that the banks were computed from the
matrices `run.py` rebuilt.** They are *producer-banked*: the digest each
carries is the one the engine recorded against the file it enumerated, and
it equals the digest this certificate pins. A self-declared digest is
metadata, not a computation. `--full` is what would close that gap — and at
this order it has not been run (see *Runtimes*).

## Why `|T4|` is an invariant

For a 4-subset of rows, `T4 = Σ_c H[i][c]H[j][c]H[k][c]H[l][c]`. Row
negation contributes one sign to `T4`, so `|T4|` is fixed; column negation
contributes `d_c⁴ = 1`; row and column permutations relabel the multiset
without changing it. Hence the multiset `{|T4|}` over all `C(n,4)` row
4-subsets is constant on Hadamard-equivalence classes
(`note/NOTE-B.md` §3.1, invariant I5), and two matrices whose profiles
differ in a single bin are inequivalent. Transpose is *not* in the group,
which is why the row-side caveat is a real restriction and not a formality —
here as at 1676 when cert 20 was written, since no transposed profile at
1772 is banked here at all.

Two warnings from §3.1 this certificate honours mechanically: the **signed**
`T4` histogram is not an invariant (every banked profile declares its
folding in its own `folded` field), and **`dim V` is not an invariant** —
only `dim W = V + ⟨1⟩` is. On the 1772 objects `dim V` reads 1770 / 1771 /
1770 and is worthless; `dim W = 1771` on all three. `run.py` measures both
(control C6) precisely so the trap is on the record and never as evidence.

## The evidence chain

**[0]** six file pins. **[1]** rebuild `H` (hypotheses H0–H4, D1/D3/D5, Σ̄,
compression lemma), verify, pin; check the layout `n = |G| = 442`, `s = 1`,
`N = 4(n+s)`; form `H″`, count the negated cells, check `S·H″·S = H_alt`,
verify, pin; re-derive the `ψ`-twist seeds and rebuild `H'`, verify, pin;
three distinct digests. **[1b]** control C6 — the dim-`V` trap on the real
objects. **[2]** audit six banks in exact integers; declared headers; matrix
identity per bank against an in-process digest; `blas == bits` ×3. **[3]**
the three separations in both arithmetics, each with identical support, the
pinned differing-bin count, the differences summing to zero, a non-zero
first moment, the tail that does *not* separate, the agreeing bins below the
top divergent one, and the first eight divergent bins printed. **[4]** the
two theorems derived in code from the counts, and the row-side caveat
asserted. **[5]** controls — C0/C1 five small Hadamard matrices profiled by
straight enumeration and by the pair-vector route `--full` takes (Sylvester
`H(8)`, `H(16)` against their **forced** profiles; Paley I `H(20)`; GS
`H(28)`, `H(36)`, whose three- and four-bin profiles keep the agreement from
being vacuous); C2 the orientation switch applied to GS `H(28)` and `H(36)`,
still Hadamard by brute force, moving exactly `12v²` cells and satisfying
the same sign-pattern identity; C3 the comparator in the null direction; C4
a total-preserving corruption of a banked profile that only the
second-moment identity can catch, required to be caught; C5 the dim-`V` trap
on Sylvester `H(16)` under a seeded signed row negation. **[6]** `--full`:
`certs/06-668-separation/full_recompute.py` imported (not copied),
smoke-tested on the forced profile of Sylvester `H(128)`, then run on the
verified rows — offered, priced, **not run here**.

## The separations

Same 89 bins, `≡ 4 (mod 8)`, on all three matrices — the support separates
nothing. The support runs `4, 12, …, 612` unbroken, then `628–644`,
`660–692`, and the isolated `708`, `740`, `772`, `1764`. Neither does the
extreme tail separate: every bin above `|T4| = 476` agrees on both
comparisons with `H` (29 bins, up to 1764), and every bin above
`|T4| = 452` on `H'` vs `H″` (32 bins). The **bulk** separates: of the 60
bins up to `|T4| = 476`, only `444`, `460`, `468` agree for `H` vs `H'` and
only `460`, `468` for `H` vs `H″`; of the 57 bins up to `|T4| = 452`, only
`380`, `404`, `428`, `436` agree for `H'` vs `H″`. The first eight bins of
each pair, as `run.py` prints them:

**`H` vs `H'` (the Lemma-T twist)** — 57 of 89 bins differ; largest
`|Δ| = 2 560 842` at `|T4| = 4`, i.e. `4.13·10⁻⁵` of that bin;
`Σ |T4|·Δ = +55 798 800`.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H` | 61 934 029 130 | 59 717 426 100 | 55 520 798 861 | 49 778 171 281 | 43 034 414 080 | 35 886 330 792 | 28 860 579 439 | 22 388 208 089 |
| `H'` | 61 931 468 288 | 59 716 764 308 | 55 522 617 277 | 49 780 050 525 | 43 034 542 248 | 35 884 341 470 | 28 861 416 579 | 22 388 284 353 |
| Δ | **−2 560 842** | −661 792 | +1 818 416 | +1 879 244 | +128 168 | −1 989 322 | +837 140 | +76 264 |

**`H` vs `H″` (the orientation switch)** — 58 of 89 bins differ; largest
`|Δ| = 3 067 068` at `|T4| = 44`, i.e. `8.55·10⁻⁵` of that bin;
`Σ |T4|·Δ = +72 493 792`.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H` | 61 934 029 130 | 59 717 426 100 | 55 520 798 861 | 49 778 171 281 | 43 034 414 080 | 35 886 330 792 | 28 860 579 439 | 22 388 208 089 |
| `H″` | 61 931 541 520 | 59 716 317 614 | 55 523 091 925 | 49 779 015 109 | 43 035 240 878 | 35 883 263 724 | 28 862 137 199 | 22 388 244 191 |
| Δ | −2 487 610 | −1 108 486 | +2 293 064 | +843 828 | +826 798 | **−3 067 068** | +1 557 760 | +36 102 |

**`H'` vs `H″`** — 53 of 89 bins differ; largest `|Δ| = 1 077 746` at
`|T4| = 44`, i.e. `3.00·10⁻⁵` of that bin; `Σ |T4|·Δ = +16 694 992`.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H'` | 61 931 468 288 | 59 716 764 308 | 55 522 617 277 | 49 780 050 525 | 43 034 542 248 | 35 884 341 470 | 28 861 416 579 | 22 388 284 353 |
| `H″` | 61 931 541 520 | 59 716 317 614 | 55 523 091 925 | 49 779 015 109 | 43 035 240 878 | 35 883 263 724 | 28 862 137 199 | 22 388 244 191 |
| Δ | +73 232 | −446 694 | +474 648 | −1 035 416 | +698 630 | **−1 077 746** | +720 620 | −40 162 |

`run.py` prints these eight rows of each pair; the full 57-, 58- and 53-bin
lists are the banked JSONs. All three difference vectors sum to zero, and
all three largest discrepancies are of order `10⁻⁵` of their bin — invisible
to any sample of practical size. As at 668, 716 and 1676, nothing cheaper
than the exact profile could have found this, and no sampled comparison of
these matrices is banked here because none would have added anything.

## Why the twist and the switch leave the class

`G = ℤ₄₄₂` and `ρ = 441 = −1`, so the only character with `ψ² = 1` has
`ψ(ρ) = −1`. `note/NOTE-B.md` §1.4's proposition — a character twist with
`ψ(ρ) = +1` is a diagonal conjugation `S H S`, and therefore manufactures
nothing — **does not apply**. The `i = 2` border is a genuinely different
border, and the exact profile confirms that the construction really does
produce a new equivalence class here, as it does at 668, 716 and 1676.
Whether it always does at `ψ(ρ) = −1` is still **not claimed**: there are now
**four** proven instances and no general argument.

The same §1.4 remark identifies the orientation switch as the *other* thing
the `ψ(ρ) = −1` conjugation produces — the array with its twelve
off-diagonal blocks negated, which is not of rank-one sign pattern and so is
not a diagonal conjugation either. At 1772 the two moves land in different
classes, as at 668, 716 and 1676. Four orders, four verdicts, one pattern,
and **no theorem**.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports) —
all three rebuilt in clause [1] of every run:

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `H` decoded `(1,1)` record | `1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2` | yes (cert 01's pin) |
| `H'` Lemma-T `i = 2` rebuild | `82484769a28ac93201f208ca3256bfd491f8edc8bb5e0309764c4f609a113378` | yes (cert 02's pin; seeds re-derived) |
| `H″` orientation switch | `7f1fae050def5b9b7bdc491c05b24551465cbea8d3d9482a9cd23c98ba607e53` | yes (from `H`) |

**Banked files** (SHA-256 of the file bytes, compared in `run.py`) — all six
are this certificate's own:

| file | SHA-256 |
| --- | --- |
| `data/sep1772-decoded-exact-blas.json` | `5985d5f9e1a7ceb54d12dc65e7d5179412eaffb78d94ac3ecf8366db5edbc0d4` |
| `data/sep1772-decoded-exact-bits.json` | `f4b2522d2b8d8ff06ed0195f72051ef28e1f47d7aa691d0afd425a1c84cf98e6` |
| `data/sep1772-twisted-exact-blas.json` | `3b010bf5406916977f060a033990427633bebff0891895b72a0e943afde0f76e` |
| `data/sep1772-twisted-exact-bits.json` | `481ba9f78be1c1d0d2072cfad67598d30c74423ec12e12702350dbba6fb72b35` |
| `data/sep1772-orient-exact-blas.json` | `f58bb4d7db1106950e6506899ffc1329787740cc35b4e0c48574b1d6f06bbea7` |
| `data/sep1772-orient-exact-bits.json` | `57ea8694738d3a6bffe5fa56bf7bab5c6e0d354425b7678f57f6b0d691472f5d` |

`data/payload-records.json` and `data/twisted-i2-records.json` are not
file-pinned here on purpose: they are shared with certs 01 and 02, and the
binding pin on each is the canonical digest of the matrix it produces, which
is checked in clause [1] — reinforced, for the twisted record, by the
`ψ`-twist re-derivation that binds it to `payload-records.json` outright.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 6.2 s** (exit 0, 132 checks; measured 2026-09-03 on the desk) |
| rebuild + full hypothesis re-check + `verify.py`, per matrix | 1.7 s |
| the orientation switch, the cell count and the `S·H″·S` identity at `N = 1772` | < 1 s |
| `dim V` / `dim W` on all three 1772 matrices | < 1 s |
| producing the banked profiles, `blas`, 16 threads | 3 391.6 s (`H`), 3 398.6 s (`H'`), 3 411.7 s (`H″`); peak RSS ≈ 788 MB |
| producing the banked profiles, `bits`, 16 threads | 8 544.6 s (`H`), 8 509.1 s (`H'`), 8 586.5 s (`H″`); peak RSS ≈ 109 MB |
| `run.py --full` at 1772 | **not run in this repository** — see below |

The six banked profiles were produced on the rented `c2d-highcpu-16`
(`prof42-1`), 16 threads; the seconds and the peak resident sets are the
values those runs recorded in the JSONs themselves.

**Why no `--full` leg was run here.** `--full` is offered and wired exactly
as in certs 11, 13, 14, 15, 19, 20, 21 and 22 —
`certs/06-668-separation/full_recompute.py` imported by a `sys.path` insert
rather than copied, BLAS threads capped at three before numpy loads,
smoke-tested against the forced profile of Sylvester `H(128)` first — and
`--matrix` and `--impl` select the leg. It is **priced, not run**. The
campaign's own price at this order is in the table above: **just under an
hour per `blas` leg and about two and a half hours per `bits` leg, on
sixteen rented threads.** This desk has three threads and a different
engine, so the honest desk price is a scaling: one 1772 leg is about **68×**
the 716 leg the same module took in this repository (cert 14, 400.3 s), i.e.
roughly **7–8 hours** for a single `blas` matrix. *Where the 68× comes
from:* it is the source laboratory's **measured** sub-`n⁵` scaling — its
desk-measured 716→2060 ratio of 137 (`experiments/pr0042/REGISTRATION.md`,
Amendment 1) implies the exponent `4.66`, and `(1772/716)^4.66 = 68`. On the
`Θ(n⁵)` law quoted everywhere else here the same leg is `(1772/716)⁵ = 93×`,
i.e. ≈ **10.3 h** — the figure cert 14's own notes already carry for 1772.
Both are estimates and both say *hours*; the smaller, measured one is quoted
so that the decision not to run is not defended with an inflated price. The
`blas` route also materialises a `C(n,2) × n` pair matrix —
`1 569 106 × 1 772`, which is 2.78 GB as `int8` and **11.1 GB** as the
`float32` copy `_profile_blas` makes — past this desk's memory; the `bits`
route (a `1 569 106 × 28` `uint64` packing) is the tractable one at this
order. So this certificate's default verdict is *banked exact computation
audited*, and it says so in its own output. Certs 06, 08, 11, 13, 14, 15 and
19 each have an in-repo `--full` `blas` leg on the record; certs 20, 21, 22
and this one do not, and the word *replayed* is not used of them.

## What is NOT claimed

* The default run recomputes nothing; only `--full` would bind a bank to a
  matrix by computation, and **no `--full` leg has been run in this
  repository at 1772**.
* **Row-side only.** The transposed profiles at 1772 **are not banked
  here** — they are a separate leg of the same campaign, complete in `blas`
  but not in `bits` when this certificate was written, and this repository
  certifies a profile only once `blas` and `bits` agree bin for bin — so
  **nothing here is claimed under the transpose-extended relation**, unlike
  668 (certs 08, 15, 19), 716 (cert 15) and 1676 (cert 21). Refuting
  `A ≈ B` needs both `A ≁ B` and `A ≁ Bᵀ`, and only the first is in hand at
  this order. Nor is anything said about `H` versus `Hᵀ` at 1772: `Hᵀ` was
  never profiled (cert 19 decides that question at 668 and at 668 only).
* **Nothing about order 2092**, and nothing about any order not computed.
  Order 2060's own three-class statement is cert 22's, and is likewise
  row-side.
* **No general statement about `ψ(ρ) = −1`.** That the Lemma-T twist leaves
  the class is now proven at 668, 716, 1676 and 1772. Four orders are four
  orders; the general claim stays unclaimed.
* **No general statement about orientation.** That the switch leaves both
  known classes is now proven at four bordered orders and, unbordered, at
  2060. Same posture.
* No novelty or priority claim of any kind at 1772: this counts classes
  among the artifacts banked here, and says nothing about who first
  exhibited a Hadamard matrix of this order. The decoded record is a public
  artifact (`PROVENANCE.md`).
* **No claim that three is the number of classes** at 1772.
* **Matching invariants prove nothing.** `dim W` agrees on all three, as
  every cheap invariant did at 668, 716 and 1676. Read every "agrees" as
  "did not separate", never as "the same".

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/23-1772-three-classes/run.py
```

The default path is standard library only, no network, and is the whole
certificate. The `--full` paths below import numpy (finder-side only, never
in the trust chain) and cap BLAS threads at three; at this order each is
hours, and the `blas` variants want memory this desk does not have. **None
of them has been run here.**

```
python certs/23-1772-three-classes/run.py --full --impl bits --matrix orient
python certs/23-1772-three-classes/run.py --full
```

Exit code 0 iff every check passed.
