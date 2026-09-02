# cert 20 — order 1676 carries at least three Hadamard equivalence classes

**Label: PROVEN.** Default run: `python certs/20-1676-three-classes/run.py`
from the repository root. Standard library only, about six seconds, exit 0,
99 checks. That run **audits a banked exact computation**; `--full` is
offered and priced below and **has not been run in this repository at this
order**. The trust boundary is the one certs 11, 13 and 14 draw, and is set
out below. **Row-side only**: the transposed profiles at 1676 are a separate
leg of the same campaign, not yet complete in both arithmetics when this
certificate was written and therefore not banked here, so nothing at all is
claimed under the transpose-extended relation at this order.

---

## The theorems

> **Theorem 1 (the twist, a third instance).** Let `H` be the decoded
> `(s, i) = (1, 1)` bordered Goethals–Seidel record at order 1676
> (`data/payload-records.json`; cert 01) and let `H'` be its Lemma-T
> `i = 2` rebuild (`data/twisted-i2-records.json`, order 1676; cert 02).
> Then `H` and `H'` are **not Hadamard-equivalent**: there is no
> `H' = D_r P_r H P_c D_c` with `P` permutation matrices and `D` diagonal
> `±1`.

> **Theorem 2 (three classes).** Let `H″` be `H` with its **twelve
> off-diagonal core blocks negated** and the `4×4` border, the four row
> strips and the four column strips unchanged. Then `H″` is a Hadamard
> matrix, and it is **not Hadamard-equivalent** to either of `H` and `H'`.
> With Theorem 1, **order 1676 carries at least three Hadamard equivalence
> classes.**

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all
`C(1676,4) = 327 588 749 775` row 4-subsets is a Hadamard-equivalence
invariant (`note/NOTE-B.md` §3.1, invariant I5). All three profiles populate
the same 142 bins; `H` differs from `H'` in **68** of them, `H″` from `H` in
**70**, and `H″` from `H'` in **66**. An invariant that differs is a
separation. ∎

**Consequence.** The Lemma-T construction at `ψ(ρ) = −1` provably leaves the
equivalence class at a **third** order — 668 (cert 06) and 716 (cert 11)
were the first two. Three instances; `note/NOTE-B.md` §1.4's question is
**still not answered in general**, and no general theorem is claimed here.

## What `H″` is, and what the theorem says about orientation

`note/NOTE-B.md` §1.0 fixes a *standard* Goethals–Seidel orientation and
calls the other one (the six transposed blocks negated) a convention. Put
`S = diag(I₄, diag(1,−1,−1,−1) ⊗ I_n)`, `n = 418`. Then `S·H″·S` is exactly
the same seeds and border assembled in the alternate orientation with the
border strips signed by superblock — `P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`
— an identity of sign patterns that `run.py` checks cell by cell (clause
[1]). So the theorem says: **at 1676, as at 668 and 716, the GS orientation
is not a gauge for Hadamard equivalence** — the two orientations of one
record are two classes — and every bordered GS record found at any order
carries a second candidate class for free, by negating twelve blocks.

**Corollary (twist versus orientation, at three orders).** The `ψ(ρ) = −1`
Lemma-T twist and the orientation switch land in **different** classes at
1676 (`H″ ≁ H'`, 66 bins) exactly as they do at 716 (cert 14, 25 bins) and
at 668 (cert 13, 27 bins). **No structural account of why is claimed here**:
at 1676 the statement rests on the profiles alone. Cert 13 records such an
account at 668 — the two moves agree on the core and differ on the border —
from a transport identity checked in the source laboratory; nothing
corresponding was computed or checked at 1676, so nothing corresponding is
asserted. That is now three orders with the same verdict, and **no general
statement is made**: nothing here says the twist and the orientation switch
differ at every order, or at any order not computed.

## The trust boundary — what a default run does and does not establish

The three `C(1676,4)` enumerations **were not run inside this repository.**
They ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py`, the engine unchanged since the
2060 registration, numpy, **16 threads on a rented `c2d-highcpu-16`
(`prof42-1`, `us-east1-b`), 2026-09-02 10:36Z–19:55Z** — under the
pre-registration `experiments/pr0042/REGISTRATION.md`, **flushed 10:17 UTC,
before any matrix it governs was built**, which fixed the decision rules and
the kill criteria in advance (`blas ≠ bits` in any bin: hard stop, no claim;
a builder digest mismatch: nothing uploaded). Its output is banked in
`data/sep1676-{decoded,twisted,orient}-exact-{blas,bits}.json`.

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
  cross-check), passes `verify/verify.py`, and carries cert 01's canonical
  digest `8e919c2b…1cdb99` — computed in-process;
* `H″` is formed from those rows by negating exactly `12·n² = 2 096 688`
  cells, passes `verify/verify.py`, and carries the canonical digest
  `16d1617c…6a84346`; the alternate-orientation identity holds cell by cell;
* `H'` rebuilds too, with its seeds **re-derived here** as the `ψ`-twist of
  the decoded seeds (`ψ(g) = (−1)^g` on `ℤ₄₁₈`; `ρ = 417` is odd, so
  `ψ(ρ) = −1`), compared character for character against the bank, so cert
  02's shared record is bound to `payload-records.json` by computation
  rather than by a file pin;
* each of the six banks declares `matrix_canonical_sha256`, and for **all
  three** matrices it is compared against the in-process digest of the
  matrix rebuilt in this same run;
* every banked profile satisfies the forced identities: 142 bins, all
  `≡ 4 (mod 8)` (`1676 ≡ 4 mod 8`), total
  `C(1676,4) = 327 588 749 775`, second moment
  `n³(n−1)(n−2)/24 = 550 023 273 154 800` — recomputed here and also
  compared against the `second_moment`, `total`, `n` and `C_n_4` fields the
  bank declares;
* `blas == bits` bin for bin on each of the three matrices;
* the three separations — 68, 70 and 66 bins — with their difference vectors
  summing to zero.

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
here more than anywhere else in this repository, since no transposed profile
at 1676 is banked here at all.

Two warnings from §3.1 this certificate honours mechanically: the **signed**
`T4` histogram is not an invariant (every banked profile declares its
folding in its own `folded` field), and **`dim V` is not an invariant** —
only `dim W = V + ⟨1⟩` is. On the 1676 objects `dim V` reads 1674 / 1675 /
1674 and is worthless; `dim W = 1675` on all three. `run.py` measures both
(control C4) precisely so the trap is on the record and never as evidence.

## The evidence chain

**[0]** six file pins. **[1]** rebuild `H` (hypotheses H0–H4, D1/D3/D5, Σ̄,
compression lemma), verify, pin; check the layout `n = |G| = 418`, `s = 1`,
`N = 4(n+s)`; form `H″`, count the negated cells, check `S·H″·S = H_alt`,
verify, pin; re-derive the `ψ`-twist seeds and rebuild `H'`, verify, pin.
**[1b]** control C4 — the dim-`V` trap on the real objects. **[2]** audit
six banks in exact integers; matrix identity per bank against an in-process
digest; `blas == bits` ×3. **[3]** the three separations, each with
identical support, the pinned differing-bin count, the differences summing
to zero, a non-zero first moment, the tail that does *not* separate, and the
first eight divergent bins printed. **[4]** controls — C0/C1 five small
Hadamard matrices profiled by straight enumeration and by the pair-vector
route `--full` takes (Sylvester `H(8)`, `H(16)` against their **forced**
profiles; Paley `H(20)`; GS `H(28)`, `H(36)`, whose three- and four-bin
profiles keep the agreement from being vacuous); the orientation switch
applied to GS `H(28)`, still Hadamard and moving exactly `12·7²` cells; C2 a
total-preserving corruption of a banked profile that only the second-moment
identity can catch, required to be caught; C3 the dim-`V` trap on Sylvester
`H(16)` under a seeded signed row negation. **[5]** `--full`:
`certs/06-668-separation/full_recompute.py` imported (not copied),
smoke-tested on the forced profile of Sylvester `H(128)`, then run on the
verified rows — offered, priced, **not run here**.

## The separations

Same 142 bins, `≡ 4 (mod 8)`, on all three matrices — the support separates
nothing, and neither does the extreme tail: every bin above `|T4| = 564`
agrees on both comparisons with `H`, and every bin above `|T4| = 540` on
`H'` vs `H″`. The **bulk** separates. The first eight bins of each pair, as
`run.py` prints them:

**`H` vs `H'` (the Lemma-T twist)** — 68 of 142 bins differ; largest
`|Δ| = 2 809 810` at `|T4| = 20`, i.e. `5.95·10⁻⁵` of that bin;
`Σ |T4|·Δ = −20 505 840`.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H` | 55 147 198 687 | 52 329 871 149 | 47 215 439 663 | 40 653 198 021 | 33 547 376 838 | 26 657 993 721 | 20 481 183 809 | 15 276 024 237 |
| `H'` | 55 146 545 415 | 52 329 311 849 | 47 218 249 473 | 40 653 475 393 | 33 547 530 716 | 26 657 747 145 | 20 479 582 049 | 15 275 249 189 |
| Δ | −653 272 | −559 300 | **+2 809 810** | +277 372 | +153 878 | −246 576 | −1 601 760 | −775 048 |

**`H` vs `H″` (the orientation switch)** — 70 of 142 bins differ; largest
`|Δ| = 3 095 240` at `|T4| = 20`, i.e. `6.56·10⁻⁵` of that bin;
`Σ |T4|·Δ = +11 699 520`.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H` | 55 147 198 687 | 52 329 871 149 | 47 215 439 663 | 40 653 198 021 | 33 547 376 838 | 26 657 993 721 | 20 481 183 809 | 15 276 024 237 |
| `H″` | 55 145 104 333 | 52 329 711 391 | 47 218 534 903 | 40 653 673 809 | 33 547 436 380 | 26 658 420 825 | 20 479 041 633 | 15 275 675 349 |
| Δ | −2 094 354 | −159 758 | **+3 095 240** | +475 788 | +59 542 | +427 104 | −2 142 176 | −348 888 |

**`H'` vs `H″`** — 66 of 142 bins differ; largest `|Δ| = 1 441 082` at
`|T4| = 4`, i.e. `2.61·10⁻⁵` of that bin; `Σ |T4|·Δ = +32 205 360`.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H'` | 55 146 545 415 | 52 329 311 849 | 47 218 249 473 | 40 653 475 393 | 33 547 530 716 | 26 657 747 145 | 20 479 582 049 | 15 275 249 189 |
| `H″` | 55 145 104 333 | 52 329 711 391 | 47 218 534 903 | 40 653 673 809 | 33 547 436 380 | 26 658 420 825 | 20 479 041 633 | 15 275 675 349 |
| Δ | **−1 441 082** | +399 542 | +285 430 | +198 416 | −94 336 | +673 680 | −540 416 | +426 160 |

`run.py` prints these eight rows of each pair; the full 68-, 70- and 66-bin
lists are the banked JSONs. All three difference vectors sum to zero, and
all three largest discrepancies are of order `10⁻⁵` of their bin — invisible
to any sample of practical size. As at 668 and 716, nothing cheaper than the
exact profile could have found this, and no sampled comparison of these
matrices is banked here because none would have added anything.

## Why the twist and the switch leave the class

`G = ℤ₄₁₈` and `ρ = 417 = −1`, so the only character with `ψ² = 1` has
`ψ(ρ) = −1`. `note/NOTE-B.md` §1.4's proposition — a character twist with
`ψ(ρ) = +1` is a diagonal conjugation `S H S`, and therefore manufactures
nothing — **does not apply**. The `i = 2` border is a genuinely different
border, and the exact profile confirms that the construction really does
produce a new equivalence class here, as it does at 668 and at 716. Whether
it always does at `ψ(ρ) = −1` is still **not claimed**: there are now
**three** proven instances and no general argument.

The same §1.4 remark identifies the orientation switch as the *other* thing
the `ψ(ρ) = −1` conjugation produces — the array with its twelve
off-diagonal blocks negated, which is not of rank-one sign pattern and so is
not a diagonal conjugation either. At 1676 the two moves land in different
classes, as at 668 and 716. Three orders, three verdicts, one pattern, and
**no theorem**.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports) —
all three rebuilt in clause [1] of every run:

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `H` decoded `(1,1)` record | `8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99` | yes (cert 01's pin) |
| `H'` Lemma-T `i = 2` rebuild | `6a4938371ddbe4ad8bd35f21d7e61dad683b15f8f2ec1c88e88ce579c4907405` | yes (cert 02's pin; seeds re-derived) |
| `H″` orientation switch | `16d1617cc62532b26c010f3b174c741f0b9388089516759834030d9056a84346` | yes (from `H`) |

**Banked files** (SHA-256 of the file bytes, compared in `run.py`) — all six
are this certificate's own:

| file | SHA-256 |
| --- | --- |
| `data/sep1676-decoded-exact-blas.json` | `57b9a43caf5246de779ad3205a45642c98f7a211be47e4ed12d718fe098781c9` |
| `data/sep1676-decoded-exact-bits.json` | `469e0b0382d479a6d917316246807cadfa1f113bb2bfcd1429ec1712622e7b94` |
| `data/sep1676-twisted-exact-blas.json` | `328ed05c9614a223d95bd35583c83433a8700249c0b50872fbfc8d846e9b5a49` |
| `data/sep1676-twisted-exact-bits.json` | `311ef88606d0543967e2b0cf46aad4f3fb3f1353cc59b018c5369e447c0c2bb1` |
| `data/sep1676-orient-exact-blas.json` | `a83b239695a3bd820de222e829e65a10a5dd66a432858af57cc950eb4ff40be2` |
| `data/sep1676-orient-exact-bits.json` | `af198e51aecd165e8a2a22ee5ece8dfa73d8ddedf314fa94684ee367db14e9d5` |

`data/payload-records.json` and `data/twisted-i2-records.json` are not
file-pinned here on purpose: they are shared with certs 01 and 02, and the
binding pin on each is the canonical digest of the matrix it produces, which
is checked in clause [1] — reinforced, for the twisted record, by the
`ψ`-twist re-derivation that binds it to `payload-records.json` outright.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path, end to end** | **≈ 5.8 s** (exit 0, 99 checks; measured 2026-09-02 on the desk) |
| rebuild + full hypothesis re-check + `verify.py`, per matrix | 1.4 s |
| the orientation switch, the cell count and the `S·H″·S` identity at `N = 1676` | < 1 s |
| `dim V` / `dim W` on all three 1676 matrices | < 1 s |
| producing the banked profiles, `blas`, 16 threads | 2 638.7 s (`H`), 2 652.4 s (`H'`), 2 649.4 s (`H″`); peak RSS ≈ 762 MB |
| producing the banked profiles, `bits`, 16 threads | 6 628.6 s (`H`), 6 642.9 s (`H'`), 6 631.5 s (`H″`); peak RSS ≈ 104 MB |
| `run.py --full` at 1676 | **not run in this repository** — see below |

The six banked profiles were produced on the rented `c2d-highcpu-16`
(`prof42-1`), 16 threads; the seconds and the peak resident sets are the
values those runs recorded in the JSONs themselves.

**Why no `--full` leg was run here.** `--full` is offered and wired exactly
as in certs 11, 13, 14, 15 and 19 —
`certs/06-668-separation/full_recompute.py` imported by
a `sys.path` insert rather than copied, BLAS threads capped at three before
numpy loads, smoke-tested against the forced profile of Sylvester `H(128)`
first — and `--matrix` and `--impl` select the leg. It is **priced, not
run**: one 1676 leg is about **52×** the 716 leg the same module took in
this repository (cert 14, 400.3 s), i.e. roughly **6–7 hours** for a single
`blas` matrix at three threads on this desk. *Where the 52× comes from,
because it is not this repository's usual law:* it is the source
laboratory's **measured** sub-`n⁵` scaling — its desk-measured 716→2060
ratio of 137 (`experiments/pr0042/REGISTRATION.md`, Amendment 1) implies the
exponent `4.66`, and `(1676/716)^4.66 = 52`. On the `Θ(n⁵)` law quoted
everywhere else here (certs 06, 08, 13 price 668→1676 at `(1676/668)⁵ = 99`)
the same leg is `(1676/716)⁵ = 70×`, i.e. ≈ **7.8 h** — the figure cert 14's
own notes carry for this comparison. Both are estimates and both say
*hours*; the smaller, measured one is quoted so that the decision not to run
is not defended with an inflated price. The `blas` route also
materialises a `C(n,2) × n` pair matrix — `1 403 650 × 1 676`, which is
2.35 GB as `int8` and **9.4 GB** as the `float32` copy `_profile_blas` makes
— past this desk's memory; the `bits` route (a `1 403 650 × 27` `uint64`
packing) is the tractable one at this order. So this certificate's default
verdict is *banked exact computation audited*, and it says so in its own
output. Certs 06, 08, 11, 13, 14, 15 and 19 each have an in-repo `--full`
`blas` leg on the record; this one does not, and the word *replayed* is not
used of it.

## What is NOT claimed

* The default run recomputes nothing; only `--full` would bind a bank to a
  matrix by computation, and **no `--full` leg has been run in this
  repository at 1676**.
* **Row-side only.** The transposed profiles at 1676 **are not banked
  here** — they are a separate leg of the same campaign and had not
  completed in both arithmetics when this certificate was written, and this
  repository certifies a profile only once `blas` and `bits` agree bin for
  bin — so **nothing here is claimed under the transpose-extended
  relation**, unlike 668 (certs 08, 15, 19) and 716 (cert 15). Refuting
  `A ≈ B` needs both `A ≁ B` and `A ≁ Bᵀ`, and only the first is in hand at
  this order.
* **Nothing about orders 1772 or 2060.** The same three constructions exist
  at 1772 and the orientation switch exists at 2060; those legs are separate
  computations under the same registration. Not made here, so not said here.
* **No general statement about `ψ(ρ) = −1`.** That the Lemma-T twist leaves
  the class is now proven at 668, 716 and 1676. Three orders are three
  orders; the general claim stays unclaimed.
* **No general statement about orientation.** That the switch leaves both
  known classes is now proven at three orders too. Same posture.
* No novelty or priority claim of any kind at 1676: this counts classes
  among the artifacts banked here, and says nothing about who first
  exhibited a Hadamard matrix of this order. The decoded record is a public
  artifact (`PROVENANCE.md`).
* **No claim that three is the number of classes** at 1676.
* **Matching invariants prove nothing.** `dim W` agrees on all three, as
  every cheap invariant did at 668 and 716. Read every "agrees" as "did not
  separate", never as "the same".

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/20-1676-three-classes/run.py
```

The default path is standard library only, no network, and is the whole
certificate. The `--full` paths below import numpy (finder-side only, never
in the trust chain) and cap BLAS threads at three; at this order each is
hours, and the `blas` variants want memory this desk does not have. **None
of them has been run here.**

```
python certs/20-1676-three-classes/run.py --full --impl bits --matrix orient
python certs/20-1676-three-classes/run.py --full
```

Exit code 0 iff every check passed.
