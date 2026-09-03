# cert 13 — order 668 carries at least four Hadamard equivalence classes

**Label: PROVEN.** Default run: `python certs/13-668-orientation/run.py` from
the repository root. Standard library only, under a second, exit 0. That run
**audits a banked exact computation**; `--full` recomputes the new profile
here. The trust boundary is the one certs 06, 08 and 11 draw, and is set out
below.

---

## The theorem

> **Theorem.** Let `H` be the decoded `(s, i) = (1, 1)` bordered
> Goethals–Seidel record at order 668 (`data/payload-records.json`; certs 01,
> 06), and let `H″` be `H` with its **twelve off-diagonal core blocks
> negated** and the `4×4` border, the four row strips and the four column
> strips unchanged. Then `H″` is a Hadamard matrix, and it is **not
> Hadamard-equivalent** to any of `H` (decoded), `H'` (the Lemma-T `i = 2`
> rebuild, cert 02/06) or `H★` (the Hall-switch matrix, cert 08).

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all `C(668,4) = 8 222 179 035`
row 4-subsets is a Hadamard-equivalence invariant (`note/NOTE-B.md` §3.1).
All four profiles populate the same 80 bins; `H″`'s differs from `H`'s in
**27** bins, from `H'`'s in **27**, from `H★`'s in **26**. An invariant that
differs is a separation. ∎

**Consequence.** With certs 06 and 08 (`H`, `H'`, `H★` pairwise inequivalent),
**order 668 carries at least four Hadamard equivalence classes.** Row-side
only in this certificate: the transposed profile of `H″` was not computed
here. It *has* since been computed — cert 15, 2026-09-02 — and separates from
`H`, `H'` and `H★` in 50, 50 and 49 of the 80 bins, so the fourth class holds
under the transpose-extended relation too; nothing here rests on that.

## What `H″` is, and what the theorem says about orientation

`note/NOTE-B.md` §1.0 fixes a *standard* Goethals–Seidel orientation and
calls the other one (the six transposed blocks negated) a convention. Put
`S = diag(I₄, diag(1,−1,−1,−1) ⊗ I_n)`. Then `S·H″·S` is exactly the same
seeds and border assembled in the alternate orientation with the border
strips signed by superblock — `P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]` — an
identity of sign patterns that `run.py` checks cell by cell (clause [1]).
So the theorem says: **at 668 the GS orientation is not a gauge for
Hadamard equivalence** — the two orientations of one record are two
classes — and every bordered GS record found at any order carries a second
candidate class for free.

The Lemma-T twist is close by. For `D = diag(I₄, I₄ ⊗ diag(ψ(g)))`,
`ψ(g) = (−1)^g`, the (D-e) transport instance `H_t` of NOTE-B §1.5 on the
twisted seeds satisfies `D·H_t·D = H″` (checked in the source laboratory;
recorded in `intel/fleet-2026-09-01/` there, and re-derived at the desk on
2026-09-02). The banked `H'` (cert 02) is *not* `H_t`: its border is a
Hall-type strip switch of the decoded border. That is why `H″ ≁ H'` is a
theorem and not a contradiction — the twist and the orientation switch
agree on the core and differ on the border, and the border is enough to
change the class.

## The trust boundary — what a default run does and does not establish

The `C(668,4)` enumeration of `H″` **was not run inside this repository.**
It ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py`, numpy, three threads, on
**2026-09-02**, under the pre-registration
`experiments/inequiv/REGISTRATION-668-orientation.md`, **flushed before the
matrix was built**, which fixed the decision rule and the kill criteria in
advance (`blas ≠ bits` in any bin: hard stop, no claim). Its output is banked
in `data/sep668-orient-exact-{blas,bits}.json`. The three comparison profiles
are the banks certs 06 and 08 already audit and pin.

**A default `run.py` establishes:**

* the eight bank files are byte-for-byte the ones pinned in `run.py`;
* `H` rebuilds from the banked record through the full master-theorem
  hypothesis re-check, passes `verify/verify.py`, and carries cert 06's
  canonical digest — computed in-process;
* `H″` is formed from those rows by negating exactly `12·n² = 330 672` cells,
  passes `verify/verify.py`, and carries the pinned canonical digest
  `af1c285c…2953c7`; the alternate-orientation identity holds;
* each bank declares `matrix_canonical_sha256`; for `H` and `H″` it is
  compared against the in-process digest, for `H'` and `H★` against the
  digests certs 06 and 08 pin (those certificates rebuild them);
* every banked profile satisfies the forced identities: bins `≡ 4 (mod 8)`,
  total `C(668,4)`, second moment `n³(n−1)(n−2)/24 = 5 517 193 410 096`;
* `blas == bits` bin for bin on each of the four matrices;
* the three separations of `H″` (27/27/26 bins) and the three prior ones.

**A default run does not establish that the `H″` banks were computed from
the matrix `run.py` rebuilt.** They are *producer-banked*: the digest each
carries is the one the engine recorded against the file it enumerated, and it
equals the digest this certificate pins. `--full` closes the gap: it
recomputes the `H″` profile here, from the rows clause [1] verified, and
compares to both banks bin for bin (see *Runtimes*).

## The evidence chain

**[0]** eight file pins. **[1]** rebuild `H` (hypotheses H0–H4, D1/D3/D5,
Σ̄, compression lemma), verify, pin; form `H″`, count the negated cells,
check `S·H″·S = H_alt`, verify, pin. **[2]** audit eight banks in exact
integers; matrix identity per bank; `blas == bits` ×4. **[3]** the three
separations, differences summing to zero, the divergent bins printed; the
three prior separations re-affirmed from the same banks. **[4]** controls —
C0/C1 five small Hadamard matrices profiled by straight enumeration and by
the pair-vector route `--full` takes (Sylvester `H(8)`, `H(16)` against their
forced profiles; Paley `H(20)`; GS `H(28)`, `H(36)`); the orientation switch
applied to GS `H(28)` is still Hadamard; C2 a total-preserving corruption of
a banked profile that only the second-moment identity can catch, required
to be caught. **[5]** `--full`: `certs/06-668-separation/full_recompute.py`
imported (not copied), smoke-tested on the forced profile of Sylvester
`H(128)`, then run on the verified `H″` rows.

## The separations

Same 80 bins, `≡ 4 (mod 8)`, on all four matrices — the support separates
nothing. The bulk does:

| pair | differing bins | largest `\|Δ\|` | at `\|T4\|` | as a fraction of that bin |
| --- | --- | --- | --- | --- |
| `H″` vs `H` | 27 | 66 814 | 52 | `2.6·10⁻⁴` |
| `H″` vs `H'` | 27 | 53 466 | 4 | `2.6·10⁻⁵` |
| `H″` vs `H★` | 26 | 81 556 | 36 | `1.1·10⁻⁴` |

`run.py` prints the first five divergent bins of each pair; the full lists
are the banked JSONs. As at 668's earlier pairs and at 716, nothing cheaper
than the exact profile could have found this.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports):

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `H` decoded `(1,1)` record | `bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0` | yes (cert 06's pin) |
| `H″` orientation switch | `af1c285cbe2def88427381ab3002a267321b282a9fa78ca37e72830b602953c7` | yes (from `H`) |
| `H'` Lemma-T rebuild | `600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3` | no — cert 06 |
| `H★` Hall switch | `7f6af1d9e9c80fae52c6f178e298144d209c01f6c103db119e6c12425aa1a722` | no — cert 08 |

**Banked files** (SHA-256 of the file bytes):

| file | SHA-256 |
| --- | --- |
| `data/sep668-orient-exact-blas.json` | `5ad283be1baea5c191d39ff3d9219e744ae209ff1df10fea5d152cebbb06c6fd` |
| `data/sep668-orient-exact-bits.json` | `f0eeeb6c89415d09fd8c3f08a13e72666ed0ee36c0066a4d9529e4a88745e726` |
| `data/sep668-exact-blas-decoded.json` | `370fffe6c2f5dc53c09d3b74f8c09dd2bc2a39a1ac2b27fb5167ab4d3559387b` |
| `data/sep668-exact-bits-decoded.json` | `7bace61441f17b5e95fff433bdc5939da212e2b8735e8738d7ed3078fae456b7` |
| `data/sep668-exact-blas-twisted.json` | `8526b3cfa7938a9af334e23f722b1c215ffd1e318c0c713ecc3da1b91f5b3afe` |
| `data/sep668-exact-bits-twisted.json` | `f40bbb8c3906d6fc7374e3e04c2b68eaf29393e50b5662f69eee2426ed3f1e9a` |
| `data/sep668-hall-exact-blas.json` | `35e716ecb43bb6190d5dd6f4160e0bc2bed4f61a3aacf07a36ff9d190810c154` |
| `data/sep668-hall-exact-bits.json` | `a6f703b499d98995f6446a1aed671284c47e99cfe869f3ce8dc8b5fd9394accb` |

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **≈ 0.7 s** (exit 0) |
| producing the banked `H″` profiles upstream | `blas` 82.3 s (peak 411 MB), `bits` 148.9 s (peak 71 MB), 3 threads |
| **`run.py --full --impl blas`** | **285.1 s** (exit 0; measured 2026-09-02 ~04:40 UTC): the fresh `blas` profile of `H″`, computed here from the verified rows by cert 06's `full_recompute.py` (the `U Uᵀ`-triangle route, a third arithmetic), matched **both** banked implementations bin for bin — 80 bins, 284 s of enumeration after the Sylvester `H(128)` smoke test |
| `run.py --full` (both paths) | `bits` leg not yet run in this repository (upstream 148.9 s) |

## What is NOT claimed

* The default run recomputes nothing; only `--full` binds the `H″` bank to
  the matrix by computation.
* Nothing under the transpose-extended relation for `H″`.
* Nothing about the orientation switch at 1772 — the same construction
  exists there and the same computation costs about 130× this one. (It has
  since been made: cert 23, 2026-09-03, row-side only. So have the 716
  instance, ≈1.4× this one, and the 1676 instance, ≈98× this one — cert 14
  and cert 20, both 2026-09-02. Nothing in this certificate depends on any
  of them.)
* No claim that four is the number of classes at 668, and no priority claim
  of any kind: order 668 is settled by the publicly posted matrix; this
  certificate counts classes among artifacts banked here.
* Matching invariants prove nothing; every "agrees" is "did not separate".

## How to re-run

```
python verify/verify.py --selftest
python certs/13-668-orientation/run.py
python certs/13-668-orientation/run.py --full --impl blas
python certs/13-668-orientation/run.py --full
```
