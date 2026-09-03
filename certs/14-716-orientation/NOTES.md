# cert 14 — order 716 carries at least three Hadamard equivalence classes

**Label: PROVEN.** Default run: `python certs/14-716-orientation/run.py` from
the repository root. Standard library only, under two seconds, exit 0, 91
checks. That run **audits a banked exact computation**; `--full` recomputes
the new profile here. The trust boundary is the one certs 11 and 13 draw, and
is set out below.

---

## The theorem

> **Theorem.** Let `H` be the decoded `(s, i) = (1, 1)` bordered
> Goethals–Seidel record at order 716 (`data/payload-records.json`; certs 01,
> 11), and let `H″` be `H` with its **twelve off-diagonal core blocks
> negated** and the `4×4` border, the four row strips and the four column
> strips unchanged. Then `H″` is a Hadamard matrix, and it is **not
> Hadamard-equivalent** to either of `H` (decoded) or `H'` (the Lemma-T
> `i = 2` rebuild, cert 02/11).

*Proof.* The multiset `{|T4(i,j,k,l)|}` over all `C(716,4) = 10 859 143 295`
row 4-subsets is a Hadamard-equivalence invariant (`note/NOTE-B.md` §3.1).
All three profiles populate the same 87 bins; `H″`'s differs from `H`'s in
**27** bins and from `H'`'s in **25**. An invariant that differs is a
separation. ∎

**Consequence.** With cert 11 (`H ≁ H'`, 27 of 87 bins), **order 716 carries
at least three Hadamard equivalence classes.** Row-side only *in this
certificate*: no transposed profile at 716 had been computed when it was
written. The three *have* since been computed and banked — cert 15,
2026-09-02 — and the three-class statement holds under the
transpose-extended relation too; nothing here rests on that.

## What `H″` is, and what the theorem says about orientation

`note/NOTE-B.md` §1.0 fixes a *standard* Goethals–Seidel orientation and
calls the other one (the six transposed blocks negated) a convention. Put
`S = diag(I₄, diag(1,−1,−1,−1) ⊗ I_n)`, `n = 178`. Then `S·H″·S` is exactly
the same seeds and border assembled in the alternate orientation with the
border strips signed by superblock — `P[a][J]·(−1)^[J≠0]`, `Q[I]·(−1)^[I≠0]`
— an identity of sign patterns that `run.py` checks cell by cell (clause
[1]). So the theorem says: **at 716, as at 668, the GS orientation is not a
gauge for Hadamard equivalence** — the two orientations of one record are two
classes — and every bordered GS record found at any order carries a second
candidate class for free, by negating twelve blocks.

**Corollary (twist versus orientation, at two orders).** The `ψ(ρ) = −1`
Lemma-T twist and the orientation switch land in **different** classes at 716
(`H″ ≁ H'`, 25 bins) exactly as they do at 668 (cert 13, `H″ ≁ H'`, 27 bins).
**No structural account of why is claimed here**: at 716 the statement rests on
the profiles alone. Cert 13 records such an account at 668 — the two moves
agree on the core and differ on the border — from a transport identity checked
in the source laboratory; nothing corresponding was computed or checked at 716,
so nothing corresponding is asserted. That is now two orders with
the same verdict, and **no general statement is made**: nothing here says the
twist and the orientation switch differ at every order, or at any order not
computed.

## The trust boundary — what a default run does and does not establish

The `C(716,4)` enumeration of `H″` **was not run inside this repository.**
It ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py`, numpy, three threads, on
**2026-09-02**, under the pre-registration
`experiments/inequiv/REGISTRATION-716-orientation.md`, **flushed before the
matrix was built**, which fixed the decision rule and the kill criteria in
advance (`blas ≠ bits` in any bin: hard stop, no claim). Its output is banked
in `data/sep716-orient-exact-{blas,bits}.json`. The two comparison profiles
are the banks cert 11 already audits and pins.

**A default `run.py` establishes:**

* the six bank files are byte-for-byte the ones pinned in `run.py`;
* `H` rebuilds from the banked record through the full master-theorem
  hypothesis re-check, passes `verify/verify.py`, and carries cert 01/11's
  canonical digest — computed in-process;
* `H″` is formed from those rows by negating exactly `12·n² = 380 208` cells,
  passes `verify/verify.py`, and carries the pinned canonical digest
  `a6b4f56e…885fcd`; the alternate-orientation identity holds;
* `H'` rebuilds too, with its seeds **re-derived here** as the `ψ`-twist of
  the decoded seeds (`ψ(g) = (−1)^g` on `Z₁₇₈`; `ρ = 177` is odd, so
  `ψ(ρ) = −1`), so cert 02's shared record is bound to `payload-records.json`
  by computation rather than by a file pin;
* each of the six banks declares `matrix_canonical_sha256`, and for **all
  three** matrices it is compared against the in-process digest of the matrix
  rebuilt in this same run;
* every banked profile satisfies the forced identities: bins `≡ 4 (mod 8)`,
  total `C(716,4) = 10 859 143 295`, second moment
  `n³(n−1)(n−2)/24 = 7 807 861 101 040` — recomputed and also compared
  against the field the bank declares;
* `blas == bits` bin for bin on each of the three matrices;
* the two separations of `H″` (27 and 25 bins) and cert 11's prior one (27).

**A default run does not establish that the `H″` banks were computed from
the matrix `run.py` rebuilt.** They are *producer-banked*: the digest each
carries is the one the engine recorded against the file it enumerated, and it
equals the digest this certificate pins. `--full` closes the gap: it
recomputes the `H″` profile here, from the rows clause [1] verified, and
compares to both banks bin for bin (see *Runtimes*).

## Why `|T4|` is an invariant

For a 4-subset of rows, `T4 = Σ_c H[i][c]H[j][c]H[k][c]H[l][c]`. Row
negation contributes one sign to `T4`, so `|T4|` is fixed; column negation
contributes `d_c⁴ = 1`; row and column permutations relabel the multiset
without changing it. Hence the multiset `{|T4|}` over all `C(n,4)` row
4-subsets is constant on Hadamard-equivalence classes
(`note/NOTE-B.md` §3.1, invariant I5), and two matrices whose profiles differ
in a single bin are inequivalent. Transpose is *not* in the group, which is
why the row-side caveat above is a real restriction and not a formality.

## The evidence chain

**[0]** six file pins. **[1]** rebuild `H` (hypotheses H0–H4, D1/D3/D5, Σ̄,
compression lemma), verify, pin; form `H″`, count the negated cells, check
`S·H″·S = H_alt`, verify, pin; re-derive the `ψ`-twist seeds and rebuild `H'`,
verify, pin. **[1b]** control C4 — the dim-`V` trap on the real objects:
`dim W = 715` on all three (invariant), `dim V` is 714 / 715 / 714 and is
worthless. **[2]** audit six banks in exact integers; matrix identity per
bank against an in-process digest; `blas == bits` ×3. **[3]** the two
separations, differences summing to zero, the first eight divergent bins
printed; cert 11's separation re-affirmed from the same banks. **[4]**
controls — C0/C1 five small Hadamard matrices profiled by straight
enumeration and by the pair-vector route `--full` takes (Sylvester `H(8)`,
`H(16)` against their forced profiles; Paley `H(20)`; GS `H(28)`, `H(36)`);
the orientation switch applied to GS `H(28)` is still Hadamard; C2 a
total-preserving corruption of a banked profile that only the second-moment
identity can catch, required to be caught; C3 the dim-`V` trap on Sylvester
`H(16)` under a seeded signed row negation. **[5]** `--full`:
`certs/06-668-separation/full_recompute.py` imported (not copied),
smoke-tested on the forced profile of Sylvester `H(128)`, then run on the
verified `H″` rows.

## The separations

Same 87 bins, `≡ 4 (mod 8)`, on all three matrices — the support separates
nothing. The bulk does. The first eight bins of each comparison, as `run.py`
prints them:

**`H″` vs `H` (decoded)** — 27 of 87 bins differ; largest `|Δ| = 163 492` at
`|T4| = 36`, i.e. `1.64·10⁻⁴` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H` | 2 650 505 561 | 2 383 887 265 | 1 944 278 842 | 1 450 375 604 | 996 815 228 | 634 922 458 | 377 261 905 | 210 018 619 |
| `H″` | 2 650 471 211 | 2 384 005 045 | 1 944 280 348 | 1 450 401 372 | 996 651 736 | 634 970 434 | 377 270 413 | 210 006 111 |
| Δ | −34 350 | +117 780 | +1 506 | +25 768 | **−163 492** | +47 976 | +8 508 | −12 508 |

**`H″` vs `H'` (Lemma-T twist)** — 25 of 87 bins differ; largest
`|Δ| = 53 806` at `|T4| = 20`, i.e. `2.77·10⁻⁵` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 | 52 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `H'` | 2 650 421 375 | 2 383 982 899 | 1 944 334 154 | 1 450 404 118 | 996 701 290 | 634 954 430 | 377 258 529 | 209 998 559 |
| `H″` | 2 650 471 211 | 2 384 005 045 | 1 944 280 348 | 1 450 401 372 | 996 651 736 | 634 970 434 | 377 270 413 | 210 006 111 |
| Δ | +49 836 | +22 146 | **−53 806** | −2 746 | −49 554 | +16 004 | +11 884 | +7 552 |

`run.py` prints these eight rows of each pair; the full 27- and 25-bin lists
are the banked JSONs. Both difference vectors sum to zero, and both largest
discrepancies are of order `10⁻⁴`–`10⁻⁵` of their bin — invisible to any
sample of practical size. As at 668, nothing cheaper than the exact profile
could have found this.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports) — all
three rebuilt in clause [1] of every run:

| matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- |
| `H` decoded `(1,1)` record | `3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6` | yes (cert 01/11's pin) |
| `H″` orientation switch | `a6b4f56ec98004e736f0ad74af52826aece4b4ab92750e4706e44486c1885fcd` | yes (from `H`) |
| `H'` Lemma-T `i = 2` rebuild | `6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7` | yes (cert 02/11's pin; seeds re-derived) |

**Banked files** (SHA-256 of the file bytes) — the first two are this
certificate's own, the other four are cert 11's, reused verbatim:

| file | SHA-256 |
| --- | --- |
| `data/sep716-orient-exact-blas.json` | `b1c6b0adf393288303f780efebef7ca40bf5d21611660d70928115ff16951cc4` |
| `data/sep716-orient-exact-bits.json` | `bd1c3c23fad9b29a169de1815a4ca44d111b89bb173ef6ac29d55905e077e69c` |
| `data/sep716-exact-blas-decoded.json` | `80ee1e151ec1f759d7213d500623603716b9afa6fc382a385ce6970efac35a6b` |
| `data/sep716-exact-bits-decoded.json` | `a0d5b3a65b83c39c905ec2a1d3b25ca1c58e0106b76aaa6eb3b2feee3748aeed` |
| `data/sep716-exact-blas-twisted.json` | `e2076eb890557775edaccda3c9dcbab7e585f5181e0da1bb5914913bf0749b46` |
| `data/sep716-exact-bits-twisted.json` | `c385773a7a2bf4506b94752406be787bac3d85b8a9ee18bb23668080c6afe7bc` |

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **≈ 1.4 s** (exit 0, 91 checks; measured 2026-09-02) |
| producing the banked `H″` profiles upstream | `blas` 114.7 s (peak 422 MB), `bits` 224.1 s (peak 73 MB), 3 threads |
| **`run.py --full --impl blas`** | **400.3 s** (exit 0; measured 2026-09-02 ~10:35 UTC by the independent auditing lane, 399 s of enumeration): the fresh `blas` profile of `H″`, computed here from the verified rows by cert 06's `full_recompute.py` (the `U Uᵀ`-triangle route, a third arithmetic), matched **both** banked implementations bin for bin |
| `run.py --full` (both paths) | not yet run in this repository (upstream 224.1 s for the `bits` leg) |

## What is NOT claimed

* The default run recomputes nothing; only `--full` binds the `H″` bank to
  the matrix by computation.
* **Row-side only in this certificate.** No transposed profile at 716 had
  been computed when this certificate was written, so **nothing is claimed
  here under the transpose-extended relation** — unlike 668, where cert 08
  makes that statement for the first three classes. The three 716 transposes
  *have* since been computed and banked — cert 15, 2026-09-02 — and the
  three-class statement holds under the transpose-extended relation too;
  nothing in this certificate rests on that.
* Nothing about the orientation switch at 1772 — the same construction
  exists there and the same computation costs about 93× this one (the `n⁵`
  scaling cert 13 uses). It has since been made — cert 23, 2026-09-03,
  row-side only, and it quotes that same 93× price — as has the 1676
  instance, priced here at about 70× on the same law (cert 20,
  2026-09-02); nothing in this certificate rests on either.
* No general statement about `ψ`. That the twist and the orientation switch
  land in different classes is now proved at 668, at 716, at 1676 (cert 20)
  and at 1772 (cert 23); four orders are four orders, and the general claim
  stays unclaimed.
* No novelty or priority claim of any kind at 716: this counts classes among
  the artifacts banked here, and says nothing about who first exhibited a
  Hadamard matrix of this order.
* No claim that three is the number of classes at 716.
* Matching invariants prove nothing; every "agrees" is "did not separate".

## How to re-run

```
python verify/verify.py --selftest
python certs/14-716-orientation/run.py
python certs/14-716-orientation/run.py --full --impl blas
python certs/14-716-orientation/run.py --full
```
