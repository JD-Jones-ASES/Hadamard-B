# cert 15 — 716 carries three classes, and 668 four, with the transpose in the group

**Label: PROVEN.** Default run:
`python certs/15-transpose-extended-668-716/run.py` from the repository root.
Standard library only, **≈ 3.1 s**, exit 0, **244 checks**. That run **audits a
banked exact computation**; `--full` recomputes the one new 668 profile here.
The trust boundary is the one certs 06, 08, 11, 13 and 14 draw, and is set out
below.

This certificate adds no new matrix and no new construction. It adds **four
exact 4-profiles** — the transposes of the three order-716 matrices and of the
order-668 orientation switch — and with them the statement every earlier
certificate at these orders explicitly withheld.

---

## The relation, and why each pair needs two refutations

`A ~ B` is Hadamard equivalence: `B = D_r P_r A P_c D_c`. The
**transpose-extended** relation is

> `A ≈ B`  iff  `A ~ B`  or  `A ~ Bᵀ`.

Since `A ~ Bᵀ ⟺ Aᵀ ~ B`, refuting `A ≈ B` takes **two** refutations:
`profile(A) ≠ profile(B)`, and `profile(A) ≠ profile(Bᵀ)` — or, reaching the
same statement from the other side, `profile(Aᵀ) ≠ profile(B)`. The multiset
`{|T4|}` over all `C(n,4)` row 4-subsets is a Hadamard-equivalence invariant
(`note/NOTE-B.md` §3.1, invariant **I5**), and the transpose of a Hadamard
matrix is Hadamard, so `profile(Bᵀ)` is that same invariant computed on that
other matrix. Nothing else is needed, and nothing weaker will do: transpose is
**not** in the equivalence group, so a row-side separation says nothing about
the transpose-extended relation on its own.

## The theorems

> **Theorem (716).** Let `H` be the decoded `(s,i) = (1,1)` bordered
> Goethals–Seidel record at order 716 (certs 01/11), `H'` its Lemma-T `i = 2`
> rebuild (certs 02/11), and `H″` the orientation switch of `H` — the twelve
> off-diagonal core blocks negated, border unchanged (cert 14). Then `H`, `H'`
> and `H″` are pairwise inequivalent **under the transpose-extended
> relation**. **Order 716 therefore carries at least three Hadamard
> equivalence classes with the transpose in the group.**

> **Theorem (668).** Let `H`, `H'`, `H★` be the decoded record, the Lemma-T
> rebuild and the Hall switch (certs 06/08), and `H″` the orientation switch
> (cert 13). Then all six pairs are separated under the transpose-extended
> relation. **Order 668 therefore carries at least four Hadamard equivalence
> classes with the transpose in the group.**

*Proof.* Every leg is an exact `|T4|` 4-profile comparison over all row
4-subsets — `C(716,4) = 10 859 143 295`, `C(668,4) = 8 222 179 035` — in two
arithmetics that agree bin for bin. At 716 all three transposes were computed,
so each of the three pairs gets **both** refutation routes (`A` vs `Bᵀ` and
`Aᵀ` vs `B`) and both separate. At 668 the one new profile is `(H″)ᵀ`; with
cert 13's `H″ ≁ X` it supplies the second refutation for each `X ∈ {H, H',
H★}`, and cert 08 already holds the other three pairs. An invariant that
differs is a separation. ∎

**What is new here and what is re-affirmed.** New: the twelve transpose-related
comparisons at 716, and the six involving `(H″)ᵀ` at 668. Re-affirmed from the
same banks, in the same run: 716's three row-side separations (certs 11, 14)
and 668's nine (certs 06, 08, 13). Cert 08 already had the 668 three-class
theorem in transpose-extended form; the **fourth** class is what this
certificate carries across.

## The trust boundary — what a default run does and does not establish

The four new `C(n,4)` enumerations **were not run inside this repository.**
They ran in the source laboratory — `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py` (the same engine as certs 11, 13 and
14, unchanged), numpy, **16 threads on a rented `c2d-highcpu-16`**
(`us-east1-b`), **2026-09-02 ~10:25–10:36 UTC** — under the pre-registration
`experiments/pr0042/REGISTRATION.md`, **flushed 10:17 UTC, before any matrix it
governs was built**. That registration fixes the objects (§2), the decision
rules (§4: "differs in any bin ⟹ inequivalent"; "equal in every bin ⟹
MEASURED, nothing proved") and the kill criteria (§5: `blas ≠ bits` in any bin
is a hard stop) in advance.

The division of labour matters and is stated: **the matrices were built and
verified at the desk** (`experiments/pr0042/build_matrices.py`, every digest in
its `manifest.json`, each through this repository's `verify/verify.py`) and
only then uploaded; the rented machine **enumerated and did nothing else**.

**A default `run.py` establishes:**

* the twenty-five bank files are byte-for-byte the ones pinned in `run.py`;
* at 716, `H` rebuilds from the banked record through the full master-theorem
  hypothesis re-check; `H″` is formed by negating exactly `12·n² = 380 208`
  cells and satisfies the alternate-orientation identity; `H'` rebuilds with
  its seeds **re-derived here** as the `ψ`-twist of the decoded seeds
  (`ψ(g) = (−1)^g` on `Z₁₇₈`, `ρ = 177` odd, so `ψ(ρ) = −1`); all three are
  **transposed in-process**, all six matrices pass `verify/verify.py`, and all
  six canonical digests match the pins;
* at 668 the same for `H`, `H″` (`12·n² = 330 672` cells) and `(H″)ᵀ`;
  `H'`, `H★`, `(H')ᵀ` and `(H★)ᵀ` are certs 06/08's, are **not** rebuilt here,
  and their banks are bound by the digests those certificates pin — each check
  label says which of the two bindings it is using;
* transposition is checked as an involution and, on every object here, as a
  genuine move: `(Mᵀ)ᵀ = M` cell for cell, and `Mᵀ ≠ M` in every row;
* every banked profile satisfies the forced identities — bins `≡ n (mod 8)`,
  total `C(n,4)`, second moment `n³(n−1)(n−2)/24` = `7 807 861 101 040` (716)
  and `5 517 193 410 096` (668) — recomputed **and** compared against the field
  the bank declares where it declares one;
* `blas == bits` bin for bin on each of the **thirteen** matrices, so every leg
  of both theorems, the transpose legs included, rests on two independent
  implementations;
* all thirty pair comparisons, with their counts, in **both** arithmetics; and
  the transpose-extended verdicts derived in code from those counts.

**A default run does not establish that the banked histograms were computed
from the matrices `run.py` rebuilt.** They are *producer-banked*: the digest
each carries is the one the engine recorded against the file it enumerated, and
it equals the digest this certificate pins and, for the nine matrices rebuilt
here, the digest computed in-process in this run. `--full` closes the gap for
the leg the 668 fourth class rests on: it recomputes the `(H″)ᵀ` profile at 668
here, from the rows clause [2] verified, and compares to both banks bin for
bin.

## The evidence chain

**[0]** twenty-five file pins. **[1]** 716: rebuild `H` (hypotheses H0–H4,
D1/D3/D5, Σ̄, compression lemma), form `H″` and check `S·H″·S = H_alt`
cell by cell, re-derive the `ψ`-twist seeds and rebuild `H'`; transpose all
three; six `verify/verify.py` calls, six digests pinned. **[2]** 668: the same
for `H`, `H″`, `(H″)ᵀ`. **[2b]** control C1 — the dim-`V` trap on the real
objects. **[3]** twenty-six banked profiles audited in exact integers; matrix
identity per bank, against an in-process digest where the matrix was rebuilt
here and against certs 06/08's pins otherwise; `blas == bits` ×13; the
populated-bin count of every matrix pinned (87/86 at 716, 80/79 at 668).
**[4]** the thirty pair comparisons, each asserted to the exact bin count, with
the union support size and the first moment printed, and the first six
divergent bins of the two new legs the theorems rest on. **[5]** the
transpose-extended verdicts, derived from those counts, with both refutation
routes shown wherever both transposes are banked. **[6]** controls C2–C5.
**[7]** `--full`: `certs/06-668-separation/full_recompute.py` imported (not
copied), smoke-tested on the forced profile of Sylvester `H(128)`, then run on
the verified `(H″)ᵀ` rows at 668.

**Controls.**

| | control | result |
| --- | --- | --- |
| **C1** | the dim-`V` trap on the real objects | `dim W` = **715** on all six 716 objects — an invariant that separates **none** of the fifteen pairs clause [5] proves inequivalent; `dim V` reads 714/714/714/714/715/715 and is worthless (Trap 1, §3.1). At 668, `dim W = 667` on `H`, `H″`, `(H″)ᵀ` |
| **C2** | the transposed-profile route, on matrices small enough for straight `O(C(n,4))` enumeration | Sylvester `H(8)` and `H(16)` are **symmetric**, so `profile(Mᵀ) = profile(M)` is *forced* — and holds, and both match the forced Sylvester profile; Paley I `H(20)` is not symmetric and is **MEASURED**, never asserted (its profiles happen to agree) |
| **C3** | the transposes are genuinely different objects | every transpose banked here populates **one bin fewer** than its original, and it is the same bin at each order: `\|T4\| = 684` at 716, `\|T4\| = 644` at 668 |
| **C4** | the comparator in the null direction | every banked profile against itself: 0 differing bins, 26 times |
| **C5** | *negative* control: a **total-preserving** corruption of the new `(H″)ᵀ` bank (one count moved from `\|T4\| = 4` to `\|T4\| = 660`), which only the second-moment identity can catch | the assert fires; the corrupted profile still totals `C(668,4)` |

## The separations

### Order 716 — fifteen pairs, `≡ 4 (mod 8)` bins throughout

The three originals populate the same **87** bins; each transpose populates
**86**. Counts as `run.py` prints them (identical in both arithmetics):

| pair | differing | union | source |
| --- | --- | --- | --- |
| `H` vs `H'` | 27 | 87 | cert 11 |
| `H` vs `H″` | 27 | 87 | cert 14 |
| `H'` vs `H″` | 25 | 87 | cert 14 |
| `H` vs `Hᵀ` | **57** | 87 | new |
| `H` vs `(H')ᵀ` | **57** | 87 | new |
| `H` vs `(H″)ᵀ` | **57** | 87 | new |
| `Hᵀ` vs `H'` | **56** | 87 | new |
| `Hᵀ` vs `H″` | **56** | 87 | new |
| `Hᵀ` vs `(H')ᵀ` | **28** | 86 | new |
| `Hᵀ` vs `(H″)ᵀ` | **28** | 86 | new |
| `H'` vs `(H')ᵀ` | **57** | 87 | new |
| `H'` vs `(H″)ᵀ` | **57** | 87 | new |
| `(H')ᵀ` vs `H″` | **57** | 87 | new |
| `(H')ᵀ` vs `(H″)ᵀ` | **26** | 86 | new |
| `H″` vs `(H″)ᵀ` | **57** | 87 | new |

The verdicts, derived from those counts:

| pair | `A` vs `B` | `A` vs `Bᵀ` | `Aᵀ` vs `B` | verdict |
| --- | --- | --- | --- | --- |
| `H` vs `H'` | 27 | 57 | 56 | **SEPARATED** |
| `H` vs `H″` | 27 | 57 | 56 | **SEPARATED** |
| `H'` vs `H″` | 25 | 57 | 57 | **SEPARATED** |

Both routes are shown because at 716 both transposes of every pair are banked;
either column alone would carry the theorem, and the third column is the
redundant crossing that must hold if the second does.

### Order 668 — fifteen pairs

The four originals populate the same **80** bins; each transpose populates
**79**.

| pair | differing | union | source |
| --- | --- | --- | --- |
| `H` vs `H'` | 26 | 80 | cert 06 |
| `H` vs `H★` | 27 | 80 | cert 08 |
| `H'` vs `H★` | 27 | 80 | cert 08 |
| `H` vs `(H')ᵀ` | 50 | 80 | cert 08 |
| `H` vs `(H★)ᵀ` | 49 | 80 | cert 08 |
| `H'` vs `(H★)ᵀ` | 50 | 80 | cert 08 |
| `H` vs `H″` | 27 | 80 | cert 13 |
| `H'` vs `H″` | 27 | 80 | cert 13 |
| `H★` vs `H″` | 26 | 80 | cert 13 |
| `H` vs `(H″)ᵀ` | **50** | 80 | new |
| `H'` vs `(H″)ᵀ` | **50** | 80 | new |
| `H★` vs `(H″)ᵀ` | **49** | 80 | new |
| `H″` vs `(H″)ᵀ` | **50** | 80 | new |
| `(H★)ᵀ` vs `(H″)ᵀ` | **25** | 79 | new |
| `(H')ᵀ` vs `(H″)ᵀ` | **24** | 79 | new |

| pair | `A` vs `B` | `A` vs `Bᵀ` | `Aᵀ` vs `B` | verdict |
| --- | --- | --- | --- | --- |
| `H` vs `H'` | 26 | 50 | n/a | **SEPARATED** |
| `H` vs `H★` | 27 | 49 | n/a | **SEPARATED** |
| `H` vs `H″` | 27 | 50 | n/a | **SEPARATED** |
| `H'` vs `H★` | 27 | 50 | 49 | **SEPARATED** |
| `H'` vs `H″` | 27 | 50 | 50 | **SEPARATED** |
| `H★` vs `H″` | 26 | 49 | 50 | **SEPARATED** |

`Hᵀ` was not computed at 668 and is not needed: every pair there has `A = H` or
a banked `Aᵀ`, so the "n/a" rows are refuted by their middle column alone.

### The two new legs the theorems rest on

`run.py` prints the first six divergent bins of each; the full lists are the
banked JSONs.

**716, `H` vs `(H″)ᵀ`** — 57 of 87 bins differ; largest `|Δ| = 154 128` at
`|T4| = 28`, i.e. `1.06·10⁻⁴` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H` | 2 650 505 561 | 2 383 887 265 | 1 944 278 842 | 1 450 375 604 | 996 815 228 | 634 922 458 |
| `(H″)ᵀ` | 2 650 443 777 | 2 383 856 094 | 1 944 313 602 | 1 450 529 732 | 996 771 316 | 634 926 754 |
| Δ | −61 784 | −31 171 | +34 760 | **+154 128** | −43 912 | +4 296 |

**668, `H` vs `(H″)ᵀ`** — 50 of 80 bins differ; largest `|Δ| = 114 365` at
`|T4| = 4`, i.e. `5.52·10⁻⁵` of that bin.

| `\|T4\|` | 4 | 12 | 20 | 28 | 36 | 44 |
| --- | --- | --- | --- | --- | --- | --- |
| `H` | 2 073 064 058 | 1 852 054 148 | 1 491 070 735 | 1 091 442 371 | 732 009 734 | 452 971 620 |
| `(H″)ᵀ` | 2 072 949 693 | 1 852 022 538 | 1 491 143 856 | 1 091 518 719 | 731 972 800 | 452 991 744 |
| Δ | **−114 365** | −31 610 | +73 121 | +76 348 | −36 934 | +20 124 |

Every difference vector sums to zero, as it must; the first moment, which
nothing forces, does not. As everywhere else in this repository, the largest
discrepancies are of order `10⁻⁴`–`10⁻⁵` of their bins — invisible to any
sample of practical size.

## The six-matrix remark, and what it does not mean

Under **plain** Hadamard equivalence the six order-716 objects
`H, H', H″, Hᵀ, (H')ᵀ, (H″)ᵀ` are **pairwise inequivalent**: all fifteen
profile comparisons separate, the least separated pair by 25 bins. Six classes
are therefore exhibited at 716 by three constructions.

**The house counts three.** The transpose-extended relation is the one under
which a matrix and its transpose are the same object, and three is the count
that survives *either* convention — so it is the number stated in the theorem,
the note and the README. The remark is recorded because it is true and because
it is what makes the transpose-extended statement non-trivial: at 716 a matrix
and its transpose are *never* equivalent among the objects banked here, so the
transpose legs were genuinely open until they were computed. Nothing further is
read into it: no claim that the six are inequivalent to anything else, no claim
about how many classes order 716 has, and no analogue asserted at 668, where
`Hᵀ` was not computed.

## Pinned digests

**Matrices** (canonical SHA-256, the digest `verify/verify.py` reports):

| order | matrix | canonical SHA-256 | rebuilt here |
| --- | --- | --- | --- |
| 716 | `H` decoded `(1,1)` record | `3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6` | yes (certs 01/11's pin) |
| 716 | `Hᵀ` | `e1c4a6fa1f3cc24f0506eddac5bbb94bcbcc5eeab8ef1881c27ad9b8a60be278` | yes (new) |
| 716 | `H'` Lemma-T `i = 2` rebuild | `6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7` | yes (certs 02/11's pin; seeds re-derived) |
| 716 | `(H')ᵀ` | `41fe458af7fe215e59cd98985d4c6835f2364ad13896a0777f260bbccc21ea72` | yes (new) |
| 716 | `H″` orientation switch | `a6b4f56ec98004e736f0ad74af52826aece4b4ab92750e4706e44486c1885fcd` | yes (cert 14's pin) |
| 716 | `(H″)ᵀ` | `7445c760ccaa45d1012828f845acafec6b753596798cba93529bf5f6119de3ef` | yes (new) |
| 668 | `H` decoded `(1,1)` record | `bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0` | yes (cert 06's pin) |
| 668 | `H″` orientation switch | `af1c285cbe2def88427381ab3002a267321b282a9fa78ca37e72830b602953c7` | yes (cert 13's pin) |
| 668 | `(H″)ᵀ` | `49f97ecfb6bdc05c16df3f46aa360202ed29a88e0618b7e6bbbe690f958538d9` | yes (new) |
| 668 | `H'` Lemma-T rebuild | `600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3` | no — cert 06 |
| 668 | `H★` Hall switch | `7f6af1d9e9c80fae52c6f178e298144d209c01f6c103db119e6c12425aa1a722` | no — cert 08 |
| 668 | `(H★)ᵀ` | `565a9ca5a9db739f74215474364202a4d35fd691542b18e8ace8bcbb3c190c65` | no — cert 08 |
| 668 | `(H')ᵀ` | `32afde351e1f44aa4236cb3c406fbbcd59c5ab2cc3e02c8380e08680db2d19d7` | no — cert 08 |

The four transposed digests marked *new* are the desk's own, recorded in
`experiments/pr0042/manifest.json` when the matrices were built and verified
there, and re-derived in-process by every run of this certificate.

**Banked files, this certificate's own** (SHA-256 of the file bytes; all eight
produced 2026-09-02, banked for cert 15):

| file | SHA-256 |
| --- | --- |
| `data/sep716-decoded-T-exact-blas.json` | `b2bd98a8c5a3403273408e009f88030b66a3994d6822be0e9f008db04bc512e5` |
| `data/sep716-decoded-T-exact-bits.json` | `854e73fc748cb4b28113f267fc9266ede5a468105a3c32138978cafc4baaf040` |
| `data/sep716-twisted-T-exact-blas.json` | `135189d8fe4dd619d72ff8ea0ff3cc4b94c07334c1cde4e227f4e08c4cb91fcb` |
| `data/sep716-twisted-T-exact-bits.json` | `2c13fb2e980782b8298e0bf9a37d7bbde569f9a832c5886aa2ec7e4198c1a87b` |
| `data/sep716-orient-T-exact-blas.json` | `b287dfcb5a26e5cf47f8ac8ec445725253cc1a4a51da67eaa219431a8d2da062` |
| `data/sep716-orient-T-exact-bits.json` | `be270c974ce237f537cb39a550819169849805f691120976e71ecfc673b30704` |
| `data/sep668-orient-T-exact-blas.json` | `536e0c136d16b8271a3c7916529fd57c86f6738a8695c1d691f4fbf98455ec2d` |
| `data/sep668-orient-T-exact-bits.json` | `81fae9e0a4f23116067769bffc80451f2c1d834101f9b65f90a2cccf4aab70cf` |

The other seventeen banks are certs 06, 08, 11, 13 and 14's, reused verbatim
and re-pinned by digest in `run.py`: `sep716-exact-{blas,bits}-{decoded,
twisted}.json`, `sep716-orient-exact-{blas,bits}.json`,
`sep668-exact-{blas,bits}-{decoded,twisted}.json`,
`sep668-hall-exact-{blas,bits}.json`, `sep668-hall-T-exact-{blas,bits}.json`,
`sep668-twisted-T-exact.json` (both implementations in one file), and
`sep668-orient-exact-{blas,bits}.json`.

Each of the eight new banks carries a `banked_note` naming the header fields
added at banking time — `schema`, `description`, `matrix`,
`matrix_canonical_sha256`, `producer_filename`, `arithmetic`, `banked_note`.
**Every numeric field is the producer's own output, unaltered.**

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **≈ 3.1 s** (exit 0, 244 checks; measured 2026-09-02, three runs, 3.12/3.14/3.13 s) |
| producing `sep716-decoded-T-exact-{blas,bits}.json` | 63.5 s (peak 508 MB) / 129.5 s (peak 76 MB) |
| producing `sep716-twisted-T-exact-{blas,bits}.json` | 63.9 s (peak 506 MB) / 130.8 s (peak 76 MB) |
| producing `sep716-orient-T-exact-{blas,bits}.json` | 63.5 s (peak 508 MB) / 130.2 s (peak 76 MB) |
| producing `sep668-orient-T-exact-{blas,bits}.json` | 47.8 s (peak 493 MB) / 97.0 s (peak 73 MB) |
| **`run.py --full --impl blas`** | **287 s** (exit 0, 248 checks; measured 2026-09-02 by the independent auditing lane, 283 s of enumeration): the fresh `blas` profile of `(H″)ᵀ` at 668, computed here from the verified rows by cert 06's `full_recompute.py`, matched **both** banked implementations bin for bin — and cert 13's comparable 668 leg took 285 s here |
| `run.py --full` (both paths) | not yet run in this repository |

The eight producer runs are the values those runs recorded in the JSONs
themselves: one rented `c2d-highcpu-16`, 16 threads, `us-east1-b`, 2026-09-02
~10:25–10:36 UTC, under `experiments/pr0042/REGISTRATION.md`. The producers are
numpy and therefore **finder-side only** — never in the trust chain, which is
`verify/verify.py` and stdlib integers.

## What is NOT claimed

* **The default run recomputes nothing.** Only `--full` binds the new 668 bank
  to the matrix by computation, and it has not yet been run here; the seven
  other new profiles are audited, not recomputed, in this repository at all.
* **Nothing at 2060, 1676 or 1772.** The same registration governs orientation
  switches and transposed profiles at those orders; those legs were still
  running when this certificate was written. They will carry their own
  certificates or they will not be claimed.
* **No general theorem.** That the twist, the Hall switch and the orientation
  switch land in different classes — and stay apart with the transpose in the
  group — is now known at two orders. Two orders are two orders.
* **No novelty or priority claim of any kind**, at either order. Order 668's
  "≥ 2 classes" statement is the anonymous preprint's (cert 08, `note/NOTE-B.md`
  §3.4); existence at both orders is the announcing team's. This certificate
  counts classes among the artifacts banked here.
* **No claim that three and four are the counts.** They are lower bounds
  exhibited by the matrices in hand.
* **Matching invariants prove nothing.** `dim W` agrees on all six 716 objects
  that the 4-profile proves pairwise inequivalent, and Paley `H(20)`'s profile
  agrees with its transpose's. Read every "agrees" in this repository as "did
  not separate", never as "the same".

## How to re-run

```
python verify/verify.py --selftest
python certs/15-transpose-extended-668-716/run.py
python certs/15-transpose-extended-668-716/run.py --full --impl blas
python certs/15-transpose-extended-668-716/run.py --full
```
