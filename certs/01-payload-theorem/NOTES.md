# cert 01 — the twelve public records against the theorem

## The claim

> The twelve parameter records banked in `data/payload-records.json`
> satisfy **every hypothesis** of Theorem A in its house form
> (Theorem B) — H0 (shape), H1 (the two-tier PAF profile), H2 (the
> corner/row-table orthogonality budget), H3 (the coset-invariant
> column-table Gram), H4 (the coupling to the Goethals–Seidel array of
> the coset sums), these four numbered as `tools/bordered_gs.py`
> numbers them (see its `LABEL MAPPING` block; the note's (H1)–(H4)
> are a permutation of the same four conditions) — and **every
> derived law** of Theorem C — D1, D2,
> D3, D5, the Σ̄ law, and the compression-lemma cross-check — and they
> **assemble to Hadamard matrices whose canonical SHA-256 digests are
> the twelve pinned below**. The twelve land on exactly the cells
> `(0,1)×5`, `(1,1)×4`, `(3,4)`, `(5,6)`, `(7,8)` of the Theorem C
> classification.

## Honesty labels

| part of the claim | label |
| --- | --- |
| the records satisfy H0–H4 and the derived laws; the cell census | **MEASURED** (a machine run on the stated inputs, exact integer arithmetic) |
| the assembled matrices are Hadamard and match the pinned digests | **PROVEN-BY-CERTIFICATE** (explicit artifacts plus a green run of `verify/verify.py`) |

## Priority

**No priority claim of any kind is made here, on the records or on
anything derived from them.** The twelve parameter records were
*decoded* from Hadamard matrices that were posted publicly by another
team; they are that team's mathematical content. The dated provenance
chain is [`../../PROVENANCE.md`](../../PROVENANCE.md); the record bank
itself carries a `provenance` block per order pointing into the decoded
tape. What this cert adds is the verification, and nothing else. No
novelty of existence is claimed at any of the twelve orders — every one
of them is long settled in the literature.

## Data provenance

- **Records:** `data/payload-records.json` (12 records, key `orders`),
  decoded from the publicly posted matrices; see
  [`../../PROVENANCE.md`](../../PROVENANCE.md) and the per-record
  `provenance` field (which names the decoded tape byte range and the
  source of the border tables).
- **Pinned digests:** the twelve canonical SHA-256 values below are the
  digests carried, in agreement, by **two independent lab runs** — the
  bordered-GS theorem check (`experiments/bordered_gs_theorem/`
  `theorem_check_report.json`, second implementation) and the earlier
  decode replay (`experiments/pr0023/replay_report.json`, first
  implementation), both in the source lab repository. The two reports
  were compared digest-for-digest before pinning; they agree on all
  twelve.
- **Checker:** `tools/bordered_gs.py`, the lab's independent second
  implementation of the theorem check, imported unmodified.

## How to re-run

From the repository root, on bare Python 3.9 or newer, stdlib only:

```
python verify/verify.py --selftest
python certs/01-payload-theorem/run.py
```

Runtime ≈ 12 s. The twelve matrices are assembled into a temporary
directory outside the repository, handed to `verify/verify.py` one at a
time, digest-compared, and deleted. Nothing is written inside the
repository and no matrix is committed.

Exit code 0 iff all twelve passed every check and every digest matched.

## Pinned digests

Canonical SHA-256 of the assembled matrix, as reported by
`verify/verify.py` (the digest of the `+`/`−` serialization, one row per
line, newline-terminated).

| order | (s, i) | w | canonical SHA-256 |
| --- | --- | --- | --- |
| 668 | (1, 1) | 166 | `bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0` |
| 716 | (1, 1) | 178 | `3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6` |
| 892 | (0, 1) | 223 | `e77fc79ab287f5f5ba5bbdc10191bdc7593839052fe1015c1fb6a2e974ab54de` |
| 1132 | (0, 1) | 283 | `7d1c1e892149e90330d58bb0cf9ef2c888078df1b35fb55f8724d580ebf7b743` |
| 1244 | (0, 1) | 311 | `4cb747cf511eba1f203582b5121bdf6ab02671133e45579c1d023add8b2da143` |
| 1388 | (5, 6) | 57 | `a6b92584eb803b87026709d64fe892dec8f7182a120e13de9edd3065cf05bf0b` |
| 1436 | (7, 8) | 44 | `e4d745a4d44f39a5671f9cd86f5c1d0aef93504dcfb2e253451cadc9e4086728` |
| 1676 | (1, 1) | 418 | `8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99` |
| 1772 | (1, 1) | 442 | `1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2` |
| 1916 | (3, 4) | 119 | `be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9` |
| 1948 | (0, 1) | 487 | `fddc841ebf951f6e17e939551d058ea5df046251ea065b5f6e7ee2fd8d0f62ce` |
| 1964 | (0, 1) | 491 | `740b907cd442f1b7fd40dcc31f2b3aae9794842da6dc579f98dac1d0d9e1493d` |

The comparison against this table is **coded** into `run.py` (the
`PINNED` dictionary); a mismatch is a hard failure with a non-zero exit
code, not a warning.

## What this cert does not say

- It does not claim any of these matrices, records, seeds or borders as
  this repository's own work.
- It does not claim novelty of existence at any of the twelve orders.
- It says nothing about Hadamard equivalence among these matrices. The
  one equivalence statement this repository proves is at order 668 and
  lives in its own cert.
