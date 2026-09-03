# cert 02 — the four twisted `i = 2` instances (Lemma T)

## The claim

> At each of the orders 668, 716, 1676, 1772 the ψ-twist of the decoded
> `s = 1, i = 1` record is a valid **`s = 1, i = 2`** bordered
> Goethals–Seidel instance: the twisted seeds satisfy the two-tier PAF
> profile with `K` the index-2 subgroup, the banked border `(E, P, Q)`
> satisfies H2, H3 and H4 (tools numbering; see
> `tools/bordered_gs.py`'s `LABEL MAPPING`), the border has exactly
> the shape Theorem D
> forces (`Q[2I+1] = −Q[2I]`, `P[r][2J+1] = −P[r][2J]`, `E` and `p`
> Hadamard), and the assembly is a Hadamard matrix whose canonical
> SHA-256 is the digest pinned below. **Each is a different artifact
> from the decoded record at the same order** — different file,
> different digest.

Here `ψ` is the character of `G = ℤ_n` whose kernel is the index-2
subgroup, i.e. `ψ(g) = (−1)^g`, and `κ(g) = g mod 2`
(`coset_divisors = [2]`). This instantiates **Lemma T** (NOTE-B.md
§1.4): `PAF_{ψx}(t) = ψ(t)·PAF_x(t)`, so the `s=1, i=1` profile
(`Σ PAF(t) = −4` for `t ≠ 0`) twists into the `s=1, i=2` profile (`−4`
on `K∖0`, `+4` off `K`).

At all four orders `ρ` is odd, so **`ψ(ρ) = −1`** — the case in which
§1.4's conjugation proposition does *not* apply and the twist is not a
diagonal conjugation of the decoded instance.

## Honesty labels

| part of the claim | label |
| --- | --- |
| the twisted seeds satisfy the `i = 2` profile; the border satisfies H2/H3/H4 and Theorem D's forced shape | **MEASURED** (exact integer arithmetic) |
| the four assembled matrices are Hadamard and match the pinned digests | **PROVEN-BY-CERTIFICATE** |
| the twisted artifact differs from the decoded artifact at the same order | **PROVEN-BY-CERTIFICATE** (distinct digests — a statement about *files*) |

## What is and is not claimed about equivalence

- At **668**, the twisted matrix and the decoded matrix are proven
  Hadamard-**inequivalent**. That is a separate result with its own
  certificate (the exact 4-profile separation, NOTE-B.md §3.4); *this*
  cert does not prove it and does not depend on it.
- At **716**, likewise: the twisted and decoded matrices are proven
  Hadamard-**inequivalent** — the exact 4-profile over all
  `C(716,4)` row 4-subsets, 27 of 87 bins differing (cert 11,
  NOTE-B.md §3.6). Again a separate result with its own certificate;
  *this* cert does not prove it and does not depend on it.
- At **1676**, likewise: the twisted and decoded matrices are proven
  Hadamard-**inequivalent** — the exact 4-profile over all
  `C(1676,4)` row 4-subsets, 68 of 142 bins differing (cert 20,
  2026-09-02, NOTE-B.md §3.7) and, since cert 21, under the
  transpose-extended relation too. A separate result with
  its own certificate; *this* cert does not prove it and does not
  depend on it, even though cert 20 rebuilds from this record.
- At **1772**, likewise: the twisted and decoded matrices are proven
  Hadamard-**inequivalent** — the exact 4-profile over all
  `C(1772,4)` row 4-subsets, 57 of 89 bins differing (cert 23,
  2026-09-03, NOTE-B.md §3.8), row-side only, the transposed 1772
  profiles being a pending leg. A separate result with its own
  certificate; *this* cert does not prove it and does not depend on
  it, even though cert 23 rebuilds from this record. With that, all
  four orders in this table are decided.
- No novelty of existence is claimed at any of these four orders. They
  are long settled; the artifacts exist to instantiate Lemma T.

## Data provenance

- `data/twisted-i2-records.json` — one complete record per order, in the
  same format as `data/payload-records.json`.
  - **Seeds:** derived deterministically as `x'_q = ψ·x_q` from the
    decoded seeds in `data/payload-records.json`, and spelled out in
    full in the JSON. `run.py` **re-derives them from
    `data/payload-records.json` and compares character for character**,
    so the bank cannot drift from its stated derivation.
  - **Group, reflection shift:** copied from the decoded record at the
    same order (also re-checked by `run.py`).
  - **Border `(E, P, Q)`:** taken from the source lab repository's
    `experiments/bordered_gs_theorem/twist_report.json`, which found it
    by exhaustive enumeration of the width-4 `i = 2` border systems over
    the twisted coset sums. `col_table[r][k] = Q[k][r]`, matching the
    reconstruction `Q[iI+c][r] = col_table[r][iI+c]` in
    `tools/bordered_gs.py`.
  - The decoded records themselves are the posting team's mathematical
    content; see [`../../PROVENANCE.md`](../../PROVENANCE.md) and cert
    01. **No priority claim of any kind is made on them or on anything
    derived from them.**
- **Checker:** `tools/bordered_gs.py`, imported unmodified.

## How to re-run

From the repository root, on bare Python 3.9 or newer, stdlib only:

```
python verify/verify.py --selftest
python certs/02-twisted-i2/run.py
```

Runtime ≈ 4 s. Matrices are assembled into a temporary directory outside
the repository, verified one at a time, digest-compared, and deleted.

Exit code 0 iff all four passed every check and every digest matched.

## Pinned digests

| order | `n` | `i` | `w` | `ψ(ρ)` | canonical SHA-256 of the twisted `i = 2` matrix |
| --- | --- | --- | --- | --- | --- |
| 668 | 166 | 2 | 83 | −1 | `600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3` |
| 716 | 178 | 2 | 89 | −1 | `6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7` |
| 1676 | 418 | 2 | 209 | −1 | `6a4938371ddbe4ad8bd35f21d7e61dad683b15f8f2ec1c88e88ce579c4907405` |
| 1772 | 442 | 2 | 221 | −1 | `82484769a28ac93201f208ca3256bfd491f8edc8bb5e0309764c4f609a113378` |

For reference, the digests of the **decoded** records at the same orders
(cert 01), which `run.py` asserts are different:

| order | canonical SHA-256 of the decoded matrix |
| --- | --- |
| 668 | `bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0` |
| 716 | `3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6` |
| 1676 | `8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99` |
| 1772 | `1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2` |

Both comparisons are **coded** into `run.py`; a mismatch is a hard
failure with a non-zero exit code.
