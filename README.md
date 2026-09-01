# Hadamard-B

Bordered Goethals–Seidel arrays: an exact characterisation, a parameter
classification, certified instances, and a proven separation at order
668. This repository states and proves an *if and only if* theorem for
extending a Goethals–Seidel array over a finite abelian group `G` to a
Hadamard matrix through a **coset border** of width `4s` — border
strips constant on the cosets of a subgroup `K ≤ G` of index `i` —
classifies the surviving parameters (`i = s+1`, `s` odd, beyond the
two classical cells), resolves the `s = 1` border system completely,
and verifies twelve publicly posted matrices through the theorem's
hypotheses. Its third movement proves that order 668 carries at least
**three** Hadamard-equivalence classes, pairwise separated by an
exactly computed invariant. Every claim's certificate replays from
`certs/`; the few measurements made in the source laboratory rather
than replayed here are labelled as such where they appear.

## The three movements

**I. The theorems** ([note/NOTE-B.md](note/NOTE-B.md) §1).
Theorem A: the bordered array is Hadamard **iff** four conditions hold
— a coset-invariant Gram on the column table, a forced two-tier PAF
profile on the seeds, an orthogonality budget on the corner, and a
coupling of the border to the Goethals–Seidel array of the *coset
sums over `G/K`* (the compression lemma). Theorem C: under the
non-degeneracy hypothesis `w > 2s`, the only surviving parameter
cells are the classical `(s,i) = (0,1)` and `(1,1)` and the family
`i = s+1` with `s` odd, where the corner is forced to be a `4s×4s`
Hadamard matrix. Theorem D: at `s = 1` the `i = 2` border system
collapses onto the `i = 1` system (the Gram is forced, both tables
are doubled 4×4 Hadamard matrices, and the coupling is the same 4×4
equation with the twisted coset-sum vector in place of the row sums);
with the character-twist lemma this makes the `i = 2` **seed
problem** a character-twist reparametrization of the `i = 1` seed
problem over any group with a unique index-2 subgroup — a bijection
of search problems, not of matrices (at order 668 the two assembled
matrices are provably inequivalent).

**II. The instances** (§2). The twelve public records at orders 668,
716, 892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948, 1964
satisfy every hypothesis literally and land exactly on the
classification (cert 01, 12/12 green, digests machine-compared
against the banked replay). Constructed here: four matrices in the
twisted `i = 2` frame at orders 668, 716, 1676, 1772 (cert 02); a
from-scratch `H(52)` on a non-cyclic group gating Theorem D
(cert 03); an `H(76)` on a non-cyclic group with a non-scalar
multiplier subgroup (cert 04); and two `H(20)` instances at the
hypothesis boundary `w = 2s` (cert 05). No novelty of existence is
claimed at any of these orders.

**III. Existence plus separation** (§3). **Order 668 carries at
least three Hadamard-equivalence classes** — the decoded record, the
Lemma-T rebuild constructed here, and the Hall-switch matrix of a
public preprint (rebuilt and verified firsthand) — pairwise
separated by the exact 4-profile over all 8 222 179 035 row
4-subsets, computed by independent implementations that agree bin
for bin and hit the closed-form second moment to the unit (cert 06).
The preprint was first to publish an inequivalence at this order;
the priority statement is NOTE-B.md §3.4. The corresponding pair at
the former frontier order 2060 carries **computational evidence of
inequivalence — not a proof** (sampled profiles; an exhausted
block-affine family), stated at exactly that strength.

## Priority posture

The twelve parameter records were decoded from public seed data
(the sign-stream posted 2026-08-12; expanded matrices for ten of
the twelve on GitHub from 2026-08-13 — see
[PROVENANCE.md](PROVENANCE.md)). **No priority claim of any kind is
made on the records themselves, on the decode, or on existence at
those twelve orders**; the decode is provenance, and this repository
leads with the theorems the instances instantiate. The `s ≤ 1` layers of the construction are classical
(Goethals–Seidel 1970; Wallis–Whiteman 1972; Spence 1975 — credits in
NOTE-B.md §4). The repository's single hedged novelty statement, with
the exact list of sources it is bounded by, is NOTE-B.md §4.

## Replay

Everything runs from the repository root on bare Python 3.9 or newer.
Standard library only, no network; cert 06's optional `--full`
recomputation is the one flag that uses numpy (finder-side only,
never in the trust chain).

```
python verify/verify.py --selftest
python certs/01-payload-theorem/run.py
python certs/02-twisted-i2/run.py
python certs/03-h52-gate/run.py
python certs/04-h76-nonscalar/run.py
python certs/05-h20-boundary/run.py
python certs/06-668-separation/run.py
python certs/07-2060-evidence/run.py
python certs/08-hall-switch-three-classes/run.py
python certs/09-1916-conjugation/run.py
```

`verify/verify.py` is the trust chain. It accepts a matrix file only
if the matrix is square, has every entry in {+1, −1}, and satisfies
H·Hᵀ = n·I, by exact integer arithmetic with no floating point
anywhere. Each cert rebuilds its matrices from the banked data in
`data/`, hands them to the trust chain, compares the canonical
SHA-256 against the digest pinned in its `NOTES.md`, and deletes
them. Generated matrices are never committed.

## Layout

| path | contents |
| --- | --- |
| `note/NOTE-B.md` | the mathematics: the theorems with proofs, the instance map, the separation, prior art and credit |
| `verify/verify.py` | the trust chain: Hadamard check, exact arithmetic, stdlib only |
| `tools/` | the bordered-GS assembler/checker shared by the certs |
| `data/` | banked parameter records, border tables, profiles, provenance pins |
| `certs/` | one directory per claim, each with `run.py` and `NOTES.md` |

## Credits, license, provenance

The mathematical credit chain is NOTE-B.md §4; the dated provenance
chain for the public artifacts is [PROVENANCE.md](PROVENANCE.md);
authorship is [DISCLOSURE.md](DISCLOSURE.md). Code is MIT
([LICENSE](LICENSE)); the note and prose documentation are
CC BY-SA 4.0 ([LICENSE-DOCS](LICENSE-DOCS.md)); the banked records
are mathematical data over which no license is claimed.

Companion repositories:
[Hadamard-T](https://github.com/JD-Jones-ASES/Hadamard-T) (T-matrix
witnesses and the Hadamard orders they close),
[Hadamard-M](https://github.com/JD-Jones-ASES/Hadamard-M) (the
Miyamoto erratum at order 515 and `H(7796)`), and
[Hadamard-formal](https://github.com/JD-Jones-ASES/Hadamard-formal)
(Lean 4 / Mathlib formalizations; PALOMAR-2026-08-31-000001).
