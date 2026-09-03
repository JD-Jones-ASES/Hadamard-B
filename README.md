# Hadamard-B

Bordered Goethals–Seidel arrays: an exact characterisation, a
parameter classification inside the house-Gram branch, certified
instances, and a proven separation at order 668. This repository
states and proves an *if and only if* theorem for extending a
Goethals–Seidel array over a finite abelian group `G` to a Hadamard
matrix through a **coset border** of width `4s` — border strips
constant on the cosets of a subgroup `K ≤ G` of index `i` —
classifies the surviving parameters within that branch under
`w > 2s` (with `n ≥ 3`: the classical `(0,1)` and `(1,1)`, and
`i = s+1` with `s` odd), resolves the `s = 1` border system
completely, and verifies twelve publicly posted matrices through the
theorem's hypotheses. Its third movement proves that order 668
carries at least **four** Hadamard-equivalence classes and orders
716, 1676, 1772 and 2060 at least **three** each, pairwise separated
by an exactly computed invariant; the 668, 716 and 1676 statements
hold with the transpose added to the group, while the 1772 and 2060
statements are row-side only and say so. Every computational
claim carries a certificate in `certs/`; the theorems are
paper-grade proofs in the note, labelled as such. A default run
rebuilds and verifies the
matrices and **audits** the banked exact profiles against the
identities and digests that bind them; the profile **recomputations**
are the `--full` paths. The few measurements made in the source
laboratory rather than replayed here are labelled as such where they
appear.

## The three movements

**I. The theorems** ([note/NOTE-B.md](note/NOTE-B.md) §1).
Theorem A: the bordered array is Hadamard **iff** four conditions hold
— a coset-invariant Gram on the column table, a forced two-tier PAF
profile on the seeds, an orthogonality budget on the corner, and a
coupling of the border to the Goethals–Seidel array of the *coset
sums over `G/K`* (the compression lemma). Theorem C: within the
house-Gram branch, under the non-degeneracy hypothesis `w > 2s`, the
parameter bounds are forced; with `n ≥ 3` the only surviving cells
are the classical `(s,i) = (0,1)` and `(1,1)` and the family
`i = s+1` with `s` odd, where the corner is forced to be a `4s×4s`
Hadamard matrix. Theorem E: the house Gram is not itself forced, but
everything short of that is — every admissible Gram is `4i` times an
orthogonal projector, and at `i = s+1` it is exactly one of the
`|Ĝ̄[2]|` real-character twists of the house form, with no others
(NOTE-B.md §1.2.1; the complete `(3,4)` classification and the
general-branch cells are cert 12). Theorem E′ sharpens the
hypothesis from `w > 2s` to `w > s`, classifies the `w = s`
boundary, and exhibits a seed quadruple at `(s,i,w) = (3,3,3)` that
escapes it — so the sharpened hypothesis is best possible
(NOTE-B.md §1.7; cert 16). Theorem F is the structure theorem for
the **border**: at `i = s+1` under `w > 2s` the corner is forced to
`E = −(1/4i)PĈᵀQ`, (H3) becomes an identity, and (H3)+(H4) reduce
to the single condition `E ∈ {±1}`; a border kit therefore depends
only on the quotient, the reflection class and the S-part of the
coset sums, not on the order. At `(s,i) = (3,4)` the census is
exhaustive — 2048 admissible S-parts per group, every one of which
admits a kit at every `κ(ρ)` — so **that cell is one-layer**: at
order 2092 a seed quadruple *is* the matrix, and the border is
never the obstruction (NOTE-B.md §1.8; cert 17).
Theorem D: at `s = 1` the
`i = 2` border system collapses onto the `i = 1` system (the Gram is
forced, both tables are doubled 4×4 Hadamard matrices, and the
coupling is the same 4×4 equation with the twisted coset-sum vector
in place of the row sums); the degenerate alternative is proved to be
the classical `i = 1` construction rewritten; cert 10 replays the
`768²` corner census of (D-e) and the compressed-block identities
that reduction and (D-d) rest on; with the
character-twist lemma this makes the `i = 2` **seed
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
multiplier subgroup (cert 04); two `H(20)` instances at the
hypothesis boundary `w = 2s` (cert 05); and four matrices in the
**even-`s`** branch of the rigidity theorem — `H(88)` twice and
`H(56)` twice, at the cell `(s,i) = (2,4)` with `Ḡ = ℤ₄`,
`S = {χ, χ³}` and Gram `(8,0,−8,0)`, outside the house form none of
the others leaves (cert 18). None of the four is equivalent to a
Kronecker product with `H(2)`, and none collapses to a smaller
index; at that cell the border layer is proved never to obstruct,
by an exhaustive census of all 215 040 `(S-part, κ(ρ), Q′)` classes.
No novelty of existence is claimed at any of these orders.

**III. Existence plus separation** (§3). **Order 668 carries at
least four Hadamard-equivalence classes** — the decoded record, the
Lemma-T rebuild constructed here, the Hall-switch matrix of a
public preprint (rebuilt and verified firsthand), and the decoded
record with its twelve off-diagonal core blocks negated (the other
Goethals–Seidel orientation; cert 13, 2026-09-02) — pairwise
separated by the exact 4-profile over all 8 222 179 035 row
4-subsets, computed by independent implementations that agree bin
for bin and hit the closed-form second moment to the unit (certs 06,
08, 13), and all four classes stay apart when the **transpose** is
added to the group (cert 15, 2026-09-02). The decoded record is
moreover **not equivalent to its own transpose** — 49 of the 80 bins
differ — so those four matrices and their four transposes are
pairwise inequivalent, all 28 comparisons separating: **eight**
classes are exhibited at 668 under plain Hadamard equivalence, while
the transpose-extended count stays **four**, a matrix and its
transpose being one class there by definition (cert 19, 2026-09-02).
As far as this laboratory's search located, the preprint was the first published
statement of an inequivalence at this order; the bounded priority
statement is NOTE-B.md §3.4. **Order 716 carries at least three** —
the decoded record, its Lemma-T rebuild, and the decoded record with
its twelve off-diagonal core blocks negated (the other
Goethals–Seidel orientation; cert 14, 2026-09-02) — pairwise
separated in 27, 27 and 25 of the 87 bins of the same exact
invariant over all 10 859 143 295 row 4-subsets (certs 11, 14), so
the twist at `ψ(ρ) = −1` provably leaves the equivalence class at a
second order, and the orientation switch leaves both classes at a
second order too; the three transposed profiles have since been
computed, so the three-class statement holds under the
transpose-extended relation as well (cert 15). **Order 1676 carries
at least three** as well — the same three constructions, the decoded
record, its Lemma-T rebuild and the orientation switch — pairwise
separated in 68, 70 and 66 of the 142 bins of the same exact
invariant over all 327 588 749 775 row 4-subsets, in two arithmetics
that agree bin for bin (cert 20, 2026-09-02); so the twist at
`ψ(ρ) = −1` provably leaves the equivalence class at a **third**
order, and the orientation switch leaves both classes at a third
order too. The transposes of the rebuild and of the orientation
switch have since been computed in both arithmetics, each differing
from every original in 139 of 144 bins, so the three-class statement
at 1676 holds under the **transpose-extended** relation as well
(cert 21, 2026-09-02); `Hᵀ` at 1676 was not computed and nothing is
claimed about it. **Order 1772 carries at least three** — the
decoded `(1,1)` record, its Lemma-T `i = 2` rebuild and the
orientation switch, separated in 57, 58 and 53 of the 89 bins of the
same exact invariant over all 409 422 905 815 row 4-subsets, in two
arithmetics that agree bin for bin (cert 23, 2026-09-03); so the
twist leaves the class at a **fourth** order — every decoded `(1,1)`
order in this repository — and the orientation switch leaves both
classes there too. That statement is **row-side only**: the
transposed 1772 profiles are a pending leg of the same campaign.
**Order 2060 carries at
least three** — the publicly posted matrix, the plain GS-array
realisation of the same decoded seed, separated in 146 of the 147
bins of the same exact invariant over all 748 155 697 135 row
4-subsets, with two independent arithmetics agreeing bin for bin
(cert 07, exact mode, banked 2026-09-02), and that realisation in the
**other Goethals–Seidel orientation** — its twelve off-diagonal
515-blocks negated, unbordered here since 2060 sits at `s = 0` —
which differs from the plain array in 107 of the 145 bins the two
share and from the posted matrix in 146 of 147 (cert 22,
2026-09-03). So the orientation switch is a class of its own at a
**fifth** order. The 2060 statement is **row-side only**, as the
1772 one is, and those two are the only ones here that are: at each,
the transposed profiles are pending legs of the same campaign, and
nothing at either order is claimed under the transpose-extended
relation. The sampled
statistic that first suggested the pair stays in the record at its own
label, computational evidence.

## Priority posture

The twelve parameter records were decoded from public seed data
(the sign-stream posted 2026-08-12; expanded matrices for ten of
the twelve in a third-party GitHub repository created and pushed
2026-08-12 UTC, per its own repository metadata — see
[PROVENANCE.md](PROVENANCE.md)). **No priority claim of any kind is
made on the records themselves, on the decode, or on existence at
those twelve orders**; the decode is provenance, and this repository
leads with the theorems the instances instantiate. The `s ≤ 1` layers of the construction are classical
(Goethals–Seidel 1970; Wallis–Whiteman 1972; Spence 1975 — credits in
NOTE-B.md §4). The repository's single hedged novelty statement, with
the exact list of sources it is bounded by, is NOTE-B.md §4.

## Replay

Everything runs from the repository root on bare Python 3.9 or newer.
Standard library only, no network. The eleven `--full` flags of certs 06,
08, 11, 13, 14, 15, 19, 20, 21, 22 and 23 are the only paths anywhere in the
repository that use numpy, and they recompute the exact 4-profiles from the matrices
rebuilt in the run rather than auditing the banks (finder-side only,
never in the trust chain). Cert 16's `--wide` and cert 17's `--full`
are wider runs of the same standard-library code, not a different
arithmetic; cert 17's `--full` (the 16 384-class kit census) **has** been
run in this repository — 2026-09-02, 5 min 28 s, 16 384 / 16 384, census
digest matching the pin — and its `NOTES.md` records the run.
The in-repo `--full` runs so far — certs 08, 06, 11, 13, 14, 15
and 19's — covered the `blas` (float32) leg only; the `bits` recomputation
has not been run in this repository for any of them. Cert 14's leg has been
run here once (400 s), cert 15's once (287 s) and cert 19's once (282 s),
each matching both banked implementations bin for bin. **Certs 20's, 21's,
22's and 23's `--full` have not been run at all**: at order 1676 one leg is of order 6–7 hours
(52× cert 14's 716 leg on the source laboratory's measured sub-`n⁵`
scaling; 70× and ≈ 7.8 h on the `Θ(n⁵)` law used elsewhere here) and
the `blas` route wants about 9.4 GB; at order 1772 one leg is of order 7–8
hours (68× that same 716 leg on the measured scaling; 93× and ≈ 10.3 h on
the `Θ(n⁵)` law) with a `blas` route wanting about 11.1 GB; and at order
2060 one leg is of order
15 hours (137× that same 716 leg on the measured scaling; 197× and ≈ 22 h on
the `Θ(n⁵)` law) with a `blas` route wanting about 17.5 GB — so the flag is
offered and priced and those certificates' verdicts are audits. Each cert's `NOTES.md` records
which legs were run and when.

```
python verify/verify.py --selftest
python certs/01-payload-theorem/run.py
python certs/02-twisted-i2/run.py
python certs/03-h52-gate/run.py
python certs/04-h76-nonscalar/run.py
python certs/05-h20-boundary/run.py
python certs/06-668-separation/run.py
python certs/07-2060-evidence/run.py
python certs/07-2060-evidence/run.py --selftest   # exact-mode acceptance layer
python certs/08-hall-switch-three-classes/run.py
python certs/09-1916-conjugation/run.py
python certs/10-theorem-d-census/run.py
python certs/11-716-separation/run.py
python certs/12-gram-rigidity/run.py
python certs/13-668-orientation/run.py
python certs/14-716-orientation/run.py
python certs/15-transpose-extended-668-716/run.py
python certs/16-theorem-eprime-boundary/run.py
python certs/16-theorem-eprime-boundary/run.py --wide
python certs/17-border-kits-34/run.py
python certs/17-border-kits-34/run.py --full
python certs/18-cell24-instances/run.py
python certs/19-668-transpose-eight-classes/run.py
python certs/20-1676-three-classes/run.py
python certs/21-transpose-extended-1676/run.py
python certs/22-2060-three-classes/run.py
python certs/23-1772-three-classes/run.py
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
authorship is [DISCLOSURE.md](DISCLOSURE.md). The entire repository
— code, note, and prose documentation — is MIT ([LICENSE](LICENSE));
the banked records are mathematical data over which no license is
claimed.

Companion repositories:
[Hadamard-T](https://github.com/JD-Jones-ASES/Hadamard-T) (T-matrix
witnesses and the Hadamard orders they close),
[Hadamard-M](https://github.com/JD-Jones-ASES/Hadamard-M) (the
Miyamoto erratum at order 515 and `H(7796)`),
[Hadamard-formal](https://github.com/JD-Jones-ASES/Hadamard-formal)
(Lean 4 / Mathlib formalizations; PALOMAR-2026-08-31-000001), and
[Hadamard-B-Formal](https://github.com/JD-Jones-ASES/Hadamard-B-Formal)
(the machine-checked tranche of this repository's §1 —
Theorem A, the `s ≤ 1` layer, Theorem D; PALOMAR-2026-09-01-000006).
