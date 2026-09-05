# Disclosure

The results are a theorem suite for bordered Goethals–Seidel arrays
over finite abelian groups — an exact (iff) characterisation of the
coset-border construction; the classification, within the house-Gram
branch and under `w > 2s`, of its surviving parameters (with `n ≥ 3`
the cells are `(0,1)`, `(1,1)`, and `i = s+1` with `s` odd); and the
complete resolution of the `s = 1` border system —
together with the theorems that order 668 carries at least four
Hadamard-equivalence classes and that orders 716, 1676, 1772 and
2060 each carry at least three, pairwise separated by an exactly
computed invariant — and, since 2026-09-03, every one of those counts
holding with the transpose added to the group, so that no separation
statement here is row-side any longer;
plus the machine verification of twelve publicly posted matrices
through the theorems' hypotheses and twelve certified instance
matrices. Added 2026-09-02: the rigidity theorem sharpened to
`w > s` with its boundary classified; the border-kit structure
theorem, whose exhaustive `(3,4)` census makes that cell one-layer
— at order 2092 a seed quadruple *is* the matrix, the border never
the obstruction; the first four certified instances in the
even-`s` branch, at the cell `(s,i) = (2,4)`, where a second
exhaustive census shows the border never obstructs either; and the
decoded 668 record proved inequivalent to its own transpose (49 of
80 bins), so those four matrices and their four transposes are
pairwise inequivalent — eight classes exhibited at 668 under plain
Hadamard equivalence, the transpose-extended count unchanged at
four, which is the count this disclosure and the abstract quote; and
the order-1676 three-class theorem carried across to the
transpose-extended relation, each of its six transpose-related
comparisons separating in 139 of 144 bins, so that order 2060 was
then the only separation statement here that was row-side only (1772
joined it the next day, below, and both were settled the day after).
Added
2026-09-03: order 2060 raised from two classes to three, by the
orientation switch of the plain Goethals–Seidel realisation — the
twelve off-diagonal `515`-blocks negated, with no border to leave
alone at that order — which differs from the plain array in 107 of
the 145 bins they share and from the publicly posted matrix in 146 of
147; and no
priority or novelty of any kind is claimed at order 2060, which was
settled by the publicly posted matrix. Also 2026-09-03: order 1772
raised from one class to three — the decoded `(1,1)` record, its
Lemma-T `i = 2` rebuild and the orientation switch, separated in 57,
58 and 53 of the 89 bins of the exact 4-profile over all
409 422 905 815 row 4-subsets in both arithmetics — so the
`ψ(ρ) = −1` twist provably leaves the equivalence class at a fourth
order, which is every decoded `(1,1)` order here, and the
orientation switch is a class of its own at a fifth; and no priority
or novelty of any kind is claimed at order 1772, which is long
settled. Both of those statements were row-side when first made and
are not now: also 2026-09-03, the four remaining transposed profiles
were computed in both arithmetics and both orders were carried across
to the transpose-extended relation — at 2060 every pair carrying both
refutations, the two involving the posted matrix differing in *every*
bin of their union (134 of 134), and at 1772 each of the six
transpose-related comparisons separating in 91 of 92 bins. **No
separation statement in this repository is row-side any longer.**
Nothing is claimed at any order about a matrix versus its *own*
transpose beyond 668, where it was decided: `Hᵀ` at 1676 and at 1772,
and `Gᵀ` at 2060, were never profiled.
Added 2026-09-05, in two lanes. **Theorem T**: the coset-border
family is closed under transposition — branch-free, at every cell
and every `w` — with remark R fixing exactly when the border-kept
orientation switch is again Hadamard (**PROVEN** +
**PROVEN-BY-CERTIFICATE**, cert 26). **Sylow-2 rigidity**: the
extension group of a hypothetical *cocyclic* `H(2092)` is cut from
twenty-four `(E, e*)` pairs to two, conditional throughout on the
cited de Launey–Flannery–Horadam equivalence (**PROVEN** +
**PROVEN-BY-CERTIFICATE**, cert 27). **Lemma F** and the records
theorem: an automorphism of order `p = 523` would force a 16-block
circulant array on an `H(2092)`, and no automorphism of odd prime
order acts on `H(668)`, `H′(668)`, `H″(668)` or `H(716)`, whose
automorphism groups are therefore 2-groups (**PROVEN** +
**PROVEN-BY-CERTIFICATE**, cert 29). The `(2,4)` **existence
theorem**: at that cell the border is never the obstruction, at
every `w` — proved by classification where cert 18's census proved
it by enumeration (**PROVEN** + **PROVEN-BY-CERTIFICATE**,
cert 28). **Theorem F in the general branch**, with the
row-factorisation lemma and its dual form, which make a border-side
negative complete in the row table; the first border kill in the
note — the `(7,12)/S3` tier admits no column table at any order
`4(12w+7)`, and its `χ₆` twin dies with it — and a banked
`(7,12)/S2` column table (**PROVEN** + **PROVEN-BY-CERTIFICATE**,
certs 30 and 31). And the **good-matrix product theorem**
`a_{2k} = −(b₀c₀d₀)·a_k b_k c_k d_k` with its closure sign law
(**PROVEN**, with small-order censuses and a Williamson control as
**COMPUTATIONAL-EVIDENCE**, cert 32). Every statement here about
`H(2092)` is **conditional** — a constraint on the shape such a
matrix would have to have, or on one layer of one construction —
and **none asserts existence**, at 2092 or anywhere else: a dead
border tier is a dead border tier, and no good quadruple is known
at `n = 523`.

AI-generated results with a human managing the workflow. Produced by
Claude Code (Fable 5, Anthropic); external reviews at earlier stages
by GPT 5.6 (OpenAI) and by Grok (xAI), intaken and adjudicated.

## What the AI stations did

Everything mathematical. The decode of the publicly posted matrices
into parameter records; the statement and proofs of the theorems; the
adversarial skeptic passes and their adjudication (including the two
documented invariant traps in `note/NOTE-B.md` §3); the exact
4-profile computations and their independent implementations; the
verifier and every certificate; the firsthand source reads that fix
the credits and the novelty statement in `note/NOTE-B.md` §4. Where
an external report and a primary source disagreed, the source
governed.

## What the human owner did

Granted the sessions, paid for the compute, obtained papers this
laboratory could not otherwise read, relayed material to and from the
outside reviewers, and ruled on publication, licensing, and scope.
No mathematical contribution, and none is claimed. The owner's name
appears here and in the copyright line; it appears in no derivation.

## Verification

The trust chain is `verify/verify.py`. It accepts a matrix only if
the matrix is square, has every entry in {+1, −1}, and satisfies
H·Hᵀ = n·I. The arithmetic is exact: rows are packed into integers
and orthogonality of a pair is a popcount identity, so no floating
point enters anywhere. `python verify/verify.py --selftest`
exercises it against known Hadamard matrices, against
Hadamard-preserving row and column operations, and against
corruptions it must reject.

Every computational claim in the note carries a certificate; the
theorems are proved on paper and labelled **PROVEN (paper-grade)**.
Each directory under
`certs/` rebuilds its objects from the small banked data in `data/`,
re-establishes the defining identities with its own exact-integer
loops, hands the resulting matrices to the trust chain, and checks
the canonical SHA-256 in the verdict against the digest pinned in
that certificate. The large matrices are regenerated, not committed.
Where a claim rests on an exact 4-profile too large to recompute on
every run, the default path **audits** the banked profile — the
file digest, the matrix binding where the bank declares one, the
forced congruence, the total, the second moment, and agreement
between two independent implementations — and the word *replay*
belongs to the optional `--full` paths of certs 06, 08, 11, 13, 14,
15, 19, 20, 21, 22, 23, 24 and 25, which recompute those profiles from
the matrices rebuilt in the same run.
The default path of every certificate uses nothing outside the
standard library and nothing on the network; those thirteen `--full`
recomputations use numpy on the finder side only, and they are the
only numpy anywhere in the repository. Certs 20's, 21's, 22's, 23's,
24's and 25's
have not been run: at order 1676 one leg is of order 6–7 hours — 52× cert 14's
716 leg on the source laboratory's measured sub-`n⁵` scaling, ≈ 7.8 h
on the `Θ(n⁵)` law used elsewhere here — and the `blas` route wants
about 9.4 GB; at order 1772 one leg is of order 7–8 hours (68× that
same 716 leg, ≈ 10.3 h on the `Θ(n⁵)` law) and the `blas` route
wants about 11.1 GB; at order 2060 one leg is of order 15 hours (137× that
same 716 leg, ≈ 22 h on the `Θ(n⁵)` law) and the `blas` route wants
about 17.5 GB, so those certificates' verdicts are audits and say so. Cert 17's `--full` is
standard-library like the rest of it, and has been run here: the
16 384-class census reproduced on 2026-09-02 at the pinned digest,
as its `NOTES.md` records.

Reading the certificates is one route. Rebuilding from the
definitions in the note is another, and it requires trusting none of
this. Where a cert's notes cite "the source laboratory" or its
`experiments/` paths, that laboratory's repository is **private**
and is not part of the trust chain: every matrix digest quoted from
it is re-derived here from `data/` at certificate run time, and the
exact 4-profiles it produced are banked in `data/`, audited on the
default path and recomputed by `--full`. Nothing in this repository
requires access to it. Independent verification is
invited; questions and verification reports are welcome via GitHub
issues.

## Credit for external mathematics

The full credit chain, with the exact sources read firsthand and the
one bounded novelty statement, is `note/NOTE-B.md` §4; the dated
provenance chain for the public artifacts is `PROVENANCE.md`. In
brief: the twelve verified matrices are the announcing team's (seed
data posted publicly 2026-08-12); the order-2060 artifact is
Schneider's; the classical spine is Goethals–Seidel (1970),
Wallis–Whiteman (1972), and Spence (1975); the compression device is
Đoković–Kotsireas's. No priority claim of any kind is made on the
public records themselves, on the decode, or on existence at those
orders.
