# Disclosure

The results are a theorem suite for bordered Goethals–Seidel arrays
over finite abelian groups — an exact (iff) characterisation of the
coset-border construction, the classification of its surviving
parameters, and the complete resolution of the `s = 1` border system —
together with the theorem that order 668 carries at least three
Hadamard-equivalence classes, pairwise separated by an exactly
computed invariant; plus the machine
verification of twelve publicly posted matrices through the theorems'
hypotheses and eight certified instance matrices.

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

Every claim in the note is replayable. Each directory under `certs/`
rebuilds its objects from the small banked data in `data/`,
re-establishes the defining identities with its own exact-integer
loops, hands the resulting matrices to the trust chain, and checks
the canonical SHA-256 in the verdict against the digest pinned in
that certificate. The large matrices are regenerated, not committed.
The default path of every certificate uses nothing outside the
standard library and nothing on the network; cert 06's optional
`--full` recomputation uses numpy on the finder side only.

Reading the certificates is one route. Rebuilding from the
definitions in the note is another, and it requires trusting none of
this. Where a cert's notes cite "the source laboratory" or its
`experiments/` paths, that laboratory's repository is **private**
and is not part of the trust chain: every pin quoted from it is
re-derived here from `data/` at certificate run time, and nothing in
this repository requires access to it. Independent verification is
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
