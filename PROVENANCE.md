# PROVENANCE — the dated chain for everything sourced from outside

Three bodies of external material enter this repository: the twelve
publicly posted Hadamard matrices verified in movement II (the
announcement's sign-stream tape and its decoder, the third-party
GitHub decode, and the relayed zip of twelve CSVs); the publicly
posted matrix of order 2060 used in movement III; and the anonymous
order-668 inequivalence preprint at `hadamard-668.vercel.app`, whose
published Hall-switch data is banked in `data/sep668-hall-switch.json`
and rebuilt by cert 08 — its chain is under *Independent public
corroboration* below. Each chain is given with its dates, its
retrieval, and the digests that pin it. Nothing else here comes from
outside; what is this laboratory's own is stated at the end.

**Dates in this document are UTC**, retrieval dates included; a
timestamp quoted from a source keeps that source's own zone and is
marked where that zone is not UTC. The SHA-256 digest of every file in
`data/` is recorded at *Frozen data digests (2026-09-01)*, the last
section of this document.

## The twelve public records

### The announcement

On **2026-08-12** Levent Alpöge (`@__alpoge__`) posted, on X, seed data
for Hadamard matrices at the twelve orders 668, 716, 892, 1132, 1244,
1388, 1436, 1676, 1772, 1916, 1948 and 1964. The announcement is three
posts and no paper:

| post | what it is |
| --- | --- |
| `https://x.com/__alpoge__/status/2087504785952182273` | the payload. `created_at` **Wed Aug 12 11:41:39 +0000 2026**; body is pure `+`/`-`, 23,828 characters, `display_text_range = [0, 23828]`, no prose at all. Read firsthand through the fxtwitter API mirror on **2026-08-16** and again through `og:description` on **2026-08-31**. |
| `https://x.com/__alpoge__/status/2087504788938510427` | the decoder: an obfuscated POSIX-shell program, 6,130 bytes, that reads the first post as an input tape. The program itself was extracted and read firsthand (SHA-256 `cf0ac77c3c6329d7d904466c18e653229521913ed19f2535322a9c5abc986f05`); the post itself was read firsthand through the fxtwitter API mirror on **2026-09-01** — author `__alpoge__`, `created_at` **Wed Aug 12 11:41:40 +0000 2026**, its text opening with the sed-substitution decoder that matches the extracted script. X's own oEmbed endpoint returned HTTP 402 to this laboratory, so no oEmbed retrieval is claimed; the pinned full-script digest remains the stronger content bind. |
| `https://x.com/__alpoge__/status/2087504790435840207` | the credits, verbatim and entire: *"weekend fun w @tehwalris, Saul Reynolds-Haertle, and of course claude:)) i only claim bad suggestions!!"* (`@tehwalris` = Philippe Voinov). Read firsthand through `og:description` on **2026-08-31**. |

**What was disclosed: data and a decoding grammar. No method.** No
construction is described in any of the three posts, no search is
described, and no preprint accompanies them (see *Independent public
corroboration* below). The public record of the event, as read on
2026-08-16, is secondary: the Wikipedia article *Hadamard matrix*
(edited 2026-08-12) credits Philippe Voinov, Saul Reynolds-Haertle,
Claude and Alpöge and states that no further details of the
construction were provided; OEIS A007299's 2026 history line lists the
same twelve orders; and an early public prose description of what the
second post does — that it expands the tape through a Goethals–Seidel
construction, and which twelve orders it covers — is a third-party
post, `https://x.com/adi_baradwaj/status/2087620459185819903`
(2026-08-12 19:21 UTC). It is *not* claimed to be the earliest public
reading: the third-party GitHub decode below, which already held the
decoder and its expanded outputs, was created and pushed the same day
at 17:28–17:35 UTC.

### The public artifacts, and how they are pinned

| artifact | retrieved | pin |
| --- | --- | --- |
| **the tape** — the 23,828-character sign stream of the payload post | fetched through the fxtwitter API mirror **2026-08-16**, re-fetched and banked **2026-08-28** | normalised sign stream, SHA-256 `5b5fe8fa42f0d6a8b4e4c9926726d82a6aab8e1070c1ae4d1b430c1277e58db4`, length 23,828. This digest is carried in `data/payload-records.json` (`tape_sha256`). |
| **the expanded matrices** — `github.com/foocker/Hadamard668`, a third-party public decode created and pushed **2026-08-12 UTC**, hours after the announcement (repository metadata read through the GitHub API on 2026-09-01: `created_at` `2026-08-12T17:28:39Z`, `pushed_at` `2026-08-12T17:35:13Z`, still exactly two commits, its head `f951552a…` dated `2026-08-12T17:35:01Z`, subject "Add Hadamard matrix, decoder, and visualizations". Those are creation and push times; the repository's own metadata does not establish when it became visible to any particular reader, and no such claim is made) | cloned read-only **2026-08-28**; metadata re-read **2026-09-01** | its `answer.md` normalises to the same tape digest. It held ten complete matrices; its H(1948) stopped at **603 complete rows plus a 305-entry partial row — 1,174,949 of 3,794,704 entries, 30.96 % of that matrix** — and it held no H(1964) at all. Its decoder script was **never executed**; its outputs were used as data only. |
| **a zip of twelve CSVs**, circulated as the announcement matrices | received by relay and verified **2026-08-31** | 4,873,392 bytes, SHA-256 `503a4352f282a3d0141a1814b4f89036961524934670a5603ec33a95ea6f179d`; twelve entries, all `open_hadamard/hadamard_<N>.csv`, internal mtimes all `2026-08-12 05:58:58`. That it was posted by Sumeet Motwani (`@sumeetrm`) on Google Drive as Alpöge's matrices is **REPORTED-FROM-SECONDARY-SOURCE**; no source URL is recorded here. Nothing in this repository rests on that attribution. |

The three copies agree. All twelve matrices of the zip verify green in
the trust chain and are **entry-identical** to the matrices this
laboratory rebuilt from tape bytes — including the parts of H(1948)
and all of H(1964) that are absent from the GitHub copy.

### What this repository redistributes

**`data/payload-records.json`, and no matrices.** The file is the
decoded parameter records: for each of the twelve orders, the abelian
group `G`, the border width `s`, the reflection and shift convention,
the Goethals–Seidel block variant, the four `±1` seed sequences, the
coset divisors, and — where `s ≥ 1` — the corner and the row/column
border tables. Schema `gs-frame-params/1`, generated 2026-08-28, 29,138
bytes, SHA-256
`9afbf60efed82e97ffe229e0059f06dfa7849a2ac2434c52a6ada09559738afb`; it
carries the tape digest and length in its own header.

The seeds are **byte-exact tape extractions**, recorded with their tape
spans. The four `s = 1` border tables are **not on the tape**: they were
measured from the public matrices at orders 668, 716, 1676 and 1772, one
table shared by all four, and were later confirmed against the zip at
all four orders.

The decode was performed here, from the public artifacts, through this
laboratory's own standard-library code; the segmentation grammar is the
announcement's own second post, read firsthand as above. **No
third-party method disclosure was used, and none was located in the
sources enumerated here as of 2026-09-01** — the re-sweep of that date
found no bordered or coset-border Goethals–Seidel characterisation
anywhere public. The
correctness of the decode is not asserted, it is replayed: each cert
rebuilds its matrices from these records and hands them to
`verify/verify.py`. No generated matrix is committed; the certificates
rebuild them and delete them.

### Priority posture

**The twelve records are the announcing team's mathematical content.**
The decode is provenance, not an achievement, and **no priority claim of
any kind is made on the records themselves, on the decode, or on
existence at the twelve orders.** (What this laboratory built *from*
the records — the theorems, the constructed instances, the
separations — is claimed exactly as NOTE-B.md §4.2 states, and no
further.) No novelty of existence is claimed at any of the twelve
orders. That those
orders were open before 2026-08-12 is not proven here either: it is
REPORTED-FROM-AUDITED-TABLE, from Table 4 of Cati–Pasechnik,
arXiv:[2411.18897](https://arxiv.org/abs/2411.18897)v2 (2025-08-30) —
re-checked 2026-09-01 and still at v2, no v3.
What this repository adds at these orders is verification against the
theorems' hypotheses (NOTE-B.md §2.1, cert 01).

## The order-2060 artifact chain (cert 07)

The matrix of order 2060 used in movement III is public and is not this
laboratory's. It comes from the gist

> `https://gist.github.com/schneiderlo/b866a2ff2fcd93934f0db54cfa4069d0`

gist id `b866a2ff2fcd93934f0db54cfa4069d0`, user `schneiderlo`,
description "Hadamard matrix of order 2060", created
**2026-08-23T02:33:38Z**, retrieved **2026-08-29**. No gist revision sha
was recorded at retrieval; the payload digest is what binds the file. Re-read
through the GitHub API on 2026-08-31: unchanged since 2026-08-23T02:33:42Z.
Re-read again **2026-09-01**: the gist history still holds a single revision
(`00784a83…`, committed 2026-08-23T02:33:38Z) and `updated_at` is still
2026-08-23T02:33:42Z, so the bytes are the bytes that were digested and the
payload pin below stands unre-fetched.

| quantity | value |
| --- | --- |
| payload SHA-256 of the gist file `H2060.txt` | `c7a145d86210740dd3f8ea21ca896a54d6916007a042638f17c8c47f097200f7` |
| canonical SHA-256 of the plain Goethals–Seidel array over the same decoded seed | `510f89b7b423c85da0c7cada52cb0f62e0415d736d2040d6701f66c4e524cf6a` |

The gist file is already in `verify/verify.py`'s canonical
serialisation, so its payload digest and its canonical digest coincide.
Cert 07 rebuilds **both** order-2060 matrices from the banked seed — the
posted one (the `×104`-twisted array) and the plain array over the same
seed — verifies each in the trust chain, and matches both digests.

**Announcement, cited separately from the artifact.** The gist is paired
with the post `https://x.com/modkin_mp/status/2091352950039785791`
(Loïc Schneider, `@modkin_mp`, 2026-08-23 02:32 UTC), which states that a
Hadamard matrix of order 2060 exists and shows a 2060 × 2060 pixel image.
Neither the gist nor the post discloses a construction method. The
artifact was located and downloaded directly from the gist URL and
verified firsthand; nothing in this repository depends on the post.

**Credit for the order-2060 matrix is Schneider's: he is the public
source and publisher of the artifact, and his own repository later
gave a producer-side construction for it.** No universal firstness and
no sole mathematical authorship is asserted here.
`github.com/schneiderlo/hadamard-2060` (created
**2026-08-25T23:38:42Z**, last push 2026-08-27T13:40:24Z) gives the
Cooper–Wallis route `T(103) ⊗ W(5) → H(2060)` from the producer's own
side; its explicit matrix file (SHA-256
`8c81090db76e2b503561ba82598b514e0c2e39940a3703783c5b1b2c6b80a37a`)
parses to the canonical digest `c7a145d8…` above — measured 2026-08-31,
the repo and the gist are one object. Re-read **2026-09-01**: `main` is still
the banked commit `9d4fdce…` (last push still 2026-08-27T13:40:24Z) and the
matrix file's git blob sha is identical at `HEAD` and at that commit, so the
`8c81090d…` pin holds. This laboratory decoded that
matrix's structure independently from the banked artifact alone, and was
not first to state it publicly; no firstness is claimed. The plain array
is a genuinely different matrix, but a derivative of the public seed, not
an independent find. The separation of the two is stated at
COMPUTATIONAL-EVIDENCE and nowhere stronger (NOTE-B.md §3.5).

## Independent public corroboration

Dated observations of the public record, first made **2026-08-31** by a
sweep lane of this laboratory and **re-run in full at the public flip on
2026-09-01** (arXiv API author sweep; GitHub API and raw sources; the
Palomar registry's CC0 feed). Each bullet carries its own read dates.

- **No preprint.** Re-run **2026-09-01**, against the arXiv API: author
  sweeps for Alpöge, Voinov and Reynolds-Haertle; `all:"Hadamard matrix"`
  against the orders 668, 1948 and 1964; a `ti:"Hadamard matrix"` sweep;
  and the newest `cat:math.CO` items. The author sweeps return **zero
  Hadamard papers**; the order sweep returns two entries, neither the
  announcement; the title and category sweeps surface nothing on these
  orders after 2026-08-26. **There is still no arXiv preprint of the
  announcement**, and its only public statement is the X posts. This is
  corroborated from a second direction: the provenance files of both
  downstream Lean repositories still **cite no source but the X posts for
  the matrix data**, `type: "web post"`. (Ramos's file additionally lists
  Hadamard 1893 as `uses-as-background`; that is background, not matrix
  data, and no other file cites anything further.)
- **`schneiderlo/hadamard-2060`** (GitHub, created 2026-08-25) addresses
  **order 2060 only**. It says nothing about the twelve orders.
- **Two independent public Lean formalizations of a single H(668)**,
  each a machine check of one supplied matrix, neither a construction
  theorem and neither reaching a second order:
  `Paul-Lez/hadamard-668-comparator` (created 2026-08-14, last push
  2026-08-17; PALOMAR-2026-08-17-000002), whose `Challenge.lean` exhibits
  the bordered structure at 668 independently — circulants on `Fin 166`,
  a width-4 border, per-block-constant strips, and byte-for-byte the same
  `4 × 4` corner these records carry; and `Arthur742Ramos/hadamard-668-lean`
  (created 2026-08-28; PALOMAR-2026-08-29-000009), whose README describes
  its input as "compact bordered block-circulant data" publicly posted on
  2026-08-12. (The quoted fragment is pinned in this laboratory's record;
  the surrounding wording is reported from the 2026-08-31 read.)
  Both re-read **2026-09-01** and both unchanged since the banked reads:
  Paul-Lez created 2026-08-14T10:50:28Z, last push still 2026-08-17T12:38:03Z;
  Ramos created 2026-08-28T14:03:24Z, last push still 2026-08-28T14:11:22Z.
  Both public, Apache-2.0, neither archived.
- **A third formalization is asserted but not publicly verifiable.**
  Ramos's `formalization.yaml` carries a *related formalizations* entry —
  `github.com/Arthur742Ramos/isabelle-afp-monorepo`, `relationship:
  "independent"` — whose note describes an Isabelle/HOL project holding "an
  earlier kernel-checked formalization of the same supplied certificate".
  **That repository returns HTTP 404 to an unauthenticated read on
  2026-09-01** and does not appear in that account's public repository
  listing. The existence, the content and the word *earlier* are therefore
  **asserted by that file and not publicly verifiable**; nothing here rests
  on them. The statement above stands as written — the claimed work is not
  public, and it is Isabelle, not Lean.
- **An anonymous inequivalence preprint at order 668.** The page
  `hadamard-668.vercel.app` hosts "Two H-inequivalent Hadamard matrices
  of order 668" (the title as the page renders it) — no author on the
  page, an empty PDF `/Author` field,
  and a PDF `CreationDate` of 2026-08-13 04:14 UTC+05:30 (2026-08-12
  22:44 UTC), which is a compilation timestamp and not by itself
  evidence of when the page became public; it was public when retrieved
  and verified firsthand on 2026-08-31, and still up when re-retrieved
  **2026-09-01**, title and separation statement unchanged. Its `H` is byte-identical to
  the decoded 668 record banked here (an independent decode of the same
  public data, border included); its second matrix is a 1,328-entry
  paired Hall switch of `H`, rebuilt and verified here (both its
  published digests reproduce). The priority statement it carries, and
  what remains this laboratory's, are in NOTE-B.md §3.4.
- **The `s = 1` layer is being read publicly.** On **2026-08-30** a
  third party (Tavis Rudd, `@tavisrudd`) posted the `s = 1` bordered
  Goethals–Seidel lemma with the row sums and the `Σ PAF ≡ −4` profile at
  the 668 instance (`x.com/tavisrudd/status/2093926416492732579`, read
  firsthand via `og:description` on **2026-08-31** and again through the
  fxtwitter API mirror on **2026-09-01**, which returns author
  `tavisrudd` and `created_at` **Sun Aug 30 04:58:55 +0000 2026** —
  matching the publication time `2026-08-30T04:58:55Z` decoded
  independently from the status id's snowflake timestamp). The payload
  read on 2026-09-01 shows the post as a **reply within his own thread**
  (`replying_to_status` `2093925739250978920`) and carries **no quote
  field**; what the thread root does or does not quote was not verified,
  and no relation to the announcement post is asserted here. That lemma
  is classical — Wallis–Whiteman 1972, with Spence 1975 the even sibling
  (NOTE-B.md §4) — and the post is recorded here as independent
  corroboration of the `s = 1` reading of the public payload, not as a
  source.
- **A fourth third-party public redistribution of the twelve orders.**
  `github.com/bbeartheancient/hoa64` (public; `created_at`
  `2026-08-13T03:39:26Z`, last push `2026-08-16T01:42:05Z`, located by the
  **2026-09-01** re-sweep) is a Hadamard construction and search toolchain
  whose README records verifying H(668) as a genuine Hadamard matrix and
  ingesting the twelve orders as gzipped CSVs. It states **no bordered
  variant, no characterisation, and no new-order or inequivalence claim** —
  it lists Goethals–Seidel only as a difference-family method. It is
  recorded here as further corroboration that the payload was broadly
  decoded within days of the announcement, and as no scoop of anything
  claimed in this repository.

## What is this laboratory's own

The theorems and their proofs (NOTE-B.md §1: the exact characterisation,
the house form, the parameter classification, the twist lemmas, the
`s = 1, i = 2` resolution with its index-2 collapse, the Gram rigidity
theorem and its `w > s` sharpening with the `w = s` boundary classified,
and the border-kit structure theorem with the exhaustive `(3,4)`
census); `verify/verify.py` and the bordered Goethals–Seidel assembler
and checkers in `tools/`; the matrices constructed here (NOTE-B.md §2.2
and §2.3), which instantiate the theorems at orders where existence has
long been settled and where no novelty is claimed; and the separation computation of NOTE-B.md §3 — the
exact 4-profile over all 8 222 179 035 row 4-subsets at order 668,
computed by two independent implementations that agree bin for bin.
Authorship and the division of labour are in
[DISCLOSURE.md](DISCLOSURE.md); the mathematical credit chain, including
the single hedged novelty statement and the exact list of sources that
bounds it, is NOTE-B.md §4.

## Release re-check (2026-09-01 — executed at the public flip)

Executed on the day the repository was made public. Each item carries the
result of its own run.

**The flip itself.** Both repositories were verified **publicly reachable by
unauthenticated HTTP GET on 2026-09-01**: anonymous reads of
`api.github.com/repos/JD-Jones-ASES/Hadamard-B` and
`.../Hadamard-B-Formal` return their metadata, and the pinned commit
`35665af5bc97bf3811bf971e39051727fcc99d10` resolves anonymously. The
server-side visibility field reads `"public"` (`private: false`) on both.
This document's edits ride on the commit made after those reads.

- [x] **Present-tense internal state.** 2026-09-01: the 2060 exact wording in
      `certs/07-2060-evidence/NOTES.md` is restated as dated fact — the
      registered `blas` pair completed and separates. A grep of both NOTES
      for in-progress wording leaves one hit,
      `certs/06-668-separation/NOTES.md`'s "in-flight artifact of another
      lane's run", which is a standing reason for citing rather than banking
      that third engine, not a stale state; left as written. The `--full`
      run and the third engine were already restated as done on 2026-08-31.
      The file-digest line now points at the digest section below.
- [x] **No-preprint statement.** 2026-09-01: arXiv API author sweeps
      (Alpöge, Voinov, Reynolds-Haertle), `all:"Hadamard matrix"` against
      668/1948/1964, a `ti:"Hadamard matrix"` sweep, and the newest
      `cat:math.CO` items — **no preprint**; the statement is re-dated, not
      retracted. Both Lean provenance files still cite the X posts as their
      only matrix-data source; the "and nothing else" phrasing was too
      strong and has been narrowed.
- [x] **The gist.** 2026-09-01: gist `b866a2ff…` re-read through the GitHub
      API — one revision only, `created_at` 2026-08-23T02:33:38Z,
      `updated_at` 2026-08-23T02:33:42Z, both unchanged, so the payload pin
      `c7a145d8…00f7` stands. `schneiderlo/hadamard-2060`'s last push is
      still 2026-08-27T13:40:24Z, `main` is still `9d4fdce…`, and the matrix
      file's blob sha is identical at `HEAD` and at the banked commit.
- [x] **Palomar registry.** 2026-09-01: `browse/index.json` reads 144
      entries over 159 versions and 18 days. Title-and-repository scan of
      2026-08-27 → 2026-09-01, plus a **full-record keyword scan of every
      one of the 26 entries dated 2026-08-30 → 2026-09-01** (each entry's
      own `entries/PALOMAR-…-v1.json` fetched and scanned for *hadamard,
      goethals, seidel, 668, 2060, 2092, sign stream, bordered*). The
      full-record scan of the earlier days was run 2026-08-31. **One hit,
      this laboratory's own**: PALOMAR-2026-08-31-000001,
      `JD-Jones-ASES/Hadamard-formal`. No Hadamard entry exists beyond the
      three already known — PALOMAR-2026-08-17-000002
      (`Paul-Lez/hadamard-668-comparator`), PALOMAR-2026-08-29-000009
      (`Arthur742Ramos/hadamard-668-lean`), PALOMAR-2026-08-31-000001.
      Hadamard-B itself is not registered.
- [x] **Corroboration items.** 2026-09-01: both Lean repositories re-dated
      and unchanged since the banked reads; the vercel preprint still up,
      its title corrected to the page's own casing; the Tavis Rudd post
      re-read and the unsupported "quote-tweeting" relation removed; two
      bullets added — the asserted, non-public Isabelle formalization, and
      `bbeartheancient/hoa64`, a fourth third-party public redistribution.
      No new public decode carrying a method, and no new public statement
      of the bordered form.
- [x] **Table-relative claims.** 2026-09-01: arXiv:2411.18897 re-checked —
      still v2 (`[v1]` 2024-11-28, `[v2]` 2025-08-30), no v3. The citation
      holds as written.
- [x] **Digests.** 2026-09-01: every file in `data/` digested at the frozen
      tree and recorded in *Frozen data digests* below; the placeholder
      line at the top of this document now points there.
- [x] **Cert 07's exact-file pins.** 2026-09-02: the bits confirmation pair
      completed under `experiments/inequiv/REGISTRATION-2060-exact.md` and
      reproduced the separating `blas` profiles bin for bin; `bank_exact.py`
      wrote the four-file bank and its digests are pinned in
      `EXACT_FILE_PINS`; cert 07 runs in EXACT mode. Post-flip data files
      (not in the frozen table below, which is the flip tree):
      `sep2060-exact-blas-plain.json` `5428aeac…0037d8`,
      `sep2060-exact-blas-gist.json` `a20b9a63…e3aef0`,
      `sep2060-exact-bits-plain.json` `e6c3af94…92ca49`,
      `sep2060-exact-bits-gist.json` `9d8cc4b5…7f8a93` (full digests in
      `certs/07-2060-evidence/NOTES.md`).
- [x] **Cert 13 banks.** 2026-09-02: the orientation-switch profiles at 668,
      producer-banked from `Hadamard-2060` under
      `experiments/inequiv/REGISTRATION-668-orientation.md`:
      `sep668-orient-exact-blas.json` `5ad283be…06c6fd`,
      `sep668-orient-exact-bits.json` `f0eeeb6c…45e726` (full digests pinned in
      `certs/13-668-orientation/run.py`). Matrix `H″` canonical digest
      `af1c285c…2953c7`, rebuilt in-run from `payload-records.json`.
- [x] **Cert 14 banks.** 2026-09-02: the orientation-switch profiles at 716,
      producer-banked from `Hadamard-2060` under
      `experiments/inequiv/REGISTRATION-716-orientation.md`:
      `sep716-orient-exact-blas.json` `b1c6b0ad…951cc4`,
      `sep716-orient-exact-bits.json` `bd1c3c23…77e69c` (full digests pinned
      in `certs/14-716-orientation/run.py` and tabled in its `NOTES.md`).
      Matrix `H″₇₁₆` canonical digest `a6b4f56e…885fcd`, rebuilt in-run from
      `payload-records.json`; the four comparison banks are cert 11's, reused
      unchanged and re-pinned by digest.
- [x] **Cert 15 banks.** 2026-09-02: the eight transposed-matrix profiles at
      716 and 668, producer-banked from `Hadamard-2060` under
      `experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC, before any
      matrix it governs was built; enumerated on a rented `c2d-highcpu-16`,
      16 threads, ~10:25–10:36 UTC). The matrices were built AND verified at
      the desk by `experiments/pr0042/build_matrices.py` — every digest in
      its `manifest.json`, each through this repository's `verify/verify.py`
      — and only then uploaded; the rented machine enumerated and nothing
      else. Transposed-matrix canonical digests, all four re-derived in-run
      by `certs/15-transpose-extended-668-716/run.py`: `H₇₁₆ᵀ`
      `e1c4a6fa…0be278`, `(H'₇₁₆)ᵀ` `41fe458a…21ea72`, `(H″₇₁₆)ᵀ`
      `7445c760…9de3ef`, `(H″₆₆₈)ᵀ` `49f97ecf…8538d9`. Bank files
      (full digests pinned in that `run.py` and tabled in its `NOTES.md`):

      | file | sha256 |
      | --- | --- |
      | `sep716-decoded-T-exact-blas.json` | `b2bd98a8c5a3403273408e009f88030b66a3994d6822be0e9f008db04bc512e5` |
      | `sep716-decoded-T-exact-bits.json` | `854e73fc748cb4b28113f267fc9266ede5a468105a3c32138978cafc4baaf040` |
      | `sep716-twisted-T-exact-blas.json` | `135189d8fe4dd619d72ff8ea0ff3cc4b94c07334c1cde4e227f4e08c4cb91fcb` |
      | `sep716-twisted-T-exact-bits.json` | `2c13fb2e980782b8298e0bf9a37d7bbde569f9a832c5886aa2ec7e4198c1a87b` |
      | `sep716-orient-T-exact-blas.json` | `b287dfcb5a26e5cf47f8ac8ec445725253cc1a4a51da67eaa219431a8d2da062` |
      | `sep716-orient-T-exact-bits.json` | `be270c974ce237f537cb39a550819169849805f691120976e71ecfc673b30704` |
      | `sep668-orient-T-exact-blas.json` | `536e0c136d16b8271a3c7916529fd57c86f6738a8695c1d691f4fbf98455ec2d` |
      | `sep668-orient-T-exact-bits.json` | `81fae9e0a4f23116067769bffc80451f2c1d834101f9b65f90a2cccf4aab70cf` |

      The seventeen comparison banks are certs 06, 08, 11, 13 and 14's,
      reused unchanged and re-pinned by digest.
- [x] **Cert 19 banks.** 2026-09-02: the exact 4-profile of `Hᵀ` at 668 — the
      transpose of the decoded record — produced **at the desk** (three
      threads; 90.2 s `blas`, 169.4 s `bits`) under
      `experiments/pr0042/REGISTRATION.md`, **Amendment 2**, flushed ~18:17
      UTC (commit `48e5e47`), before the matrix it governs was built; its
      result paragraph is dated 18:27Z. The matrix was built by transposing
      the cert-01 record rebuilt through `tools/bordered_gs.py`, and passed
      this repository's `verify/verify.py` before the run. Matrix `Hᵀ` at
      668, canonical digest `7a25f3f7…62edb25`, re-derived in-run by
      `certs/19-668-transpose-eight-classes/run.py`. Bank files (full digests
      pinned in that `run.py` and tabled in its `NOTES.md`):

      | file | sha256 |
      | --- | --- |
      | `sep668-decoded-T-exact-blas.json` | `fb6917010cbd9147312c3908a407bcd86252d22bbf343ae3434f182b71e302f5` |
      | `sep668-decoded-T-exact-bits.json` | `eff70d552abfe6a85327243a6cbb75d81d36589ab87806c890e5545ea9ef57a5` |

      The thirteen comparison banks are certs 06, 08, 13 and 15's, reused
      unchanged and re-pinned by digest.
- [x] **Cert 20 banks.** 2026-09-02: the six exact 4-profiles at order 1676
      — the decoded `(1,1)` record `H`, its Lemma-T `i = 2` rebuild `H′`,
      and the orientation switch `H″`, each in both arithmetics —
      producer-banked from `Hadamard-2060` under
      `experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC, before any
      matrix it governs was built; Amendment 1 ~11:05 UTC re-priced the
      campaign). Enumerated by the unchanged engine
      `experiments/inequiv/exact_profile_big.py` on a rented
      `c2d-highcpu-16` (`prof42-1`, `us-east1-b`), 16 threads,
      2026-09-02 10:36Z–19:55Z; the `blas` legs took 2 638.7 s / 2 652.4 s
      / 2 649.4 s at ≈ 762 MB peak and the `bits` legs 6 628.6 s /
      6 642.9 s / 6 631.5 s at ≈ 104 MB. As for cert 15, the matrices were
      built AND verified at the desk by
      `experiments/pr0042/build_matrices.py` — every digest in its
      `manifest.json`, each through this repository's `verify/verify.py`,
      the two pinned ones reproducing their cert-01 and cert-02 pins — and
      only then uploaded; the rented machine enumerated and nothing else.
      Matrix canonical digests, all three re-derived in-run by
      `certs/20-1676-three-classes/run.py`: `H₁₆₇₆` `8e919c2b…1cdb99`
      (cert 01's pin), `H′₁₆₇₆` `6a493837…4907405` (cert 02's pin, seeds
      re-derived here as the `ψ`-twist of the decoded seeds), `H″₁₆₇₆`
      `16d1617c…6a84346` (formed in-run from `H`). Bank files (full digests
      pinned in that `run.py` and tabled in its `NOTES.md`):

      | file | sha256 |
      | --- | --- |
      | `sep1676-decoded-exact-blas.json` | `57b9a43caf5246de779ad3205a45642c98f7a211be47e4ed12d718fe098781c9` |
      | `sep1676-decoded-exact-bits.json` | `469e0b0382d479a6d917316246807cadfa1f113bb2bfcd1429ec1712622e7b94` |
      | `sep1676-twisted-exact-blas.json` | `328ed05c9614a223d95bd35583c83433a8700249c0b50872fbfc8d846e9b5a49` |
      | `sep1676-twisted-exact-bits.json` | `311ef88606d0543967e2b0cf46aad4f3fb3f1353cc59b018c5369e447c0c2bb1` |
      | `sep1676-orient-exact-blas.json` | `a83b239695a3bd820de222e829e65a10a5dd66a432858af57cc950eb4ff40be2` |
      | `sep1676-orient-exact-bits.json` | `af198e51aecd165e8a2a22ee5ece8dfa73d8ddedf314fa94684ee367db14e9d5` |

      All six are cert 20's own; it reuses no comparison bank from another
      certificate, and the two shared parameter records
      (`payload-records.json`, `twisted-i2-records.json`) are bound to it
      by the canonical digests above rather than by a file pin. Cert 20's
      `--full` has **not** been run in this repository — at this order one
      leg is of order 6–7 h and the `blas` route wants ≈ 9.4 GB — so its
      verdict is an audit of these banks and says so.
- [x] **Cert 21 banks.** 2026-09-02: the two remaining exact 4-profiles at
      order 1676 — the transposes `(H′)ᵀ` and `(H″)ᵀ` of the Lemma-T rebuild
      and of the orientation switch, each in both arithmetics —
      producer-banked from `Hadamard-2060` under the same
      `experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC, before any
      matrix it governs was built; §2 lists `H_1676-twisted-T` and
      `H_1676-orient-T` and no third 1676 transpose, so **`Hᵀ` at 1676 was
      never enumerated** and nothing is claimed about it). Enumerated by the
      unchanged engine `experiments/inequiv/exact_profile_big.py` on the same
      rented `c2d-highcpu-16` (`prof42-1`, `us-east1-b`), 16 threads,
      2026-09-02: the `blas` legs took 2 650.8 s (`twisted-T`) and 2 665.8 s
      (`orient-T`) at ≈ 762 MB peak, the `bits` legs 6 661.6 s and 6 660.6 s
      at ≈ 104 MB. As for certs 15 and 20, the matrices were built AND
      verified at the desk by `experiments/pr0042/build_matrices.py` — digests
      in its `manifest.json`, each through this repository's
      `verify/verify.py` — and only then uploaded; the rented machine
      enumerated and nothing else. Matrix canonical digests, both re-derived
      in-run by transposing the rows this repository rebuilds:
      `(H′)ᵀ` `46f29432…4af261ce3`, `(H″)ᵀ` `ad870e1c…ea85303d`. Bank files
      (full digests pinned in that `run.py` and tabled in its `NOTES.md`):

      | file | sha256 |
      | --- | --- |
      | `sep1676-twisted-T-exact-blas.json` | `d951a78a0af94cb3979eec3878583d750ecb2c417946ab4c6159a3a7bfaaa040` |
      | `sep1676-twisted-T-exact-bits.json` | `2a249d15e2d8c0aaa71d0fbeb503da43175aed7b9466c32dfaa77db1401bad6c` |
      | `sep1676-orient-T-exact-blas.json` | `4d7506d93e35125d236432d92457d60321c1e1c4f014febc45e5ded0b06b820a` |
      | `sep1676-orient-T-exact-bits.json` | `46750b71d6400a51ddfde19ec574cc6a8d92c345aabc521f5cd866ca3eb84e41` |

      The six comparison banks are cert 20's, reused unchanged and re-pinned
      by digest. Banking was by the laboratory's `experiments/pr0042/bank.py`
      under `--cert 21`, which adds only the seven header fields each bank's
      `banked_note` names and nothing numeric. Cert 21's `--full` has **not**
      been run in this repository either, for cert 20's reasons and at cert
      20's price, so its verdict is an audit of these banks and says so.
- [x] **Cert 22 banks.** 2026-09-03: the exact 4-profile of `H_2060-orient`
      — the plain Goethals–Seidel realisation at order 2060 with its twelve
      off-diagonal `515×515` blocks negated, the orientation switch in its
      **unbordered** form (2060 is the `s = 0` layer: no border, four
      circulant blocks, `N = 4·515`) — in both arithmetics, producer-banked
      from `Hadamard-2060` under the same
      `experiments/pr0042/REGISTRATION.md` (flushed 10:17 UTC 2026-09-02,
      before any matrix it governs was built; §2 registers this object in as
      many words, and §4 fixes the decision rule "`H″ ≁ plain` and
      `H″ ≁ gist` ⟹ order 2060 carries at least THREE classes (row-side)"
      in advance. Amendment 1, ~11:05 UTC, added the second instance for the
      2060 legs and re-priced the campaign; both predate the leg).
      Enumerated by the unchanged engine
      `experiments/inequiv/exact_profile_big.py` on a second rented
      `c2d-highcpu-16` (`prof42-2`, `us-east1-b`, created 18:14 UTC
      2026-09-02), 16 threads: the `blas` leg landed 20:10 UTC 2026-09-02 in
      6 631.0 s at 860.7 MB peak, the `bits` leg 01:00 UTC 2026-09-03 in
      17 654.7 s at 131.9 MB, and they agree bin for bin. As for certs 15,
      20 and 21, the matrix was built AND verified at the desk by
      `experiments/pr0042/build_matrices.py` — digest in its
      `manifest.json`, through this repository's `verify/verify.py`, with
      the plain source re-verified against cert 07's pin first — and only
      then uploaded; the rented machine enumerated and nothing else. Matrix
      canonical digest, re-derived in-run by negating the twelve blocks of
      the plain array this repository rebuilds: `H″₂₀₆₀`
      `4e1891b0…0870b801`. Bank files (full digests pinned in that `run.py`
      and tabled in its `NOTES.md`):

      | file | sha256 |
      | --- | --- |
      | `sep2060-orient-exact-blas.json` | `38135aef205b4428760dc0439b29196b6776f215397fb34880156c66ee283f00` |
      | `sep2060-orient-exact-bits.json` | `18dcf4c2e5d603324182eeef45c13e89ef80d5c0c0a5add3a7f90c333e4e5e87` |

      The four comparison banks are cert 07's, reused unchanged and re-pinned
      by digest at exactly the values cert 07's own `EXACT_FILE_PINS` carries,
      so the two certificates cannot drift apart on the shared evidence.
      Banking was by the laboratory's `experiments/pr0042/bank.py` under
      `--cert 22 --date 2026-09-03` (2026-09-03 UTC, the date the banks' own
      `banked_note` carries), which adds only the seven header fields that
      note names and nothing numeric; the `description` field names the base
      as the plain GS realisation (cert 07 pin) with its twelve off-diagonal
      515×515 blocks negated, no border. Cert 22's `--full` has **not**
      been run in this repository — at this order one leg is of order 15 h
      and the `blas` route wants ≈ 17.5 GB — so its verdict is an audit of
      these banks and says so. The 2060 transposes (`H_2060-orient-T`,
      `H_2060-plain-T`) are pending legs of the same campaign and are **not**
      banked, so that order's statement stays row-side.
- [x] **`CITATION.cff`.** 2026-09-01: `date-released` set to 2026-09-01, the
      flip date. Version left at 0.1.0.

## Frozen data digests (2026-09-01)

Every file in `data/`, SHA-256, computed with the standard library over a
clean working tree at commit `35665af5bc97bf3811bf971e39051727fcc99d10` —
the commit that was pushed at the public flip. These are the twenty files
this repository redistributes; no matrix is among them.

| file | sha256 |
| --- | --- |
| `h20-boundary.json` | `716610543b79ab9e1c9f1adb142c114544e70e58d370500130b83c17a18cf254` |
| `h52-gate.json` | `ef60c4ff9f245eec5ba7f035e5968152836207fcd9235a0b1851150d2fb1d170` |
| `h76-nonscalar.json` | `500ba1d22787407183d91cc303a80cfc14d250a5fe2f2d5c0d665553ffc7b8bf` |
| `n1916-twist.json` | `1a3f92228074f69a7ead11d66371d18dfb39aeb2d17155f3b7fc9782b7b8d51b` |
| `payload-records.json` | `9afbf60efed82e97ffe229e0059f06dfa7849a2ac2434c52a6ada09559738afb` |
| `sep2060-records.json` | `1c9742fe485f5cc7232c2e876d322ca0094270bc40e537d16da422c0a13202bc` |
| `sep2060-sampled-histograms.json` | `586d9fe51a4c24448dd44efbe1e3e60a3ec0c83167e263fb008f538efe2ecdc8` |
| `sep668-exact-bits-decoded.json` | `7bace61441f17b5e95fff433bdc5939da212e2b8735e8738d7ed3078fae456b7` |
| `sep668-exact-bits-twisted.json` | `f40bbb8c3906d6fc7374e3e04c2b68eaf29393e50b5662f69eee2426ed3f1e9a` |
| `sep668-exact-blas-decoded.json` | `370fffe6c2f5dc53c09d3b74f8c09dd2bc2a39a1ac2b27fb5167ab4d3559387b` |
| `sep668-exact-blas-twisted.json` | `8526b3cfa7938a9af334e23f722b1c215ffd1e318c0c713ecc3da1b91f5b3afe` |
| `sep668-hall-T-exact-bits.json` | `48fdb26f8b1ee5135ed278ec866e204c1ab47df168c043fabff8699c0f4fd8bb` |
| `sep668-hall-T-exact-blas.json` | `151fb5d6e70cf56d6a1c2aa124a597a837bca0ecf5d64958b43a34c05383e0db` |
| `sep668-hall-exact-bits.json` | `a6f703b499d98995f6446a1aed671284c47e99cfe869f3ce8dc8b5fd9394accb` |
| `sep668-hall-exact-blas.json` | `35e716ecb43bb6190d5dd6f4160e0bc2bed4f61a3aacf07a36ff9d190810c154` |
| `sep668-hall-switch.json` | `13efd2402b8394c62c901af4f7cfbec7b2e474832dd3055c6b9e9e220b351c85` |
| `sep668-sampled-histograms.json` | `a21973871c0ad80a8a3b95e057066ae72710194a3638ce3228c647ce502804a6` |
| `sep668-twisted-T-exact.json` | `38355274ec61d33fcd96e24255e4a7b02874150cd914fdfb928d28cee751fc4a` |
| `sep668-twisted-record.json` | `fe8154179ba2ebfe097c82e468368cdc8a070548555bb10140949af0560611fb` |
| `twisted-i2-records.json` | `aafa83e070d2dc59da80aec1bcb6457b6bfcb8d7bcc758f20d9344047bfcb079` |

`payload-records.json`'s digest is the same `9afbf60e…` quoted in *What this
repository redistributes* above, at the same 29,138 bytes.

**Added after that frozen tree.** Four exact 4-profiles at order 716, banked
2026-09-01 for cert 11 (producer: `Hadamard-2060`,
`experiments/inequiv/exact_profile_big.py`, under the pre-registration
`experiments/inequiv/REGISTRATION-716-exact.md`). The table above is the
state at commit `35665af5…` and is not rewritten; these four are listed here
so the ledger stays complete. Each is also SHA-256 pinned in
`certs/11-716-separation/run.py`.

| file | sha256 |
| --- | --- |
| `sep716-exact-bits-decoded.json` | `a0d5b3a65b83c39c905ec2a1d3b25ca1c58e0106b76aaa6eb3b2feee3748aeed` |
| `sep716-exact-bits-twisted.json` | `c385773a7a2bf4506b94752406be787bac3d85b8a9ee18bb23668080c6afe7bc` |
| `sep716-exact-blas-decoded.json` | `80ee1e151ec1f759d7213d500623603716b9afa6fc382a385ce6970efac35a6b` |
| `sep716-exact-blas-twisted.json` | `e2076eb890557775edaccda3c9dcbab7e585f5181e0da1bb5914913bf0749b46` |

And one more, banked 2026-09-01 for cert 12: the six `(H1)`-realizability
witnesses at `(s,i) = (3,4)`, six `±1` tables `Q ∈ {±1}^{16×12}`, derived and
verified at the source laboratory 2026-09-01 by an exhaustive search with two
verdict-exact symmetry reductions. Nothing about that search is trusted in
this repository: a witness certifies itself, and
`certs/12-gram-rigidity/run.py` re-checks every inner product of `QQᵀ`
against `I₄ ⊗ M`. The file is SHA-256 pinned in that `run.py`.

| file | sha256 |
| --- | --- |
| `gram34-witnesses.json` | `abaf4728e8ba5cd737024b9ab319640c8c634e497884a444b683fb5ee4b93307` |

And one more, banked 2026-09-02 for cert 18: `data/cell24-records.json`,
the four even-`s` `(2,4)` coset-border records — seeds, coset sums, S-part,
the `8×8` halves `Q′` and `P′` of the column and row tables, the reflection
`ρ` and its image in `Ḡ`, the corner, and the canonical digest of the
assembled matrix. Everything in it except the seeds, `Q′`, `P′` and `ρ` is
redundant and is **checked rather than trusted**: cert 18 recomputes the
coset sums and the S-part from the seeds, recomputes the corner as
`E = −(1/16)PĈᵀQ` from the compression lemma, assembles the matrix with its
own assembler, and compares the digest three ways. The file's schema is
documented in `certs/18-cell24-instances/NOTES.md`, and its SHA-256 is
pinned in that certificate's `run.py`.

| file | sha256 |
| --- | --- |
| `cell24-records.json` | `9727b392940d416d3f25dca5d51d2db71cd499bc73c3b8dc4efd22801180f179` |

## Certificates 16, 17 and 18 — where the code came from (2026-09-02)

Certificates `16-theorem-eprime-boundary`, `17-border-kits-34` and
`18-cell24-instances` were first written in the source laboratory
(`Hadamard-2060`, private) by a **Cursor cloud lane running Fable 5.1**, on
**2026-09-02**, as the adversarial pass over that laboratory's fleet
adjudication of 2026-09-01 — the lane's job was to try to break the claims
before they entered this note, and it filed the three certificates as the
part that survived. The desk replayed all three (default paths, and cert
17's `--full`) before adopting them, and they were then **ported here** in
this repository's house form: rebuilt against `verify/verify.py` and
`data/`, with every read of a source-laboratory path removed, so that
nothing in `certs/16`, `certs/17` or `certs/18` imports or opens anything
outside this repository.

Cert 17's full-census digest `d60b83d8…7f4d` first travelled as a pin — the
one the source laboratory's run of the identical code recorded — and the
default path of the port reproduced the *sample* digest `47d5d44e…e36c35`
bit for bit. The 16 384-class run has since been made **in this repository**
(2026-09-02, Python 3.14, one worker, 5 min 28 s; 2048 / 2048 on each of the
eight `(group, κ(ρ))` sweeps, 16 384 / 16 384, the digest matching the pin),
so that number is now a reproduction rather than a quotation. The four
canonical matrix digests of cert 18 were first computed upstream; they are
re-derived here from `data/cell24-records.json` at certificate run time, by
`verify/verify.py`, on every run. As everywhere else in this repository,
the source laboratory is **not** part of the trust chain.
