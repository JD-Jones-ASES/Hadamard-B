# PROVENANCE — the dated chain for everything sourced from outside

Two bodies of external material enter this repository: the twelve
publicly posted Hadamard matrices verified in movement II, and the
publicly posted matrix of order 2060 used in movement III. Both chains
are given below with their dates, their retrieval, and the digests that
pin them. Nothing else here comes from outside; what is this
laboratory's own is stated at the end.

File digests are recorded at freeze.

## The twelve public records

### The announcement

On **2026-08-12** Levent Alpöge (`@__alpoge__`) posted, on X, seed data
for Hadamard matrices at the twelve orders 668, 716, 892, 1132, 1244,
1388, 1436, 1676, 1772, 1916, 1948 and 1964. The announcement is three
posts and no paper:

| post | what it is |
| --- | --- |
| `https://x.com/__alpoge__/status/2087504785952182273` | the payload. `created_at` **Wed Aug 12 11:41:39 +0000 2026**; body is pure `+`/`-`, 23,828 characters, `display_text_range = [0, 23828]`, no prose at all. Read firsthand through the fxtwitter API mirror on **2026-08-16** and again through `og:description` on **2026-08-31**. |
| `https://x.com/__alpoge__/status/2087504788938510427` | the decoder: an obfuscated POSIX-shell program, 6,130 bytes, that reads the first post as an input tape. The program itself was extracted and read firsthand (SHA-256 `cf0ac77c3c6329d7d904466c18e653229521913ed19f2535322a9c5abc986f05`); **the status URL is recorded from a relayed station's script, not fetched from X by this laboratory** — the 2026-08-16 sweep could not enumerate the author's timeline. |
| `https://x.com/__alpoge__/status/2087504790435840207` | the credits, verbatim and entire: *"weekend fun w @tehwalris, Saul Reynolds-Haertle, and of course claude:)) i only claim bad suggestions!!"* (`@tehwalris` = Philippe Voinov). Read firsthand through `og:description` on **2026-08-31**. |

**What was disclosed: data and a decoding grammar. No method.** No
construction is described in any of the three posts, no search is
described, and no preprint accompanies them (see *Independent public
corroboration* below). The public record of the event, as read on
2026-08-16, is secondary: the Wikipedia article *Hadamard matrix*
(edited 2026-08-12) credits Philippe Voinov, Saul Reynolds-Haertle,
Claude and Alpöge and states that no further details of the
construction were provided; OEIS A007299's 2026 history line lists the
same twelve orders; and the earliest public reading of the second post
— that it expands the tape through a Goethals–Seidel construction, and
which twelve orders it covers — is a third-party post,
`https://x.com/adi_baradwaj/status/2087620459185819903`
(2026-08-12 19:21 UTC).

### The public artifacts, and how they are pinned

| artifact | retrieved | pin |
| --- | --- | --- |
| **the tape** — the 23,828-character sign stream of the payload post | fetched through the fxtwitter API mirror **2026-08-16**, re-fetched and banked **2026-08-28** | normalised sign stream, SHA-256 `5b5fe8fa42f0d6a8b4e4c9926726d82a6aab8e1070c1ae4d1b430c1277e58db4`, length 23,828. This digest is carried in `data/payload-records.json` (`tape_sha256`). |
| **the expanded matrices** — `github.com/foocker/Hadamard668`, a third-party public decode created **2026-08-13**, the day after the announcement | cloned read-only **2026-08-28** | its `answer.md` normalises to the same tape digest. It held ten complete matrices; its H(1948) stopped at **603 complete rows plus a 305-entry partial row — 1,174,949 of 3,794,704 entries, 30.96 % of that matrix** — and it held no H(1964) at all. Its decoder script was **never executed**; its outputs were used as data only. |
| **a zip of twelve CSVs**, circulated as the announcement matrices | received by relay and verified **2026-08-31** | 4,873,392 bytes, SHA-256 `503a4352f282a3d0141a1814b4f89036961524934670a5603ec33a95ea6f179d`; twelve entries, all `open_hadamard/hadamard_<N>.csv`, internal mtimes all `2026-08-12 05:58:58`. That it was posted by Sumeet Motwani (`@sumeetrm`) on Google Drive as Alpöge's matrices is **REPORTED-FROM-SECONDARY-SOURCE**; no source URL is recorded here. Nothing in this repository rests on that attribution. |

The three copies agree. All twelve matrices of the zip verify green in
the trust chain and are **entry-identical** to the matrices this
laboratory rebuilt from tape bytes — including the 2,619,755 entries of
H(1948) that were never public and all 3,857,296 entries of H(1964),
which had no public copy anywhere.

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
third-party method disclosure was used — none exists.** The
correctness of the decode is not asserted, it is replayed: each cert
rebuilds its matrices from these records and hands them to
`verify/verify.py`. No generated matrix is committed; the certificates
rebuild them and delete them.

### Priority posture

**The twelve records are the announcing team's mathematical content.**
The decode is provenance, not an achievement, and **no priority claim of
any kind is made on the records or on anything derived from them.** No
novelty of existence is claimed at any of the twelve orders. That those
orders were open before 2026-08-12 is not proven here either: it is
REPORTED-FROM-AUDITED-TABLE, from Table 4 of Cati–Pasechnik,
arXiv:[2411.18897](https://arxiv.org/abs/2411.18897)v2 (2025-08-30).
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

**Credit for the order-2060 matrix is Schneider's, and he published its
construction first.** `github.com/schneiderlo/hadamard-2060` (created
**2026-08-25T23:38:42Z**, last push 2026-08-27T13:40:24Z) gives the
Cooper–Wallis route `T(103) ⊗ W(5) → H(2060)` from the producer's own
side; its explicit matrix file (SHA-256
`8c81090db76e2b503561ba82598b514e0c2e39940a3703783c5b1b2c6b80a37a`)
parses to the canonical digest `c7a145d8…` above — measured 2026-08-31,
the repo and the gist are one object. This laboratory decoded that
matrix's structure independently from the banked artifact alone, and was
not first to state it publicly; no firstness is claimed. The plain array
is a genuinely different matrix, but a derivative of the public seed, not
an independent find. The separation of the two is stated at
COMPUTATIONAL-EVIDENCE and nowhere stronger (NOTE-B.md §3.5).

## Independent public corroboration

Dated observations of the public record, all made **2026-08-31** by a
sweep lane of this laboratory (arXiv API author sweep; GitHub API and raw
sources; the Palomar registry's CC0 feed). They are reported as of that
date and are re-checked at release.

- **No preprint.** An arXiv author sweep for Alpöge, Voinov and
  Reynolds-Haertle returns zero Hadamard papers. As of 2026-08-31 there
  is no arXiv preprint of the announcement, and its only public statement
  is the X posts. This is corroborated from a second direction: the
  provenance files of both downstream Lean repositories cite the X posts
  as `type: "web post"` and nothing else.
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
- **An anonymous inequivalence preprint at order 668.** The page
  `hadamard-668.vercel.app` hosts "Two H-inequivalent Hadamard Matrices
  of Order 668" — no author on the page, an empty PDF `/Author` field,
  and a PDF `CreationDate` of 2026-08-13 04:14 UTC+05:30 (2026-08-12
  22:44 UTC), which is a compilation timestamp and not by itself
  evidence of when the page became public; it was public when retrieved
  and verified firsthand on 2026-08-31. Its `H` is byte-identical to
  the decoded 668 record banked here (an independent decode of the same
  public data, border included); its second matrix is a 1,328-entry
  paired Hall switch of `H`, rebuilt and verified here (both its
  published digests reproduce). The priority statement it carries, and
  what remains this laboratory's, are in NOTE-B.md §3.4.
- **The `s = 1` layer is being read publicly.** On 2026-08-31 a third
  party (Tavis Rudd, `@tavisrudd`) posted, quote-tweeting the
  announcement, the `s = 1` bordered Goethals–Seidel lemma with the row
  sums and the `Σ PAF ≡ −4` profile at the 668 instance
  (`x.com/tavisrudd/status/2093926416492732579`, read firsthand via
  `og:description`). That lemma is classical — Wallis–Whiteman 1972, with
  Spence 1975 the even sibling (NOTE-B.md §4) — and the post is recorded
  here as independent corroboration of the `s = 1` reading of the public
  payload, not as a source.

## What is this laboratory's own

The theorems and their proofs (NOTE-B.md §1: the exact characterisation,
the house form, the parameter classification, the twist lemmas, and the
`s = 1, i = 2` resolution with its index-2 collapse); `verify/verify.py`
and the bordered Goethals–Seidel assembler and checkers in `tools/`; the
matrices constructed here (NOTE-B.md §2.2), which instantiate the
theorems at orders where existence has long been settled and where no
novelty is claimed; and the separation computation of NOTE-B.md §3 — the
exact 4-profile over all 8 222 179 035 row 4-subsets at order 668,
computed by two independent implementations that agree bin for bin.
Authorship and the division of labour are in
[DISCLOSURE.md](DISCLOSURE.md); the mathematical credit chain, including
the single hedged novelty statement and the exact list of sources that
bounds it, is NOTE-B.md §4.

## Release re-check (____-__-__)

To be executed at the public flip, and the date filled in above.

- [ ] **No-preprint statement.** Re-run the arXiv author sweep (Alpöge,
      Voinov, Reynolds-Haertle; also a title/abstract sweep for the
      twelve orders). Re-date the statement, or retract it if a preprint
      has appeared, and re-check the two Lean repositories' provenance
      files still cite the X posts as their only source.
- [ ] **The gist.** Re-read
      `gist.github.com/schneiderlo/b866a2ff2fcd93934f0db54cfa4069d0`
      through the GitHub API: created/updated timestamps, and the payload
      digest still `c7a145d8…00f7`. Re-check the last push of
      `schneiderlo/hadamard-2060` and that its matrix still parses to the
      same canonical digest.
- [ ] **Palomar registry.** Re-sweep for Hadamard entries
      (`browse/index.json` → per-day pages → each entry record; keyword-scan
      the full record text, not titles). Send a browser user-agent — the
      default `urllib` agent gets HTTP 403 and the sweep then silently
      returns nothing. Re-state the entry count and the entries.
- [ ] **Corroboration items.** Re-date each bullet above: creation and
      last-push dates of the two Lean repositories and their scope lines;
      any new public decode of the payload; any new public statement of
      the bordered form.
- [ ] **Table-relative claims.** If any claim citing Cati–Pasechnik
      ships, re-check arXiv:2411.18897 is still at v2 (2025-08-30), the
      version cited here.
- [ ] **Digests.** Record the `data/` file digests and pin them in this
      document, replacing the placeholder line at the top.
