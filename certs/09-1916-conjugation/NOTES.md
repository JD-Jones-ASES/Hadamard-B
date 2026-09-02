# cert 09 — order 1916: a non-house Gram that is only a diagonal conjugation

**Label: PROVEN-BY-CERTIFICATE.** Replay:
`python certs/09-1916-conjugation/run.py` from the repository root. Standard
library only, about twelve seconds, 69 checks, exit 0.

---

## What this certifies

Two statements in `note/NOTE-B.md` that were previously carried on upstream
reports are now replayable **inside this repository**, from the banked record
and nothing else:

* **§1.2's non-house witness.** "Any admissible Gram `M` is house" is
  **false as a literal statement**: the character twist of the order-1916
  record satisfies every hypothesis of the master theorem with a Gram that is
  not `(4s+4)I_i − 4J_i`.
* **§3.3's conjugation.** That witness is nevertheless `S · H · S` for a
  diagonal `S`, so it **manufactures no new matrix**. This is the
  `ψ(ρ) = 1` case of §1.4's proposition, doing exactly what the proposition
  promises.

Both halves hold together, and that is the point: the counterexample to
Gram-forcing is real, and it is worthless as a source of matrices.

## The theorem

> **Theorem.** Let `H` be the decoded `(s, i) = (3, 4)` bordered
> Goethals–Seidel record at order 1916 (`data/payload-records.json`), with
> `G = Z₂ × Z₂ × Z₁₁₉`, `K` of index `i = 4`, `Ḡ = G/K = Z₂ × Z₂`, `w = 119`,
> and `ρ = 0`. Let `ψ = ψ̄ ∘ κ` for a character `ψ̄` of `Ḡ`, and let
>
> ```
> x'_q = ψ · x_q ,   P' = P D̄ ,   Q' = D̄ Q ,   E' = E
> ```
>
> be the character twist (`data/n1916-twist.json`). Then
>
> **(a)** the twisted instance satisfies H1–H4 in their **general** form,
> with Gram `M'(ē) = ψ̄(ē) · M(ē)`, and for `ψ̄` nontrivial that `M'` is
> **not** the house form; and
>
> **(b)** the assembled matrix is `H' = S · H · S` exactly, with
> `S = diag(I_{4s}, I₄ ⊗ diag(ψ(g)))` — 952 rows and the **same** 952
> columns negated, **no permutation at all**.

Both halves are checked in exact integer arithmetic, (b) over all
`1916² = 3 671 056` cells, for **all four** characters of `Ḡ`, not just one.

## Why `ψ(ρ) = 1` here, and why that matters

`r_shift = [0,0,0]`, so `ρ` is the identity of `G` and `κ(ρ) = 0`; hence
`ψ(ρ) = ψ̄(0) = +1` for **every** character. `Ḡ = Z₂ × Z₂` is elementary
abelian, so every `ψ̄` is real and `ψ² = 1` and `ψ|_K = 1` come for free.
§1.4's proposition therefore applies to all four, and none of them can leave
the equivalence class.

This is the exact opposite of the situation at order 668, where `G = Z₁₆₆`,
`ρ = −1`, and the unique `ψ` has `ψ(ρ) = −1`: there the proposition does
**not** apply, and the twist provably **does** leave the class (certs 06
and 08). The two certificates are the two sides of §1.4.

## The evidence chain

`run.py` runs all of the following before it prints a verdict.

**[0] The bank is pinned.** `data/n1916-twist.json`, SHA-256 compared in code.

**[1] The house record is rebuilt, not assumed.** It goes through
`tools/bordered_gs.py`'s `check_record`, which re-checks **every hypothesis of
the master theorem** — H0 shape, H1 the two-tier PAF profile, H2 the
corner/row-table budget, H3 the column-table Gram, H4 the coupling, the
derived D1/D3/D5 and the Σ̄ law, and the compression-lemma cross-check — then
goes to `verify/verify.py` (all 1 834 570 row pairs) and its canonical digest
is compared against the pin.

**[2] The character group is derived, not transcribed.** From the record's own
`coset_divisors` `run.py` builds `Ḡ`, all four characters
`χ_t(c) = ∏_j (−1)^{t_j c_j}` in the record's own coset indexing, and checks
for each that `ψ̄² = 1`, `ψ̄|_K = 1` and `ψ̄(κ(ρ)) = +1`. It then checks that
the banked `psibar` is the first-coordinate character — the one the upstream
witness used.

**[3] Each twist is checked in the GENERAL form.** `tools/bordered_gs.py`'s
`check_record` hard-codes the **house** two-tier PAF profile, so it cannot be
used on a non-house instance. `run.py` therefore re-implements H1–H4 with an
**arbitrary** Gram — labelled as `tools/bordered_gs.py` labels them (see its
`LABEL MAPPING` block; `note/NOTE-B.md`'s Theorem A permutes these labels):

| | general form checked here |
| --- | --- |
| **H1** | `Σ_q PAF_q(t) = 4n` at `t = 0` and `−M(κ(t))` elsewhere. In the house form `M(0) = 4s` and `M(ē) = −4`, which is exactly the two-tier profile `−4s` on `K∖{0}` and `+4` off `K` — so this is a strict generalisation, not a different check. |
| **H2** | `E Eᵀ + w P Pᵀ = N·I_{4s}` |
| **H3** | `Q Qᵀ = I₄ ⊗ M`, with `M` a function on `Ḡ` — the off-diagonal superblocks vanish and each diagonal superblock's `(a,b)` entry depends only on `a − b`. **This is where the Gram lives**, and the only place the house form is normally imposed. |
| **H4** | `E Qᵀ + P Ĉᵀ = 0`, `Ĉ` the `G/K` Goethals–Seidel array of the coset sums |

Then, per twist: the compression lemma is cross-checked against the *twisted*
`Ĉ`, the instance is assembled, handed to `verify/verify.py`, and its
canonical digest compared against the banked value.

**[4] The Gram is measured, not asserted.** For each twist `run.py` reports
`M`, checks `M'(ē) = ψ̄(ē)·M(ē)` exactly, and checks that the trivial
character reproduces the house Gram and the house matrix while the three
nontrivial ones do not.

| `ψ̄` | `M(0), M(1), M(2), M(3)` | house? | canonical SHA-256 |
| --- | --- | --- | --- |
| `χ₀ = [+,+,+,+]` | `12, −4, −4, −4` | **yes** | `be2073ee…577ae9` (= the record) |
| `χ₁ = [+,−,+,−]` | `12, +4, −4, +4` | no | `9a746db1676429edd6099d1f…` |
| `χ₂ = [+,+,−,−]` | `12, −4, +4, +4` | no | `05d411fa…97175c` **(the banked witness)** |
| `χ₃ = [+,−,−,+]` | `12, +4, +4, −4` | no | `a5db39097f566a73f17d047c…` |

The house form at `(s, i) = (3, 4)` is `(4s+4)I₄ − 4J₄`, i.e.
`M(0) = 12`, `M(ē≠0) = −4`. All three nontrivial twists leave it.

**The eigenvalue window does not save it.** `run.py` computes the spectrum of
each `M` over the (real) characters of `Ḡ`. All four give the **same
multiset** `{0, 16, 16, 16}` — rank exactly `s = 3`, trace `4si = 48`, every
nonzero eigenvalue equal to `4i = 16`, which §1.2.1's Theorem E says is the
only value available at any cell. `M' = ψ̄·M` merely relabels the Fourier
coefficients, so **no eigenvalue argument can exclude these instances**; a
forcedness proof at `s ≥ 2` must quotient by the twist, exactly as §1.2 says.

**[5] The conjugation is verified cell by cell.** `S` is built from the banked
rule `sgn[j] = +1` for `j < 4s`, `sgn[4s + I·n + g] = ψ̄(κ(g))`. For each of
the four twists, `run.py` checks
`H'[r][c] = S[r]·S[c]·H[r][c]` at **every one of the 3 671 056 cells** —
`4 × 3 671 056 = 14 684 224` cell comparisons in total — and confirms the row
permutation and the column permutation are both the identity. Each nontrivial
twist negates 952 rows and the same 952 columns and differs from `H` in
exactly `2·952·964 = 1 835 456` cells, so the map is diagonal but far from
trivial.

**[6] Controls.**

| | control | result |
| --- | --- | --- |
| **C1** | the house checker must **reject** the twisted instance: the twisted record is fed back to `tools/bordered_gs.py`'s `check_record` | rejected, and for the right reasons — `H1 PAF profile`, `D5 row-sum law`, `H3`, `Sigma-bar law`. All four are house-form bookkeeping. The assembler is orientation-agnostic, so the rows it returns are still the twisted matrix (digest checked) — the *checker* is what says "not a house instance". |
| **C2** | *negative* control: one cell of the twisted matrix is flipped and the cell-by-cell checker must fire | fires, at exactly `(7, 11)`, and at no other cell |
| **C3** | the conjugation is an involution: `S (S H S) S` | returns `be2073ee…`, the house digest |
| **C4** | the twist is not a relabelling in disguise: exactly `4·(n/2) = 952` of the 1916 rows carry `ψ(g) = −1`, and the same 952 columns; both permutations are the identity | holds |

## Pinned digests

| object | canonical SHA-256 |
| --- | --- |
| house `(3,4)` record | `be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9` |
| twisted sibling (`χ₂`) | `05d411faed301863d1e068651976f2e0f8e200b495af265ec20bc1bd6597175c` |
| `data/n1916-twist.json` | `1a3f92228074f69a7ead11d66371d18dfb39aeb2d17155f3b7fc9782b7b8d51b` |

The two matrix digests are pinned by earlier, independent work: the house one
by the payload replay (cert 01 here; `theorem_check_report.json` upstream),
the twisted one by the upstream skeptic lane's attack 5. `run.py` re-derives
both from `data/payload-records.json` plus the banked twist rule and nothing
else. The `χ₁` and `χ₃` digests are new here and are banked in
`data/n1916-twist.json`.

`data/payload-records.json` is not file-pinned, on purpose: it is shared with
certs 01, 06 and 08, and the binding pin on it is the canonical digest of the
matrix it produces.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, end to end** | **≈ 12.4 s** (69 checks, exit 0) |
| rebuild + full hypothesis re-check + `verify.py`, house | ≈ 1.9 s |
| per twist: general H1–H4 + compression + assemble + `verify.py` | ≈ 1.5 s |
| per twist: the cell-by-cell conjugation check, `1916²` cells | ≈ 0.6 s |

No numpy anywhere; no flag adds any. The heaviest single step is the general
H1 check, which recomputes all four periodic autocorrelations over `G` from
the definition.

## What is NOT claimed

* **Nothing about `ψ(ρ) = −1` in general.** The proposition's proof breaks
  there, and at orders 668, 716 and 1676 that case provably **does** leave
  the equivalence class (certs 06 and 08; cert 11; cert 20). Whether it
  always does is not claimed; 668, 716 and 1676 are the proven cases, and
  1772 stays unclaimed.
* **This certificate does not, by itself, show there are no other non-house
  Grams.** What it does is exhibit the non-house instances that any
  forcedness statement must quotient out. That those are the *only* ones is
  Theorem E and its Corollary E1 (`note/NOTE-B.md` §1.2.1), with the
  complete `(s,i) = (3,4)` classification measured in cert 12. Read
  together, the four Grams swept here **are** the complete admissible list
  at this cell over `Ḡ = Z₂ × Z₂`, so this sweep is an exhaustive
  classification at order 1916 and not merely a sweep of one construction.
  The `s ≥ 2` Gram-forcing question is closed — up to the twist, and the
  twist family is now proved complete rather than assumed to be.
* **No new matrix at 1916 is claimed.** The opposite is proved: all four
  twists give the house matrix up to a diagonal conjugation. Nothing here
  says anything about how many Hadamard equivalence classes order 1916
  carries.
* **No claim about the twist at other orders in this repository's bank.** The
  sweep here is over the character group of `Ḡ` at 1916 only.
