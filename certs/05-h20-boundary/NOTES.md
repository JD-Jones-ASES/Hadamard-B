# cert 05 — the hypothesis boundary `w ≤ 2s`, at order 20

## What the two instances are

Both live on `G = ℤ₂×ℤ₂` with `K = ⟨(1,1)⟩`, the **diagonal** subgroup;
`s = 1`, `i = 2`, `w = 2`, `N = 20`. Three things put them outside
everything the twelve decoded records exercise:

1. **`K` is not a coordinate-divisor kernel.** `tools/bordered_gs.py`'s
   `coset_map` only builds kernels of mixed-radix reductions, so it
   cannot even *express* this record. This cert therefore carries its own
   stdlib implementation of the ansatz, written against the displayed
   definitions of NOTE-B.md §1.0–§1.1, for an arbitrary `K` given by
   generators.
2. **`w = 2 = 2s` exactly**, so Theorem C's non-degeneracy hypothesis
   `w > 2s` **fails**.
3. **T2 uses the transpose-negated orientation** (the six transposed
   blocks of the Goethals–Seidel array negated), which the house checker
   refuses outright (`gs_variant != "standard"`).

## The claim

> Both matrices are Hadamard (`H Hᵀ = 20·I`), matching the digests
> pinned below, and both satisfy Theorem A's (H1)–(H4) with the house
> Gram `M = 8I₂ − 4J₂` — with an **arbitrary** subgroup `K` and in
> **either** orientation. So no unstated regularity of the decoded
> records is doing secret work: the stated hypotheses suffice on their
> own.
>
> At `w = 2s` exactly, D3's **forcing argument lapses**. D3 argues:
> off-diagonal, `(E Eᵀ)[a][b] = −w·(P Pᵀ)[a][b]`; entries of `P Pᵀ` are
> sums of `4i` signs, hence even; entries of `E Eᵀ` are bounded by `4s`;
> so `w > 2s` gives `|w·(even)| ≤ 4s < 2w` and forces the even number to
> be `0`. At `w = 2s` the bound reads `4s = 2w` and the step fails —
> `(P Pᵀ)[a][b] = ∓2` with `(E Eᵀ)[a][b] = ±4s` is no longer excluded by
> the inequality.
>
> **Measured**, on both banked instances: `E Eᵀ = 4·I₄` and
> `P Pᵀ = 8·I₄ = 4i·I₄`. **D3's conclusion holds anyway.**
>
> And that is not an accident of the search that found them. A
> **complete finite check**, run inside this cert, shows that at
> `(s, i, w) = (1, 2, 2)` the (H3) system `E Eᵀ + w·P Pᵀ = N·I₄` has
> **no solution at all** with a non-Hadamard corner. At these parameters
> `w > 2s` is therefore **sufficient but not necessary** for D3's
> conclusion: crossing the boundary costs the *proof*, not — here — the
> *statement*.

**This is the documented edge of Theorem C, and nothing more.** No
counterexample to Theorem A, to Theorem C, or to D3 is claimed; none was
found; the skeptic lane that built these instances found none either.
Order 20 is long settled: **no novelty of existence and no priority of
any kind is claimed.**

## A correction to a natural misreading

These two instances are sometimes described as showing that "the corner
is not Hadamard" or that "D3's conclusion fails" at the boundary. **They
do not.** Both corners *are* Hadamard, `E Eᵀ = 4I₄` on both, and the
exhaustive check above says no admissible corner at these parameters
could have been otherwise. What fails at `w = 2s` is D3's *inequality*,
i.e. its proof. `run.py` asserts the measured statement programmatically
(a non-Hadamard corner here would be a **hard failure**, not a silent
change), so this NOTES.md cannot drift away from what the code finds.

## How the exhaustive boundary check is exhaustive

(H3) sees `E` only through `E Eᵀ` and `P` only through `P Pᵀ`, so
negating a column of `E`, negating a column of `P`, and permuting the
`4s` rows of `(E, P)` jointly are all free. A non-Hadamard `4×4` sign
corner has some off-diagonal `(E Eᵀ)[a][b] ≠ 0`; for `±1` rows of length
4 the only nonzero values are `±4`, i.e. `row_b = ±row_a`. Relabel that
pair to rows 0 and 1, normalise row 0 of `E` and row 0 of `P` to
all-plus, and enumerate `e₁ ∈ {+e₀, −e₀}`, every `p₁`, and every
completion to four rows. That is what `run.py` does: 112 normalised
branches survive at `(1,2,2)`, none completes to a non-Hadamard corner.

The control at `(s, i, w) = (1, 2, 3)` — the same cell with `w > 2s`,
where D3 proves it — is infeasible at the very first step (0 branches):
`w·(P Pᵀ)[0][1] = ∓4` has no integer solution when `w = 3`. That is D3's
forcing in miniature, and it shows the boundary case is the
non-vacuous one.

## Honesty labels

| part of the claim | label |
| --- | --- |
| both matrices are Hadamard and match their pinned digests | **PROVEN-BY-CERTIFICATE** |
| both satisfy Theorem A (H1)–(H4) with an arbitrary `K` and in either orientation | **MEASURED** (exact integer arithmetic) |
| `w = 2s` exactly, so D3's hypothesis fails; `E Eᵀ = 4I₄`, `P Pᵀ = 8I₄` | **MEASURED** |
| at `(s,i,w) = (1,2,2)` no (H3) solution has a non-Hadamard corner | **PROVEN** by a complete finite check (exhaustive up to the exact symmetries of (H3), as described above) |
| whether a non-Hadamard corner is possible at any *other* `w ≤ 2s` parameters | **not claimed** — the check above covers `(1,2,2)` and nothing else |

## Data provenance

- `data/h20-boundary.json` — one record per instance:
  `group`, `K_generators`, `orientation`, `r_shift`, `seeds`, `corner`,
  `row_table`, `col_rows`, `pinned_sha256`. Note that `col_rows` is the
  `4i × 4s` matrix `Q` written out **row by row**, *not* the transposed
  `col_table` of the `tools/bordered_gs.py` record format — that format
  cannot express this instance at all.
- **Source:** the source lab repository's skeptic lane,
  `experiments/bordered_gs_theorem/skeptic/attack3_outside.py`, tests T1
  and T2, whose artifacts are `H20_T1-diagK-w2s.txt` and
  `H20_T2-flipped.txt` and whose write-up is the attack-3 section of
  `skeptic/LANE-D-SKEPTIC-REPORT.md`. The search was pre-registered in
  that script's docstring: family, exhaustive enumeration over the 1 536
  seed quadruples meeting the profile, 400 admissible `Q`, all
  `(E-row, P-row)` pairs by exact hash join; kill criterion exhaustion.
- **How the parameters were banked:** the T1 and T2 searches were re-run
  in a scratch copy of that lane to recover the generating parameters,
  and the rebuild from those parameters was confirmed **byte-identical**
  to the lab's two `H20_*.txt` artifacts, with digests equal to the two
  values quoted in the skeptic report. `run.py` rebuilds from the banked
  parameters alone.

## How to re-run

From the repository root, on bare Python 3.9 or newer, stdlib only:

```
python verify/verify.py --selftest
python certs/05-h20-boundary/run.py
```

Runtime ≈ 0.8 s (of which the exhaustive boundary check is ≈ 0.7 s).
Matrices are assembled into a temporary directory outside the
repository, verified, digest-compared, and deleted.

Exit code 0 iff both instances rebuilt, verified, matched their digests,
returned the measured Grams above, and the exhaustive boundary check
returned "none exists".

## Pinned digests

| instance | orientation | `K` | `(n, s, i, w)` | `E Eᵀ` | `P Pᵀ` | canonical SHA-256 |
| --- | --- | --- | --- | --- | --- | --- |
| `T1-diagK-w2s` | standard | `⟨(1,1)⟩` | (4, 1, 2, 2) | `4·I₄` | `8·I₄` | `50eecc761e12b76944b301b7aaeb03a61cb6b88cfc52c67caaacf20eef0e6c9b` |
| `T2-flipped` | transpose-negated | `⟨(1,1)⟩` | (4, 1, 2, 2) | `4·I₄` | `8·I₄` | `f279337f4c2376dfaff2f5ec82cec21a055a1888008c1972763c0f6721066caf` |

The comparison is **coded** into `run.py` (`PINNED`), and is also
cross-checked against each record's `pinned_sha256`; a mismatch on
either is a hard failure with a non-zero exit code.
