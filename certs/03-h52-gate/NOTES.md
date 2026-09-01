# cert 03 — the Theorem-D gate: `H(52)` on `G = ℤ₂×ℤ₂×ℤ₃`

## The claim

> The banked parameter record is a valid `s = 1, i = 2` bordered
> Goethals–Seidel instance on the **non-cyclic** group
> `G = ℤ₂×ℤ₂×ℤ₃` (`n = 12`, `N = 4(n+1) = 52`), with
> `K = ker(t ↦ t₀)` of index 2 (`w = 6 > 2s`) and `ρ = (0,0,0)`, hence
> `κ(ρ) = 0` and **`ε = +1`**. It satisfies every hypothesis of
> Theorem A/B and every derived law of Theorem C; it satisfies every
> clause of **Theorem D** — (D-a) the Gram is the house form
> `M = 8I₂ − 4J₂` (the genuine `i = 2` branch, not the degenerate
> `M = 4J₂`), (D-b) `Q` is the ±-doubling of a 4×4 Hadamard `U`,
> (D-c) `P[r][2J+1] = −P[r][2J]` with `p` and `E` both 4×4 Hadamard,
> (D-d) `E = −¼·p·Λ(d)ᵀ·U` exactly, with the forced Parseval value
> `Σ_q δ_q² = 4`; and it assembles to a Hadamard matrix of order 52
> whose canonical SHA-256 is the digest pinned below.

Measured values: `δ = (2, 0, 0, 0)`, `d = (2, 0, 0, 0)`,
`Σ_q δ_q² = 4`.

## Why this artifact exists — it is a gate, not a new order

**Order 52 is long settled. No novelty of existence is claimed, and no
priority of any kind is claimed.** The artifact's entire role is to be a
*gate*: an instance that is **new** — not a decoded record, not a twist
of one — exercising the branch of Theorem D that the four decoded
`s = 1` instances do not reach. Those four are all cyclic with `ρ` odd,
i.e. `ε = −1` and `G` with a unique index-2 subgroup. This gate is
non-cyclic (so `G` has three subgroups of index 2, and `K` is a choice)
and `ε = +1`. The whole point of building it was to run the derived
`s = 1, i = 2` border system end to end on something the derivation had
not already been fitted to.

## Pre-registration

The gate was **pre-registered before the search ran**, under the lab's
No-Noise Law: family, group, subgroup, reflection, orientation, symmetry
class (none — the raw `2⁴⁸` seed space), enumeration method, budget,
kill criterion, the three acceptance conditions, and the exact meaning
of a negative outcome were all fixed in advance.

- Registration document: the source lab repository's
  `experiments/pr0038/GATE_REGISTRATION.md`, **flushed to disk
  2026-08-31T20:53:28Z**, before `i2_border.py gate` was executed for
  the first time.
- Result document: `experiments/pr0038/gate.json`, written
  **2026-08-31T21:09:30Z** — sixteen minutes later.

The registration states the negative outcome in advance: exhaustion with
no quadruple would have been a **BOUNDED-NEGATIVE-SEARCH** closing
exactly the cell `(s=1, i=2, G = ℤ₂×ℤ₂×ℤ₃, K = ker(t↦t₀))` and nothing
else. The outcome was positive.

## Honesty labels

| part of the claim | label |
| --- | --- |
| the record satisfies H0–H4, Theorem C's derived laws, and Theorem D's clauses (D-a)–(D-d) | **MEASURED** (exact integer arithmetic) |
| the assembled `H(52)` is Hadamard and matches its pinned digest | **PROVEN-BY-CERTIFICATE** |
| Theorem D itself | **PROVEN** (paper-grade; NOTE-B.md §1.5). This cert is a gate on the theorem, not its proof. |

## Data provenance

- `data/h52-gate.json` — the complete record, copied verbatim from the
  `params` block of the source lab repository's
  `experiments/pr0038/gate.json`, plus a `pinned_sha256` field and a
  `provenance` block naming the registration document and the search.
- **Search:** exhaustive meet-in-the-middle over the raw `2⁴⁸` seed
  space of `G`, shell-restricted by the condition `Σ_q δ_q² = 4` that
  Theorem D **forces** (exactly one `δ_q = ±2`, three zero; the sign is
  a gauge, so only `+2` was enumerated). Border solved as
  `E = −¼·p·Λ(d)ᵀ·U` over all 768×768 ordered pairs of 4×4 Hadamard
  `(p, U)`, first `±1` corner taken.
- **Checker:** `tools/bordered_gs.py`, imported unmodified. Theorem D's
  clauses are re-derived in this cert's `run.py`, not read off any
  report.

## How to re-run

From the repository root, on bare Python 3.9 or newer, stdlib only:

```
python verify/verify.py --selftest
python certs/03-h52-gate/run.py
```

Runtime ≈ 0.2 s. The matrix is assembled into a temporary directory
outside the repository, verified, digest-compared, and deleted.

Exit code 0 iff every check passed and the digest matched.

## Pinned digest

| artifact | order | `(n, s, i, w)` | `ε` | canonical SHA-256 |
| --- | --- | --- | --- | --- |
| `H(52)` gate, `G = ℤ₂×ℤ₂×ℤ₃`, `K = ker(t↦t₀)` | 52 | (12, 1, 2, 6) | +1 | `e2c3e48b0fc65f5283e833096824b4fec651d8c57694ae45b3842c23c87ad7ca` |

The comparison is **coded** into `run.py` (`PINNED_SHA`), and is also
cross-checked against the `pinned_sha256` field of the banked record; a
mismatch on either is a hard failure with a non-zero exit code.
