# cert 17 — the `(3,4)` border-kit census, and Theorem F on the banked records

**Label: PROVEN** (Theorem F is a paper-grade proof, `note/NOTE-B.md` §1.8)
**+ PROVEN-BY-CERTIFICATE** (the S-part census, the shell arithmetic, the kit
census, the from-scratch instances) **+ MEASURED** (Theorem F literally on
7 / 7 banked coset-border records). Default run:
`python certs/17-border-kits-34/run.py` from the repository root. Standard
library only, exact integers only. **8.6 s, 46 checks, exit 0.** The default
path censuses a fixed deterministic **256-class sample**; `--full` censuses
all **16 384** classes.

---

## The statements

> **Theorem F (structure).** At `i = s+1`, house Gram `M = (4s+4)I − 4J`,
> `w > 2s`: `QᵀQ = 4i·I_{4s}`; `P(I₄⊗P_{S^c}) = 0` and `PᵀP = QQᵀ = I₄⊗M`;
> `ĈᵀĈ = ĈĈᵀ = I₄⊗Σ̄`; `E = −(1/4i)·P·Ĉᵀ·Q`; and given `P, Q` with those
> properties, `EEᵀ = 4s·I` is automatic and (H3) is an identity — so
> **(H3)+(H4) reduce to the single condition `E ∈ {±1}`**. (H4) depends on
> the seeds only through the **S-part** of the coset-sum table, whose kernel
> is the per-seed constants; and a kit transports to any `(G′,K′,ρ′,x′)` with
> the same `Ḡ`, the same `κ(ρ)` and the same S-part, **at any `w`**.

> **Theorem F (the `(3,4)` census).** `Ḡ ∈ {ℤ₄, ℤ₂²}`, `w` even. (i) exactly
> one seed is coset-balanced and the other three have
> `|σ̂_q(χ)| = 2` at every nontrivial `χ` — the *silent-seed lemma*;
> (ii) exactly **2048** admissible S-parts per group; (iii) **every** S-part
> admits a kit at **every** `κ(ρ)`; (iv) at `N = 2092` the coset-sum tables
> form a shell of `4192 × 64 = 2048 × 131 = 268 288` tables per group.

> **Corollary (the `(3,4)` cell is one-layer).** Every house-profile seed
> quadruple on an abelian group of order `4w`, `w > 6`, with an index-4
> subgroup extends to a Hadamard matrix of order `16w + 12`. At `N = 2092`
> (`w = 130`): a seed quadruple **is** `H(2092)`.

The proofs are in `note/NOTE-B.md` §1.8. This certificate carries the finite
parts.

## What `run.py` checks (exit 0 iff every check passes)

**[0]** `data/payload-records.json`, SHA-256 pinned in `run.py`.

**[A] the S-part census, the silent-seed lemma, the shell.** All spectra at
the three nontrivial characters with `Σ_q |σ̂_q(χ)|² = 4s = 12` (Corollary E2)
and entries in `2ℤ` are enumerated by **brute force**; those for which
`σ = (r + T)/4` is integral with the right parity for some integer row-sum
vector are the S-parts. **2048 per group**, both groups, and every one has
the silent-seed shape. The `σ`-shell at `N = 2092` — ordered even `r` with
`Σ r_q² = 2092` — has **4192** members (`= r₄(523)`, Jacobi), each with
exactly one coordinate `≡ 0 (mod 8)` and three `≡ 2 (mod 4)`; each S-part is
compatible with exactly **131** of them and each shell vector with exactly
**64** S-parts, giving **268 288** tables per group; and the silent seed is
exactly the seed whose row sum is `≡ 0 (mod 8)`.

**[B] the kit census.** For every S-part in the sweep and **every** `κ(ρ)` —
all four per group, no symmetry reduction — `kitlib.py` finds `(E,P,Q)`, and
every kit is re-verified exactly with a **full `σ`** (r-part included) at
`w = 130`: (H1) `QQᵀ = I₄⊗M`, `PPᵀ = 16I`, `E ∈ H(12)`, (H3)
`EEᵀ + wPPᵀ = NI`, (H4) `EQᵀ + PĈᵀ = 0`. Default: every 64th S-part, both
groups, all four `κ(ρ)` — **256 / 256**, in 6.0 s, at most **60** of the 301
candidate `Q`s tried before success. Census digests pinned in `pins.json`;
the digest covers `(group, κ(ρ), S-part, found)`, not the kit itself.

**[C] controls, through `verify/verify.py`.** From-scratch `(3,4)` instances
at `n = 8` (`w = 2`, below `2s`): the house-profile quadruples are found by a
**full meet-in-the-middle over all 256 sequences** on `ℤ₈` (`ρ = 0` and
`ρ = 7`) and on `ℤ₂×ℤ₄` (both index-4 quotients), the kits by this engine,
the matrices by this certificate's own assembler — **four `H(44)`, all
green**. `ℤ₂³` with `|K| = 2` carries **no** house quadruple, and the same
exhaustive search says so. Then `H(124)` on `ℤ₂×ℤ₂×ℤ₇` (`n = 28`, `w = 7`)
**twice**: once with the order-1916 record's kit transported verbatim
(Theorem F(g) — same `Ḡ`, same `κ(ρ)`, same S-part, a different `w`), at the
pinned digest `dc3fe6b3…c9e7`; once with this engine's own kit, a different
artifact, also green.

**[D] Theorem F on the banked records.** Every record in
`data/payload-records.json` with `s ≥ 1` and `w > 2s` — **seven**: 668, 716,
1676, 1772 at `(1,1)`; 1916 at `(3,4)`; 1388 at `(5,6)`; 1436 at `(7,8)` —
has `QᵀQ = 4i·I`, `PᵀP = QQᵀ = I₄⊗M`, `ĈᵀĈ = ĈĈᵀ = I₄⊗Σ̄`,
`E = −(1/4i)·P·Ĉᵀ·Q` and (H4) re-multiplied entry by entry; and at `i = s+1`,
`P` is superblock-balanced and the (H4)-kernel is exhibited — shifting `σ` by
per-seed constants leaves `PĈᵀ` unchanged. At `(1,1)` the kernel is
`0`-dimensional and the row sums enter (H4) directly (Theorem D's (D-e)), so
Theorem F(f) is vacuous there and the run says so.

## The engine (`kitlib.py`), and why it is a second implementation

The source laboratory's engine backtracks over 12-cliques of `P`-rows,
intersecting compatibility masks, then looks up a `Q`-clique. This engine
inverts the order: it **fixes `Q`** from a list of 301 candidate 12-cliques
(the Kronecker table `H₄ ⊗ [b₁ b₂ b₃]` first, then random cliques from a
fixed seed), computes `V = (4Ĉ⁰)ᵀQ` once, filters the 648 block-balanced
sign representatives `p` by `p·V ∈ {±64}¹²`, and finds a 12-clique among the
survivors. Both rest on the same theorem: under Theorem E at `i = s+1`,
(H1) ⟺ `Q` has 12 orthogonal block-balanced columns, `P` likewise,
`E = −(1/16)PĈ⁰ᵀQ`, `EEᵀ = 12I` automatic, and the only open condition is
`E ∈ {±1}`. The GS array over the quotient, the group labelling, the S-part
parametrisation and the assembler are all written in this certificate;
`kitlib.py` imports nothing but `itertools` and `random`.

## Pinned digests

| what | SHA-256 |
| --- | --- |
| `data/payload-records.json` (read by clause [D]) | `9afbf60efed82e97ffe229e0059f06dfa7849a2ac2434c52a6ada09559738afb` |
| `H(124)` with the transported order-1916 kit (canonical) | `dc3fe6b38d9007db9261b5739da13baaffd7a8f7dc316918e203abf6a169c9e7` |
| sample census (256 classes) | `47d5d44e8e7087eeba199b720dfba84cc97749f2326d2e219f5f1f44cbe36c35` |
| full census (16 384 classes) | `d60b83d8a121cb5926c67503463356edb733f64e5deb4c34b13acba566497f4d` |

The `H(124)` built with this engine's own kit has canonical digest
`133b3f96d744e202…` in the run above; it is deliberately **not** pinned — a
different candidate-`Q` list would produce a different (equally valid)
border, and pinning it would bind the certificate to the search order rather
than to the statement.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **8.6 s** (46 checks, exit 0; measured here 2026-09-02, Python 3.14, one worker) |
| clause [A], the two S-part censuses | 0.9 s |
| clause [B], the 256-class sample | 6.0 s |
| clause [C] + [D] | 1.7 s |
| **`run.py --full`** | **5 min 28 s** (46 checks, exit 0; measured here 2026-09-02, Python 3.14, one worker). The 16 384-class census itself is 325.8 s of that — the eight `(group, κ(ρ))` sweeps at 39.8–41.8 s each — and its digest matched the pin. The source laboratory measured 4 min 2 s for the same census on 2026-09-02. |

## What is NOT claimed

* Nothing that rests on the pin alone. The `--full` census **has** been run
  in this repository (2026-09-02, Python 3.14, one worker, 5 min 28 s):
  2048 / 2048 for each of the eight `(group, κ(ρ))` sweeps, 16 384 / 16 384,
  digest `d60b83d8…7f4d` matching the pin, 46 checks, exit 0. The pin is
  therefore no longer a travelling number — it is a reproduction. Anyone who
  runs `--full` and gets a different digest has found a real discrepancy and
  should say so.
* Nothing about the **seed layer**. The census says the border never
  obstructs at `(3,4)`; a house-profile seed quadruple on a group of order
  520 is still the whole problem at 2092, and none is known.
* Nothing about the **number or structure** of kits per S-part — existence
  only. A different candidate-`Q` list finds different kits and the same
  statement.
* Nothing at `w ≤ 2s` beyond the `H(44)` controls: those show the same kit
  shape *works* at `w = 2` (sufficiency), not that every kit at `w ≤ 2s` has
  this shape. Theorem F(b) is a necessity statement and uses `w > 2s`.
* No novelty is claimed for the `(3,4)` mechanism itself; the source
  laboratory computed these kits first. This certificate is the second
  implementation the house's two-implementation rule requires before the
  claim enters `note/NOTE-B.md`.
* The `H(124)` parameters are **embedded in `run.py`**, not banked in
  `data/`: they are not a public record but a control this certificate
  constructs, and the transport from the order-1916 record is exactly what
  clause [C] is testing.

## How to re-run

```
python verify/verify.py --selftest
python certs/17-border-kits-34/run.py
python certs/17-border-kits-34/run.py --full
```
