# cert 32 — the good-matrix product theorem: controls and the closure sign

**Label: COMPUTATIONAL-EVIDENCE** beside the **PROVEN** theorem of
`note/NOTE-B.md` §1.13 — *plus* **PROVEN-BY-CERTIFICATE** for the
enumerations themselves, which are complete, exact and unprefiltered.
**This certificate proves no theorem.** The proof is in the note; what is
here is the finite evidence around it, and four controls that can fail.

Default run: `python certs/32-good-product-theorem/run.py` from the
repository root. Standard library only, exact integers only, no floats, no
network, nothing imported or opened outside this repository. **1.2 s, 73
checks, exit 0.** No `--full`: everything is run.

---

## The theorem

Let `n` be odd and let `(A, B, C, D)` be a Goethals–Seidel quadruple on
`ℤ_n` — `Σ_q PAF_q(t) = 0` for every `t ≠ 0` — with `A` of **skew type**
(`a₀ = +1`, `a_{−k} = −a_k`) and `B, C, D` **symmetric**. Then

> **`a_{2k} = −(b₀c₀d₀)·a_k b_k c_k d_k`**  for every `k ≠ 0`.

Under the normalisation `b₀ = c₀ = d₀ = +1` this reads
`a_{2k} = −a_k b_k c_k d_k`. **Without the prefactor the identity is false**,
and this certificate exhibits exactly how false: each unprefixed variant holds
on precisely the half of the solution set where `b₀c₀d₀` has the matching
sign.

## What `run.py` checks (exit 0 iff every check passes)

**[A] the complete small-order enumeration.** Every good quadruple at
`n = 7, 9, 11, 13, 15` is found by meet-in-the-middle over **all**
`2^{(n−1)/2}` skew seeds and **all** `2^{(n+1)/2}` symmetric seeds with `b₀`
**free** — no normalisation, no coset prefilter, nothing assumed about the
shape of a solution.

| `n` | skew seeds | symmetric seeds | good quadruples | with `b₀ = c₀ = d₀ = +1` |
| --- | --- | --- | --- | --- |
| 7 | 8 | 16 | **528** | 66 |
| 9 | 16 | 32 | **288** | 36 |
| 11 | 32 | 64 | **1 440** | 180 |
| 13 | 64 | 128 | **3 456** | 432 |
| 15 | 128 | 256 | **4 224** | 528 |
| | | | **9 936** | 1 242 |

The corrected identity holds in **9 936 / 9 936**. Each unprefixed variant —
`a_{2k} = −a_kb_kc_kd_k` and `a_{2k} = +a_kb_kc_kd_k` — holds in **exactly
half** at every `n`; under the normalisation the `−` form holds in all and
the `+` form in none.

**[A′] the trust chain.** One quadruple per order is assembled into the
standard-orientation Goethals–Seidel array over `ℤ_n` (§1.0, at `κ(ρ) = 0`)
and handed to `verify/verify.py`, which must return exit 0 on a Hadamard
matrix of order `4n` — orders 28, 36, 44, 52, 60. The matrices are written,
verified and **deleted**; nothing is committed.

**[B] the group-ring identity, instantiated.** The proof's central identity
over `𝔽₂[ℤ_n]`,

```
   T^[2] + U = e + ε·J ,   T = U + V + W + Z ,
   ε = 1 + |V| + |W| + |Z|  (mod 2),
```

with `U, V, W, Z` the `{0,1}` indicators of the `−1` positions of `A, B, C, D`
and `X^[2]` the squaring map `(X^[2])_{2a} = x_a`, is re-checked coefficient
by coefficient on **every one of the 9 936**. Its `x^{2k}` coefficient *is*
the theorem. **Control:** on 400 quadruples per order that are **not** good —
random skew/symmetric quadruples with a nonzero aggregate PAF — the identity
must fail, and it does in 341–395 of each 400. It is a consequence of
goodness, not of the seed shapes; a control that could have been vacuous and
is not.

**[C] the Williamson control — it could fail, and it does, where it should.**
With **four symmetric** seeds every `U_q* = U_q`, so the same lemma gives
`T^[2] = e + (σ+1)J` and the conclusion is the **constant-product
corollary**

> `a_k b_k c_k d_k = −(a₀b₀c₀d₀)` for every `k ≠ 0`,

and **not** the doubling relation. Complete meet-in-the-middle over all
symmetric seeds with `b₀` free:

| `n` | 7 | 9 | 11 | 13 | 15 | total |
| --- | --- | --- | --- | --- | --- | --- |
| Williamson quadruples | 960 | 2 112 | 1 920 | 5 184 | 4 608 | **14 784** |
| constant-product corollary | 960 | 2 112 | 1 920 | 5 184 | 4 608 | **14 784** |
| doubling relation | 0 | 0 | 0 | 0 | 0 | **0** |

**Skewness of `A` is load-bearing.** (Recorded so no one is surprised by it:
the *unprefixed* doubling variant `a_{2k} = −a_kb_kc_kd_k` does hold in **48**
of the 960 at `n = 7`, and in none at any larger `n` — a small-`n`
coincidence, and not the theorem's form.)

**[D] the closure sign, computed and never assumed.** Doubling permutes
`ℤ_n ∖ {0}`; along one cycle the theorem determines `A` from a single entry,
and consistency around the cycle is exactly **one** parity condition. In the
multiplier setting — `M₀ ≤ ℤ_n^*` of odd order with the seeds constant on the
orbits of `⟨M₀, −1⟩`, and `2` reaching `⟨M₀, −1⟩` in `L` steps — that
condition reads

```
   Π_j π_j = closure · (−(b₀c₀d₀))^L ,   π_j = b_{2^j} c_{2^j} d_{2^j},
   closure = +1  if 2^L ∈ M₀ ,   −1  if 2^L ∈ −M₀ .
```

Computed here, from the definition of `ℤ_n^*` and nothing else:

* **`n = 7`, `M₀ = {1}`.** `⟨M₀,−1⟩` has 2 elements and 3 orbits on
  `ℤ₇ ∖ {0}`; `L = 3` and `2³ = 1 ∈ M₀`, so `closure = +1` and
  `Π_j π_j = −(b₀c₀d₀)` — `−1` under the normalisation. **That prediction is
  then checked against all 528 enumerated good quadruples at `n = 7`**, with
  the general prefactor: 528 / 528. It is a control, not a restatement.
* **`n = 523`, `M₀` of order 3.** `523` is prime, `|ℤ_523^*| = 522 = 2·3²·29`,
  and there is exactly one subgroup of order 3, `M₀ = {1, 60, 462}`.
  `⟨M₀,−1⟩` has order 6 and **87** orbits on `ℤ_523 ∖ {0}`, all of size 6;
  `2` has order exactly `L = 87` in `ℤ_523^*/⟨M₀,−1⟩`, so doubling is a
  **single 87-cycle** on those orbits; and `2^87 = 463`, with
  `−463 = 60 ∈ M₀`, so `2^87 ∈ −M₀` and `closure = −1`. Hence
  `Π_j π_j = (−1)·(−1)^87 = +1` under the normalisation: **one** parity
  condition, not 87.

## What is **not** claimed

* **The theorem is not machine-checked.** §1.13's proof is a paper proof.
  This certificate instantiates the identity it turns on and exhausts the
  identity's consequences at five small orders; it does not verify the proof.
* **Nothing about `n = 523` beyond the group arithmetic of [D].** No seed
  quadruple is known there, the closure condition is *necessary* for a
  hypothetical one, and no existence statement of any kind follows. In
  particular nothing is claimed about `H(2092)`.
* **Nothing about good matrices at `n > 15`**, and nothing about how many
  there are at any order not enumerated here.
* Nothing about the literature. A product-type relationship for good matrices
  is reported in print (arXiv:1811.05094, whose abstract announces "a new
  relationship between the entries of good matrices" without giving its
  form); the theorem stated and proved in §1.13 is derived in this
  laboratory's chain and is not claimed to be new. See `PROVENANCE.md`.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`** | **1.18 s** (73 checks, exit 0; measured here 2026-09-05, Python 3.14, one worker) |
| [A], the five complete enumerations | 0.2 s |
| [A′], five assemblies through `verify/verify.py` | 0.5 s |
| [B] + [C] + [D] | 0.5 s |

## How to re-run

```
python verify/verify.py --selftest
python certs/32-good-product-theorem/run.py
```
