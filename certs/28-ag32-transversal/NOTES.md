# cert 28 — the `(2,4)` existence theorem (the AG(3,2) transversal census)

**Label: PROVEN** (the `W(8,4)` reformulation, the Walsh classification, the
orthogonality rule and the `2^{3−rk D}` count; `note/NOTE-B.md` §2.5) **+
PROVEN-BY-CERTIFICATE** (the two finite inputs the proof takes from here, and
the transversal census by two constructions of `W`). Default run:
`python certs/28-ag32-transversal/run.py` from the repository root. Standard
library only, **≈ 0.7 s**, exit 0, **14 checks**. **No `--full`**: the census is
exhaustive on the default path — 448 of 448 types, 105 of 105 matchings, 256 of
256 sign functions — and there is no heavier leg to offer.

`note/NOTE-B.md` §2.4 already carried the `(2,4)` border proposition (a),(b)
and the census (c), and cert 18 [4] established (c) by an exhaustive
215 040-class enumeration — a **count**. §2.5 replaces the existence half by a
**theorem**, and this is its certificate. The four `(2,4)` matrices themselves
remain cert 18's: **cited, not re-banked.**

---

## The theorem

> **Theorem ((2,4) existence, `note/NOTE-B.md` §2.5).** Let `(s,i) = (2,4)`,
> `Ḡ = ℤ₄`, `S = {χ, χ³}`, `M = (8,0,−8,0)`. For every `Q′ ∈ H(8)`, every
> admissible S-part (`τ_q = (σ_q(0)−σ_q(2), σ_q(1)−σ_q(3))`, even entries,
> `Σ_q |τ_q|² = 8` — 112 of them) and every `κ(ρ)`, an anti-periodic border kit
> exists. Equivalently: (H4) holds **iff** there is a weighing matrix `W(8,4)`
> whose rows are odd-signed affine planes of `Q′`'s `AG(3,2)` structure, each
> **transversal** to the perfect matching `Π` read off `C*`. Every perfect
> matching of the eight labels admits `2^{3−rk D}` transversal parallel
> classes, `D` its set of four pair-differences. **Hence the border is never
> the obstruction at `(2,4)`, at every `w`.**

## Which inputs are the proof's and which are this certificate's

This split is why the section's label is `PROVEN + PROVEN-BY-CERTIFICATE`
rather than one or the other.

**Proved on paper (§2.5):** the reduction of (H4) to `E Q′ᵀ + P′ C*ᵀ = 0` on
the anti-periodic doubling; `W := (1/4) P′ C*ᵀ` being an integer `W(8,4)` with
`E = −(1/2) W Q′` and `P′ = (1/2) W C*`, and the converse; the Walsh
classification of the 112 admissible signed rows (support an affine plane, an
odd number of minus signs); the orthogonality rule for two signed planes; and
the `2^{3−rk D}` count, with the argument that no odd subset of `D` sums to
zero.

**Carried by this certificate, and not re-derived on paper:** the **448-type
table of the true `Ĉ`**, from which `C* C*ᵀ = 8I` and the two `±2` per row and
column are read; and the step **"`Π` is always a translation"**, verified type
by type.

## Only the rank-1 matchings occur

Of the 105 perfect matchings of the eight labels, only the **7 rank-1** ones —
the translations `k ↦ k ⊕ v`, `v ≠ 0` — actually arise as `Π`. The rank-2 and
rank-3 orbits (42 and 56 matchings, 4 and 2 transversal planes) are **surplus
generality**: this certificate settles them because a uniform statement over all
105 is cheaper to certify than a case distinction, not because the border needs
them. The theorem uses the `rk D = 1`, eight-plane case and nothing else.

## The evidence chain

**[1],[4] `C*` from the actual `Ĉ`** for all `112 × 4 = 448` types:
`C* C*ᵀ = 8I`; exactly two `±2` per row and per column; the column supports form
a perfect matching; and that matching is a **translation** of the label space in
every one of the 448 types — **7 distinct matchings** in all.

**[3] The 256-function Walsh exhaustion.** Over all 256 sign functions on
`𝔽₂³`, the Walsh spectrum is supported on exactly 4 points with values `±4`
**iff** the support is one of the 14 affine planes and the signing has an odd
number of minus signs: `112 = 14 × 8`.

**[6] The transversal census, two ways (D-008).** `|AGL(3,2)| = 1344`; its
orbits on the 105 perfect matchings have sizes **7 / 42 / 56**, matching ranks
**1 / 2 / 3**. *Route A* runs a clique search for a `W(8,4)` per orbit,
independently. *Route B* applies the uniform closed-form `W` — a transversal
plane, its complement, all four sign origins `u` — to **all 105** matchings,
asserting `#classes = 2^{3−rk D}` matching by matching. The two routes agree
everywhere; the closed form for translation by `e₃` is written out.

**[7] End to end, through the trust chain.** For **cert 18's** two `H(88)` seed
quadruples (`n = 20`, `w = 5`, `ρ̄ = 0` and `1`) a `Q′` is drawn with a random
affine structure and random signs; this certificate's own construction produces
`(E, P′)`; the anti-periodic doubling is checked to satisfy (H1), (H3) at
`w = 5` and (H4) **against the true `Ĉ`**; and the matrix is assembled by this
certificate's own block-explicit assembler and handed to `verify/verify.py`,
which returns `HADAMARD`.

**The digests clause [7] produces are this certificate's own, and are not cert
18's pins** — the borders differ, so the matrices differ. Cert 18's four
matrices are the certified `(2,4)` instances and are cited, never re-banked
here.

## The relationship to cert 18

Cert 18 [4]'s exhaustive 215 040-class census is an **independent second proof**
of the same existence statement, by enumeration where this one argues by
structure. It is not re-run here. Cert 18 also carries the four matrices, the
not-a-Kronecker-product test, and Theorem 3's small-case exhaustion; none of
that is duplicated.

C6's repair stands and is stated in §2.4: anti-periodicity of `P` is **forced**
only for `w > 4`. Existence of an anti-periodic kit is proved for **every** `w`,
which is the statement this certificate carries.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **≈ 0.7 s** (exit 0, 14 checks; measured 2026-09-05 on the desk) |
| the 448-type `C*` table | milliseconds |
| the 105-matching census, both routes | milliseconds |
| clause [7], two `H(88)` through `verify/verify.py` | the bulk of the run |

There is no `--full`. Everything the certificate asserts is exhausted on the
default path, and the run says so.

## What is NOT claimed

* **Nothing about `H(2092)`.** The `(2,4)` cell does not land at 2092.
* **Anti-periodicity of `P` is not shown necessary at `w ≤ 4`.** The two
  `H(56)` instances of §2.3 have `w = 3`, where the anti-periodic kit is what
  they use and is sufficient — but is not shown forced.
* **Nothing about the seed layer** at this cell. The theorem says the border
  never obstructs; it says nothing about which seed quadruples exist.
* **The "Theorem AP" of the source laboratory — the same statement for all
  anti-periodic cells — remains a SKETCH** and is not ported.
* **No equivalence statement** among cert 18's four `(2,4)` matrices, or
  against other known `H(56)` and `H(88)`. Cert 26's census profiles three of
  them against their own transposes and nothing else.
* **No novelty of existence** at orders 56 or 88, which were never open.

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/28-ag32-transversal/run.py
```

Standard library only, exact integers only, no network and no numpy. The run
writes its two assembled `H(88)` into a private temporary directory, hands them
to `verify/verify.py`, and deletes them; generated matrices are never
committed. Exit code 0 iff every check passed.

## Provenance

The theorem and this certificate's code are the source laboratory's
(`Hadamard-2060`, `certs/0027-ag32-transversal`, filed by a Cursor cloud lane
under skeptic pass III and adopted as decision **D-068** after a desk replay and
an independent audit that reproduced the orbit-and-rank dictionary — 7 / 42 / 56
and `2^{3−rk D}` — with independent code, and re-assembled the `H(88)`
instances with fresh borders). The mathematics is carried across unchanged; what
changed in transit is the house framing and the scratch directory, which is now
a private temporary directory rather than a path inside the repository. Credit
is to stations, as everywhere in this repository.
