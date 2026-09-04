# cert 27 — Sylow-2 rigidity of a cocyclic `H(4t)`, `t` an odd prime `> 7`

**Label: PROVEN** (the derivation; `note/NOTE-B.md` §1.10) **+
PROVEN-BY-CERTIFICATE** (the rule against existence, two implementations; the
`t = 3` census, two predicates) **+ CITED** (de Launey–Flannery–Horadam 2000,
for "RDS ⟺ cocyclic Hadamard matrix", which is **not** re-derived). Default
run: `python certs/27-sylow2-rigidity/run.py` from the repository root.
Standard library only, **≈ 12.6 s**, exit 0, **14 checks**. `--full` widens the
equation box and re-runs the existence sweep with no early exit: **≈ 13.9 s**,
**run in this repository**.

This is the first statement in `note/NOTE-B.md` about the **cocyclic** shape of
a hypothetical `H(2092)`. It shuts ten of the twelve doors and says nothing
about the two it leaves open.

---

## The theorem

> **Theorem (Sylow-2 rigidity, `note/NOTE-B.md` §1.10).** Let `t` be an odd
> prime `> 7` and let `E` be a group of order `8t` with a central involution
> `e*`, carrying a `(4t, 2, 4t, 2t)` relative difference set relative to
> `⟨e*⟩`. Then the Sylow-2 subgroup `P` of `E` is **not** `ℤ₈` and **not**
> `D₄`. `P = ℤ₂³`, or `ℤ₄×ℤ₂` with `e*` a non-square, forces `t` a **square**.
> `ℤ₄×ℤ₂` with `e*` the square forces `t` a **sum of two squares**. At
> `t = 523`: `P = Q₈`, and `E ∈ {Q₈ × ℤ₅₂₃, Q₄₁₈₄}`.

**The bookkeeping, repaired.** The extension is **by `ℤ₂`**, so the extension
group has order `2N = 8t`; at `N = 2092` that is **4184**, not `4·2092 = 8368`.
There are **twelve** groups of order `8t` of this shape — five direct products
and seven semidirect, the latter from the Aut-orbits of index-2 subgroups of
the five groups of order 8 (orbit counts `1, 2, 1, 2, 1`) — and **twenty-four**
`(E, e*)` pairs, counting central involutions of `P` in the kernel of the
action: `1+1+3+1+3+7+3+1+1+1+1+1 = 24`. Two survive at `t = 523`, one central
involution each.

**Why the Sylow structure is forced at all.** `n_t | 8` and `n_t ≡ 1 (mod t)`
force the Sylow-`t` subgroup normal, so `E = ℤ_t ⋊ P`; the action
`P → Aut(ℤ_t) = ℤ_{t−1}` has image of order dividing `gcd(8, t−1)`, and at
`t = 523`, `gcd(8, 522) = 2`, so the action is through `ℤ₂` at most.

## The proof's engine, and what this certificate carries

The proof (in the note) projects the row-sum identity `g g⁽⁻¹⁾ = 8t(1 − e*)`
along `ℤ_t` into `ℤ[P]`. With `f(q) = Σ_{h ∈ ℤ_t} g(h,q)` — odd-valued, since
it is a sum of `t` signs — one gets `f(q e*) = −f(q)` and, in `ℤ[P]`,

```
Σ_q f(q)² = 8t,        Σ_q f(q) f(m⁻¹q) = 0   for m ∉ ⟨e*⟩,
```

i.e. `a² + b² + c² + d² = 4t` on coset representatives of `⟨e*⟩`, together with
cross conditions read off the **algebraic role** of the elements — labelling-free,
which is what makes the five cases a classification rather than a computation.
That derivation is **paper-grade and lives in the note**. This certificate
carries the finite content around it.

## The evidence chain

**[A] Group facts.** `523` is prime, `≡ 3 (mod 4)`, neither a square nor a sum
of two squares; `n_523 | 8` and `n_523 ≡ 1 (mod 523)` force `n_523 = 1`;
`gcd(8, 522) = 2`; `|E| = 2N = 4184`; the five groups of order 8; the Aut-orbits
of their index-2 subgroups; **12** groups of order `8t`; **24** `(E, e*)` pairs;
exactly two with Sylow-2 `= Q₈`.

**[B] Two implementations of the projected system (D-008).** *Route A* expands
`Σ_q f(q) f(m⁻¹q)` brutally in the group table; *route B* evaluates the five
hand-derived cross conditions. They agree on **all 4096** odd 4-tuples in
`[−7,7]⁴` (default) or `[−9,9]⁴` (`--full`), for **every** `(P, e*)` pair. The
two routes share no code: one knows only the multiplication table, the other
only the five algebraic conditions.

**[C] The rule against existence.** Existence of an odd solution equals the
theorem's rule at **every odd `t ≤ 201`** and at `t = 523` — **1326 cells**. At
`t = 523` only `Q₈` admits solutions: **8384** of them, exactly the ordered odd
four-square representations of 2092.

**[D] The `t = 3` census — a control, not an instance.** All **12** groups of
order 24 with a normal `ℤ₃`, every central involution, all `2¹²` transversals,
under **two predicates** (D-008): the difference-count definition of an RDS, and
the group-ring identity `g g⁽⁻¹⁾ = 8t(1 − e*)` in `ℤ[E]`. They agree everywhere.
RDS exist only in `Q₈ × ℤ₃` (**192**) and the dicyclic `Q₂₄` (**576**), and each
develops into a `12 × 12` Hadamard matrix `H[x,y] = g(xy⁻¹)`, checked here.

**`t = 3` is outside the theorem's range** (`t` an odd prime `> 7`). It is not
an instance of the theorem and is not offered as one. It is a **consistency
control on the machinery**: it exhibits the objects the projected system is
about, and it shows the two predicates agree on a case where the objects
actually exist — the projected system alone would be a screen with nothing to
screen.

## The cross conditions, case by case

| `P` | `e*` | cross conditions | obstruction |
| --- | --- | --- | --- |
| `ℤ₈` | the unique central involution | `ab + bc + cd − da = 0` | impossible for odd values (`≡ 2 mod 4`) |
| `D₄` | the unique central involution | `ac = bd`, `ad = −bc` | `cd(a²+b²) = 0` |
| `ℤ₂³` | any of the 7 | `ab = −cd`, `ac = −bd`, `ad = −bc` | `a² = b² = c² = d²`, so `t` a square |
| `ℤ₄×ℤ₂` | the square | `ac + bd = 0` | `a + bi = λ i (c + di)` ⇒ `t` a sum of two squares |
| `ℤ₄×ℤ₂` | a non-square | `(a₀+a₂)(a₁+a₃) = 0`, `a₀a₂ + a₁a₃ = 0` | `4t = (a₁+a₃)²` on one branch and `4t = (a₀+a₂)²` on the other — `t` a square either way |
| `Q₈` | the unique central involution | identities | only `a²+b²+c²+d² = 4t` |

At `t = 523` — prime, `≡ 3 (mod 4)`, not a square, not a sum of two squares —
only `Q₈` survives.

## What is CITED versus what is re-derived

**CITED, not re-derived:** the equivalence *"a `(4t,2,4t,2t)`-RDS in a central
extension `E` of `ℤ₂` by a group of order `4t`, relative to the central `ℤ₂`, is
the same thing as a cocyclic Hadamard matrix of order `4t` with extension group
`E`"* — de Launey–Flannery–Horadam 2000, with Flannery 1997. This certificate
constrains the **extension group** *on the assumption of* that equivalence, and
would say nothing without it.

**Re-derived here:** the projected row-sum system, the five cases, the group
census, the 12-group / 24-pair bookkeeping, the `|E| = 2N` correction, and the
`t = 3` control.

## Runtimes

| step | cost |
| --- | --- |
| **`run.py`, default path** | **≈ 12.6 s** (exit 0, 14 checks; measured 2026-09-05 on the desk) |
| **`run.py --full`** | **≈ 13.9 s** (**run here**, 2026-09-05) |
| [D], the `2¹²` transversal census over 24 `(E, e*)` pairs | the bulk of the run |

`--full` is a **wider run of the same standard-library code**, not a different
arithmetic — the same relation cert 16's `--wide` has to its default path. It
widens [B]'s equation box from `[−7,7]⁴` to `[−9,9]⁴` and re-runs [C] with no
early exit, counting every solution at every odd `t ≤ 51` as well as at 523.

## What is NOT claimed

* **Nothing about the existence of a cocyclic `H(2092)`.** The two surviving
  groups, `Q₈ × ℤ₅₂₃` and `Q₄₁₈₄`, are open doors. This certificate shuts the
  other ten.
* **The DLFH equivalence is not re-derived** — see above. Every statement here
  is conditional on it.
* **Nothing about the Goethals–Seidel array.** The cocyclic question and the
  coset-border construction of §1 are different objects; §1.11's Lemma F
  records that the GS array is not of the 16-block circulant shape a `p = 523`
  automorphism would force, which is a separate matter.
* **Nothing about Ito's constructions, Williamson-type seeds, or the
  Butson/quasi-LP material** that survived the same pass in the source
  laboratory. Those are not theorems of this note.

## How to re-run

From the repository root, on bare Python 3.9 or newer:

```
python verify/verify.py --selftest
python certs/27-sylow2-rigidity/run.py
python certs/27-sylow2-rigidity/run.py --full
```

Standard library only, exact integers only, no network and no numpy on any
path. Exit code 0 iff every check passed.

## Provenance

The theorem, the projected system and this certificate's code are the source
laboratory's (`Hadamard-2060`, `certs/0026-sylow2-rigidity`, filed by a Cursor
cloud lane under skeptic pass III and adopted as decision **D-068** after a desk
replay and an independent audit that re-derived the case analysis and the
12-group / 24-pair bookkeeping by hand). The mathematics is carried across
unchanged; what changed in transit is the house framing — the docstring, the
check count, the verdict block and this file. Credit is to stations, as
everywhere in this repository.
