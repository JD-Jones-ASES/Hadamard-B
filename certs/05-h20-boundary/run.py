#!/usr/bin/env python3
"""cert 05 -- the two order-20 instances at the hypothesis boundary w = 2s.

Run from the repository root:

    python certs/05-h20-boundary/run.py

WHAT THESE TWO INSTANCES ARE

Both live on `G = Z2 x Z2` with `K` the DIAGONAL subgroup <(1,1)>,
`s = 1`, `i = 2`, `w = 2`, `N = 20`.  Two things about that are outside
everything the twelve decoded records exercise:

  * `K` is not the kernel of a coordinate-divisor map, so
    `tools/bordered_gs.py` cannot even express the record -- its
    `coset_map` only builds kernels of mixed-radix reductions.  This
    cert therefore carries its own stdlib implementation of the ansatz,
    written against the displayed definitions, for arbitrary `K` given
    by generators.
  * `w = 2 = 2s` EXACTLY, so the non-degeneracy hypothesis `w > 2s` of
    Theorem C FAILS.  D3's forcing argument -- "entries of `P P^T` are
    even and `|E E^T|_off <= 4s`, so `w > 2s` makes
    `|w * (even)| <= 4s < 2w` force that even number to be 0" --
    degenerates: at `w = 2s` the bound reads `4s = 2w`, and an
    off-diagonal `(P P^T) = -+2` with `(E E^T) = +-4s` is no longer
    excluded by the inequality.

  T2 additionally uses the TRANSPOSE-NEGATED orientation (the six
  transposed blocks of the Goethals-Seidel array negated), which the
  house checker refuses outright.

WHAT THIS CERT DEMONSTRATES

  1. Both matrices are Hadamard (`H H^T = 20 I`), by the trust chain,
     matching the digests pinned below.
  2. Both satisfy Theorem A's (H1)-(H4) with the house Gram, with an
     ARBITRARY subgroup K and in EITHER orientation.  So no unstated
     regularity of the decoded records is doing secret work: the stated
     hypotheses suffice on their own.
  3. `w = 2s` exactly, so D3's HYPOTHESIS fails -- and yet the MEASURED
     corner Grams are `E E^T = 4 I_4` and `P P^T = 8 I_4 = 4i I_4`:
     D3's CONCLUSION holds at these two instances anyway.
  4. That is not an accident of the search.  A complete finite check,
     run here, shows that at `(s, i, w) = (1, 2, 2)` the (H3) system
     `E E^T + w P P^T = N I_4` has NO solution at all with a
     non-Hadamard corner.  So at these parameters `w > 2s` is
     SUFFICIENT but not NECESSARY for D3's conclusion; the boundary
     costs the PROOF, not (here) the statement.

     The check is exhaustive up to the exact symmetries of (H3): (H3)
     sees E only through `E E^T` and P only through `P P^T`, so column
     sign flips of E and of P are free, and the 4s rows may be permuted.
     A non-Hadamard 4x4 sign corner has some off-diagonal
     `(E E^T)[a][b] != 0`, which for +-1 rows of length 4 forces
     `|(E E^T)[a][b]| = 4`, i.e. `row_b = +- row_a`.  Relabel that pair
     to rows 0 and 1 and normalise row 0 of E and of P to all-plus; the
     cert then enumerates `e_1 in {+e_0, -e_0}`, every `p_1`, and every
     completion.

     A control at `(s, i, w) = (1, 2, 3)` -- the same cell with
     `w > 2s`, where D3 proves it -- must also return "none".

Stdlib only.  Exact integers only.  No network.  Nothing is written
inside the repository.
"""

import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

DATA = os.path.join("data", "h20-boundary.json")
VERIFY = os.path.join("verify", "verify.py")

PINNED = {
    "T1-diagK-w2s":
        "50eecc761e12b76944b301b7aaeb03a61cb6b88cfc52c67caaacf20eef0e6c9b",
    "T2-flipped":
        "f279337f4c2376dfaff2f5ec82cec21a055a1888008c1972763c0f6721066caf",
}

# The displayed standard Goethals-Seidel array
#     [  A     BR     CR     DR   ]
#     [ -BR     A    D^T R  -C^T R]
#     [ -CR   -D^T R   A     B^T R]
#     [ -DR    C^T R -B^T R   A   ]
# as (sign, seed index, mode) per block; mode 'd' = dev, 'r' = X R,
# 't' = X^T R.  The transpose-negated orientation negates every 't' block.
STD = {
    (0, 0): (1, 0, 'd'), (1, 1): (1, 0, 'd'),
    (2, 2): (1, 0, 'd'), (3, 3): (1, 0, 'd'),
    (0, 1): (1, 1, 'r'), (0, 2): (1, 2, 'r'), (0, 3): (1, 3, 'r'),
    (1, 0): (-1, 1, 'r'), (2, 0): (-1, 2, 'r'), (3, 0): (-1, 3, 'r'),
    (1, 2): (1, 3, 't'), (1, 3): (-1, 2, 't'),
    (2, 1): (-1, 3, 't'), (2, 3): (1, 1, 't'),
    (3, 1): (1, 2, 't'), (3, 2): (-1, 1, 't'),
}


def die(msg):
    print("FAIL: %s" % msg)
    sys.exit(1)


class Grp(object):
    """Z_{f0} x ... x Z_{fk}, elements as tuples in mixed-radix row-major
    order (the same order the repository's record format uses)."""

    def __init__(self, factors):
        self.f = tuple(factors)
        self.elts = [tuple(t) for t in
                     itertools.product(*[range(m) for m in self.f])]
        self.pos = {g: k for k, g in enumerate(self.elts)}
        self.n = len(self.elts)

    def add(self, a, b):
        return tuple((x + y) % m for x, y, m in zip(a, b, self.f))

    def sub(self, a, b):
        return tuple((x - y) % m for x, y, m in zip(a, b, self.f))

    def zero(self):
        return tuple(0 for _ in self.f)


def subgroup(G, gens):
    K = {G.zero()}
    front = [G.zero()]
    while front:
        a = front.pop()
        for g in gens:
            b = G.add(a, tuple(g))
            if b not in K:
                K.add(b)
                front.append(b)
    return K


def cosets(G, K):
    seen, reps = {}, []
    for g in G.elts:
        if g in seen:
            continue
        c = len(reps)
        reps.append(g)
        for k in K:
            seen[G.add(g, k)] = c
    m = len(reps)
    qadd = [[seen[G.add(reps[a], reps[b])] for b in range(m)]
            for a in range(m)]
    qsub = [[seen[G.sub(reps[a], reps[b])] for b in range(m)]
            for a in range(m)]
    return seen, reps, qadd, qsub


def signs(s):
    if set(s) - {"+", "-"}:
        raise ValueError("string outside {+,-}: %r" % (s,))
    return [1 if ch == "+" else -1 for ch in s]


def mmT(A, B):
    return [[sum(x * y for x, y in zip(ra, rb)) for rb in B] for ra in A]


def core(G, xs, rho, flip):
    n = G.n
    out = [[0] * (4 * n) for _ in range(4 * n)]
    for I in range(4):
        for J in range(4):
            sg, sd, mode = STD[(I, J)]
            if flip and mode == 't':
                sg = -sg
            x = xs[sd]
            for gi, g in enumerate(G.elts):
                row = out[I * n + gi]
                for hi, h in enumerate(G.elts):
                    if mode == 'd':
                        v = x[G.sub(h, g)]
                    elif mode == 'r':
                        v = x[G.sub(rho, G.add(g, h))]
                    else:
                        v = x[G.sub(G.add(g, h), rho)]
                    row[J * n + hi] = sg * v
    return out


def gs_quotient(m, qadd, qsub, sig, rho_bar, flip):
    out = [[0] * (4 * m) for _ in range(4 * m)]
    for I in range(4):
        for J in range(4):
            sg, sd, mode = STD[(I, J)]
            if flip and mode == 't':
                sg = -sg
            x = sig[sd]
            for g in range(m):
                for h in range(m):
                    if mode == 'd':
                        v = x[qsub[h][g]]
                    elif mode == 'r':
                        v = x[qsub[rho_bar][qadd[g][h]]]
                    else:
                        v = x[qsub[qadd[g][h]][rho_bar]]
                    out[I * m + g][J * m + h] = sg * v
    return out


def build(inst):
    """Rebuild one instance: returns (details, rows-as-+-strings)."""
    G = Grp(inst["group"])
    n = G.n
    s = int(inst["s"])
    N = int(inst["order"])
    flip = (inst["orientation"] == "transpose-negated")
    K = subgroup(G, inst["K_generators"])
    kappa, reps, qadd, qsub = cosets(G, K)
    i, w = len(reps), n // len(reps)
    rho = tuple(inst["r_shift"])
    seeds = [signs(x) for x in inst["seeds"]]
    xs = [dict(zip(G.elts, x)) for x in seeds]
    E = [signs(r) for r in inst["corner"]]
    P = [signs(r) for r in inst["row_table"]]
    Q = [signs(r) for r in inst["col_rows"]]

    det = {"n": n, "s": s, "i": i, "w": w, "N": N, "flip": flip,
           "K": sorted(K), "rho": rho}
    fail = []
    if N != 4 * (n + s):
        fail.append("order != 4(|G|+s)")
    if (i, w) != (int(inst["i"]), int(inst["w"])):
        fail.append("declared (i,w)=(%s,%s) but computed (%d,%d)"
                    % (inst["i"], inst["w"], i, w))
    if len(E) != 4 * s or any(len(r) != 4 * s for r in E):
        fail.append("corner shape")
    if len(P) != 4 * s or any(len(r) != 4 * i for r in P):
        fail.append("row_table shape")
    if len(Q) != 4 * i or any(len(r) != 4 * s for r in Q):
        fail.append("col_rows shape")

    # ---- Theorem A (H1): Q Q^T = I_4 (x) M, M Gbar-invariant symmetric ----
    QQ = mmT(Q, Q)
    h1 = True
    for I in range(4):
        for J in range(4):
            if I != J and any(QQ[I * i + a][J * i + b]
                              for a in range(i) for b in range(i)):
                h1 = False
    Mfun = {}
    for I in range(4):
        for a in range(i):
            for b in range(i):
                e = qsub[a][b]
                v = QQ[I * i + a][I * i + b]
                if e in Mfun and Mfun[e] != v:
                    h1 = False
                Mfun.setdefault(e, v)
    for e in range(i):
        if Mfun.get(e) != Mfun.get(qsub[0][e]):
            h1 = False
    det["M"] = [Mfun[e] for e in range(i)]
    det["M_is_house"] = (Mfun == {e: (4 * s + 4 if e == 0 else 0) - 4
                                  for e in range(i)})
    det["H1"] = h1
    if not h1:
        fail.append("H1")

    # ---- (H2): the profile, against the measured M ----
    h2 = True
    for t in G.elts:
        tot = sum(sum(x[g] * x[G.add(g, t)] for g in G.elts) for x in xs)
        want = 4 * n if t == G.zero() else -Mfun[kappa[t]]
        if tot != want:
            h2 = False
    det["H2"] = h2
    if not h2:
        fail.append("H2")

    # ---- (H3): E E^T + w P P^T = N I ----
    EE, PP = mmT(E, E), mmT(P, P)
    h3 = all(EE[a][b] + w * PP[a][b] == (N if a == b else 0)
             for a in range(4 * s) for b in range(4 * s))
    det["H3"] = h3
    det["EEt"] = EE
    det["PPt"] = PP
    det["E_is_Hadamard"] = all(EE[a][b] == (4 * s if a == b else 0)
                               for a in range(4 * s) for b in range(4 * s))
    det["PPt_eq_4i_I"] = all(PP[a][b] == (4 * i if a == b else 0)
                             for a in range(4 * s) for b in range(4 * s))
    if not h3:
        fail.append("H3")

    # ---- (H4): E Q^T + P Chat^T = 0 ----
    sig = [[sum(x[g] for g in G.elts if kappa[g] == c) for c in range(i)]
           for x in xs]
    Chat = gs_quotient(i, qadd, qsub, sig, kappa[rho], flip)
    h4 = all(sum(E[r][j] * Q[k][j] for j in range(4 * s)) +
             sum(P[r][c] * Chat[k][c] for c in range(4 * i)) == 0
             for r in range(4 * s) for k in range(4 * i))
    det["H4"] = h4
    det["sigma"] = sig
    if not h4:
        fail.append("H4")

    # ---- assemble ----
    C = core(G, xs, rho, flip)
    S = 4 * s
    ch = {1: "+", -1: "-"}
    rows = []
    for r in range(S):
        line = [ch[v] for v in E[r]]
        for J in range(4):
            line += [ch[P[r][i * J + kappa[g]]] for g in G.elts]
        rows.append("".join(line))
    for I in range(4):
        for gi, g in enumerate(G.elts):
            line = [ch[Q[i * I + kappa[g]][c]] for c in range(S)]
            line += [ch[v] for v in C[I * n + gi]]
            rows.append("".join(line))
    if len(rows) != N or any(len(r) != N for r in rows):
        fail.append("assembled shape")
    det["failures"] = fail
    return det, rows


# ------------------------------------------------- the boundary statement

def h3_nonhadamard_corner_exists(s, i, w):
    """Is there ANY solution of E E^T + w P P^T = N I_{4s} whose corner is
    not Hadamard?  Exhaustive up to the symmetries of (H3); see the module
    docstring.  Returns (found, witness_or_None, branches)."""
    S, W = 4 * s, 4 * i
    N = 4 * (w * i + s)

    def vecs(k):
        return [[1 if (b >> j) & 1 else -1 for j in range(k)]
                for b in range(1 << k)]

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    Ev, Pv = vecs(S), vecs(W)
    allp = [(e, p) for e in Ev for p in Pv
            if dot(e, e) + w * dot(p, p) == N]

    def orth(r1, r2):
        return dot(r1[0], r2[0]) + w * dot(r1[1], r2[1]) == 0

    e0, p0 = [1] * S, [1] * W
    r0 = (e0, p0)
    branches = 0
    for e1 in ([1] * S, [-1] * S):
        for p1 in Pv:
            r1 = (e1, p1)
            if not orth(r0, r1):
                continue
            branches += 1
            L = [r for r in allp if orth(r, r0) and orth(r, r1)]
            for a in range(len(L)):
                for b in range(a + 1, len(L)):
                    if orth(L[a], L[b]):
                        E = [e0, e1, L[a][0], L[b][0]]
                        EE = [[dot(x, y) for y in E] for x in E]
                        return True, (E, EE), branches
    return False, None, branches


def main():
    for p in (DATA, VERIFY):
        if not os.path.isfile(p):
            die("missing %s -- run this from the repository root: "
                "python certs/05-h20-boundary/run.py" % p)
    with open(DATA, "r", encoding="ascii") as fh:
        insts = json.load(fh)["instances"]
    if sorted(x["label"] for x in insts) != sorted(PINNED):
        die("banked instances != pinned instances")

    bad = []
    t0 = time.time()
    tmp = tempfile.mkdtemp(prefix="cert05-")
    try:
        for inst in insts:
            label = inst["label"]
            det, rows = build(inst)
            if det["failures"]:
                bad.append("%s: %s" % (label, det["failures"]))
                continue
            if det["w"] > 2 * det["s"]:
                bad.append("%s: w = %d > 2s = %d -- this instance is NOT at "
                           "the boundary" % (label, det["w"], 2 * det["s"]))
            if det["w"] != 2 * det["s"]:
                bad.append("%s: w = %d, expected exactly 2s = %d"
                           % (label, det["w"], 2 * det["s"]))
            if len(det["K"]) != det["w"] or det["K"] != sorted(
                    [(0, 0), (1, 1)]):
                bad.append("%s: K = %s is not the diagonal subgroup"
                           % (label, det["K"]))
            path = os.path.join(tmp, "H20_%s.txt" % label)
            with open(path, "w", encoding="ascii", newline="\n") as fh:
                fh.write("\n".join(rows) + "\n")
            proc = subprocess.run([sys.executable, VERIFY, path],
                                  capture_output=True, text=True)
            line = ((proc.stdout or proc.stderr).strip().splitlines()
                    or ["(no output)"])[-1]
            os.remove(path)
            if proc.returncode != 0 or "canonical_sha256=" not in line:
                bad.append("%s: verify.py rc=%d: %s"
                           % (label, proc.returncode, line))
                continue
            got = line.rsplit("canonical_sha256=", 1)[1].strip()
            if got != PINNED[label]:
                bad.append("%s: DIGEST MISMATCH got=%s pinned=%s"
                           % (label, got, PINNED[label]))
            if got != inst.get("pinned_sha256"):
                bad.append("%s: banked pinned_sha256 disagrees with the "
                           "observed digest" % label)
            print("  %-14s orientation=%-18s K=<(1,1)>  n=%d s=%d i=%d w=%d "
                  " (w > 2s ? %s)"
                  % (label, inst["orientation"], det["n"], det["s"],
                     det["i"], det["w"], det["w"] > 2 * det["s"]))
            print("                 Theorem A: H1=%s H2=%s H3=%s H4=%s   "
                  "M = %s (house: %s)"
                  % (det["H1"], det["H2"], det["H3"], det["H4"],
                     det["M"], det["M_is_house"]))
            print("                 E E^T = %s" % (det["EEt"],))
            print("                 P P^T = %s" % (det["PPt"],))
            print("                 corner Hadamard = %s ; P P^T = 4i I = %s"
                  % (det["E_is_Hadamard"], det["PPt_eq_4i_I"]))
            print("                 %s" % line)
            if not det["E_is_Hadamard"]:
                bad.append("%s: the corner is NOT Hadamard -- this cert's "
                           "measured statement has changed" % label)
            if not det["PPt_eq_4i_I"]:
                bad.append("%s: P P^T != 4i I -- this cert's measured "
                           "statement has changed" % label)
            sys.stdout.flush()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # --- the boundary statement, exhaustively ---
    found_b, wit_b, br_b = h3_nonhadamard_corner_exists(1, 2, 2)
    found_c, wit_c, br_c = h3_nonhadamard_corner_exists(1, 2, 3)
    print("  boundary check (s,i,w) = (1,2,2), N=20, w = 2s:  exhaustive "
          "search for an (H3) solution with a non-Hadamard corner over %d "
          "normalised branches -> %s"
          % (br_b, "FOUND" if found_b else "none exists"))
    print("  control        (s,i,w) = (1,2,3), N=28, w > 2s:  same search "
          "-> %s, and the two-row system is already infeasible (%d "
          "branches survive vs %d at the boundary) -- D3's forcing in "
          "miniature: w*(P P^T)[0][1] = -+4 has no solution when w = 3"
          % ("FOUND" if found_c else "none exists", br_c, br_b))
    if found_b:
        bad.append("(1,2,2): a non-Hadamard corner DOES solve (H3); this "
                   "cert's exhaustive statement has changed: %s" % (wit_b,))
    if found_c:
        bad.append("(1,2,3): a non-Hadamard corner solves (H3) even though "
                   "w > 2s -- this contradicts D3: %s" % (wit_c,))

    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 05: FAIL (%d problems)" % len(bad))
        return 1
    print("CERT 05: PASS -- two Hadamard matrices of order 20 satisfying "
          "Theorem A with an arbitrary subgroup K and in either "
          "orientation, at w = 2s exactly, where D3's forcing argument "
          "lapses; both corners are Hadamard anyway, and no non-Hadamard "
          "corner solves (H3) at these parameters (%.1fs)"
          % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
