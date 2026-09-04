#!/usr/bin/env python3
"""cert 28 -- the (2,4) EXISTENCE theorem: a border kit exists for every
Q', every admissible S-part and every kappa(rho), by an AG(3,2)
transversal argument.

  WHAT THIS IS THE CERTIFICATE OF.  note/NOTE-B.md S2.4 carries the
  (2,4) border proposition (a),(b) and the census (c).  Cert 18 [4]
  established (c) by an exhaustive 215 040-class census -- a count.
  S2.5 replaces the existence half by a THEOREM, and this certificate
  carries the two finite inputs that theorem does not re-derive on
  paper, plus a second construction of the object it produces.

  THE THEOREM -- note/NOTE-B.md S2.5.  Let (s,i) = (2,4), Gbar = Z_4,
  S = {chi, chi^3}, M = (8,0,-8,0).  For every Q' in H(8), every
  admissible S-part (tau_q = (sigma_q(0)-sigma_q(2),
  sigma_q(1)-sigma_q(3)), even entries, sum_q |tau_q|^2 = 8 -- 112 of
  them) and every kappa(rho), an anti-periodic border kit exists.
  Equivalently: (H4) holds iff there is a weighing matrix W(8,4) whose
  rows are odd-signed affine planes of Q's AG(3,2) structure, each
  TRANSVERSAL to the perfect matching Pi read off C*.  Every perfect
  matching of the eight labels admits 2^(3 - rk D) transversal
  parallel classes, D the set of four pair-differences.  HENCE THE
  BORDER IS NEVER THE OBSTRUCTION AT (2,4), AT EVERY w.

  WHICH INPUTS ARE THE CERTIFICATE'S AND WHICH THE PROOF'S.  The
  W(8,4) reformulation, the Walsh classification of the 112 admissible
  signings, the orthogonality rule for two signed planes and the
  2^(3 - rk D) count are PROVEN on paper (S2.5).  Two inputs are
  finite computations and live HERE, not there: the 448-TYPE TABLE OF
  THE TRUE Chat -- from which C* C*^T = 8I and the two +-2 per row and
  column are read -- and the step "Pi IS ALWAYS A TRANSLATION",
  verified type by type.  That is why the section's label is
  PROVEN + PROVEN-BY-CERTIFICATE.

  ONLY THE 7 RANK-1 MATCHINGS -- the translations k -> k xor v, v
  nonzero -- actually occur as Pi.  The rank-2 and rank-3 orbits (42
  and 56 matchings, 4 and 2 transversal planes) are surplus
  generality: this certificate settles them because a uniform
  statement over all 105 matchings is cheaper to certify than a case
  distinction, not because the border needs them.

  THE FOUR (2,4) MATRICES ARE CERT 18's -- cited, not re-banked.  This
  certificate re-uses two of cert 18's seed quadruples (the two H(88)
  at n = 20) to build FRESH borders of its own construction and put
  the assembled matrices through verify/verify.py; the digests it
  produces are its own, and are NOT cert 18's pins, because the
  borders differ.

WHAT THIS SCRIPT DOES  (standard library only, exact integers only)

  [1],[4] C* FROM THE ACTUAL Chat for all 112 S-parts x 4 kappa(rho) =
      448 types: C* C*^T = 8I; exactly two +-2 per row and per column;
      the column supports form a perfect matching; that matching is a
      TRANSLATION of the label space in every one of the 448 types --
      7 distinct matchings in all.
  [3] THE 256-FUNCTION WALSH EXHAUSTION: over all 256 sign functions on
      F_2^3, the spectrum is supported on exactly 4 points with values
      +-4 iff the support is one of the 14 affine planes and the
      signing has an odd number of minus signs -- 112 = 14 x 8.
  [6] THE TRANSVERSAL CENSUS, TWO WAYS (D-008).  |AGL(3,2)| = 1344;
      its orbits on the 105 perfect matchings have sizes 7 / 42 / 56,
      matching ranks 1 / 2 / 3.  ROUTE A: a clique search finds a
      W(8,4) per orbit independently.  ROUTE B: the uniform closed-form
      W (a transversal plane, its complement, all four sign origins u)
      on ALL 105 matchings, with #classes = 2^(3 - rk D) asserted
      matching by matching.  The two routes agree everywhere; the
      closed form for translation by e_3 is written out.
  [7] END TO END, THROUGH THE TRUST CHAIN.  For cert 18's two H(88)
      seed quadruples (n = 20, w = 5, rho-bar 0 and 1), a random Q'
      (random affine structure and random signs) is drawn, this
      certificate's own construction produces (E, P'), the
      anti-periodic doubling is checked to satisfy (H1) and (H3) and
      (H4) against the TRUE Chat, and the matrix is assembled by this
      certificate's own block-explicit assembler and handed to
      verify/verify.py.

  NO --full.  The census is exhaustive on the default path: 448 of 448
  types, 105 of 105 matchings, 256 of 256 functions.  There is no
  heavier leg to offer.  Measured here: 0.7 s.

  Cert 18 [4]'s 215 040-class census is an INDEPENDENT SECOND PROOF of
  the same existence statement, by enumeration rather than by
  structure; it is not re-run here.

  PORTED from the source laboratory (Hadamard-2060 certs/0027-ag32-
  transversal, D-068), unchanged in its mathematics; see NOTES.md.

Usage:
  python certs/28-ag32-transversal/run.py
"""

import hashlib
import itertools
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
VERIFY = os.path.join(ROOT, "verify", "verify.py")
SCRATCH = None   # a private temporary directory, made in main()

T0 = time.time()
OK = True
NCHECK = [0]


def pm(s):
    return [1 if ch == "+" else -1 for ch in s]


def to_pm(row):
    return "".join("+" if v == 1 else "-" for v in row)


def canonical_sha256(rows):
    text = "\n".join(to_pm(r) for r in rows) + "\n"
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def gs_blocks(seqs, m, rho):
    def dev(x):
        return [[x[(h - g) % m] for h in range(m)] for g in range(m)]

    def xr(x):
        return [[x[(rho - g - h) % m] for h in range(m)] for g in range(m)]

    def xtr(x):
        return [[x[(g + h - rho) % m] for h in range(m)] for g in range(m)]

    def ng(X):
        return [[-v for v in r] for r in X]
    A = dev(seqs[0])
    BR, CR, DR = xr(seqs[1]), xr(seqs[2]), xr(seqs[3])
    BtR, CtR, DtR = xtr(seqs[1]), xtr(seqs[2]), xtr(seqs[3])
    return [[A, BR, CR, DR],
            [ng(BR), A, DtR, ng(CtR)],
            [ng(CR), ng(DtR), A, BtR],
            [ng(DR), CtR, ng(BtR), A]]


def assemble_bordered(xs, n, rho, E, P, Q):
    """House-orientation BGS: [E P~; Q~ C], strips constant on Z_4-cosets."""
    N = 4 * (n + 2)
    core = gs_blocks(xs, n, rho)
    rows = []
    for r in range(8):
        line = list(E[r])
        for J in range(4):
            line += [P[r][4 * J + (h % 4)] for h in range(n)]
        rows.append(line)
    for I in range(4):
        for g in range(n):
            line = list(Q[4 * I + (g % 4)])
            for J in range(4):
                line += core[I][J][g]
            rows.append(line)
    assert len(rows) == N and all(len(r) == N for r in rows)
    return rows


def run_verify(rows, scratch_dir, name):
    os.makedirs(scratch_dir, exist_ok=True)
    p = os.path.join(scratch_dir, name + ".txt")
    with open(p, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(to_pm(r) for r in rows) + "\n")
    out = subprocess.run([sys.executable, VERIFY, p], capture_output=True, text=True,
                         env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    v = [l for l in (out.stdout + out.stderr).splitlines() if l.startswith("VERDICT")]
    os.unlink(p)
    return out.returncode, (v[-1] if v else "")


def say(label, cond, extra=""):
    global OK
    NCHECK[0] += 1
    OK &= bool(cond)
    print("  [%s] %s%s" % ("PASS" if cond else "FAIL", label, ("  -- " + str(extra)) if extra else ""), flush=True)


def matmul(A, B):
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) for col in Bt] for row in A]


def transpose(A):
    return [list(c) for c in zip(*A)]


# ---------------------------------------------------------------- Chat over Z_4 and C*
GS_TAB = {(0, 1): (1, "R", 1), (0, 2): (2, "R", 1), (0, 3): (3, "R", 1),
          (1, 0): (1, "R", -1), (1, 2): (3, "TR", 1), (1, 3): (2, "TR", -1),
          (2, 0): (2, "R", -1), (2, 1): (3, "TR", -1), (2, 3): (1, "TR", 1),
          (3, 0): (3, "R", -1), (3, 1): (2, "TR", 1), (3, 2): (1, "TR", -1)}


def chat(sigma, rho):
    m = 4
    out = []
    for I in range(4):
        for c in range(m):
            row = []
            for J in range(4):
                for c2 in range(m):
                    if I == J:
                        row.append(sigma[0][(c2 - c) % m])
                    else:
                        q, form, sg = GS_TAB[(I, J)]
                        row.append(sg * (sigma[q][(rho - c - c2) % m] if form == "R" else sigma[q][(c + c2 - rho) % m]))
            out.append(row)
    return out


def cstar_from_chat(C):
    """C*[(I,c),(J,c'')] for c, c'' in {0,1}, and the anti-periodicity check in the row index."""
    Cs = [[0] * 8 for _ in range(8)]
    for I in range(4):
        for c in range(2):
            for J in range(4):
                for c2 in range(2):
                    Cs[2 * I + c][2 * J + c2] = sum((-1) ** (c1 // 2) * C[4 * I + c][4 * J + c1] for c1 in range(4) if c1 % 2 == c2)
    anti = True
    for I in range(4):
        for c in range(2):
            for J in range(4):
                for c2 in range(2):
                    v = sum((-1) ** (c1 // 2) * C[4 * I + c + 2][4 * J + c1] for c1 in range(4) if c1 % 2 == c2)
                    if v != -Cs[2 * I + c][2 * J + c2]:
                        anti = False
    return Cs, anti


def sigma_from_tau(tau, rpart=None, s2=None):
    """sigma_q on Z_4 with sigma_hat_q(chi) = tau_q[0] + i tau_q[1] (chi(c) = i^c), row sum rpart[q],
    chi^2-component s2[q]; integrality needs the parities to match (see below)."""
    out = []
    for q in range(4):
        a, b = tau[q]
        r = 0 if rpart is None else rpart[q]
        t2 = 0 if s2 is None else s2[q]
        # sigma(c) = (1/4)[r + t2 (-1)^c + 2 Re((a + ib) i^{-c})]
        row = []
        for c in range(4):
            re = [a, b, -a, -b][c]           # Re((a+ib) i^{-c}) = a cos - ... : c=0: a; c=1: b; c=2: -a; c=3: -b
            v = r + t2 * (-1) ** c + 2 * re
            assert v % 4 == 0, (tau, rpart, s2)
            row.append(v // 4)
        out.append(row)
    return out


def sylvester8():
    return [[(-1) ** bin(x & y).count("1") for y in range(8)] for x in range(8)]


H0 = sylvester8()
PLANES = []
for a in range(1, 8):
    for b in (0, 1):
        PLANES.append(frozenset(y for y in range(8) if bin(a & y).count("1") % 2 == b))
assert len(PLANES) == 14


def structures():
    """The 30 affine structures on the labels {0..7}: label k carries point pi[k]; planes as label sets."""
    out = {}
    for pi in itertools.permutations(range(8)):
        inv = {pi[k]: k for k in range(8)}
        key = frozenset(frozenset(inv[y] for y in pl) for pl in PLANES)
        out.setdefault(key, pi)
    return out


def is_hadamard(H):
    n = len(H)
    return all(sum(H[i][c] * H[j][c] for c in range(n)) == (n if i == j else 0) for i in range(n) for j in range(n))


def matching_from_cstar(Cs):
    pairs = []
    for b in range(8):
        nz = [k for k in range(8) if Cs[k][b] != 0]
        pairs.append(frozenset(nz))
    return pairs


def agl32():
    """AGL(3,2) as permutations of {0..7} = F_2^3 (order 1344)."""
    mats = []
    for cand in itertools.product((0, 1), repeat=9):
        M = [cand[0:3], cand[3:6], cand[6:9]]
        # invertible over F_2 ?
        rows = [sum(M[r][c] << c for c in range(3)) for r in range(3)]
        span = {0}
        for r in rows:
            span |= {s ^ r for s in span}
        if len(span) == 8:
            mats.append(M)
    perms = []
    for M in mats:
        for t in range(8):
            perm = []
            for x in range(8):
                bits = [(x >> c) & 1 for c in range(3)]
                y = sum(((sum(M[r][c] * bits[c] for c in range(3)) % 2) << r) for r in range(3)) ^ t
                perm.append(y)
            perms.append(tuple(perm))
    return perms


def matching_orbits(group):
    """Orbits of the group on perfect matchings of {0..7}."""
    def matchings(pts):
        if not pts:
            yield frozenset()
            return
        a = pts[0]
        for b in pts[1:]:
            rest = [p for p in pts if p not in (a, b)]
            for m in matchings(rest):
                yield m | {frozenset((a, b))}
    allm = list(matchings(list(range(8))))
    assert len(allm) == 105
    seen, orbits = set(), []
    for m in allm:
        if m in seen:
            continue
        orb = {frozenset(frozenset(g[x] for x in pair) for pair in m) for g in group}
        seen |= orb
        orbits.append(orb)
    return orbits


def find_W(Pi):
    """8 mutually orthogonal signed planes of the STANDARD structure, each a transversal of the matching Pi.
    Signs d = 1 (existence is d-free).  Returns W or None."""
    cands = []
    for pl in PLANES:
        if all(len(pl & pair) == 1 for pair in Pi):
            for u in pl:
                w = [0] * 8
                for k in pl:
                    w[k] = -1 if k == u else 1
                cands.append(tuple(w))
    n = len(cands)
    orth = [[sum(a * b for a, b in zip(cands[i], cands[j])) == 0 for j in range(n)] for i in range(n)]
    sol = []

    def rec(pool):
        if len(sol) == 8:
            return True
        for i in pool:
            sol.append(i)
            if rec([j for j in pool if j > i and orth[i][j]]):
                return True
            sol.pop()
        return False
    return [list(cands[i]) for i in sol] if rec(list(range(n))) else None


def explicit_W_parallel(v):
    """(6): the closed-form W for Pi = translation by v, in coordinates where v = e_3 (only v = 4 needed)."""
    assert v == 4
    M = [[1, 1], [0, 1]]
    W = []
    for alpha in range(4):
        a0, a1 = alpha & 1, (alpha >> 1) & 1
        m0 = (M[0][0] * a0 + M[0][1] * a1) % 2
        m1 = (M[1][0] * a0 + M[1][1] * a1) % 2
        xu = m0 | (m1 << 1)
        for c in (0, 1):
            w = [0] * 8
            for x in range(4):
                x0, x1 = x & 1, (x >> 1) & 1
                z = (a0 * x0 + a1 * x1 + c) % 2
                pt = x | (z << 2)
                w[pt] = -1 if x == xu else 1
            W.append(w)
    return W


def main():
    global SCRATCH
    SCRATCH = tempfile.mkdtemp(prefix="cert28-")
    print("(1),(4) C* from the actual Chat: all 112 S-parts (two loud slots, both signs) x 4 kappa(rho)")
    slots = [(q, ab) for q in range(4) for ab in (0, 1)]
    types = []
    for (s1, s2) in itertools.combinations(slots, 2):
        for e1 in (2, -2):
            for e2 in (2, -2):
                tau = [[0, 0] for _ in range(4)]
                tau[s1[0]][s1[1]] += e1
                tau[s2[0]][s2[1]] += e2
                for rho in range(4):
                    types.append((tuple(map(tuple, tau)), rho))
    say("448 types = 112 S-parts x 4 kappa(rho) (C6's count, both signs of tau)", len(types) == 448, len(types))
    ok_struct = True
    matchings_seen = {}
    for tau, rho in types:
        # sigma = (r + t2 (-1)^c + 2 Re(tau_q i^-c))/4 is integral with r = t2 = 0 since the tau entries are even;
        # (H4) sees only the S-part (C1(f)), so this choice loses nothing.
        sigma = sigma_from_tau(tau)
        C = chat(sigma, rho)
        Cs, anti = cstar_from_chat(C)
        gram = matmul(Cs, transpose(Cs))
        two = all(sorted(abs(v) for v in row) == [0, 0, 0, 0, 0, 0, 2, 2] for row in Cs) and \
            all(sorted(abs(Cs[k][b]) for k in range(8)) == [0, 0, 0, 0, 0, 0, 2, 2] for b in range(8))
        if not (anti and two and all(gram[a][b] == (8 if a == b else 0) for a in range(8) for b in range(8))):
            ok_struct = False
        pairs = matching_from_cstar(Cs)
        Pi = frozenset(pairs)
        if len(Pi) != 4 or any(len(p) != 2 for p in Pi) or len(set().union(*Pi)) != 8:
            ok_struct = False
        # translation in label space?
        vs = {min(p) ^ max(p) for p in Pi}
        matchings_seen.setdefault(Pi, []).append((tau, rho))
        if len(vs) != 1:
            ok_struct = False
    say("every type: C* anti-periodic in the row index, C* C*^T = 8I, exactly two +-2 per row and per column; the column "
        "supports form a perfect matching (each pair twice) that is a translation k -> k xor v of the label space",
        ok_struct)
    say("the 448 types realise %d distinct matchings (the lane: 7 pairings)" % len(matchings_seen), len(matchings_seen) == 7)

    print("\n(3) admissible rows: the 256 +-1 functions on F_2^3 with a 4-point Walsh support of values +-4")
    good = 0
    supports = set()
    for g in itertools.product((1, -1), repeat=8):
        W = [sum(g[x] * H0[y][x] for x in range(8)) for y in range(8)]
        nz = [y for y in range(8) if W[y] != 0]
        if len(nz) == 4 and all(abs(W[y]) == 4 for y in nz):
            good += 1
            supports.add(frozenset(nz))
            signs = [W[y] // 4 for y in nz]
            if signs.count(-1) % 2 == 0:
                good = -10 ** 6
    say("exactly 112 = 14 planes x 8 odd signings; every support is an affine plane", good == 112 and supports == set(PLANES),
        (good, len(supports)))

    print("\n(6) classification: AGL(3,2)-orbits on the 105 perfect matchings, and a W(8,4) of transversal signed planes per orbit")
    G = agl32()
    say("|AGL(3,2)| = 1344", len(G) == 1344, len(G))
    orbits = matching_orbits(G)
    print("   orbit sizes: %s" % sorted(len(o) for o in orbits))
    all_have_W = True
    reps = {}
    for orb in orbits:
        Pi = next(iter(orb))
        W = find_W(Pi)
        ntrans = sum(1 for pl in PLANES if all(len(pl & pair) == 1 for pair in Pi))
        reps[Pi] = W
        print("   orbit of size %3d, representative %s: transversal planes %d, W found: %s"
              % (len(orb), sorted(sorted(p) for p in Pi), ntrans, W is not None))
        if W is None:
            all_have_W = False
        else:
            WWt = matmul(W, transpose(W))
            assert all(WWt[a][b] == (4 if a == b else 0) for a in range(8) for b in range(8))
            used = {}
            for w in W:
                used.setdefault(frozenset(k for k in range(8) if w[k]), []).append([k for k in range(8) if w[k] == -1][0])
            print("      W rows: %s" % "  ".join("".join({0: ".", 1: "+", -1: "-"}[v] for v in w) for w in W))
            print("      planes used (support -> list of u): %s" % "; ".join("%s -> %s" % (sorted(p), us) for p, us in used.items()))
    say("every orbit admits a W(8,4) of admissible transversal signed planes", all_have_W)
    # (6) the transversal argument on all 105 matchings: rank of D, 2^{3-rank} transversal classes, uniform W
    def f2rank(vecs):
        basis = []
        for v in vecs:
            for b in basis:
                v = min(v, v ^ b)
            if v:
                basis.append(v)
        return len(basis)
    allm = set().union(*orbits)
    rank_ok, uni_ok = True, True
    rank_hist = {}
    for Pi in allm:
        D = [min(p) ^ max(p) for p in Pi]
        rk = f2rank(D)
        trans = [pl for pl in PLANES if all(len(pl & pair) == 1 for pair in Pi)]
        # classes: functionals l (1..7) with l(d) = 1 for all d in D
        classes = [l for l in range(1, 8) if all(bin(l & d).count("1") % 2 == 1 for d in D)]
        if len(classes) != 2 ** (3 - rk) or len(trans) != 2 * len(classes):
            rank_ok = False
        rank_hist[rk] = rank_hist.get(rk, 0) + 1
        p = trans[0]
        pc = frozenset(range(8)) - p
        W = []
        for pl in (p, pc):
            for u in sorted(pl):
                W.append([(-1 if k == u else 1) if k in pl else 0 for k in range(8)])
        WWt = matmul(W, transpose(W))
        if not all(WWt[a][b] == (4 if a == b else 0) for a in range(8) for b in range(8)):
            uni_ok = False
    say("all 105 matchings: #transversal classes = 2^(3 - rank D) (rank histogram %s = orbit sizes 7/42/56), and the "
        "uniform W (a transversal plane and its complement, all four u each) has W W^T = 4I" % rank_hist,
        rank_ok and uni_ok and rank_hist == {1: 7, 2: 42, 3: 56})
    # the explicit construction for the parallel-class orbit
    Wp = explicit_W_parallel(4)
    Pi4 = frozenset(frozenset((k, k ^ 4)) for k in range(8))
    WWt = matmul(Wp, transpose(Wp))
    say("(6) closed form for Pi = translation by e_3: W W^T = 4I, every row an admissible transversal plane",
        all(WWt[a][b] == (4 if a == b else 0) for a in range(8) for b in range(8)) and
        all(frozenset(k for k in range(8) if w[k]) in set(PLANES) and all(len(frozenset(k for k in range(8) if w[k]) & p) == 1 for p in Pi4)
            and sum(1 for k in range(8) if w[k] == -1) % 2 == 1 for w in Wp))
    # every (structure, label-translation) pair in the census lands in a classified orbit -- tautological (the
    # orbits cover all matchings), recorded as the count of orbits met by the 30 x 7 relative positions:
    structs = structures()
    say("30 affine structures on the labels", len(structs) == 30, len(structs))
    met = set()
    for key, pi in structs.items():
        for v in range(1, 8):
            Pi = frozenset(frozenset((k, k ^ v)) for k in range(8))
            # transport Pi to the standard structure through pi (label k -> point pi[k])
            Pis = frozenset(frozenset(pi[k] for k in pair) for pair in Pi)
            for idx, orb in enumerate(orbits):
                if Pis in orb:
                    met.add(idx)
    say("the 30 x 7 (structure, translation) relative positions of the census meet %d of the %d orbits" % (len(met), len(orbits)),
        len(met) == len(orbits))

    print("\n(7) end to end: borders built by this construction for cert 18's H(88) seed quadruples, verify.py")
    # seeds as banked in certs/18-cell24-instances/run.py (from cell24/common.py DESK_SEEDS_20)
    instances = [
        dict(name="H88_seed0_rho0", n=20, rho=0,
             seeds=["+--++-+-++++-+++++--", "-++----+-+-----+++-+", "++-++-+--++----+-+-+", "-----++++---+--+-+--"]),
        dict(name="H88_seed1_rho1", n=20, rho=1,
             seeds=["++-+-+-+-+++--+-++++", "----++-++-+--+----++", "+-++-----++-+++++---", "+-++---++-+-+++++++-"]),
    ]
    rng = random.Random(24)
    for inst in instances:
        n, rho = inst["n"], inst["rho"]
        x = [pm(s) for s in inst["seeds"]]
        # coset sums over Z_4 (kappa(g) = g mod 4), Chat, C*, matching
        sigma = [[sum(x[q][g] for g in range(n) if g % 4 == c) for c in range(4)] for q in range(4)]
        C = chat(sigma, rho % 4)
        Cs, anti = cstar_from_chat(C)
        Pi = frozenset(matching_from_cstar(Cs))
        # a random Q' = D P_pi H0 with a random affine structure and random signs
        pi = list(range(8))
        rng.shuffle(pi)
        d = [rng.choice((1, -1)) for _ in range(8)]
        Qp = [[d[k] * H0[pi[k]][y] for y in range(8)] for k in range(8)]
        # W: transport Pi to the standard structure, find W there, transport back, apply d
        Pis = frozenset(frozenset(pi[k] for k in pair) for pair in Pi)
        Ws = find_W(Pis)
        assert Ws is not None
        inv = {pi[k]: k for k in range(8)}
        W = [[0] * 8 for _ in range(8)]
        for r in range(8):
            for y in range(8):
                if Ws[r][y]:
                    W[r][inv[y]] = Ws[r][y] * d[inv[y]]
        E = [[-v // 2 for v in row] for row in matmul(W, Qp)]
        Pp = [[v // 2 for v in row] for row in matmul(W, Cs)]
        ok = all(abs(v) == 1 for row in E for v in row) and all(abs(v) == 1 for row in Pp for v in row)
        ok &= is_hadamard(E) and is_hadamard(Pp)
        # anti-periodic doubling and the full (H1),(H3),(H4) at w = 5 (N = 88)
        Q = [[(-1) ** (c // 2) * Qp[2 * I + (c % 2)][r] for r in range(8)] for I in range(4) for c in range(4)]
        P = [[(-1) ** (c // 2) * Pp[r][2 * J + (c % 2)] for J in range(4) for c in range(4)] for r in range(8)]
        QQt = matmul(Q, transpose(Q))
        Mv = {0: 8, 1: 0, 2: -8, 3: 0}
        h1 = all(QQt[4 * I + c][4 * J + c2] == (Mv[(c - c2) % 4] if I == J else 0) for I in range(4) for J in range(4) for c in range(4) for c2 in range(4))
        EEt, PPt = matmul(E, transpose(E)), matmul(P, transpose(P))
        h3 = all(EEt[a][b] + 5 * PPt[a][b] == (88 if a == b else 0) for a in range(8) for b in range(8))
        h4 = all(sum(E[r][j] * Q[k][j] for j in range(8)) + sum(P[r][cc] * C[k][cc] for cc in range(16)) == 0 for r in range(8) for k in range(16))
        say("%s: E, P' in H(8) from W; anti-periodic P, Q satisfy (H1), (H3) at w = 5, (H4) against the true Chat" % inst["name"],
            ok and h1 and h3 and h4)
        rows = assemble_bordered(x, n, rho, E, P, Q)
        rc, v = run_verify(rows, SCRATCH, inst["name"])
        say("%s assembled with THIS border: verify.py %s, digest %s..." % (inst["name"], v[:17] if v else "?", canonical_sha256(rows)[:12]), rc == 0)

    print("\nCONCLUSION: the existence step is a theorem: (H4) reduces to a weighing matrix W(8,4) whose rows are signed affine")
    print("planes of Q''s AG(3,2) structure, transversal to the matching read off C*; the relative position (structure, matching)")
    print("falls into finitely many AGL(3,2)-orbits, each of which admits a W (exhibited), the parallel-class orbit in closed")
    print("form.  With C6's reduction this proves 'the border is never the obstruction at (2,4)' for every Q', S-part and")
    print("kappa(rho); the census (215 040 / 215 040) is its independent second proof.")
    shutil.rmtree(SCRATCH, ignore_errors=True)
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if not OK:
        print("cert 28: SOME CHECK FAILED  (%.1fs)" % (time.time() - T0))
        return 1
    print("ALL CHECKS PASS  (%.1fs)" % (time.time() - T0))
    print("VERDICT: AT THE CELL (s,i) = (2,4) THE BORDER IS NEVER THE")
    print("         OBSTRUCTION, AT EVERY w: for every Q' in H(8), every")
    print("         one of the 112 admissible S-parts and every kappa(rho)")
    print("         an anti-periodic border kit exists.")
    print("         LABEL: PROVEN (the W(8,4) reformulation, the Walsh")
    print("         classification, the orthogonality rule and the")
    print("         2^(3 - rk D) count -- note/NOTE-B.md S2.5)")
    print("         + PROVEN-BY-CERTIFICATE (the 448-type table of the true")
    print("         Chat and the step 'Pi is always a translation', which")
    print("         the proof takes from here; and the transversal census")
    print("         by two constructions of W).")
    print("         NOT CLAIMED: anything about H(2092) -- the (2,4) cell")
    print("         does not land at 2092; that anti-periodicity of P is")
    print("         necessary at w <= 4 (it is not shown, and the two")
    print("         H(56) instances have w = 3); any statement about the")
    print("         SEED layer at this cell; and any equivalence among the")
    print("         four (2,4) matrices of cert 18, or to other known")
    print("         H(56) / H(88).")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
