#!/usr/bin/env python3
"""cert 16 -- Theorem E' (Gram rigidity under w > s) and the w = s boundary.

  THEOREM E' (note/NOTE-B.md S1.7).  Assume (H1), (H2), s >= 1 and
  w > s.  Then S = {chi : Mhat(chi) != 0} has |S| = s, is closed under
  conjugation and under the Galois action, and M = 4 sum_{chi in S} chi
  = 4i P_S, with spectrum {0^(i-s), (4i)^s}.  Theorem E (S1.2.1) needed
  w > 2s; the mod-4 lemma halves the lock.

  LEMMA (mod 4).  For every x : G -> {+-1} and every t != 0,
  PAF_x(t) == n (mod 4).  Hence under (H2), M(c) in 4Z for every c.

  THEOREM E' (boundary).  Let w = s and let M satisfy the necessary
  conditions (M(0) = 4s, M in 4Z, |M(c)| <= 4s, PSD, Mhat <= N/w,
  rank <= s).  Then exactly one of: (a) s < i and M = 4i P_S; (b) s = i
  and M = 4i I - 4C with C(c) = xi(c)[c in Gbar_0 \\ 0] for a subgroup
  Gbar_0 <= Gbar and a real character xi of it; (c) s = i+1 and
  M = 4s I; (d) s >= i+2 cannot occur.

  SHARPNESS.  Case (b) is realised by seeds at (s,i,w) = (3,3,3), so the
  hypothesis w > s cannot be weakened to w >= s.

BOTH THEOREMS ARE PAPER PROOFS (note/NOTE-B.md S1.7).  This certificate
carries the three finite things the proofs lean on or are sharpened by:

  [A] THE MOD-4 LEMMA AS A CONTROL.  The lemma is one line -- along each
      cycle of u -> u+t the number of sign changes of x is even, so the
      disagreement count d(t) is even and PAF = n - 2d(t) == n (mod 4).
      It is instantiated on eleven abelian groups and, separately, the
      cycle statement itself is checked on Z12.  A control, not a proof:
      the proof is in the note.

  [B] THE (3,3,3) ESCAPE, REBUILT FROM NOTHING.  Every +-1 quadruple on
      G = Z3 x Z3 with K = {(0,b)} whose aggregate PAF is -12 on K\\0 and
      +4 off K is found by a meet-in-the-middle over ALL 512 sequences
      per seed -- no coset-sum prefilter, so nothing is assumed about
      the shape of a solution.  A witness is re-checked lag by lag; the
      same enumeration on Z9 with K = <3> returns nothing; and a 12 x 12
      +-1 matrix Q with Q Q^T = I4 (x) (16 I3 - 4 J3) is found by
      depth-first search and re-multiplied.  So (H1)+(H2) hold at
      (s,i,w) = (3,3,3) with M = 16 I - 4 J, spectrum {4,16,16} -- not a
      projector multiple.

  [C] THE BOUNDARY CLASSIFICATION, EXHAUSTED AT SMALL PARAMETERS BY TWO
      ROUTES THAT SHARE NO CODE (the house's two-implementation rule):
        route A (sieve) -- every Gbar-invariant symmetric M : Gbar -> 4Z
          with M(0) = 4s and |M(c)| <= 4s that is PSD, has N I - w M PSD
          (the Parseval window) and rank <= s, decided by exact Fraction
          LDL^T.  No floats anywhere.
        route B (construction) -- the classification's list (a)/(b)/(c)
          built explicitly, the Galois-orbit character sums evaluated as
          Ramanujan sums so every entry is an exact integer.
      At every w = s cell the two sets must COINCIDE.  At every w > s
      cell in the same sweep every survivor must be a projector form
      (Theorem E' proper, on the range s < w <= 2s that Theorem E did
      not reach).  The w < s cells are run to show the sieve is not
      vacuous there -- non-projector survivors exist below the boundary,
      so no rigidity theorem of this type is available in that regime.

WHAT A DEFAULT RUN ESTABLISHES.  Everything above, at the cells listed
in CELLS_DEFAULT: 45 (s,i,w) triples over every abelian Gbar of the
stated order.  `--wide` adds the (7,8,w) cells on Z8.  Nothing is
audited from a bank -- this certificate reads no data file at all.

Usage:
  python certs/16-theorem-eprime-boundary/run.py
  python certs/16-theorem-eprime-boundary/run.py --wide
"""

import argparse
import itertools
import random
import sys
import time
from fractions import Fraction
from math import gcd

sys.dont_write_bytecode = True

_T0 = time.time()
FAIL = []
NCHECK = [0]


def log(msg):
    print("\n[%6.1fs] %s" % (time.time() - _T0, msg), flush=True)


def check(label, cond, extra=""):
    ok = bool(cond)
    NCHECK[0] += 1
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", label,
                           ("  -- " + str(extra)) if extra else ""), flush=True)
    if not ok:
        FAIL.append(label)
    return ok


# ======================================================================
# finite abelian groups
# ======================================================================

class Ab:
    """Finite abelian group Z_f0 x ... x Z_fk, elements as tuples."""

    def __init__(self, factors):
        self.f = tuple(factors)
        self.elts = list(itertools.product(*[range(m) for m in self.f]))
        self.idx = {e: k for k, e in enumerate(self.elts)}
        self.n = len(self.elts)
        self.zero = self.elts[0]
        self.e = 1
        for m in self.f:
            self.e = self.e * m // gcd(self.e, m)

    def add(self, a, b):
        return tuple((x + y) % m for x, y, m in zip(a, b, self.f))

    def neg(self, a):
        return tuple((-x) % m for x, m in zip(a, self.f))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def order(self, a):
        k, x = 1, a
        while x != self.zero:
            x = self.add(x, a)
            k += 1
        return k

    def chi_exp(self, a, c):
        """chi_a(c) = zeta_e ^ chi_exp(a,c)  (exact exponent mod e)."""
        return sum(x * y * (self.e // m)
                   for x, y, m in zip(a, c, self.f)) % self.e

    def subgroups(self):
        """All subgroups, as frozensets (small groups only; <= 3 generators
        suffice for every group this certificate runs)."""
        out = set()
        for k in range(0, 4):
            for gens in itertools.combinations(self.elts, k):
                H = {self.zero}
                frontier = [self.zero]
                while frontier:
                    x = frontier.pop()
                    for g in gens:
                        y = self.add(x, g)
                        if y not in H:
                            H.add(y)
                            frontier.append(y)
                out.add(frozenset(H))
        return sorted(out, key=lambda H: (len(H), sorted(H)))


def abelian_groups(n):
    """Factor tuples of all abelian groups of order n, up to isomorphism."""
    def partitions(k):
        if k == 0:
            yield ()
            return
        for first in range(k, 0, -1):
            for rest in partitions(k - first):
                if not rest or rest[0] <= first:
                    yield (first,) + rest
    fac = {}
    m, p = n, 2
    while m > 1:
        while m % p == 0:
            fac[p] = fac.get(p, 0) + 1
            m //= p
        p += 1
    per_prime = []
    for p, a in sorted(fac.items()):
        per_prime.append([tuple(p ** x for x in part) for part in partitions(a)])
    out = []
    for combo in itertools.product(*per_prime):
        out.append(tuple(sorted(sum(combo, ()), reverse=True)))
    return out or [(1,)]


# ======================================================================
# exact linear algebra -- no floats
# ======================================================================

def psd_and_rank(M):
    """Exact: is the integer symmetric matrix M PSD?  Returns (psd, rank).

    LDL^T over Fractions, symmetric pivoting on the largest remaining
    diagonal entry.  A PSD matrix has every pivot >= 0, and whenever the
    pivot is 0 the whole remaining block must vanish.
    """
    n = len(M)
    A = [[Fraction(v) for v in row] for row in M]
    active = list(range(n))
    rank = 0
    while active:
        p = max(active, key=lambda k: A[k][k])
        piv = A[p][p]
        if piv < 0:
            return False, rank
        if piv == 0:
            for k in active:
                for l in active:
                    if A[k][l] != 0:
                        return False, rank
            return True, rank
        rank += 1
        active.remove(p)
        for k in active:
            if A[k][p] == 0:
                continue
            f = A[k][p] / piv
            for l in active:
                A[k][l] -= f * A[p][l]
    return True, rank


def invariant_matrix(G, Mf):
    return [[Mf[G.sub(c2, c1)] for c2 in G.elts] for c1 in G.elts]


# ======================================================================
# clause A -- the mod-4 lemma, as a control
# ======================================================================

def paf(G, x, t):
    return sum(x[G.idx[u]] * x[G.idx[G.add(u, t)]] for u in G.elts)


def clause_a():
    log("[A] the mod-4 lemma: PAF_x(t) == n (mod 4) for every t != 0")
    rng = random.Random(20260902)
    bad = tested = 0
    for f in [(7,), (8,), (12,), (2, 2, 3), (3, 3), (2, 6), (4, 4),
              (2, 2, 2, 2), (5, 5), (9,), (2, 9)]:
        G = Ab(f)
        for _ in range(30):
            x = [rng.choice((1, -1)) for _ in G.elts]
            for t in G.elts[1:]:
                tested += 1
                if (paf(G, x, t) - G.n) % 4:
                    bad += 1
    check("mod-4 lemma holds on %d random (group, sequence, lag) triples "
          "over eleven abelian groups" % tested, bad == 0,
          "violations = %d" % bad)
    # the one-line proof, instantiated: along each cycle of u -> u+t the
    # number of sign changes is even.
    G = Ab((12,))
    x = [rng.choice((1, -1)) for _ in G.elts]
    ok = True
    for t in G.elts[1:]:
        seen = set()
        for u in G.elts:
            if u in seen:
                continue
            cyc = []
            v = u
            while v not in seen:
                seen.add(v)
                cyc.append(v)
                v = G.add(v, t)
            changes = sum(1 for a, b in zip(cyc, cyc[1:] + cyc[:1])
                          if x[G.idx[a]] != x[G.idx[b]])
            if changes % 2:
                ok = False
    check("every cycle of u -> u + t carries an even number of sign changes "
          "(the proof of the lemma, instantiated on Z12)", ok)


# ======================================================================
# clause B -- the (3,3,3) escape, rebuilt from nothing
# ======================================================================

def clause_b():
    log("[B] the (3,3,3) escape rebuilt: seeds on Z3xZ3 by full MITM, Q by DFS")
    results = {}
    for f, Kgens in (((3, 3), [(0, 1)]), ((9,), [(3,)])):
        G = Ab(f)
        K = {G.zero}
        frontier = [G.zero]
        while frontier:
            v = frontier.pop()
            for g in Kgens:
                y = G.add(v, g)
                if y not in K:
                    K.add(y)
                    frontier.append(y)
        assert len(K) == 3
        lags = G.elts[1:]
        target = tuple(-12 if t in K else 4 for t in lags)
        seqs = list(itertools.product((1, -1), repeat=G.n))       # all 512
        P = [tuple(paf(G, x, t) for t in lags) for x in seqs]
        table = {}
        for a in range(len(seqs)):
            for b in range(a, len(seqs)):
                key = tuple(u + v for u, v in zip(P[a], P[b]))
                table.setdefault(key, []).append((a, b))
        hits = 0
        witness = None
        for a in range(len(seqs)):
            for b in range(a, len(seqs)):
                need = tuple(t - u - v for t, u, v in zip(target, P[a], P[b]))
                lst = table.get(need)
                if lst:
                    hits += len(lst)
                    if witness is None:
                        c, d = lst[0]
                        witness = (seqs[a], seqs[b], seqs[c], seqs[d])
        results[f] = (hits, witness, G, K, lags)
        print("    G = Z%s: %d ordered (pair,pair) solutions of the (3,3,3) "
              "profile" % ("x".join(map(str, f)), hits), flush=True)
    h33, wit, G, K, lags = results[(3, 3)]
    check("Z3xZ3 carries seed quadruples with the (3,3,3) profile "
          "(-12 on K\\0, +4 off K)", h33 > 0 and wit is not None,
          "count = %d" % h33)
    check("Z9 with K = <3> carries none", results[(9,)][0] == 0,
          "count = %d" % results[(9,)][0])
    if wit:
        agg = {t: sum(paf(G, x, t) for x in wit) for t in lags}
        check("witness re-checked lag by lag: aggregate PAF = -12 on K\\0 "
              "and +4 off K",
              all(agg[t] == (-12 if t in K else 4) for t in lags))
        print("      witness seeds (index (a,b) -> 3a+b): %s   row sums %s"
              % (["".join("+" if v == 1 else "-" for v in x) for x in wit],
                 [sum(x) for x in wit]))
        # the witness printed in note/NOTE-B.md S1.7 must satisfy it too
        noted = ["++-++-++-", "++-++-+-+", "++-+-+-++", "++--+++-+"]
        lx = [[1 if ch == "+" else -1 for ch in s] for s in noted]
        check("the witness quoted in note/NOTE-B.md S1.7 satisfies the same "
              "profile",
              all(sum(paf(G, x, t) for x in lx) == (-12 if t in K else 4)
                  for t in lags))
    # (H1): a 12 x 12 +-1 matrix Q with Q Q^T = I4 (x) (16 I3 - 4 J3)
    def tgt(a, b):
        if a == b:
            return 12
        return -4 if a // 3 == b // 3 else 0
    allv = list(itertools.product((1, -1), repeat=12))
    rows = []

    def dfs(k):
        if k == 12:
            return True
        for v in allv:
            if k == 0 and v != allv[0]:
                continue
            if all(sum(x * y for x, y in zip(v, rows[j])) == tgt(k, j)
                   for j in range(k)):
                rows.append(v)
                if dfs(k + 1):
                    return True
                rows.pop()
        return False
    found = dfs(0)
    gram_ok = found and all(
        sum(x * y for x, y in zip(rows[a], rows[b])) == tgt(a, b)
        for a in range(12) for b in range(12))
    check("(H1) at (3,3): a 12x12 +-1 matrix Q with Q Q^T = I4 (x) "
          "(16 I - 4 J) exists (DFS) and re-multiplies exactly", gram_ok)
    M = [[16 if a == b else -4 for b in range(3)] for a in range(3)]
    psd, rk = psd_and_rank(M)
    check("M = 16 I - 4 J on Z3 is PSD of full rank 3 = s = i, spectrum "
          "{4,16,16}: NOT 4i P_S (a projector multiple has all nonzero "
          "eigenvalues equal)", psd and rk == 3)


# ======================================================================
# clause C -- the boundary classification, two routes
# ======================================================================

def ramanujan(d, j):
    """c_d(j) = sum over k in (Z/d)^* of zeta_d^{jk}  (an exact integer)."""
    if d == 1:
        return 1
    g = gcd(j % d, d)
    q = d // g

    def phi(m):
        r, x, p = m, m, 2
        while p * p <= x:
            if x % p == 0:
                while x % p == 0:
                    x //= p
                r -= r // p
            p += 1
        if x > 1:
            r -= r // x
        return r

    def mu(m):
        res, x, p = 1, m, 2
        while p * p <= x:
            if x % p == 0:
                x //= p
                if x % p == 0:
                    return 0
                res = -res
            p += 1
        if x > 1:
            res = -res
        return res
    return mu(q) * phi(d) // phi(q)


def galois_orbits(G):
    """Galois orbits of characters, indexed by a in G (generators of <a>)."""
    seen = set()
    orbits = []
    for a in G.elts:
        if a in seen:
            continue
        d = G.order(a)
        orb = frozenset(tuple((k * x) % m for x, m in zip(a, G.f))
                        for k in range(1, d + 1) if gcd(k, d) == 1)
        seen |= orb
        orbits.append((d, orb))
    return orbits


def M_S(G, orbit_list):
    """M_S(c) = 4 sum_{chi in S} chi(c), S a union of Galois orbits; exact."""
    out = {}
    for c in G.elts:
        tot = 0
        for d, orb in orbit_list:
            a = next(iter(orb))
            x = G.chi_exp(a, c)
            assert x % (G.e // d) == 0
            tot += ramanujan(d, x // (G.e // d))
        out[c] = 4 * tot
    return out


def real_characters(G, H):
    """The real characters of the subgroup H, as dicts c -> +-1."""
    out = {}
    for a in G.elts:
        vals = {}
        ok = True
        for c in H:
            x = G.chi_exp(a, c)
            if x == 0:
                vals[c] = 1
            elif 2 * x == G.e:
                vals[c] = -1
            else:
                ok = False
                break
        if ok:
            out[tuple(sorted(vals.items()))] = vals
    return list(out.values())


def route_b(G, s, i):
    """The classification's list, built explicitly, as a set of M tuples."""
    out = set()
    orbs = galois_orbits(G)
    if s < i:
        for r in range(1, len(orbs) + 1):
            for combo in itertools.combinations(orbs, r):
                if sum(len(o) for d, o in combo) != s:
                    continue
                Mf = M_S(G, list(combo))
                out.add(tuple(Mf[c] for c in G.elts))
    elif s == i:
        for H in G.subgroups():
            for xi in real_characters(G, H):
                Mf = {c: (4 * i if c == G.zero else 0)
                      - 4 * (xi[c] if (c in H and c != G.zero) else 0)
                      for c in G.elts}
                out.add(tuple(Mf[c] for c in G.elts))
    elif s == i + 1:
        out.add(tuple(4 * s if c == G.zero else 0 for c in G.elts))
    return out


def route_a(G, s, i, w):
    """The sieve: every survivor of the necessary conditions, exact."""
    N = 4 * (w * i + s)
    seen = set()
    orbits = []
    for c in G.elts[1:]:
        if c in seen:
            continue
        o = {c, G.neg(c)}
        seen |= o
        orbits.append(sorted(o))
    vals = list(range(-4 * s, 4 * s + 1, 4))
    out = set()
    count = 0
    for choice in itertools.product(vals, repeat=len(orbits)):
        count += 1
        Mf = {G.zero: 4 * s}
        for o, v in zip(orbits, choice):
            for c in o:
                Mf[c] = v
        M = invariant_matrix(G, Mf)
        psd, rk = psd_and_rank(M)
        if not psd or rk > s:
            continue
        W = [[(N if a == b else 0) - w * M[a][b] for b in range(i)]
             for a in range(i)]
        wpsd, _ = psd_and_rank(W)
        if not wpsd:
            continue
        out.add(tuple(Mf[c] for c in G.elts))
    return out, count


def is_projector_form(G, Mt, s):
    """Is M = 4 sum_{chi in S} chi for a Galois-stable S of size s?"""
    i = len(G.elts)
    if s < i:
        return Mt in route_b(G, s, i)
    if s == i:
        return Mt == tuple(4 * s if c == G.zero else 0 for c in G.elts)
    return False        # |S| = s > i is impossible: no projector form exists


CELLS_DEFAULT = [
    # (s, i, w): the w = s boundary (the classification), w > s (Theorem E'
    # proper, including the range s < w <= 2s Theorem E did not reach), and
    # w < s (below the boundary: recorded, nothing asserted).
    (1, 2, 1), (1, 2, 2), (1, 2, 3),
    (2, 2, 2), (2, 2, 1), (2, 2, 3),
    (3, 3, 3), (3, 3, 2), (3, 3, 4),
    (3, 4, 3), (3, 4, 2), (3, 4, 4), (3, 4, 6), (3, 4, 7),
    (4, 4, 4), (4, 4, 3), (4, 4, 5),
    (5, 5, 5), (5, 5, 4), (5, 5, 6),
    (5, 6, 5), (5, 6, 4), (5, 6, 6), (5, 6, 10), (5, 6, 11),
    (3, 6, 3), (3, 6, 4),
    (6, 6, 6), (6, 6, 5), (6, 6, 7),
    (3, 2, 3), (3, 2, 2), (3, 2, 4),
    (4, 3, 4), (4, 3, 3), (4, 3, 5),
    (5, 3, 5), (5, 4, 5), (6, 5, 6),
    (1, 3, 1), (1, 4, 1), (1, 4, 2), (1, 4, 3), (2, 4, 2), (2, 4, 3),
]
CELLS_WIDE = [(7, 8, 7), (7, 8, 6), (7, 8, 8), (7, 8, 15)]


def clause_c(cells, allow_groups=None):
    log("[C] the w = s boundary: route A (exact sieve) vs route B "
        "(construction)")
    summary = []
    for s, i, w in cells:
        for f in abelian_groups(i):
            if allow_groups is not None and f not in allow_groups:
                continue
            G = Ab(f)
            t0 = time.time()
            A, ncand = route_a(G, s, i, w)
            gname = "x".join("Z%d" % m for m in f)
            nonproj = [m for m in A if not is_projector_form(G, m, s)]
            regime = ("w>2s" if w > 2 * s else "s<w<=2s" if w > s else
                      "w=s" if w == s else "w<s")
            line = ("(%d,%d,%d) %-8s %-7s candidates=%-6d survivors=%-3d "
                    "non-projector=%d"
                    % (s, i, w, gname, regime, ncand, len(A), len(nonproj)))
            if w == s:
                B = route_b(G, s, i)
                agree = (A == B)
                check("%s  route A == route B [%d forms] (%.1fs)"
                      % (line, len(B), time.time() - t0), agree,
                      "" if agree else "A-B=%s B-A=%s"
                      % (sorted(A - B), sorted(B - A)))
                summary.append((s, i, w, gname, len(A), len(nonproj),
                                "A==B" if agree else "MISMATCH"))
            elif w > s:
                check("%s  Theorem E' regime: every survivor is 4i P_S "
                      "(%.1fs)" % (line, time.time() - t0),
                      len(nonproj) == 0,
                      "" if not nonproj else str(nonproj[:3]))
                summary.append((s, i, w, gname, len(A), len(nonproj), "rigid"))
            else:
                print("    [info] %s  (below the boundary; recorded, nothing "
                      "asserted)" % line, flush=True)
                summary.append((s, i, w, gname, len(A), len(nonproj), "below"))
    return summary


def clause_d(summary):
    log("[D] global checks over every cell run")
    below_nonproj = [x for x in summary if x[6] == "below" and x[5] > 0]
    check("the sieve is not vacuous: non-projector survivors appear BELOW "
          "the boundary (e.g. %s)"
          % (str(below_nonproj[0][:4]) if below_nonproj else "none"),
          len(below_nonproj) > 0)
    esc = [x for x in summary if x[2] == x[0] and x[0] == x[1] and x[5] > 0]
    esc_below_i = [x for x in summary
                   if x[2] == x[0] and x[0] < x[1] and x[5] > 0]
    check("on the boundary the escapes occur exactly at s = i (and never "
          "with s < i): s=i cells with escapes = %d, s<i boundary cells "
          "with escapes = %d" % (len(esc), len(esc_below_i)),
          len(esc) > 0 and not esc_below_i)
    nb = len([x for x in summary if x[6] == "A==B"])
    nr = len([x for x in summary if x[6] == "rigid"])
    nlo = len([x for x in summary if x[6] == "below"])
    check("cell tally: %d boundary cells (two routes agree), %d cells with "
          "w > s (rigid), %d below the boundary; %d in all"
          % (nb, nr, nlo, len(summary)), nb + nr + nlo == len(summary))
    return nb, nr, nlo


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--wide", action="store_true",
                    help="add the (7,8,w) cells on Z8 (slow: the sieve has "
                         "15^4 candidates there)")
    args = ap.parse_args(argv)

    print("=" * 78)
    print("cert 16 -- Theorem E' (w > s) and the w = s boundary, exact "
          "integers")
    print("=" * 78)

    clause_a()
    clause_b()
    summary = clause_c(CELLS_DEFAULT)
    if args.wide:
        summary += clause_c(CELLS_WIDE, allow_groups=[(8,)])
    nb, nr, nlo = clause_d(summary)

    log("SUMMARY  (s, i, w, Gbar, survivors, non-projector, status)")
    for row in summary:
        print("   ", row)

    print("\n" + "=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 16: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: the mod-4 lemma holds on every triple tested; the")
    print("         (3,3,3) escape is REBUILT FROM NOTHING -- seeds by full")
    print("         meet-in-the-middle, Q by depth-first search -- so")
    print("         Theorem E''s hypothesis w > s cannot be weakened to")
    print("         w >= s; and the boundary classification is EXHAUSTED at")
    print("         %d cells with w = s by two exact routes that share no"
          % nb)
    print("         code, with rigidity confirmed at %d cells with w > s" % nr)
    print("         and the sieve shown non-vacuous at %d cells below the"
          % nlo)
    print("         boundary.  LABEL: PROVEN (the theorems, paper-grade,")
    print("         note/NOTE-B.md S1.7) + PROVEN-BY-CERTIFICATE (the escape")
    print("         and the exhausted classification at the cells run).")
    print("         NOT claimed: anything about (H3)/(H4) at (3,3,3); that")
    print("         every listed boundary form is realised by seeds; any")
    print("         rigidity statement below the boundary; the (7,8,.)")
    print("         cells on Z2xZ4 and Z2^3, which are not swept.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
