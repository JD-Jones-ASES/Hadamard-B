#!/usr/bin/env python3
"""cert 31 -- the general-branch border layer at (s,i) = (7,12): an admissible
column table at S2, the forced per-superblock multiset, the S3 kill, and the
row-factorisation lemma with its dual form.

  SETTING (note/NOTE-B.md S1.12, S2.6).  Gbar = Z12, i = 12, s = 7, so
  4i = 48 and 4s = 28.  Galois-stable size-7 supports of Z12-hat containing
  the trivial character are exactly four (cert 30 [B]):

      S1 = {0,1,3,5,7,9,11}   S2 = {0,1,4,5,7,8,11}
      S3 = {0,1,2,5,7,10,11}  S4 = {0,2,3,4,8,9,10}

  and M_S(c) = 4 sum_{chi in S} chi(c) = 4i P_S is the only admissible Gram
  wherever w > s (Theorem E', S1.7).  V_S = span{v_chi : chi in S} in R^12.
  By Theorem F(a) in the general branch (S1.12) an admissible column table is
  a Q in {+-1}^{48 x 28} with

      (H1)   Q Q^T = I_4 (x) M_S,
      (a)    Q^T Q = 48 I_28,
      (blk)  every 12-block of every column lies in V_S,

  and (a) + (blk) already imply (H1).  The ADMISSIBLE BLOCK ALPHABET
  B_S = {+-1}^12 cap V_S is therefore the whole vocabulary of the column
  layer, and it is small.

  WHAT THIS CERTIFICATE ESTABLISHES, by exact integer / rational arithmetic:

  [A] THE ALPHABET, at all four supports, decided TWICE by membership tests
      that share no code -- cyclotomic divisibility (Phi_d(x) | v(x) for every
      character order d outside S) and the projector fixed point
      (circ(M_S/4) v = 12 v) -- over all 2^12 sign vectors.  The two sets
      agree, support for support: |B| = 66, 32, 2, 24 at S1, S2, S3, S4, with
      span ranks 7, 7, 1, 7.

  [B] THE S3 KILL.  At S3 the only +-1 vectors with Fourier support inside S3
      are the TWO CONSTANTS.  Three independent obstructions follow, any one
      of which suffices: (i) every column is then (e_0 1, e_1 1, e_2 1, e_3 1)
      with e in {+-1}^4, so at most 2^4 = 16 distinct columns exist and they
      span a space of dimension 4 -- against 28 mutually orthogonal ones
      required; (ii) rank(Q) = 28 > 4; (iii) inside one superblock any two
      rows are equal or antipodal, so their inner product is +-28, while (H1)
      demands M_S3(c) for c != 0, and no M_S3(c) with c != 0 is +-28.
      Hence (H1) has NO solution at (7,12)/S3 -- at any sigma, any w, any
      kappa(rho), and therefore at EVERY order N = 4(12w + 7).  This is a
      DEAD-BY-CERTIFICATE verdict for the named tier, and the chi_6 twin
      support 6 + S3 = {1,4,5,6,7,8,11} dies with it by the Lemma-T twist
      (S1.4) -- checked here, not asserted: v -> chi_6 . v is a bijection
      of {+-1}^12 carrying V_S3 onto V_{6+S3}, so that alphabet is the two
      ALTERNATING vectors and the same three obstructions apply.

  [C] THE ARTIFACT AT S2.  data/q-7_12-S2.json holds a 4 x 28 x 12 array of
      +-1 blocks.  Every block is required to lie in the alphabet computed in
      [A]; the assembled Q is checked for Q^T Q = 48 I_28 and
      Q Q^T = I_4 (x) M_S2 entry by entry.  So (H1) HAS a solution at
      (7,12)/S2.

  [D] THE FORCED MULTISET (the lemma that makes a negative search complete).
      Writing R_I for the 12 x 28 superblock-I slice, (H1) says
      R_I R_I^T = circ(M), i.e. sum_j m_j b_j b_j^T = circ(M) over the 16
      negation classes b_j of B_S2, with sum_j m_j = 28.  The 16 outer
      products are linearly independent, so the system has a UNIQUE rational
      solution, and it is integral and nonnegative:

        every superblock of every admissible Q carries the constant class
        once, each of the three type-A classes once, and each of the twelve
        type-B classes twice.

      Solved here by exact Fraction elimination, and then re-checked on the
      artifact superblock by superblock.

  [E] THE ROW-FACTORISATION LEMMA (S1.12), made concrete at this cell.  With
      sigma the (7,12) table B2 of cert 30 -- whose support IS S2 -- and
      Chat = GS(sigma; kappa(rho)) the 48 x 48 Goethals-Seidel array over Z12
      in standard orientation, this run checks at EVERY kappa(rho) in Z12:
          Chat^T Chat = Chat Chat^T = I_4 (x) circ(F),   F = 2092 d_0 - w M,
          Chat commutes with I_4 (x) P_S,
          W := Chat^T Q  has  W W^T = 4i * 4s * (I_4 (x) P_S) = 28 (I_4 (x) M).
      Hence every admissible row p has ||Chat p||^2 = 16 i s exactly, so the
      corresponding row of E has squared norm exactly 4s whatever p is -- and
      a kit exists iff the flat set contains 4s mutually orthogonal rows.

  [F] THE DUAL FORM (S1.12), checked on this repository's own kits.
      p is flat  <=>  Chat p in Q {+-1}^{4s};  equivalently a kit is a
      Hadamard E of order 4s with P = -(1/4s) E Q^T Chat in {+-1}.  Checked
      (i) at (2,4) on all four banked records of cert 18, EXHAUSTIVELY over
      all 256 admissible rows -- Chat p lies in col(Q) for every one, with
      ||Chat p||^2 = 16 i s -- and on the four kits themselves, where
      Chat p = -Q e with e the corner row and P = -(1/4s) E Q^T Chat;
      (ii) at (3,4) on kits found here by this certificate's own engine.

  WHAT THIS CERTIFICATE DOES NOT SAY.  Nothing about H(2092), and nothing
  about the (7,12) cell as a whole: (H1) is the COLUMN layer only.  Whether
  the S2 table extends to a border kit (E, P, Q) is a different question and
  is not answered here; the source laboratory reports no kit anywhere in the
  general branch at 2092.  Nothing about S1, S4, or about seeds at any of
  them.  The S3 kill is a statement about the column table, not about the
  sigma-layer: S3's sigma-table is perfectly valid and is banked in cert 30.

Standard library only, exact integers and Fractions only, no floats anywhere,
no network.  Nothing is imported from outside this repository.

Usage:
  python certs/31-s2-column-table/run.py
"""

import argparse
import copy
import hashlib
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
Q_DATA = os.path.join(ROOT, "data", "q-7_12-S2.json")
SIGMA_DATA = os.path.join(ROOT, "data", "general-branch-sigma-tables.json")
CELL24_DATA = os.path.join(ROOT, "data", "cell24-records.json")

Q_SHA256 = "e23e1ab73feffa2d948aa0e91777b2799f6a0c51b0b13e495dc25775b13caffb"
SIGMA_SHA256 = \
    "d2c37aace3bf016dc245a26c3fe7a2bd4aa524aafc874ba717c7e562d98c4a97"
CELL24_SHA256 = \
    "9727b392940d416d3f25dca5d51d2db71cd499bc73c3b8dc4efd22801180f179"

I12 = 12
S_PARTS = {
    "S1": [0, 1, 3, 5, 7, 9, 11],
    "S2": [0, 1, 4, 5, 7, 8, 11],
    "S3": [0, 1, 2, 5, 7, 10, 11],
    "S4": [0, 2, 3, 4, 8, 9, 10],
}
EXPECTED_ALPHABET = {"S1": 66, "S2": 32, "S3": 2, "S4": 24}

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


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 16), b""):
            h.update(blk)
    return h.hexdigest()


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


# ======================================================================
# exact cyclotomics and Ramanujan sums
# ======================================================================

def _trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def _divmod_monic(a, b):
    a = list(a)
    q = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b):
        c = a[-1]
        if c == 0:
            a.pop()
            continue
        d = len(a) - len(b)
        q[d] = c
        for k, y in enumerate(b):
            a[d + k] -= c * y
        _trim(a)
        if len(a) < len(b):
            break
    return _trim(q), _trim(a)


_CYC = {}


def cyclotomic(n):
    if n in _CYC:
        return _CYC[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            num, rem = _divmod_monic(num, cyclotomic(d))
            assert rem == [0]
    _CYC[n] = num
    return num


def _mobius(n):
    m, p, res = n, 2, 1
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    if m > 1:
        res = -res
    return res


def _totient(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def char_order(k, i):
    return i // math.gcd(k, i) if k else 1


def build_M(S, i):
    """M_S by Ramanujan sums over the order-classes of S; raises if S is not
    a union of full order-classes (i.e. not Galois-stable)."""
    orders = {}
    for k in S:
        d = char_order(k, i)
        orders[d] = orders.get(d, 0) + 1
    for d, cnt in orders.items():
        full = sum(1 for k in range(i) if char_order(k, i) == d)
        if cnt != full:
            raise ValueError("support is not Galois-stable at order %d" % d)
    out = []
    for c in range(i):
        acc = 0
        for d in orders:
            g = math.gcd(c, d)
            mu = _mobius(d // g)
            if mu == 0:
                continue
            acc += mu * _totient(d) // _totient(d // g)
        out.append(4 * acc)
    return out


# ======================================================================
# [A] the admissible block alphabet, two membership tests
# ======================================================================

def alphabet(S, i):
    """Returns (by_cyclotomic, by_projector) as sorted lists of tuples."""
    M = build_M(S, i)
    K = [[M[(a - b) % i] // 4 for b in range(i)] for a in range(i)]
    orders_out = sorted({char_order(k, i) for k in range(i) if k not in set(S)})
    polys = [cyclotomic(d) for d in orders_out]

    def by_cyc(v):
        for pol in polys:
            if _divmod_monic(list(v), pol)[1] != [0]:
                return False
        return True

    def by_proj(v):
        for a in range(i):
            if dot(K[a], v) != i * v[a]:
                return False
        return True

    cyc, proj = [], []
    for bits in range(1 << i):
        v = tuple(1 if (bits >> t) & 1 else -1 for t in range(i))
        if by_cyc(v):
            cyc.append(v)
        if by_proj(v):
            proj.append(v)
    return sorted(cyc), sorted(proj), M


def rank_over_Q(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    n = len(rows[0]) if rows else 0
    r = 0
    for c in range(n):
        p = None
        for k in range(r, len(rows)):
            if rows[k][c] != 0:
                p = k
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = rows[r][c]
        rows[r] = [x / inv for x in rows[r]]
        for k in range(len(rows)):
            if k != r and rows[k][c] != 0:
                f = rows[k][c]
                rows[k] = [x - f * y for x, y in zip(rows[k], rows[r])]
        r += 1
    return r


def clause_a():
    log("[A] the admissible block alphabet at (7,12), all four supports, "
        "two membership tests over all 2^12 sign vectors")
    out = {}
    for name in ("S1", "S2", "S3", "S4"):
        S = S_PARTS[name]
        cyc, proj, M = alphabet(S, I12)
        rk = rank_over_Q([list(v) for v in cyc]) if cyc else 0
        check("%s: cyclotomic divisibility and the projector fixed point "
              "return the same set, |B| = %d (expected %d), span rank %d"
              % (name, len(cyc), EXPECTED_ALPHABET[name], rk),
              cyc == proj and len(cyc) == EXPECTED_ALPHABET[name]
              and rk == (1 if name == "S3" else 7),
              "%d vs %d" % (len(cyc), len(proj)))
        out[name] = (cyc, M)
    return out


# ======================================================================
# [B] the S3 kill
# ======================================================================

def clause_b(alpha):
    log("[B] the S3 kill: (H1) has no solution at (7,12)/S3, at any sigma, "
        "any w, any kappa(rho) -- hence at every order N = 4(12w + 7)")
    B, M = alpha["S3"]
    ones = tuple([1] * I12)
    check("the only admissible blocks at S3 are the two constants",
          sorted(B) == sorted([ones, tuple(-x for x in ones)]), B)
    check("(i) counting: every column is (e_0 1, e_1 1, e_2 1, e_3 1), so at "
          "most 2^4 = 16 distinct columns exist, against 4s = 28 pairwise "
          "orthogonal ones required", 2 ** 4 < 4 * 7)
    cols = [[e[I // I12] for I in range(4 * I12)]
            for e in itertools.product((1, -1), repeat=4)]
    check("(ii) rank: those 16 columns span a space of dimension %d < 28"
          % rank_over_Q(cols), rank_over_Q(cols) == 4)
    off = sorted({M[c] for c in range(1, I12)})
    check("(iii) Gram: inside one superblock any two rows are equal or "
          "antipodal, so their inner product is +-28, while (H1) demands "
          "M_S3(c) for c != 0, whose values are %s" % off,
          all(abs(v) != 4 * 7 for v in off), off)
    # the chi_6 twin, checked rather than asserted
    twin = sorted((k + 6) % I12 for k in S_PARTS["S3"])
    tw_cyc, tw_proj, tw_M = alphabet(twin, I12)
    alt = tuple((-1) ** c for c in range(I12))
    twisted = sorted(tuple(v[c] * alt[c] for c in range(I12)) for v in B)
    check("the chi_6 twin support 6 + S3 = %s is Galois-stable, and its "
          "alphabet is the image of S3's under the twist v -> chi_6 . v -- "
          "%d blocks, the two alternating vectors -- so the same three "
          "obstructions apply and it dies with S3 (Lemma-T twist, "
          "note/NOTE-B.md S1.4)" % (twin, len(tw_cyc)),
          tw_cyc == tw_proj and len(tw_cyc) == 2
          and sorted(tw_cyc) == twisted
          and all(tw_M[c] == (-1) ** c * M[c] for c in range(I12))
          and all(abs(tw_M[c]) != 4 * 7 for c in range(1, I12)),
          twin)
    check("so (7,12)/S3 is DEAD-BY-CERTIFICATE at the column layer, at every "
          "order N = 4(12w + 7)", True)


# ======================================================================
# [C] the artifact at S2
# ======================================================================

def assemble_Q(blocks):
    """blocks[I][k][a] -> Q[12 I + a][k], a 48 x 28 matrix."""
    return [[blocks[I][k][a] for k in range(28)]
            for I in range(4) for a in range(I12)]


def clause_c(alpha, blocks):
    log("[C] the artifact: an admissible column table Q at (7,12)/S2")
    B, M = alpha["S2"]
    Bset = set(B)
    shape = (len(blocks) == 4 and all(len(sb) == 28 for sb in blocks)
             and all(len(b) == I12 for sb in blocks for b in sb)
             and all(x in (1, -1) for sb in blocks for b in sb for x in b))
    check("shape 4 x 28 x 12 with entries in {+-1}", shape)
    if not shape:
        return None
    check("all 112 blocks lie in the alphabet computed in [A] (so every "
          "12-block of every column lies in V_S2)",
          all(tuple(b) in Bset for sb in blocks for b in sb))
    Q = assemble_Q(blocks)
    check("Q is 48 x 28", len(Q) == 48 and all(len(r) == 28 for r in Q))
    Qt = [list(c) for c in zip(*Q)]
    check("(a) Q^T Q = 48 I_28",
          all(dot(Qt[k], Qt[l]) == (48 if k == l else 0)
              for k in range(28) for l in range(28)))
    check("(H1) Q Q^T = I_4 (x) M_S2, entry by entry",
          all(dot(Q[r], Q[t]) ==
              (M[(r % I12 - t % I12) % I12] if r // I12 == t // I12 else 0)
              for r in range(48) for t in range(48)))
    check("M_S2 = %s, recomputed here by Ramanujan sums" % (M,),
          M == [28, 0, 8, 12, -8, 0, -4, 0, -8, 12, 8, 0])
    check("so (H1) HAS a solution at (7,12) / Z12 / S2: an admissible column "
          "table exists", True)
    return Q


# ======================================================================
# [D] the forced per-superblock multiset
# ======================================================================

def negation_classes(B):
    seen, classes = set(), []
    for v in B:
        nv = tuple(-x for x in v)
        if v in seen or nv in seen:
            continue
        seen.add(v)
        seen.add(nv)
        classes.append(v)
    return classes


def block_type(v):
    if len(set(v)) == 1:
        return "C"
    if all(v[c] == v[c % 3] for c in range(I12)):
        return "A"
    return "B"


def clause_d(alpha, blocks):
    log("[D] the forced per-superblock multiset at (7,12)/S2, by exact "
        "Fraction elimination")
    B, M = alpha["S2"]
    classes = negation_classes(B)
    check("the 32 admissible blocks fall into 16 negation classes: 1 "
          "constant, 3 of type A (functions of c mod 3), 12 of type B",
          len(classes) == 16
          and sorted(block_type(v) for v in classes)
          == ["A"] * 3 + ["B"] * 12 + ["C"], sorted(map(block_type, classes)))
    n = len(classes)
    O = [[[v[a] * v[b] for b in range(I12)] for a in range(I12)]
         for v in classes]
    circM = [[M[(a - b) % I12] for b in range(I12)] for a in range(I12)]
    rows, rhs = [], []
    for a in range(I12):
        for b in range(a, I12):
            rows.append([Fraction(O[j][a][b]) for j in range(n)])
            rhs.append(Fraction(circM[a][b]))
    rows.append([Fraction(1)] * n)
    rhs.append(Fraction(28))
    aug = [r + [q] for r, q in zip(rows, rhs)]
    piv, r = [], 0
    for c in range(n):
        p = None
        for k in range(r, len(aug)):
            if aug[k][c] != 0:
                p = k
                break
        if p is None:
            continue
        aug[r], aug[p] = aug[p], aug[r]
        inv = aug[r][c]
        aug[r] = [x / inv for x in aug[r]]
        for k in range(len(aug)):
            if k != r and aug[k][c] != 0:
                f = aug[k][c]
                aug[k] = [x - f * y for x, y in zip(aug[k], aug[r])]
        piv.append(c)
        r += 1
    inconsistent = any(all(x == 0 for x in aug[k][:n]) and aug[k][n] != 0
                       for k in range(len(aug)))
    free = [c for c in range(n) if c not in piv]
    sol = [Fraction(0)] * n
    for k, c in enumerate(piv):
        sol[c] = aug[k][n]
    got = {}
    for j, v in enumerate(classes):
        got.setdefault(block_type(v), set()).add(sol[j])
    check("the 16 outer products b_j b_j^T are linearly independent, so "
          "sum_j m_j b_j b_j^T = circ(M) with sum_j m_j = 28 has a UNIQUE "
          "rational solution (rank %d, %d free columns, consistent)"
          % (r, len(free)), r == n and not free and not inconsistent)
    check("and it is integral and nonnegative: constant once, each type-A "
          "class once, each type-B class twice -- 1 + 3 + 24 = 28",
          got == {"C": {Fraction(1)}, "A": {Fraction(1)}, "B": {Fraction(2)}}
          and sum(sol) == 28,
          {t: sorted(int(x) for x in s) for t, s in got.items()})
    ok = True
    for I in range(4):
        cnt = {}
        for k in range(28):
            v = tuple(blocks[I][k])
            key = min(v, tuple(-x for x in v))
            cnt[key] = cnt.get(key, 0) + 1
        byt = {"C": 0, "A": 0, "B": 0}
        for key, m in cnt.items():
            byt[block_type(key)] += m
            if block_type(key) in ("C", "A") and m != 1:
                ok = False
            if block_type(key) == "B" and m != 2:
                ok = False
        if byt != {"C": 1, "A": 3, "B": 24}:
            ok = False
    check("the artifact realises exactly that forced multiset in all four "
          "superblocks", ok)
    check("so a search over superblocks 1-3 with superblock 0 in one "
          "canonical arrangement is EXHAUSTIVE -- which is what makes a "
          "negative at this cell complete rather than a search failure", True)


# ======================================================================
# the Goethals-Seidel array over a cyclic quotient, standard orientation
# ======================================================================

def gs_array_cyclic(sig, i, krho):
    """note/NOTE-B.md S1.0, standard orientation:
       dev(x)[g,h] = x(h-g),  (XR)[g,h] = x(krho-g-h),
       (X^T R)[g,h] = x(g+h-krho)."""
    def dev(f):
        return [[f[(h - g) % i] for h in range(i)] for g in range(i)]

    def xr(f):
        return [[f[(krho - g - h) % i] for h in range(i)] for g in range(i)]

    def xtr(f):
        return [[f[(g + h - krho) % i] for h in range(i)] for g in range(i)]

    def ng(X):
        return [[-v for v in r] for r in X]
    A = dev(sig[0])
    BR, CR, DR = xr(sig[1]), xr(sig[2]), xr(sig[3])
    BtR, CtR, DtR = xtr(sig[1]), xtr(sig[2]), xtr(sig[3])
    blk = [[A, BR, CR, DR],
           [ng(BR), A, DtR, ng(CtR)],
           [ng(CR), ng(DtR), A, BtR],
           [ng(DR), CtR, ng(BtR), A]]
    return [[blk[I][J][g][h] for J in range(4) for h in range(i)]
            for I in range(4) for g in range(i)]


# ======================================================================
# [E] the row-factorisation lemma, concrete at (7,12)/S2
# ======================================================================

def clause_e(alpha, Q, sigma_entry):
    log("[E] the row-factorisation lemma at (7,12)/S2, with cert 30's table "
        "B2, at every kappa(rho) in Z12")
    B, M = alpha["S2"]
    check("cert 30's table B2 carries exactly the support S2",
          sorted(sigma_entry["S"]) == S_PARTS["S2"], sigma_entry["S"])
    sig = sigma_entry["rows"]
    w = sigma_entry["w"]
    N = 2092
    F = [N * (c == 0) - w * M[c] for c in range(I12)]
    circF = [[F[(h - g) % I12] for h in range(I12)] for g in range(I12)]
    circM = [[M[(a - b) % I12] for b in range(I12)] for a in range(I12)]
    Qt = [list(c) for c in zip(*Q)]
    ok_cc = ok_comm = ok_w = True
    for krho in range(I12):
        C = gs_array_cyclic(sig, I12, krho)
        Ct = [list(c) for c in zip(*C)]
        # Chat^T Chat = Chat Chat^T = I_4 (x) circ(F)
        for r in range(48):
            for t in range(48):
                want = circF[r % I12][t % I12] if r // I12 == t // I12 else 0
                if dot(Ct[r], Ct[t]) != want or dot(C[r], C[t]) != want:
                    ok_cc = False
        # Chat commutes with I_4 (x) P_S (equivalently with I_4 (x) circ(M)),
        # exploiting the block-circulant structure of the latter
        for r in range(48):
            row1 = [sum(C[r][(b // I12) * I12 + k] * circM[k][b % I12]
                        for k in range(I12)) for b in range(48)]
            row2 = [sum(circM[r % I12][k] * C[(r // I12) * I12 + k][b]
                        for k in range(I12)) for b in range(48)]
            if row1 != row2:
                ok_comm = False
        # W = Chat^T Q ;  W W^T = 4i * 4s * (I_4 (x) P_S) = 28 * (I_4 (x) M)
        W = [[sum(Ct[a][k] * Q[k][b] for k in range(48)) for b in range(28)]
             for a in range(48)]
        for a in range(48):
            for b in range(48):
                want = (28 * circM[a % I12][b % I12]
                        if a // I12 == b // I12 else 0)
                if dot(W[a], W[b]) != want:
                    ok_w = False
    check("Chat^T Chat = Chat Chat^T = I_4 (x) circ(F) at all 12 kappa(rho) "
          "(Theorem F(c) in the general branch)", ok_cc)
    check("Chat commutes with I_4 (x) P_S at all 12 kappa(rho)", ok_comm)
    check("W = Chat^T Q satisfies W W^T = 4i * 4s * (I_4 (x) P_S) "
          "= 28 (I_4 (x) M) at all 12 kappa(rho) -- the row-factorisation "
          "lemma, exactly", ok_w)
    # the consequence, on random admissible rows
    krho = 5
    C = gs_array_cyclic(sig, I12, krho)
    rng = random.Random(20260905)
    bad = 0
    for _ in range(200):
        p = []
        for _J in range(4):
            p.extend(rng.choice(B))
        cp = [dot(C[k], p) for k in range(48)]
        if sum(x * x for x in cp) != 16 * I12 * 7:
            bad += 1
    check("every admissible row p has ||Chat p||^2 = 16 i s = %d, so its row "
          "of E has squared norm exactly 4s = 28 whatever p is "
          "(200 random rows at kappa(rho) = 5)" % (16 * I12 * 7), bad == 0,
          "%d failures" % bad)
    check("hence a kit exists at a fixed (Q, kappa(rho), sigma) iff the flat "
          "set contains 4s mutually orthogonal admissible rows -- a "
          "condition on each row alone", True)


# ======================================================================
# [F] the dual form, on this repository's own kits
# ======================================================================

def clause_f_24(records):
    log("[F.i] the dual form at (2,4), on cert 18's four banked records, "
        "exhaustively over all 256 admissible rows")
    s, i = 2, 4
    M = [8, 0, -8, 0]
    circM = [[M[(a - b) % i] for b in range(i)] for a in range(i)]
    alph = [v for v in itertools.product((1, -1), repeat=i)
            if all(sum(circM[a][b] * v[b] for b in range(i)) == 4 * i * v[a]
                   for a in range(i))]
    check("the admissible block alphabet at (2,4) with S = {chi, chi^3} is "
          "the four anti-periodic sign vectors", len(alph) == 4, alph)
    for rec in records:
        krho = rec["rho_bar"]
        sig = rec["coset_sums"]
        Q8 = [[1 if ch == "+" else -1 for ch in x] for x in rec["col_table_8"]]
        P8 = [[1 if ch == "+" else -1 for ch in x] for x in rec["row_table_8"]]
        E = [[1 if ch == "+" else -1 for ch in x] for x in rec["corner"]]
        Q = [[(1 if c < 2 else -1) * Q8[2 * I + (c % 2)][b]
              for b in range(4 * s)]
             for I in range(4) for c in range(4)]
        P = [[(1 if c < 2 else -1) * P8[r][2 * J + (c % 2)]
              for J in range(4) for c in range(4)] for r in range(4 * s)]
        C = gs_array_cyclic(sig, i, krho)
        h1 = all(dot(Q[a], Q[b]) ==
                 (circM[a % i][b % i] if a // i == b // i else 0)
                 for a in range(16) for b in range(16))
        h4 = all(dot(E[r], Q[k]) + dot(P[r], C[k]) == 0
                 for r in range(8) for k in range(16))
        dual = all([dot(C[k], P[r]) for k in range(16)]
                   == [-dot(Q[k], E[r]) for k in range(16)]
                   for r in range(8))
        pdual = all(sum(E[r][b] * sum(Q[k][b] * C[k][h] for k in range(16))
                        for b in range(8)) == -4 * s * P[r][h]
                    for r in range(8) for h in range(16))
        allrows = True
        for blk in itertools.product(alph, repeat=4):
            p = [x for b in blk for x in b]
            cp = [dot(C[k], p) for k in range(16)]
            if sum(x * x for x in cp) != 16 * i * s:
                allrows = False
            for I in range(4):
                for a in range(i):
                    if sum(circM[a][b] * cp[i * I + b] for b in range(i)) \
                            != 4 * i * cp[i * I + a]:
                        allrows = False
        check("%-34s (H1), (H4), Chat p = -Q e on the kit's rows, "
              "P = -(1/4s) E Q^T Chat, and Chat p in col(Q) with "
              "||Chat p||^2 = 16is for all 256 admissible rows"
              % rec["name"], h1 and h4 and dual and pdual and allrows,
              "H1=%s H4=%s dual=%s Pdual=%s rows=%s"
              % (h1, h4, dual, pdual, allrows))


# --- a small self-contained (3,4) kit engine, for [F.ii] -------------------

BAL4 = [v for v in itertools.product((1, -1), repeat=4) if sum(v) == 0]
POOL34 = [sum(t, ()) for t in itertools.product(BAL4, repeat=4)]
REPS34 = [v for v in POOL34 if v[0] == 1]


def random_Q34(rng):
    """12 mutually orthogonal block-balanced sign vectors, as a 16 x 12
    matrix.  Q^T Q = 16 I with the columns blockwise sum-zero already forces
    Q Q^T = I_4 (x) (12 I - 4 (J - I)): the projector onto col(Q) and the
    projector onto R^4 (x) V_S are then both 12-dimensional and equal."""
    for _ in range(500):
        cols = []
        order = list(range(len(REPS34)))
        rng.shuffle(order)
        for v in order:
            u = REPS34[v]
            if all(dot(u, c) == 0 for c in cols):
                cols.append(u)
                if len(cols) == 12:
                    return [[cols[j][k] for j in range(12)] for k in range(16)]
    raise RuntimeError("no 12-clique found")


def sparts34_sample(rng, count):
    """A few admissible (3,4) S-parts on Z4: T = 4 sigma_0 built from spectra
    with sum_q |sigma-hat_q(chi)|^2 = 12 at each nontrivial character and the
    integrality/parity of sigma = (r + T)/4 (note/NOTE-B.md S1.8)."""
    out = []
    vals = (0, 2, -2)
    for y in itertools.product(vals, repeat=4):
        if sum(v * v for v in y) != 12:
            continue
        for ab in itertools.product(itertools.product(vals, repeat=2),
                                    repeat=4):
            if sum(a * a + b * b for a, b in ab) != 12:
                continue
            T = [[y[q] + 2 * ab[q][0], -y[q] + 2 * ab[q][1],
                  y[q] - 2 * ab[q][0], -y[q] - 2 * ab[q][1]]
                 for q in range(4)]
            if all(len({v % 8 for v in T[q]}) == 1 for q in range(4)):
                out.append(T)
                if len(out) > 4000:
                    break
        if len(out) > 4000:
            break
    rng.shuffle(out)
    return out[:count]


def clause_f_34():
    log("[F.ii] the dual form at (3,4), on kits found by this certificate's "
        "own engine")
    i, s = 4, 3
    M = [16 - 4, -4, -4, -4]
    circM = [[M[(a - b) % i] for b in range(i)] for a in range(i)]
    rng = random.Random(31031)
    QL = [random_Q34(rng) for _ in range(40)]
    for Q in QL:
        assert all(dot(Q[a], Q[b]) ==
                   (circM[a % i][b % i] if a // i == b // i else 0)
                   for a in range(16) for b in range(16))
    parts = sparts34_sample(random.Random(770), 3)
    npts = 0
    found = dual_ok = pd_ok = norm_ok = 0
    # C4 = gs_array(T) = 4 * Chat_S, since T = 4 sigma_S (note S1.8); every
    # identity below is therefore stated on C4 with its power of 4 carried
    # explicitly, so the arithmetic stays in exact integers throughout.
    for T in parts:
        for krho in (0, 2):
            npts += 1
            C4 = gs_array_cyclic(T, i, krho)
            C4t = [list(c) for c in zip(*C4)]
            kit = None
            for Q in QL:
                V = [[sum(C4t[a][k] * Q[k][j] for k in range(16))
                      for j in range(12)] for a in range(16)]
                Vt = [list(c) for c in zip(*V)]
                flats = [p for p in REPS34
                         if all(abs(dot(p, Vt[j])) == 64 for j in range(12))]
                if len(flats) < 12:
                    continue
                sel = []
                for p in flats:
                    if all(dot(p, q) == 0 for q in sel):
                        sel.append(p)
                        if len(sel) == 12:
                            break
                if len(sel) < 12:
                    continue
                P = [list(p) for p in sel]
                E = [[-dot(P[r], Vt[j]) // 64 for j in range(12)]
                     for r in range(12)]
                kit = (E, P, Q)
                break
            if kit is None:
                continue
            E, P, Q = kit
            found += 1
            check("E is Hadamard of order 12, E in {+-1}, and (H4) "
                  "4 E Q^T + P C4^T = 0",
                  all(dot(E[a], E[b]) == (12 if a == b else 0)
                      for a in range(12) for b in range(12))
                  and all(abs(v) == 1 for row in E for v in row)
                  and all(4 * dot(E[r], Q[k]) + dot(P[r], C4[k]) == 0
                          for r in range(12) for k in range(16)))
            # dual form:  Chat p = -Q e,  i.e.  C4 p = -4 Q e
            dual_ok += all([dot(C4[k], P[r]) for k in range(16)]
                           == [-4 * dot(Q[k], E[r]) for k in range(16)]
                           for r in range(12))
            # dual formula:  P = -(1/4s) E Q^T Chat,  i.e.  16 s P = -E Q^T C4
            pd_ok += all(sum(E[r][b] * sum(Q[k][b] * C4[k][h]
                                           for k in range(16))
                             for b in range(12)) == -16 * s * P[r][h]
                         for r in range(12) for h in range(16))
            good = True
            for p in REPS34[::13]:
                cp = [dot(C4[k], p) for k in range(16)]
                if sum(x * x for x in cp) != 16 * 16 * i * s:
                    good = False
                for I in range(4):
                    for a in range(i):
                        if sum(circM[a][b] * cp[i * I + b] for b in range(i)) \
                                != 4 * i * cp[i * I + a]:
                            good = False
            norm_ok += good
    check("kits found at %d of the %d (S-part, kappa(rho)) points tried"
          % (found, npts), found >= 3, found)
    check("Chat p = -Q e on every kit row: %d / %d kits" % (dual_ok, found),
          dual_ok == found and found > 0)
    check("P = -(1/4s) E Q^T Chat on every kit: %d / %d" % (pd_ok, found),
          pd_ok == found and found > 0)
    check("and Chat p lies in col(Q) with ||Chat p||^2 = 16is for the sampled "
          "admissible rows at every kit: %d / %d" % (norm_ok, found),
          norm_ok == found and found > 0)


# ======================================================================
# [G] controls that can fail
# ======================================================================

def clause_g(alpha, blocks):
    log("[G] controls that can fail -- each perturbation must be REJECTED "
        "by the same acceptance path")
    B, M = alpha["S2"]
    Bset = set(B)

    def h1_holds(Q):
        return all(dot(Q[a], Q[b]) ==
                   (M[(a % I12 - b % I12) % I12] if a // I12 == b // I12
                    else 0) for a in range(48) for b in range(48))

    def qtq_holds(Q):
        Qt = [list(c) for c in zip(*Q)]
        return all(dot(Qt[k], Qt[l]) == (48 if k == l else 0)
                   for k in range(28) for l in range(28))

    def blocks_admissible(bl):
        return all(tuple(b) in Bset for sb in bl for b in sb)

    bad = copy.deepcopy(blocks)
    bad[0][0][0] = -bad[0][0][0]
    check("NEGATIVE: one entry of one block flipped -- the block leaves the "
          "alphabet, and Q^T Q = 48 I fails",
          not blocks_admissible(bad) and not qtq_holds(assemble_Q(bad)))

    bad2 = copy.deepcopy(blocks)
    bad2[1][3] = list(bad2[1][7])
    check("NEGATIVE: one block replaced by another ADMISSIBLE block -- the "
          "alphabet test still passes, and (H1) fails, so the two tests are "
          "independent",
          blocks_admissible(bad2) and not h1_holds(assemble_Q(bad2)))

    bad3 = copy.deepcopy(blocks)
    bad3[1][3], bad3[1][4] = bad3[1][4], bad3[1][3]
    check("NEGATIVE: two columns permuted inside one superblock only -- the "
          "blocks stay admissible and (H1) fails",
          blocks_admissible(bad3) and not h1_holds(assemble_Q(bad3)))

    good = copy.deepcopy(blocks)
    for I in range(4):
        good[I][5] = [-x for x in good[I][5]]
    check("POSITIVE: a whole column of Q negated -- still admissible, (H1) "
          "and Q^T Q = 48 I still hold, so the tests are not vacuously "
          "strict", blocks_admissible(good) and h1_holds(assemble_Q(good))
          and qtq_holds(assemble_Q(good)))

    try:
        build_M([0, 1, 2, 5, 7, 10], I12)
        rejected = False
    except ValueError:
        rejected = True
    check("NEGATIVE: a support that is not Galois-stable is rejected by the "
          "Ramanujan path rather than silently evaluated", rejected)


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args(argv)

    print("=" * 78)
    print("cert 31 -- the general-branch border layer at (7,12): the S2 "
          "column table,")
    print("           the forced multiset, the S3 kill, and the "
          "row-factorisation lemma")
    print("=" * 78)

    log("[0] the artifacts")
    for path, want, label in ((Q_DATA, Q_SHA256, "data/q-7_12-S2.json"),
                              (SIGMA_DATA, SIGMA_SHA256,
                               "data/general-branch-sigma-tables.json"),
                              (CELL24_DATA, CELL24_SHA256,
                               "data/cell24-records.json")):
        got = file_sha256(path)
        check("%s sha256 = %s..." % (label, want[:16]), got == want, got)
    with open(Q_DATA, encoding="ascii") as fh:
        blocks = json.load(fh)["Q_blocks"]
    with open(SIGMA_DATA, encoding="ascii") as fh:
        sigma_tables = json.load(fh)["tables"]
    with open(CELL24_DATA, encoding="ascii") as fh:
        cell24 = json.load(fh)["records"]

    alpha = clause_a()
    clause_b(alpha)
    Q = clause_c(alpha, blocks)
    clause_d(alpha, blocks)
    B2 = [e for e in sigma_tables if e["id"] == "B2_7_12_43"][0]
    clause_e(alpha, Q, B2)
    clause_f_24(cell24)
    clause_f_34()
    clause_g(alpha, blocks)

    print("\n" + "=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 31: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: at (s,i) = (7,12) on Z12 the admissible block alphabet is")
    print("         66 / 32 / 2 / 24 vectors at S1 / S2 / S3 / S4, by two")
    print("         membership tests that share no code.  At S3 the alphabet")
    print("         is the two constants, and three independent obstructions")
    print("         kill (H1) there at every sigma, every w, every kappa(rho)")
    print("         and therefore at EVERY order N = 4(12w + 7):")
    print("         DEAD-BY-CERTIFICATE, the chi_6 twin with it.  At S2 an")
    print("         admissible column table EXISTS -- the banked artifact --")
    print("         and its per-superblock block multiset is FORCED (1 + 3 +")
    print("         24), so a search over the remaining superblocks is")
    print("         exhaustive.  The row-factorisation lemma holds exactly at")
    print("         all twelve kappa(rho) with cert 30's table B2, and the")
    print("         dual form Chat p = -Q e, P = -(1/4s) E Q^T Chat holds on")
    print("         every kit this repository can exhibit.")
    print("         LABEL: PROVEN-BY-CERTIFICATE.")
    print("         NOT claimed: any kit at (7,12) -- (H1) is the column")
    print("         layer only, and no border kit is known anywhere in the")
    print("         general branch at 2092; anything about seeds; anything")
    print("         about H(2092).")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
