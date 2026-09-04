#!/usr/bin/env python3
"""cert 30 -- the general-branch sigma-layer at N = 2092: six coset-sum tables,
complete support enumerations, and the (3,8) collapse onto the house (3,4) cell.

  STATEMENT (note/NOTE-B.md S2.6).  Let N = 2092 = 4 * 523.  A general-branch
  cell is a quotient Gbar of order i, a fibre size w = |K| with N = 4(wi + s),
  and a Galois-stable spectral support S subset of Gbar-hat with |S| = s, so
  that -- by Theorem E' (S1.7), since w > 2s at every cell here -- the only
  admissible Gram is

        M_S(c) = 4 * sum_{chi in S} chi(c) = 4i * P_S,
        F(c)   = N * delta_0(c) - w * M_S(c).

  The SIGMA-LAYER of the cell is the set of integral coset-sum tables
  sigma = (sigma_0, ..., sigma_3) on Gbar with

    (i)   sum_q PAF_{sigma_q}(c) = F(c)  for every c in Gbar;
    (ii)  sigma_q(c) = w (mod 2) and |sigma_q(c)| <= w  (each entry is a sum
          of w signs);
    (iii) sum_{q,c} sigma_q(c)^2 = F(0)  and  sum_q (sum_c sigma_q(c))^2 = 4s.

  THIS CERTIFICATE ESTABLISHES, by exact integer arithmetic:

  [A] The sigma-layer of every advertised general branch at 2092 is NON-EMPTY,
      by explicit artifact -- six tables:
          A            (s,i,w) = (11,16,32)   Gbar = Z16    row sums (6,2,2,0)
          B1..B4       (s,i,w) = ( 7,12,43)   Gbar = Z12    row sums (4,2,2,2)
          C            (s,i,w) = ( 3, 8,65)   Gbar = F_2^3  row sums (0,-2,-2,-2)
      M_S and F are RE-DERIVED here from (S, i, w) alone; the listed M and F
      are read only to be compared against what this run computed.

  [B] The support enumerations are COMPLETE.  Exactly ONE Galois-stable
      size-11 subset of Z16-hat contains the trivial character (table A's), and
      exactly FOUR Galois-stable size-7 subsets of Z12-hat do (B1..B4's).  So
      the (7,12) cell is populated on EVERY one of its rational branches, not
      on a lucky one.

  [C] The (3,8) cell is not an independent door.  All 56 size-3 subsets of the
      dual of F_2^3 and all 12 Galois-stable size-3 subsets of the dual of
      Z4 x Z2 are real-character twists of W \\ {1} for a subgroup W of order 4,
      so every such cell coarsens to a house (3,4) profile.  Exhibited: table C
      twisted by the real character 3 and compressed along the line J = {0,4}
      has aggregate PAF [532, 520, 520, 520] = the house (3,4) profile
      [4w + 12, 4w, 4w, 4w] at w = 130, the width of N = 2092 in that cell.

  [D] The theorems these cells sit under, re-checked cell by cell: N = 4(wi+s);
      w > 2s (so Theorem E' -- indeed Theorem E -- applies and M = 4i P_S is a
      conclusion, not an ansatz); M/(4i) idempotent of rank s; and every S
      contains a real character, as Theorem 3 (S2.4) requires at N = 4 (mod 8).

  WHAT THIS CERTIFICATE DOES NOT SAY.  Nothing about SEEDS: a compressed
  sigma-table is a NECESSARY shadow of a seed quadruple on a group of order
  wi, never a sufficient one, and no seed quadruple is known at any of these
  cells.  Nothing about BORDERS (that is cert 31).  Nothing about H(2092).
  Nothing about the number of tables in any sigma-layer.

  TWO IMPLEMENTATIONS.  Every PAF and every character sum is decided twice by
  code paths that share nothing:
      PAF_sigma   direct double sum  sum_j sigma_j sigma_{j+t}
                  vs. row 0 of G G^T for the group matrix G[a][b] = sigma(b-a)
                  (or sigma(a xor b)) -- a matrix is built and multiplied, no
                  autocorrelation formula anywhere in that path;
      M_S cyclic  root-of-unity sums reduced modulo Phi_i(x)
                  vs. RAMANUJAN sums c_d(c) = mu(d/g) phi(d) / phi(d/g),
                  g = gcd(c,d), over the order-classes of S;
      M_S on F_2^3   direct character sum vs. a fast Walsh-Hadamard butterfly.
  A support that is not Galois-stable makes the first path irrational and the
  second path's order-class decomposition fail, so it is REJECTED rather than
  silently evaluated.  The Fourier budget (4s on S, N off S) is a third,
  spectral route to the same aggregate.

  FOUR CONTROLS THAT CAN FAIL are run through the identical acceptance path
  and must be rejected: one table entry moved by +2; one support element
  swapped (breaking Galois stability); a row sum falsified; w changed to 30.
  The run fails if any control is accepted.

Standard library only, exact integers only, no floats anywhere, no network.
Nothing is imported from outside this repository.

Usage:
  python certs/30-general-branch-sigma/run.py
  python certs/30-general-branch-sigma/run.py --verbose
"""

import argparse
import copy
import hashlib
import itertools
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(ROOT, "data", "general-branch-sigma-tables.json")

DATA_SHA256 = "d2c37aace3bf016dc245a26c3fe7a2bd4aa524aafc874ba717c7e562d98c4a97"

N = 2092

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


# ======================================================================
# exact cyclotomic arithmetic (integer coefficient lists, low degree first)
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
    """Phi_n(x) as an exact integer coefficient list."""
    if n in _CYC:
        return _CYC[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            num, rem = _divmod_monic(num, cyclotomic(d))
            assert rem == [0]
    _CYC[n] = num
    return num


def rational_root_sum(coeffs, i):
    """sum_j coeffs[j] * zeta_i^j as a rational integer, or None if the value
    is not rational (which is exactly the Galois-instability signal)."""
    _, rem = _divmod_monic(_trim(list(coeffs)), cyclotomic(i))
    rem = _trim(rem)
    return rem[0] if len(rem) == 1 else None


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


# ======================================================================
# implementation A -- direct
# ======================================================================

def paf_cyclic_A(x):
    n = len(x)
    return [sum(x[j] * x[(j + t) % n] for j in range(n)) for t in range(n)]


def paf_xor_A(x):
    n = len(x)
    return [sum(x[j] * x[j ^ t] for j in range(n)) for t in range(n)]


def M_cyclic_A(S, i):
    """M_S by root-of-unity sums reduced mod Phi_i(x)."""
    out = []
    for c in range(i):
        co = [0] * i
        for k in S:
            co[(k * c) % i] += 1
        v = rational_root_sum(co, i)
        if v is None:
            return None
        out.append(4 * v)
    return out


def M_xor_A(S, nbits):
    return [4 * sum((-1) ** (bin(k & c).count("1") & 1) for k in S)
            for c in range(1 << nbits)]


# ======================================================================
# implementation B -- group matrix / Ramanujan sums / Walsh butterfly
# ======================================================================

def _group_matrix(x, op):
    n = len(x)
    return [[x[op(a, b)] for b in range(n)] for a in range(n)]


def _gram_row0(G):
    """Row 0 of G G^T, by an explicit matrix product.  No autocorrelation
    formula appears anywhere on this path."""
    n = len(G)
    return [sum(G[0][k] * G[t][k] for k in range(n)) for t in range(n)]


def paf_cyclic_B(x):
    n = len(x)
    return _gram_row0(_group_matrix(x, lambda a, b: (b - a) % n))


def paf_xor_B(x):
    return _gram_row0(_group_matrix(x, lambda a, b: a ^ b))


def M_cyclic_B(S, i):
    """M_S by RAMANUJAN sums.  A Galois-stable S is a union of full
    order-classes {k : ord(zeta_i^k) = d}, and the class of order d sums to
    c_d(c) = mu(d/g) phi(d) / phi(d/g) with g = gcd(c,d).  Returns None when S
    is not such a union -- i.e. exactly when S is not Galois-stable."""
    classes = {}
    for k in range(i):
        d = i // math.gcd(k, i) if k else 1
        classes.setdefault(d, set()).add(k)
    Sset = set(S)
    orders = []
    for d, cls in classes.items():
        inter = cls & Sset
        if not inter:
            continue
        if inter != cls:
            return None
        orders.append(d)
    if sum(len(classes[d]) for d in orders) != len(Sset):
        return None
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


def M_xor_B(S, nbits):
    """M_S by a fast Walsh-Hadamard butterfly on the indicator of S."""
    n = 1 << nbits
    a = [0] * n
    for k in S:
        a[k] += 1
    h = 1
    while h < n:
        for base in range(0, n, h * 2):
            for j in range(base, base + h):
                u, v = a[j], a[j + h]
                a[j], a[j + h] = u + v, u - v
        h *= 2
    return [4 * v for v in a]


# ======================================================================
# the single acceptance path
# ======================================================================

def accept(entry, record=None):
    """The one acceptance path for a sigma-table.  Returns True iff every
    condition holds.  Nothing the caller supplies is trusted: M and F are
    recomputed from (S, i, w) and only then compared with the listed ones."""
    local = []

    def sub(label, cond, extra=""):
        if not cond:
            local.append((label, str(extra)))
        return cond

    i, w, s, S = entry["i"], entry["w"], entry["s"], entry["S"]
    rows = entry["rows"]
    kind = entry["kind"]

    if not sub("shape", len(rows) == 4 and all(len(r) == i for r in rows)):
        return False, local
    if not sub("support", len(set(S)) == len(S) == s
               and all(0 <= k < i for k in S)):
        return False, local

    if kind == "cyclic":
        MA, MB = M_cyclic_A(S, i), M_cyclic_B(S, i)
        pafA, pafB = paf_cyclic_A, paf_cyclic_B
    else:
        nb = i.bit_length() - 1
        MA, MB = M_xor_A(S, nb), M_xor_B(S, nb)
        pafA, pafB = paf_xor_A, paf_xor_B
    if not sub("M_two_implementations",
               MA is not None and MB is not None and MA == MB,
               "M_S irrational (S not Galois-stable), or the two "
               "implementations disagree"):
        return False, local
    sub("M_equals_listed", MA == entry["M"], MA)

    F = [N * (c == 0) - w * MA[c] for c in range(i)]
    sub("F_equals_listed", F == entry["F"], F)

    aggA = [0] * i
    aggB = [0] * i
    for r in rows:
        pa, pb = pafA(r), pafB(r)
        sub("PAF_two_implementations", pa == pb, "%s vs %s" % (pa, pb))
        for t in range(i):
            aggA[t] += pa[t]
            aggB[t] += pb[t]
    sub("aggregate_PAF_equals_F", aggA == F and aggB == F, aggA)
    sub("parity", all((v - w) % 2 == 0 for r in rows for v in r))
    sub("box", all(abs(v) <= w for r in rows for v in r))

    rs = [sum(r) for r in rows]
    sub("row_sums", rs == entry["row_sums"], rs)
    sub("total_norm", sum(v * v for r in rows for v in r) == F[0],
        sum(v * v for r in rows for v in r))
    sub("row_sum_shell", sum(v * v for v in rs) == 4 * s,
        sum(v * v for v in rs))

    # the Fourier budget: 4s on S, N off S -- a third, spectral route
    if kind == "cyclic":
        bud = []
        for k in range(i):
            co = [0] * i
            for d in range(i):
                co[(k * d) % i] += F[d]
            bud.append(rational_root_sum(co, i))
    else:
        bud = [sum(F[d] * (-1) ** (bin(k & d).count("1") & 1)
                   for d in range(i)) for k in range(i)]
    want = [4 * s if k in set(S) else N for k in range(i)]
    sub("fourier_budget", bud == want, bud)

    if record is not None:
        record.update(M=MA, F=F, aggregate=aggA, row_sums=rs, budget=bud)
    return (not local), local


# ======================================================================
# [A] the six tables
# ======================================================================

def clause_a(tables, verbose):
    log("[A] the six general-branch sigma-tables, every quantity re-derived "
        "from (S, i, w)")
    for entry in tables:
        rec = {}
        ok, bad = accept(entry, rec)
        check("%-12s (s,i,w) = (%2d,%2d,%3d) on %-5s: aggregate PAF = F, "
              "parity, box, norms, budget"
              % (entry["id"], entry["s"], entry["i"], entry["w"],
                 "Z%d" % entry["i"] if entry["kind"] == "cyclic" else "F_2^3"),
              ok, "; ".join("%s %s" % b for b in bad[:3]))
        if verbose and ok:
            print("        M = %s" % (rec["M"],))
            print("        F = %s" % (rec["F"],))
            print("        budget = %s" % (rec["budget"],))
    # the budgets, stated once
    for entry in tables:
        rec = {}
        accept(entry, rec)
        on = sorted({rec["budget"][k] for k in entry["S"]})
        off = sorted({rec["budget"][k] for k in range(entry["i"])
                      if k not in set(entry["S"])})
        check("%-12s Fourier budget is %s on S and %s off S"
              % (entry["id"], on, off),
              on == [4 * entry["s"]] and off == [N])


# ======================================================================
# [B] support uniqueness
# ======================================================================

def galois_orbits(i):
    units = [k for k in range(1, i) if math.gcd(k, i) == 1] or [1]
    seen, orbits = set(), []
    for k in range(i):
        if k in seen:
            continue
        o = frozenset((k * u) % i for u in units)
        seen |= o
        orbits.append(o)
    return orbits


def galois_stable_supports(i, size, must_contain_trivial=True):
    orbits = galois_orbits(i)
    out = []
    for r in range(1, len(orbits) + 1):
        for comb in itertools.combinations(orbits, r):
            if sum(len(o) for o in comb) != size:
                continue
            u = set()
            for o in comb:
                u |= set(o)
            if must_contain_trivial and 0 not in u:
                continue
            out.append(sorted(u))
    return sorted(out)


def clause_b(tables):
    log("[B] the support enumerations are complete")
    s11 = galois_stable_supports(16, 11)
    check("exactly one Galois-stable size-11 subset of Z16-hat contains the "
          "trivial character", len(s11) == 1, s11)
    banked_A = sorted(e["S"] for e in tables if e["i"] == 16)
    check("and it is table A's support", s11 == banked_A, s11)
    s7 = galois_stable_supports(12, 7)
    check("exactly four Galois-stable size-7 subsets of Z12-hat contain the "
          "trivial character", len(s7) == 4, s7)
    banked_B = sorted(sorted(e["S"]) for e in tables if e["i"] == 12)
    check("and they are exactly B1..B4's supports -- so (7,12) is populated "
          "on every rational branch", sorted(s7) == banked_B, s7)
    # the orbit structure that makes the count what it is
    sizes = sorted(len(o) for o in galois_orbits(12))
    check("Z12-hat has Galois orbits of sizes %s, and exactly eight unions "
          "total 7 (four contain chi_0, four are their chi_6 twists)" % (sizes,),
          sizes == [1, 1, 2, 2, 2, 4]
          and len(galois_stable_supports(12, 7, False)) == 8)


# ======================================================================
# [C] the (3,8) collapse and the transport
# ======================================================================

def size3_supports_collapse(mods):
    """Every Galois-stable size-3 subset of the dual of prod Z_mods is a real
    twist of W \\ {1} for some subgroup W of order 4.  Returns
    (n_stable, n_collapsing, n_subgroups_of_order_4)."""
    els = list(itertools.product(*[range(m) for m in mods]))
    exp = 1
    for m in mods:
        exp = exp * m // math.gcd(exp, m)
    units = [k for k in range(1, exp + 1) if math.gcd(k, exp) == 1]

    def smul(k, x):
        return tuple((k * xi) % m for xi, m in zip(x, mods))

    def add(x, y):
        return tuple((a + b) % m for a, b, m in zip(x, y, mods))

    ident = tuple(0 for _ in mods)
    seen, orbits = set(), []
    for x in els:
        if x in seen:
            continue
        o = frozenset(smul(k, x) for k in units)
        seen |= o
        orbits.append(o)
    stable = set()
    for r in range(1, 4):
        for comb in itertools.combinations(orbits, r):
            if sum(len(o) for o in comb) != 3:
                continue
            u = set()
            for o in comb:
                u |= set(o)
            stable.add(frozenset(u))
    reals = [x for x in els if smul(2, x) == ident]
    subs4 = set()
    for a, b in itertools.combinations(els, 2):
        cur = {ident}
        while True:
            new = set(cur)
            for u in cur:
                for v in (a, b):
                    new.add(add(u, v))
            if new == cur:
                break
            cur = new
            if len(cur) > 8:
                break
        if len(cur) == 4:
            subs4.add(frozenset(cur))
    good = 0
    for S in stable:
        for rho in reals:
            tw = frozenset(add(rho, x) for x in S)
            if any(tw == frozenset(W) - {ident} for W in subs4):
                good += 1
                break
    return len(stable), good, len(subs4)


def clause_c(tables, transport):
    log("[C] the (3,8) cell is not an independent door -- it coarsens onto "
        "the house (3,4) cell")
    n1, k1, w1 = size3_supports_collapse([2, 2, 2])
    check("all %d size-3 supports of the dual of F_2^3 are real twists of "
          "W \\ {1} for one of its %d subgroups of order 4" % (n1, w1),
          n1 == k1 == 56)
    n2, k2, w2 = size3_supports_collapse([4, 2])
    check("all %d Galois-stable size-3 supports of the dual of Z4 x Z2 are "
          "too (%d subgroups of order 4)" % (n2, w2), n2 == k2 == 12)

    C = [e for e in tables if e["id"] == "C_3_8_65"][0]
    rho = transport["twist_character"]
    J = transport["compressed_line"]
    wH = transport["house_w"]
    check("the compression line J = {0,4} is a subgroup of F_2^3 of order 2",
          sorted(J) == [0, 4] and (J[0] ^ J[1]) in J)
    tw = [[r[x] * (-1) ** (bin(rho & x).count("1") & 1) for x in range(8)]
          for r in C["rows"]]
    comp = [[tw[q][c] + tw[q][c + 4] for c in range(4)] for q in range(4)]
    check("the twisted, compressed table is the one banked", comp ==
          transport["expected_table"], comp)
    agg = [0] * 4
    for r in comp:
        pa, pb = paf_xor_A(r), paf_xor_B(r)
        check("transport PAF, two implementations", pa == pb)
        for t in range(4):
            agg[t] += pa[t]
    house = [4 * wH + 12] + [4 * wH] * 3
    check("its aggregate PAF %s is the house (3,4) compressed profile "
          "[4w+12, 4w, 4w, 4w] at w = %d (note/NOTE-B.md S1.8)" % (agg, wH),
          agg == house, house)
    check("and w = %d is the width of N = 2092 in the (3,4) cell "
          "(4(4w + 3) = 2092)" % wH, 4 * (4 * wH + 3) == N)


# ======================================================================
# [D] the theorems these cells sit under
# ======================================================================

def clause_d(tables):
    log("[D] the theorems the cells sit under, re-checked cell by cell")
    for entry in tables:
        i, w, s, S = entry["i"], entry["w"], entry["s"], entry["S"]
        rec = {}
        accept(entry, rec)
        M = rec["M"]
        check("%-12s N = 4(wi + s) = %d" % (entry["id"], 4 * (w * i + s)),
              4 * (w * i + s) == N)
        check("%-12s w = %d > 2s = %d, so Theorem E (a fortiori E') applies "
              "and M = 4i P_S is a conclusion" % (entry["id"], w, 2 * s),
              w > 2 * s)
        # M / 4i is an idempotent of rank s
        if entry["kind"] == "cyclic":
            def conv(a, b):
                return [sum(a[(c - d) % i] * b[d] for d in range(i))
                        for c in range(i)]
        else:
            def conv(a, b):
                return [sum(a[c ^ d] * b[d] for d in range(i))
                        for c in range(i)]
        MM = conv(M, M)
        check("%-12s M * M = 4i M (M/4i idempotent) and trace = 4i s"
              % entry["id"],
              MM == [4 * i * v for v in M] and i * M[0] == 4 * i * s)
        # Theorem 3 (S2.4): at N = 4 (mod 8) every realisable Gram has a live
        # real character.  chi_0 is real and lies in every S here.
        if entry["kind"] == "cyclic":
            reals = [k for k in range(i) if (2 * k) % i == 0]
        else:
            reals = list(range(i))
        check("%-12s S contains a real character, as Theorem 3 (S2.4) "
              "requires at N = 4 (mod 8)" % entry["id"],
              N % 8 == 4 and any(k in reals for k in S))


# ======================================================================
# [E] controls that must be rejected
# ======================================================================

def clause_e(tables):
    log("[E] four controls that can fail, through the identical acceptance "
        "path -- each must be REJECTED")
    base = [e for e in tables if e["id"] == "A_11_16_32"][0]
    c1 = copy.deepcopy(base)
    c1["id"] = "CTRL-entry-moved-by-2"
    c1["rows"][0][3] += 2
    c2 = copy.deepcopy(base)
    c2["id"] = "CTRL-support-swapped"
    c2["S"] = [k for k in c2["S"] if k != 4] + [2]
    c3 = copy.deepcopy(base)
    c3["id"] = "CTRL-row-sum-falsified"
    c3["row_sums"] = [6, 2, 2, 1]
    c4 = copy.deepcopy(base)
    c4["id"] = "CTRL-w-changed-to-30"
    c4["w"] = 30
    for c in (c1, c2, c3, c4):
        ok, bad = accept(c)
        check("%-24s rejected" % c["id"], not ok,
              "first failing condition: %s" % (bad[0][0] if bad else "none"))
    # and a positive control on the same path: the unperturbed table passes
    ok, _ = accept(base)
    check("the unperturbed table A still passes the same path", ok)


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verbose", action="store_true",
                    help="print each table's M, F and Fourier budget")
    args = ap.parse_args(argv)

    print("=" * 78)
    print("cert 30 -- the general-branch sigma-layer at N = 2092, exact "
          "integers")
    print("=" * 78)

    log("[0] the artifact")
    got = file_sha256(DATA)
    check("data/general-branch-sigma-tables.json sha256 = %s"
          % DATA_SHA256[:16] + "...", got == DATA_SHA256, got)
    with open(DATA, encoding="ascii") as fh:
        data = json.load(fh)
    tables = data["tables"]
    check("the artifact holds six tables", len(tables) == 6,
          [e["id"] for e in tables])
    check("and declares N = %d" % N, data["N"] == N)

    clause_a(tables, args.verbose)
    clause_b(tables)
    clause_c(tables, data["transport"])
    clause_d(tables)
    clause_e(tables)

    print("\n" + "=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 30: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: the sigma-layer of every advertised general branch at")
    print("         N = 2092 is NON-EMPTY -- six explicit integral coset-sum")
    print("         tables, every one of whose M, F, aggregate PAF, parity,")
    print("         box, norms and Fourier budgets was re-derived here from")
    print("         (S, i, w) alone, by two implementations sharing no code.")
    print("         The support enumerations are COMPLETE: one Galois-stable")
    print("         size-11 support on Z16, four size-7 on Z12, so (7,12) is")
    print("         populated on every rational branch.  The (3,8) cell is")
    print("         REDUNDANT: all 56 + 12 size-3 supports are real twists of")
    print("         W \\ {1}, and table C transports onto the house (3,4)")
    print("         profile at w = 130.  LABEL: PROVEN-BY-CERTIFICATE for")
    print("         'these tables exist with these profiles, and the support")
    print("         enumerations are complete'.")
    print("         NOT claimed: anything about SEEDS -- a compressed table is")
    print("         a necessary shadow, never sufficient, and no seed")
    print("         quadruple is known at any of these cells; anything about")
    print("         BORDERS (cert 31); anything about H(2092); the number of")
    print("         tables in any sigma-layer.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
