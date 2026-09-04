#!/usr/bin/env python3
"""cert 29 -- Lemma F, and no odd prime acts on the 668 / 716 records.

  LEMMA F -- note/NOTE-B.md S1.11.  Let H be Hadamard of ORDER N (the
  note fixes n = |G| and N = 4v; this certificate uses N throughout)
  and let (P, Q) be an automorphism -- signed permutation matrices with
  P H Q^T = H -- of odd prime order p.  Then

    (i)   the signs are removable by diagonal conjugation;
    (ii)  #fixed rows = #fixed columns, since P = H Q H^T / N;
    (iii) if (P, Q) is nontrivial then f := #fixed rows <= N/p, and if
          f <= m := (N-f)/p then f <= N/(p+1);
    (iv)  the m orbits give p x p circulant blocks, the fixed rows are
          constant on column orbits, and A A^T + p B B^T = N I_f.

  AT N = 2092, p = 523.  523 | 2092, so f = 0 (mod 523), and (iii)
  gives f <= 4: hence f = 0 and m = 4 -- a 16-BLOCK CIRCULANT ARRAY.
  That is all it is: "the all-type-1 Williamson array" is a 4-seed
  SUBFAMILY, and Lemma F forces only the block shape.  The
  Goethals-Seidel array is not of that shape (its off-diagonal blocks
  are back-circulants) and its translation automorphisms would need
  2a = 0, which no odd |G| has.  22 odd primes are admissible at 2092
  and 293 are excluded, 131 and 349 among them.

  THE RECORDS.  For row i write pi(i) for the multiset of |T4(i,j,k,l)|
  over 3-subsets of the other rows; a signed row/column permutation
  preserves pi, so the pi-classes are an invariant partition of the
  rows.  If an automorphism of odd prime order p exists, every class is
  a union of p-cycles and fixed rows, so it contributes at least
  |C| mod p fixed rows, and (iii) caps the total at N/p:

      sum_C (|C| mod p) > N/p    ==>    NO automorphism of order p.

  On H(668), H'(668), H''(668) the per-row |T4| classes are the four
  border rows as singletons and 332 tau-pairs {g, g + 83}; on H(716),
  four singletons and 356 tau-pairs {g, g + 89}.  Every class has size
  1 or 2, so the forced fixed-row count is N itself, which exceeds N/p
  for every odd p.  HENCE NO ODD PRIME ACTS ON ANY OF THE FOUR, AND
  Aut IS A 2-GROUP.  Every automorphism fixes the border rows up to
  sign and acts within tau-pairs.  (|T3| is NOT invariant under column
  signs and is not used anywhere here.)

  WHERE THE CLASS PARTITION WAS COMPUTED, AND BY WHOM.  The IMPLICATION
  above is PROVEN and is checked here from the pinned class sizes.  Its
  INPUT -- the per-row |T4| partition at the record orders -- is a
  heavy computation, and it has been made TWICE, both times in the
  SOURCE LABORATORY (Hadamard-2060):

    * first by a numpy pair-histogram finder (skeptic-pass/
      c9b_lemma_f_records.py, 2026-09-02), which produced the class
      counts banked in pins.json;
    * then by THIS SCRIPT's route B under --full -- standard library,
      exact integers, no numpy -- run on all four records as detached
      single-core processes on 2026-09-03/04: ALL CHECKS PASS at
      9 044.0 s for H(668), 9 062.2 s for H'(668), 9 049.1 s for
      H''(668) and 11 273.6 s for H(716).

  That second run is the independent implementation D-008 asks for, and
  it is why the records claim is PROVEN-BY-CERTIFICATE rather than
  COMPUTATIONAL-EVIDENCE.  IT HAS NOT BEEN RUN IN THIS REPOSITORY: the
  --full flag below is OFFERED AND PRICED, at about 2.5 h per 668
  matrix and 3.1 h at 716 on one core -- 10.7 core-hours for the four.
  Say "banked exact computation AUDITED" of a default run here.

WHAT THIS SCRIPT DOES  (standard library only, exact integers only)

  [A] LEMMA F ON SMALL MATRICES, and the arithmetic at 2092.  (i)-(iv)
      exercised on Sylvester H(8) (automorphisms of order 7 and 3) and
      Paley H(12) (orders 11, 5 and a signed Mobius map of order 3),
      each time checking the fixed-row bound, the circulant orbit
      blocks and A A^T + p B B^T = N I_f.  Then: 22 admissible / 293
      excluded odd primes at 2092, with 131 and 349 among the excluded;
      523 => f = 0; and 4 | b for a 4 x b partial Hadamard matrix
      (b = 2..12), which is what a constant-strip multi-block border
      with four fixed rows would need.

  [B] TWO |T4| IMPLEMENTATIONS ON CONTROLS (D-008).  Route A (direct
      C(N-1,3) per row) and route B (the pair-histogram) agree on
      Sylvester H(8), Paley H(12) and the GS(28) array on Williamson(7)
      seeds.  The two Aut-transitive controls return ONE class, as they
      must; GS(28) returns more than one, which is the case the
      implication is about.

  [C] THE FOUR RECORDS, REASSEMBLED HERE, AND THE IMPLICATION.  Each of
      H(668), H'(668), H''(668), H(716) is rebuilt in this run from
      data/payload-records.json and data/twisted-i2-records.json (both
      file-pinned), the twisted record's seeds re-derived as the
      psi-twist of the decoded record's and required to agree, and each
      matrix's canonical SHA-256 matched against the pin.  tau is
      exhibited as an automorphism of each, with the border sign the
      structure predicts (+1 on the (1,1) records, -1 on H').  Then
      the implication is applied to the PINNED class sizes: those sizes
      kill every odd prime.

  [D] --full [ --matrix NAME ].  Recomputes route B on the named record
      (or all four) IN THIS REPOSITORY and requires the pinned class
      sizes, the four border singletons, the tau-pairing and an empty
      surviving-prime list.  PRICE: about 2.5 h per 668 matrix and
      3.1 h at 716 on one core, 10.7 core-hours for the four.  It has
      NOT been run here; the source laboratory ran it on all four (the
      walls above).

  PORTED from the source laboratory (Hadamard-2060 certs/0028-lemma-f-
  records, D-068), unchanged in its mathematics; the lab-only inputs
  were rewired onto this repository's data/ -- see NOTES.md.

Usage:
  python certs/29-lemma-f-records/run.py
  python certs/29-lemma-f-records/run.py --full --matrix "H(668)"
"""

import argparse
import hashlib
import itertools
import json
import os
import sys
import time
from collections import Counter
from math import isqrt

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DECODE = os.path.join(ROOT, "data", "payload-records.json")
TWIST = os.path.join(ROOT, "data", "twisted-i2-records.json")
PINS = os.path.join(HERE, "pins.json")

# SHA-256 of the two banked data FILES, so a silently edited bank is a hard
# error rather than a different theorem.  Both are shared with other certs
# (payload-records.json with certs 01, 06, 08, 11, 13, 14, 15, 17, 19, 20,
# 23; twisted-i2-records.json with cert 02 and the same separation certs)
# and are re-pinned here at exactly the values those certificates carry.
FILE_PINS = {
    "data/payload-records.json": "9afbf60efed82e97ffe229e0059f06dfa7849a2ac2434c52a6ada09559738afb",
    "data/twisted-i2-records.json": "aafa83e070d2dc59da80aec1bcb6457b6bfcb8d7bcc758f20d9344047bfcb079",
}


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()

_T0 = time.time()
FAIL = []
NCHECK = [0]


def log(msg):
    print("\n[%6.1fs] %s" % (time.time() - _T0, msg), flush=True)


def check(label, cond, extra=""):
    NCHECK[0] += 1
    ok = bool(cond)
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", label,
                           ("  -- " + str(extra)) if extra else ""), flush=True)
    if not ok:
        FAIL.append(label)
    return ok


def pm(s):
    return [1 if ch == "+" else -1 for ch in s]


def to_pm(row):
    return "".join("+" if v == 1 else "-" for v in row)


def canonical_sha256(rows):
    text = "\n".join(to_pm(r) for r in rows) + "\n"
    return hashlib.sha256(text.encode("ascii")).hexdigest()


# ---------------------------------------------------------------- groups / assembler (self-contained)

class Grp:
    def __init__(self, factors):
        self.f = list(factors)
        self.n = 1
        for a in self.f:
            self.n *= a
        self.elts = list(itertools.product(*[range(a) for a in self.f]))
        self.index = {e: k for k, e in enumerate(self.elts)}

    def idx(self, t):
        return self.index[tuple(a % f for a, f in zip(t, self.f))]

    def add(self, p, q):
        return self.idx(tuple(a + b for a, b in zip(self.elts[p], self.elts[q])))

    def sub(self, p, q):
        return self.idx(tuple(a - b for a, b in zip(self.elts[p], self.elts[q])))

    def kappa(self, divisors):
        kap = []
        for t in self.elts:
            c = 0
            for a, d in zip(t, divisors):
                c = c * d + (a % d)
            kap.append(c)
        i = 1
        for d in divisors:
            i *= d
        return kap, i


GS_TAB = {(0, 1): (1, "R", 1), (0, 2): (2, "R", 1), (0, 3): (3, "R", 1),
          (1, 0): (1, "R", -1), (1, 2): (3, "TR", 1), (1, 3): (2, "TR", -1),
          (2, 0): (2, "R", -1), (2, 1): (3, "TR", -1), (2, 3): (1, "TR", 1),
          (3, 0): (3, "R", -1), (3, 1): (2, "TR", 1), (3, 2): (1, "TR", -1)}


def gs_core_rows(G, x, rho, switch=False):
    n = G.n
    sub = [[G.sub(h, g) for h in range(n)] for g in range(n)]
    rmh = [[G.sub(G.sub(rho, g), h) for h in range(n)] for g in range(n)]
    gph = [[G.sub(G.add(g, h), rho) for h in range(n)] for g in range(n)]
    rows = []
    for I in range(4):
        for g in range(n):
            row = []
            for J in range(4):
                if I == J:
                    row.extend(x[0][k] for k in sub[g])
                else:
                    q, form, sg = GS_TAB[(I, J)]
                    if switch:
                        sg = -sg
                    src = rmh[g] if form == "R" else gph[g]
                    xq = x[q]
                    row.extend(sg * xq[k] for k in src)
            rows.append(row)
    return rows


def bordered_rows(G, x, rho, s, divisors, E, P, Q, switch=False):
    n = G.n
    core = gs_core_rows(G, x, rho, switch)
    if s == 0:
        return core
    kap, i = G.kappa(divisors)
    rows = []
    for r in range(4 * s):
        row = list(E[r])
        for J in range(4):
            row.extend(P[r][i * J + kap[h]] for h in range(n))
        rows.append(row)
    for I in range(4):
        for g in range(n):
            rows.append([Q[i * I + kap[g]][c] for c in range(4 * s)] + core[I * n + g])
    return rows


def record_to_parts(rec):
    G = Grp(rec["group"])
    x = [pm(t) for t in rec["seeds"]]
    rho = G.idx(tuple(rec["r_shift"]))
    s = int(rec["s"])
    div = list(rec["coset_divisors"])
    E = [pm(t) for t in rec["corner"]]
    P = [pm(t) for t in rec["row_table"]]
    colT = [pm(t) for t in rec["col_table"]]
    _, i = G.kappa(div)
    Q = [[colT[r][k] for r in range(4 * s)] for k in range(4 * i)]
    return G, x, rho, s, div, E, P, Q


def apply_aut(H, rperm, rsign, cperm, csign):
    N = len(H)
    return [[rsign[i] * csign[j] * H[rperm[i]][cperm[j]] for j in range(N)] for i in range(N)]


# ---------------------------------------------------------------- |T4| : two implementations

def pack(rows):
    return [sum(1 << c for c, v in enumerate(r) if v == -1) for r in rows]


def row_t4_direct(rows):
    """Route A: for each row i, C(N-1,3) 4-tuple popcounts."""
    N = len(rows)
    ints = pack(rows)
    out = []
    for i in range(N):
        others = [ints[j] for j in range(N) if j != i]
        cnt = Counter()
        hi = ints[i]
        m = len(others)
        for a in range(m):
            xa = hi ^ others[a]
            for b in range(a + 1, m):
                xab = xa ^ others[b]
                for c in range(b + 1, m):
                    cnt[abs(N - 2 * (xab ^ others[c]).bit_count())] += 1
        out.append(tuple(cnt.get(s, 0) for s in range(N + 1)))
    return out


def row_t4_pairs(rows, progress=False):
    """Route B: pair-histogram (pass II's method), stdlib bit-packed ints."""
    N = len(rows)
    ints = pack(rows)
    pairs = [(k, l, ints[k] ^ ints[l]) for k in range(N) for l in range(k + 1, N)]
    tot = [[0] * (N + 1) for _ in range(N)]
    t0 = time.time()
    for i in range(N):
        hi = ints[i]
        for j in range(i + 1, N):
            hij = hi ^ ints[j]
            hist = [0] * (N + 1)
            for _k, _l, wp in pairs:
                hist[abs(N - 2 * (hij ^ wp).bit_count())] += 1
            hist[0] -= 2 * (N - 2)
            hist[N] -= 1
            for s in range(N + 1):
                tot[i][s] += hist[s]
                tot[j][s] += hist[s]
        if progress and (i % 50 == 0 or i == N - 1):
            print("      row %d / %d  (%.0fs)" % (i, N, time.time() - t0), flush=True)
    for i in range(N):
        if any(v % 3 for v in tot[i]):
            raise AssertionError("pair-histogram not divisible by 3 at row %d" % i)
        tot[i] = [v // 3 for v in tot[i]]
    return [tuple(tot[i]) for i in range(N)]


def classes_of(prof):
    d = {}
    for i, p in enumerate(prof):
        d.setdefault(p, []).append(i)
    return d


def fmin_killed(classes, N):
    odd = [p for p in range(3, N + 1) if all(p % d for d in range(2, isqrt(p) + 1))]
    killed, alive = [], []
    for p in odd:
        fmin = sum(len(v) % p for v in classes.values())
        if fmin > N / p:
            killed.append(p)
        else:
            alive.append((p, fmin))
    return killed, alive


def implication_kills_all(class_sizes, N):
    """class_sizes is a dict size -> count.  fmin = sum count*(size mod p)."""
    odd = [p for p in range(3, N + 1) if all(p % d for d in range(2, isqrt(p) + 1))]
    for p in odd:
        fmin = sum(sz % p * cnt for sz, cnt in class_sizes.items())
        if fmin <= N / p:
            return False
    return True


# ---------------------------------------------------------------- small Hadamard matrices

def sylvester(k):
    n = 1 << k
    return [[(-1) ** bin(a & b).count("1") for b in range(n)] for a in range(n)]


def paley12():
    q = 11
    chi = {x: 0 for x in range(q)}
    for x in range(1, q):
        chi[x * x % q] = 1
    for x in range(1, q):
        if chi[x] == 0:
            chi[x] = -1
    n = q + 1
    S = [[0] * n for _ in range(n)]
    for j in range(q):
        S[0][j + 1] = 1
        S[j + 1][0] = -1
    for a in range(q):
        for b in range(q):
            S[a + 1][b + 1] = chi[(b - a) % q]
    return [[S[a][b] + (1 if a == b else 0) for b in range(n)] for a in range(n)]


def check_lemma_f(H, rperm, cperm, p, label):
    N = len(H)
    assert apply_aut(H, rperm, [1] * N, cperm, [1] * N) == H
    fr = [i for i in range(N) if rperm[i] == i]
    fc = [j for j in range(N) if cperm[j] == j]
    f = len(fr)
    ok = (len(fr) == len(fc)) and (f <= N / p) and (N - f) % p == 0

    def orbits(perm):
        seen, out = set(), []
        for a in range(N):
            if a in seen:
                continue
            o = [a]
            b = perm[a]
            while b != a:
                o.append(b)
                b = perm[b]
            seen |= set(o)
            if len(o) > 1:
                out.append(o)
        return out
    ro, co = orbits(rperm), orbits(cperm)
    ok &= all(len(o) == p for o in ro + co)
    A = [[H[i][j] for j in fc] for i in fr]
    B = [[H[i][o[0]] for o in co] for i in fr]
    ok &= all(len(set(H[i][j] for j in o)) == 1 for i in fr for o in co)
    for a in range(f):
        for b in range(f):
            v = sum(A[a][k] * A[b][k] for k in range(f)) + p * sum(B[a][k] * B[b][k] for k in range(len(co)))
            ok &= (v == (N if a == b else 0))
    for o in ro:
        for o2 in co:
            vals = {}
            for k in range(p):
                for l in range(p):
                    vals.setdefault((l - k) % p, set()).add(H[o[k]][o2[l]])
            ok &= all(len(v) == 1 for v in vals.values())
    check("%s: p = %d, fixed %d = %d <= n/p; circulant orbits; A A^T + p B B^T = n I"
          % (label, p, len(fr), len(fc)), ok)
    return ok


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for p in range(2, isqrt(n) + 1):
        if s[p]:
            s[p * p::p] = bytearray(len(s[p * p::p]))
    return [p for p in range(n + 1) if s[p]]


def admissible_f(n, p):
    return [f for f in range(0, n // p + 1) if (n - f) % p == 0 and (f > 0 or n % p == 0)]


def load_matrix(name):
    with open(DECODE, encoding="ascii") as fh:
        recs = {int(r["order"]): r for r in json.load(fh)["orders"]}
    if name == "H(668)":
        N, which = 668, "H"
    elif name == "H'(668)":
        N, which = 668, "H'"
    elif name == "H''(668)":
        N, which = 668, "H''"
    elif name == "H(716)":
        N, which = 716, "H"
    else:
        raise SystemExit("unknown matrix " + name)
    G, x, rho, s, div, E, P, Q = record_to_parts(recs[N])
    n = G.n
    if which == "H":
        return bordered_rows(G, x, rho, s, div, E, P, Q)
    if which == "H''":
        return bordered_rows(G, x, rho, s, div, E, P, Q, switch=True)
    # H'(668): the psi-twisted i = 2 instance.  This repository banks it
    # whole (data/twisted-i2-records.json, cert 02's record), seeds and
    # border together, where the source laboratory kept the twisted seeds
    # in one file and the border in another.  The seeds are re-derived
    # here from the decoded record and required to agree, so the twist is
    # exercised rather than trusted.
    with open(TWIST, encoding="ascii") as fh:
        tw = {int(r["order"]): r for r in json.load(fh)["orders"]}
    t = tw[N]
    Gt, tx, trho, ts, tdiv, tE, tP, tQ = record_to_parts(t)
    psi = [(-1) ** g for g in range(n)]
    px = [[psi[g] * x[q][g] for g in range(n)] for q in range(4)]
    if tx != px:
        raise SystemExit("twisted record's seeds are not the psi-twist of the decoded record's")
    return bordered_rows(Gt, tx, trho, ts, tdiv, tE, tP, tQ)


# ---------------------------------------------------------------- stages

def stage_lemma_f():
    log("[A] Lemma F (i)-(iv) on H(8) / H(12)")
    H8 = sylvester(3)

    def matvec(M, v):
        return sum(((sum(M[r][c] * ((v >> c) & 1) for c in range(3)) % 2) << r) for r in range(3))

    def inv_t(M):
        for cand in itertools.product((0, 1), repeat=9):
            Cm = [list(cand[3 * r:3 * r + 3]) for r in range(3)]
            prod = [[sum(M[r][k] * Cm[c][k] for k in range(3)) % 2 for c in range(3)] for r in range(3)]
            if prod == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
                return Cm
    M7 = [[0, 0, 1], [1, 0, 1], [0, 1, 0]]
    M3 = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    for M, p in ((M7, 7), (M3, 3)):
        Mi = inv_t(M)
        rperm = [matvec(M, a) for a in range(8)]
        cperm = [matvec(Mi, b) for b in range(8)]
        check_lemma_f(H8, rperm, cperm, p, "Sylvester H(8), order %d" % p)

    H12 = paley12()
    q, N = 11, 12
    pts = ["inf"] + list(range(q))

    def perm_from(fun):
        idx = {p: i for i, p in enumerate(pts)}
        return [idx[fun(p)] for p in pts]
    tr = perm_from(lambda x: "inf" if x == "inf" else (x + 1) % q)
    check("Paley H(12): translation is an unsigned automorphism",
          apply_aut(H12, tr, [1] * N, tr, [1] * N) == H12)
    check_lemma_f(H12, tr, tr, 11, "Paley H(12), translation")
    m4 = perm_from(lambda x: "inf" if x == "inf" else (4 * x) % q)
    check("Paley H(12): x -> 4x is an unsigned automorphism",
          apply_aut(H12, m4, [1] * N, m4, [1] * N) == H12)
    check_lemma_f(H12, m4, m4, 5, "Paley H(12), multiplier 4")

    def mobius(a, b, c, d):
        def f(x):
            if x == "inf":
                return "inf" if c == 0 else (a * pow(c, -1, q)) % q
            den = (c * x + d) % q
            if den == 0:
                return "inf"
            return ((a * x + b) * pow(den, -1, q)) % q
        return f
    found = None
    for a, b, c, d in itertools.product(range(q), repeat=4):
        if (a * d - b * c) % q == 0:
            continue
        f = mobius(a, b, c, d)
        perm = perm_from(f)
        x = list(range(N))
        x = [perm[v] for v in x]
        x = [perm[v] for v in x]
        x = [perm[v] for v in x]
        if x != list(range(N)) or perm == list(range(N)):
            continue
        Hp = [[H12[perm[i]][perm[j]] for j in range(N)] for i in range(N)]
        cs = [Hp[0][j] * H12[0][j] for j in range(N)]
        rs = [Hp[i][0] * cs[0] * H12[i][0] for i in range(N)]
        if all(rs[i] * cs[j] * Hp[i][j] == H12[i][j] for i in range(N) for j in range(N)):
            found = (perm, rs, cs)
            break
    check("Paley H(12): a Mobius map of order 3 is a signed automorphism", found is not None)
    if found:
        perm, rs, cs = found
        fr = sum(1 for i in range(N) if perm[i] == i)
        trP = sum(rs[i] for i in range(N) if perm[i] == i)
        trQ = sum(cs[j] for j in range(N) if perm[j] == j)
        check("order-3: %d fixed points <= 4; tr P = tr Q" % fr, fr <= 4 and trP == trQ)

    log("[A] admissible odd primes; 523 => 16-block circulant array")
    for n in (2092, 668, 716):
        odd = [p for p in primes_upto(n) if p > 2]
        adm = [(p, admissible_f(n, p)) for p in odd if admissible_f(n, p)]
        print("   n = %d: %d odd primes, %d admissible" % (n, len(odd), len(adm)))
        if n == 2092:
            check("2092: 22 admissible / 293 excluded; 131, 349 excluded; 523 forces f = 0",
                  len(adm) == 22 and len(odd) - len(adm) == 293
                  and 131 not in dict(adm) and 349 not in dict(adm)
                  and dict(adm)[523] == [0])

    log("[A] constant strips: 4 | b")

    def partial_had(k, b):
        rows = [tuple([1] * b)]
        vecs = [v for v in itertools.product((1, -1), repeat=b) if v[0] == 1]

        def rec(rows):
            if len(rows) == k:
                return True
            for v in vecs:
                if all(sum(a * c for a, c in zip(v, r)) == 0 for r in rows):
                    if rec(rows + [v]):
                        return True
            return False
        return rec(rows)
    tab = {b: partial_had(4, b) for b in range(2, 13, 2)}
    check("4 x b partial Hadamard exists iff 4 | b (b = 2..12 even)",
          all(v == (b % 4 == 0) for b, v in tab.items()), tab)


def stage_controls():
    log("[B] two |T4| implementations on controls")
    for name, H in (("Sylvester H(8)", sylvester(3)), ("Paley H(12)", paley12())):
        a = row_t4_direct(H)
        b = row_t4_pairs(H)
        check("%s: route A == route B" % name, a == b)
        ca, cb = classes_of(a), classes_of(b)
        check("%s: one class (transitive Aut)" % name, len(ca) == 1 == len(cb))

    def paf(x):
        n = len(x)
        return [sum(x[u] * x[(u + s) % n] for u in range(n)) for s in range(n)]
    sym = []
    for bits in itertools.product((1, -1), repeat=4):
        sym.append([bits[0]] + list(bits[1:]) + list(bits[1:][::-1]))
    W = next(q for q in itertools.product(sym, repeat=4)
             if all(sum(paf(s)[k] for s in q) == 0 for k in range(1, 7)))
    H28 = gs_core_rows(Grp([7]), [list(s) for s in W], 0)
    a = row_t4_direct(H28)
    b = row_t4_pairs(H28)
    check("GS(28) Williamson(7): route A == route B", a == b)
    check("GS(28): more than one |T4| class (not cocyclic by this test)",
          len(classes_of(a)) >= 2, len(classes_of(a)))


def stage_records_audit():
    log("[C] file pins, the four records reassembled here, and the implication")
    for relp, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, *relp.split("/")))
        check("%s: file SHA-256 == pin" % relp, got == want, got)
    with open(PINS, encoding="ascii") as fh:
        pins = json.load(fh)
    with open(DECODE, encoding="ascii") as fh:
        recs = {int(r["order"]): r for r in json.load(fh)["orders"]}
    with open(TWIST, encoding="ascii") as fh:
        tw = {int(r["order"]): r for r in json.load(fh)["orders"]}
    check("the twisted i = 2 record at 668 carries its own pinned digest, and this "
          "run's assembly matches it",
          canonical_sha256(load_matrix("H'(668)")) == tw[668]["pinned_sha256"])

    for name, pin in pins.items():
        rows = load_matrix(name)
        sha = canonical_sha256(rows)
        check("%s: assembled digest == pin" % name, sha == pin["digest"], sha)
        check("%s: order %d" % (name, pin["N"]), len(rows) == pin["N"])
        N = pin["N"]
        n = (N - 4) // 4
        u = n // 2
        perm = list(range(4)) + [4 + I * n + (g + u) % n for I in range(4) for g in range(n)]
        if name.startswith("H'(") and "H''" not in name:
            sgn = [-1] * 4 + [1] * (4 * n)
        else:
            sgn = [1] * N
        check("%s: tau is an automorphism (predicted border sign)" % name,
              apply_aut(rows, perm, sgn, perm, sgn) == rows)
        sizes = {int(k): int(v) for k, v in pin["class_sizes"].items()}
        check("%s: pinned class sizes imply no odd prime acts" % name,
              implication_kills_all(sizes, N), sizes)
        check("%s: pinned partition is %d classes totalling %d rows"
              % (name, pin["n_classes"], N),
              (len(sizes) > 0
               and sum(k * v for k, v in sizes.items()) == N
               and sum(sizes.values()) == pin["n_classes"]
               and sizes.get(1, 0) == pin["border_singletons"]))

    # 523 => 16-block, not Williamson
    check("523 | 2092 and f <= 2092/523 = 4 with f == 0 (mod 523) forces f = 0, m = 4: "
          "a 16-block circulant array, not the Williamson array",
          2092 % 523 == 0 and admissible_f(2092, 523) == [0])


def stage_full(which):
    log("[D] --full: recompute |T4| pair-histogram on %s" % which)
    with open(PINS, encoding="ascii") as fh:
        pins = json.load(fh)
    names = list(pins) if which == "all" else [which]
    for name in names:
        pin = pins[name]
        rows = load_matrix(name)
        check("%s --full: digest" % name, canonical_sha256(rows) == pin["digest"])
        prof = row_t4_pairs(rows, progress=True)
        cls = classes_of(prof)
        sizes = Counter(len(v) for v in cls.values())
        want = {int(k): int(v) for k, v in pin["class_sizes"].items()}
        check("%s --full: class sizes %s" % (name, dict(sizes)), dict(sizes) == want)
        check("%s --full: %d classes" % (name, pin["n_classes"]), len(cls) == pin["n_classes"])
        border = [prof[i] for i in range(4)]
        check("%s --full: 4 border rows are distinct singletons" % name,
              len(set(border)) == 4 and all(len(cls[p]) == 1 for p in border))
        N = pin["N"]
        n = (N - 4) // 4
        u = n // 2
        tau_ok = all(prof[4 + I * n + g] == prof[4 + I * n + (g + u) % n]
                     for I in range(4) for g in range(n))
        check("%s --full: tau-pairs share profiles" % name, tau_ok)
        killed, alive = fmin_killed(cls, N)
        check("%s --full: no odd prime survives the fixed-row bound" % name,
              not alive, alive)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="recompute |T4| classes on the banked records (slow)")
    ap.add_argument("--matrix", default="all",
                    help="with --full: H(668), H'(668), H''(668), H(716), or all")
    args = ap.parse_args()
    print("=" * 78)
    print("cert 29 -- Lemma F, and no odd prime acts on the 668/716 records")
    print("=" * 78)
    stage_lemma_f()
    stage_controls()
    stage_records_audit()
    replayed = ""
    if args.full:
        stage_full(args.matrix)
        replayed = args.matrix
    log("SUMMARY")
    print("=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 29: FAIL (%d)" % len(FAIL))
        for f in FAIL:
            print("    FAILED:", f)
        return 1
    print("ALL CHECKS PASS")
    print("VERDICT: NO AUTOMORPHISM OF ODD PRIME ORDER ACTS ON H(668),")
    print("         H'(668), H''(668) OR H(716) -- Aut is a 2-group at")
    print("         each of the four, every automorphism fixing the four")
    print("         border rows up to sign and acting within tau-pairs.")
    print("         LABEL: PROVEN (Lemma F (i)-(iv), the prime lists, the")
    print("         16-block shape at 523, and the implication from the")
    print("         class sizes -- note/NOTE-B.md S1.11)")
    print("         + PROVEN-BY-CERTIFICATE (the class partition itself,")
    print("         computed twice in the source laboratory: a numpy")
    print("         pair-histogram finder, and this script's stdlib route B")
    print("         under --full on all four records, 2026-09-03/04).")
    if replayed:
        print("         PARTITION: RECOMPUTED in this run (%s)." % replayed)
    else:
        print("         PARTITION: banked exact computation AUDITED in this")
        print("         run, not recomputed.  --full is the in-repo replay;")
        print("         it is priced above and HAS NOT BEEN RUN HERE.")
    print("         NOT CLAIMED: the 2-PART of Aut at 668 or 716 -- ")
    print("         Aut = <-I, tau> is consistent with the row classes and")
    print("         is NOT proved; anything about H'(716) or H''(716),")
    print("         which were never run; anything about Aut of a")
    print("         hypothetical H(2092).")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
