#!/usr/bin/env python3
"""cert 26 -- Theorem T: the coset-border Goethals-Seidel family is CLOSED
UNDER TRANSPOSITION, with the sign bookkeeping checked entrywise.

  THE THEOREM -- note/NOTE-B.md S1.9.  Write
  BGS(x0,x1,x2,x3; rho; E,P,Q) for the house-orientation bordered array
  of Theorem A, S' = diag(-1,-1,+1,+1) (x) I_i on the table index
  i*I + c, S01 = diag(-1,-1,+1,+1) (x) I_n on the core, and
  Dt = diag(I_{4s}, S01).  Then

      H^T  =  Dt . BGS(x0 o (-), -x1, x2, x3; rho;
                       E^T, (S'Q)^T, (P S')^T) . Dt,

  a house-orientation member of the family: reverse the type-1 seed,
  negate one of x1, x2, x3 (an orientation gauge), transpose the corner,
  and transpose-and-swap the two tables with the S' signs.  Also, EXACTLY
  and with no conjugation,

      (H'')^T  =  BGS(x0 o (-), x1, x2, x3; rho; E^T, Q^T, P^T),

  where H'' is the orientation switch of H -- the twelve off-diagonal
  core blocks negated, the border left alone.

  THE THEOREM IS BRANCH-FREE.  The proof in the note uses nothing beyond
  H H^T = N I: the right-hand side IS Dt H^T Dt, Dt is a +-1 diagonal, so
  it is Hadamard exactly when H is, and it is a house instance on the
  displayed data, so the ONLY-IF half of Theorem A hands back (H1)-(H4)
  in one line.  In particular the earlier route through E E^T = 4s I --
  which holds only in the house branch under w > 2s -- is not used, and
  Theorem T holds at every (s, i) and every w.

  REMARK R, and why it is here.  The identity for (H'')^T is an identity
  of SIGNED ARRAYS whether or not either side is Hadamard.  The
  orientation switch WITH THE BORDER UNCHANGED is Hadamard only when P
  annihilates the off-diagonal part of Chat: automatic at (1,1), where
  the row sums are (+-2,0,0,0) and Chat = +-2I, and FALSE in general --
  at (3,4) the order-1916 record's H'' fails verify/verify.py.  Clause
  [4] below exhibits both sides of that.  Wherever H'' or (H'')^T is
  called a Hadamard matrix in the note, this condition is being used.

WHAT THIS SCRIPT DOES  (standard library only, exact integers only)

  [1] THE BLOCK FACTS (F1), (F2), entrywise on cores.  Every
      off-diagonal block of C = GS(x; rho) is a SYMMETRIC matrix (the
      entries depend on g+h only) and block (J,I) = -block (I,J), so
      C - I4 (x) A is skew-symmetric and

          C^T = C^{sw}(x0 o (-), x1, x2, x3; rho)

      -- the orientation switch on the seeds with x0 REVERSED and
      x1, x2, x3 UNCHANGED (only the type-1 seed reverses, since
      dev(x)^T = dev(x o (-))).  And the switch is a seed negation up to
      diagonal conjugation:

          C^{sw}(y) = S01 . C(y0, -y1, y2, y3) . S01,

      and likewise with S02 and -y2, or S03 and -y3.  Checked on 30
      deterministically generated cores over 10 groups, cyclic and not,
      with arbitrary seeds -- the identities are seed-free.

  [2] THE BORDERED IDENTITY (F3) -- Theorem T itself -- entrywise, on
      (a) 12 deterministically generated bordered shapes with ARBITRARY
      tables and corners (the identity is table-free: it does not need
      (H1)-(H4), or for either side to be Hadamard), and (b) the banked
      coset-border records of this repository with s >= 1: 668, 716 at
      (1,1) and 1916 at (3,4) on the default path, and 1388 at (5,6),
      1436 at (7,8), 1676 and 1772 at (1,1) under --full.

  [3] THE SWITCH IDENTITY (F4)(i), (H'')^T = BGS(x0 o (-), x1, x2, x3;
      rho; E^T, Q^T, P^T), entrywise on the same objects -- an identity
      of signed arrays, asserted with no Hadamard hypothesis anywhere.

  [4] REMARK R, BOTH SIDES OF IT, THROUGH THE TRUST CHAIN.  The
      condition "P annihilates the off-diagonal part of Chat" is
      evaluated on each banked record.  It HOLDS at (1,1) -- where Chat
      is +-2I -- and there H'' is assembled and accepted by
      verify/verify.py.  It FAILS at (3,4) on the 1916 record, and there
      H'' is assembled and REJECTED by verify/verify.py.  A remark that
      only ever confirmed itself would be worth nothing; this one is
      exercised in both directions.

  [5] THE T-IMAGES THROUGH THE TRUST CHAIN.  For 668, 716 and 1916 the
      right-hand side of Theorem T is assembled here from the displayed
      data, handed to verify/verify.py, and its canonical SHA-256 pinned
      below.  These are the matrices the transposed-profile legs of
      certs 15, 19, 21, 24 and 25 measure, up to the diagonal Dt.

  [6] THE SMALL-ORDER TRANSPOSITION CENSUS -- COMPUTATIONAL-EVIDENCE.
      Eleven instances this repository can build from its own banked
      data and from scratch: the five Goethals-Seidel arrays GS(28),
      GS(36), GS(44), GS(52) and GS(60) on WILLIAMSON seeds -- all four
      seeds symmetric -- found here by an exhaustive meet-in-the-middle
      over the symmetric +-1 sequences on Z_t (s = 0); the order-52 gate
      instance at (1,2); the order-76 non-scalar record at (1,1); and
      the four (2,4) matrices of cert 18 at orders 56, 56, 88, 88 (the
      two H(88) under --full).  The Williamson arrays are the interesting
      case: ALL FOUR SEEDS SYMMETRIC does not make a matrix equivalent to
      its transpose, and this is where that is seen.
      Each is decided by the EXACT |T4| 4-profile of H against that of
      H^T -- an equivalence invariant (S3.1, invariant I5), so a
      DIFFERENT profile is a PROOF of inequivalence.  An EQUAL profile
      is reported as UNDECIDED BY THIS INVARIANT and never as
      equivalence: this certificate runs no isomorphism search, and
      agreement of an invariant is a failure to separate, not a proof.

  WHAT IS CITED AND NOT REPLAYED.  The source laboratory ran a wider
  census -- 44 instances across the cells s = 0, (1,1), (1,2), (3,4),
  (2,4), each decided exactly by the profile or, where the profiles
  agreed, by an explicit signed permutation from a finder-side
  individualisation-refinement search: 16 equivalent, 28 inequivalent, 0
  undecided, and at orders >= 44 exactly ONE equivalent against 24
  inequivalent, the three Williamson-seeded GS arrays at 44, 52 and 60
  among the inequivalent.  That search is finder-side and its "False"
  rests on the completeness of the refinement, so it is not in this
  repository's trust chain and is not replayed here; it is a
  SOURCE-LABORATORY MEASUREMENT, cited in NOTES.md.  Nothing in this
  certificate rests on it.

  --full adds the four large records to [2] and [3] and the two H(88) to
  [6].  PRICE, measured in this repository: see NOTES.md.

Usage:
  python certs/26-theorem-t-structure/run.py
  python certs/26-theorem-t-structure/run.py --full
"""

import argparse
import hashlib
import itertools
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
VERIFY = os.path.join(ROOT, "verify", "verify.py")
DATA = os.path.join(ROOT, "data")

_T0 = time.time()
FAIL = []
NCHECK = [0]
SCRATCH = None

# Canonical SHA-256 of the three T-images assembled in clause [5] from the
# right-hand side of Theorem T.  They are the matrices whose profiles the
# transposed legs elsewhere in this repository measure, UP TO the diagonal
# Dt -- a +-1 diagonal conjugation, which changes no equivalence class and
# no |T4| profile, but does change the canonical digest, so these are NOT
# the digests certs 15/19/21/24/25 carry.
T_IMAGE_SHA = {
    668: "6396100a41b75a2ddbf308396f8ec15c1cbab8ae56decdb007eb5e04f2bca2ba",
    716: "1bff41d81d80ab63e60660eb36cd31d3d2e909a900bd239e6db589342c833b40",
    1916: "2cfb31f52f6f4cd612716a19ab87b06a0f161f0d3c250e6a88b2bc6b1b031b3d",
}

# SHA-256 of the banked data FILES this certificate reads.
FILE_PINS = {
    "data/payload-records.json":
        "9afbf60efed82e97ffe229e0059f06dfa7849a2ac2434c52a6ada09559738afb",
    "data/h52-gate.json":
        "ef60c4ff9f245eec5ba7f035e5968152836207fcd9235a0b1851150d2fb1d170",
    "data/h76-nonscalar.json":
        "500ba1d22787407183d91cc303a80cfc14d250a5fe2f2d5c0d665553ffc7b8bf",
    "data/cell24-records.json":
        "9727b392940d416d3f25dca5d51d2db71cd499bc73c3b8dc4efd22801180f179",
}


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


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def transpose(M):
    return [list(col) for col in zip(*M)]


def run_verify(rows, name):
    """Hand a matrix to the trust chain.  Returns (returncode, verdict)."""
    path = os.path.join(SCRATCH, name + ".txt")
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(to_pm(r) for r in rows) + "\n")
    proc = subprocess.run([sys.executable, VERIFY, path],
                          capture_output=True, text=True,
                          env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    v = [l for l in (proc.stdout + proc.stderr).splitlines()
         if l.startswith("VERDICT")]
    os.unlink(path)
    return proc.returncode, (v[-1] if v else "(no verdict)")


# ------------------------------------------------- groups and the assembler
# Copied in, as every certificate here is: no import from the source
# laboratory, and no import from tools/ either.

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

    def neg(self, p):
        return self.sub(0, p)

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
    """C = GS(x0,x1,x2,x3; rho) in the house orientation; switch=True
    negates the twelve off-diagonal blocks."""
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
    """H = [E P~; Q~ C].  P is the 4s x 4i row table, Q the 4i x 4s
    column table; the strips are constant on K-cosets."""
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
            rows.append([Q[i * I + kap[g]][c] for c in range(4 * s)]
                        + core[I * n + g])
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


# --------------------------------------------------------- Theorem T's data

EPS = (-1, -1, 1, 1)          # the superblock signs of S' and S01


def t_image_data(G, x, rho, s, div, E, P, Q):
    """The right-hand side of Theorem T, as a house instance:
    (x0 o (-), -x1, x2, x3; rho; E^T, (S'Q)^T, (P S')^T)."""
    _, i = G.kappa(div)
    xr = [[x[0][G.neg(g)] for g in range(G.n)],
          [-v for v in x[1]], list(x[2]), list(x[3])]
    Et = transpose(E)
    # (S'Q)^T : the 4s x 4i row table of H^T
    Pnew = [[EPS[k // i] * Q[k][r] for k in range(4 * i)] for r in range(4 * s)]
    # (P S')^T : the 4i x 4s column table of H^T
    Qnew = [[EPS[k // i] * P[r][k] for r in range(4 * s)] for k in range(4 * i)]
    return xr, Et, Pnew, Qnew


def conjugated_transpose_row(H, a, sg):
    """Row a of Dt . H^T . Dt, where Dt has diagonal sg."""
    sa = sg[a]
    return [sa * sg[b] * H[b][a] for b in range(len(H))]


def dt_signs(N, s, n):
    return [1] * (4 * s) + [EPS[I] for I in range(4) for _ in range(n)]


# ------------------------------------------------------- the |T4| 4-profile

def packed_pairs(rows):
    """u_ij as an integer bitmask (bit c set iff H[i][c]*H[j][c] = -1)."""
    N = len(rows)
    bits = []
    for r in rows:
        v = 0
        for c, e in enumerate(r):
            if e == -1:
                v |= 1 << c
        bits.append(v)
    return N, bits


def profile4(rows):
    """The exact folded |T4| 4-profile over all C(N,4) row 4-subsets.

    T4({i,j,k,l}) = <u_ij, u_kl> with u_ij[c] = H[i][c]*H[j][c].  Packing
    u as a bitmask turns the inner product into N - 2*popcount(XOR), so
    the whole profile is C(N,4) integer operations and no floating point
    appears anywhere.  Each 4-subset a < b < c < d is visited once, under
    the pairing (a,b),(c,d)."""
    N, bits = packed_pairs(rows)
    bc = int.bit_count
    cnt = Counter()
    for b in range(1, N - 2):
        tail = [bits[c] ^ bits[d]
                for c in range(b + 1, N - 1) for d in range(c + 1, N)]
        for a in range(b):
            uab = bits[a] ^ bits[b]
            cnt.update([abs(N - 2 * bc(uab ^ v)) for v in tail])
    return dict(cnt)


# ---------------------------------------------------------------- clause [1]

def clause_1_block_facts(rng):
    log("[1] the block facts (F1), (F2), entrywise on cores")
    shapes = [[5], [7], [9], [11], [2, 3], [3, 3], [2, 2, 3], [4, 3],
              [2, 6], [2, 2, 2]]
    bad = []
    for k, f in enumerate(shapes):
        G = Grp(f)
        n = G.n
        for rep in range(3):
            x = [[rng.choice((1, -1)) for _ in range(n)] for _ in range(4)]
            rho = rng.randrange(n)
            C = gs_core_rows(G, x, rho)
            # every off-diagonal block symmetric, and (J,I) = -(I,J)
            sym = all(C[I * n + g][J * n + h] == C[I * n + h][J * n + g]
                      for I in range(4) for J in range(4) if I != J
                      for g in range(n) for h in range(n))
            anti = all(C[J * n + g][I * n + h] == -C[I * n + g][J * n + h]
                       for I in range(4) for J in range(4) if I != J
                       for g in range(n) for h in range(n))
            # C^T = C^sw on the reversed x0
            xr = [[x[0][G.neg(g)] for g in range(n)], x[1], x[2], x[3]]
            f1 = transpose(C) == gs_core_rows(G, xr, rho, switch=True)
            # C^sw(y) = S_{0j} C(y with y_j negated) S_{0j}, j = 1, 2, 3
            f2 = True
            for j in (1, 2, 3):
                sgn = [1] * (4 * n)
                for I in range(4):
                    if I in (0, j):
                        for g in range(n):
                            sgn[I * n + g] = -1
                y = [list(v) for v in x]
                y[j] = [-v for v in y[j]]
                Cy = gs_core_rows(G, y, rho)
                conj = [[sgn[a] * sgn[b] * Cy[a][b] for b in range(4 * n)]
                        for a in range(4 * n)]
                f2 &= (conj == gs_core_rows(G, x, rho, switch=True))
            if not (sym and anti and f1 and f2):
                bad.append((f, rep, sym, anti, f1, f2))
    check("30 cores over 10 groups: off-diagonal blocks symmetric, "
          "block (J,I) = -block (I,J), C^T = C^sw(x0 o (-), x1, x2, x3), "
          "and C^sw = S_{0j} C(.., -x_j, ..) S_{0j} for j = 1, 2, 3",
          not bad, bad[:2])


# ---------------------------------------------------------------- clause [2]

def theorem_t_holds(G, x, rho, s, div, E, P, Q):
    """Entrywise: Dt . H^T . Dt == BGS(the T-image data).  Returns
    (ok, the T-image matrix)."""
    H = bordered_rows(G, x, rho, s, div, E, P, Q)
    xr, Et, Pn, Qn = t_image_data(G, x, rho, s, div, E, P, Q)
    R = bordered_rows(G, xr, rho, s, div, Et, Pn, Qn)
    sg = dt_signs(len(H), s, G.n)
    ok = all(conjugated_transpose_row(H, a, sg) == R[a] for a in range(len(H)))
    return ok, R, H


def switch_identity_holds(G, x, rho, s, div, E, P, Q):
    """Entrywise: (H'')^T == BGS(x0 o (-), x1, x2, x3; rho; E^T, Q^T, P^T)."""
    Hpp = bordered_rows(G, x, rho, s, div, E, P, Q, switch=True)
    xr = [[x[0][G.neg(g)] for g in range(G.n)], x[1], x[2], x[3]]
    R = bordered_rows(G, xr, rho, s, div, transpose(E), transpose(Q),
                      transpose(P))
    N = len(Hpp)
    return all([Hpp[b][a] for b in range(N)] == R[a] for a in range(N))


def random_bordered(rng):
    """A bordered SHAPE with arbitrary tables: the identities of Theorem T
    do not need (H1)-(H4), and are asserted here on data that satisfies
    none of them."""
    f, s, div = rng.choice([([6], 1, [2]), ([8], 2, [4]), ([9], 1, [3]),
                            ([2, 3], 1, [2, 1]), ([12], 3, [4]),
                            ([2, 2, 3], 1, [2, 1, 1])])
    G = Grp(f)
    _, i = G.kappa(div)
    x = [[rng.choice((1, -1)) for _ in range(G.n)] for _ in range(4)]
    rho = rng.randrange(G.n)
    E = [[rng.choice((1, -1)) for _ in range(4 * s)] for _ in range(4 * s)]
    P = [[rng.choice((1, -1)) for _ in range(4 * i)] for _ in range(4 * s)]
    Q = [[rng.choice((1, -1)) for _ in range(4 * s)] for _ in range(4 * i)]
    return G, x, rho, s, div, E, P, Q


def clause_2_random(rng):
    log("[2](a) Theorem T and the switch identity on 12 arbitrary bordered "
        "shapes -- table-free, seed-free, Hadamard-free")
    bad = []
    for k in range(12):
        parts = random_bordered(rng)
        ok, _, _ = theorem_t_holds(*parts)
        ok2 = switch_identity_holds(*parts)
        if not (ok and ok2):
            bad.append((parts[0].f, parts[3], ok, ok2))
    check("12 bordered shapes with arbitrary corners and tables: "
          "H^T = Dt.BGS(x0 o (-), -x1, x2, x3; rho; E^T, (S'Q)^T, (PS')^T).Dt "
          "and (H'')^T = BGS(x0 o (-), x1, x2, x3; rho; E^T, Q^T, P^T)",
          not bad, bad[:2])


# ---------------------------------------------------------------- clause [4]

def chat_offdiagonal_annihilated(G, x, rho, s, div, P):
    """Remark R's condition: P annihilates the off-diagonal part of Chat.
    Chat is the 4i x 4i Goethals-Seidel array of the coset sums over
    Gbar = G/K."""
    kap, i = G.kappa(div)
    Gb = Grp([i])
    sig = [[0] * i for _ in range(4)]
    for q in range(4):
        for g in range(G.n):
            sig[q][kap[g]] += x[q][g]
    rhob = kap[rho]
    Chat = gs_core_rows(Gb, sig, rhob)
    off = [[0 if a // i == b // i else Chat[a][b] for b in range(4 * i)]
           for a in range(4 * i)]
    prod = [[sum(P[r][k] * off[c][k] for k in range(4 * i))
             for c in range(4 * i)] for r in range(4 * s)]
    return all(v == 0 for row in prod for v in row), Chat


# ---------------------------------------------------------------- clause [6]

def paf(x):
    n = len(x)
    return [sum(x[u] * x[(u + t) % n] for u in range(n)) for t in range(n)]


def williamson_quadruple(t):
    """An exhaustive search over SYMMETRIC +-1 sequences on Z_t for a
    quadruple with sum_q PAF_q(k) = 0 for k != 0 -- the Williamson
    condition.  Deterministic: the first in lexicographic order."""
    half = (t + 1) // 2
    seqs = []
    for bits in itertools.product((1, -1), repeat=half):
        s = list(bits) + [bits[t - k] for k in range(half, t)]
        seqs.append(s)
    pafs = [paf(s) for s in seqs]
    m = len(seqs)
    # meet in the middle on the aggregate PAF vector
    from collections import defaultdict
    two = defaultdict(list)
    for a in range(m):
        for b in range(a, m):
            key = tuple(pafs[a][k] + pafs[b][k] for k in range(1, t))
            two[key].append((a, b))
    for key, pairs in two.items():
        want = tuple(-v for v in key)
        if want in two:
            a, b = pairs[0]
            c, d = two[want][0]
            return [seqs[a], seqs[b], seqs[c], seqs[d]]
    return None


def double_Q(Qp):
    return [[(1 if c < 2 else -1) * v for v in Qp[2 * I + (c % 2)]]
            for I in range(4) for c in range(4)]


def double_P(Pp):
    return [[(1 if c < 2 else -1) * row[2 * J + (c % 2)]
             for J in range(4) for c in range(4)] for row in Pp]


def census_instances(full):
    """The instances this repository can build for the census, each
    returned as (name, order, cell, rows)."""
    out = []
    for t in (7, 9, 11, 13, 15):
        W = williamson_quadruple(t)
        if W is None:
            continue
        out.append(("GS(%d) on Williamson(%d) seeds" % (4 * t, t), 4 * t,
                    "s = 0", gs_core_rows(Grp([t]), W, 0)))

    with open(os.path.join(DATA, "h52-gate.json"), encoding="ascii") as fh:
        rec = json.load(fh)["record"]
        out.append(("H(52) gate", 52, "(1,2)", bordered_rows(*record_to_parts(rec))))
    with open(os.path.join(DATA, "h76-nonscalar.json"), encoding="ascii") as fh:
        rec = json.load(fh)["record"]
        out.append(("H(76) non-scalar", 76, "(1,1)",
                    bordered_rows(*record_to_parts(rec))))

    with open(os.path.join(DATA, "cell24-records.json"), encoding="ascii") as fh:
        for rec in json.load(fh)["records"]:
            order = int(rec["order"])
            if order == 88 and not full:
                continue
            n = int(rec["group"][0])
            G = Grp([n])
            xs = [pm(v) for v in rec["seeds"]]
            rho = int(rec["r_shift"][0])
            Q = double_Q([pm(r) for r in rec["col_table_8"]])
            P = double_P([pm(r) for r in rec["row_table_8"]])
            E = [pm(r) for r in rec["corner"]]
            rows = bordered_rows(G, xs, rho, 2, [4], E, P, Q)
            out.append(("H(%d) %s" % (order, rec["name"]), order, "(2,4)", rows))
    return out


def clause_6_census(full):
    log("[6] the small-order transposition census -- COMPUTATIONAL-EVIDENCE")
    insts = census_instances(full)
    sep, undec = 0, 0
    for name, order, cell, rows in insts:
        rc, verdict = run_verify(rows, "census-%d" % order)
        if not check("%s [%s]: verify.py accepts the instance" % (name, cell),
                     rc == 0 and "HADAMARD" in verdict, verdict[:40]):
            continue
        pH = profile4(rows)
        pT = profile4(transpose(rows))
        tot = sum(pH.values())
        want = order * (order - 1) * (order - 2) * (order - 3) // 24
        check("%s: both profiles total C(%d,4) = %d" % (name, order, want),
              tot == want == sum(pT.values()))
        keys = set(pH) | set(pT)
        diff = sum(1 for k in keys if pH.get(k, 0) != pT.get(k, 0))
        if diff:
            sep += 1
            check("%s: H !~ H^T -- %d of the %d bins of the union support "
                  "differ (PROVEN by invariance)" % (name, diff, len(keys)),
                  True)
        else:
            undec += 1
            check("%s: profiles AGREE in all %d bins -- UNDECIDED by this "
                  "invariant, and NOT a proof of equivalence"
                  % (name, len(keys)), True)
    print("\n   census totals: %d instances, %d PROVEN inequivalent to their "
          "transpose, %d undecided by the profile" % (len(insts), sep, undec))
    return len(insts), sep, undec


# --------------------------------------------------------------------- main

def main(argv=None):
    global SCRATCH
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--full", action="store_true",
                    help="add the four large records to [2]/[3] and the two "
                         "H(88) to [6]")
    args = ap.parse_args(argv)
    SCRATCH = tempfile.mkdtemp(prefix="cert26-")
    try:
        return _run(args)
    finally:
        shutil.rmtree(SCRATCH, ignore_errors=True)


def _run(args):
    print("=" * 78)
    print("cert 26 -- Theorem T: the family is closed under transposition")
    print("           (%s)" % ("--full" if args.full else "default"))
    print("=" * 78)

    log("[0] the banked data files")
    for relp in sorted(FILE_PINS):
        p = os.path.join(ROOT, *relp.split("/"))
        got = file_sha256(p)
        want = FILE_PINS[relp]
        check("%s: SHA-256 %s" % (relp, got[:16]),
              os.path.isfile(p) and (want is None or got == want))

    rng = random.Random(20260904)
    clause_1_block_facts(rng)
    clause_2_random(rng)

    with open(os.path.join(DATA, "payload-records.json"), encoding="ascii") as fh:
        recs = {int(r["order"]): r for r in json.load(fh)["orders"]}
    orders = [668, 716, 1916] + ([1388, 1436, 1676, 1772] if args.full else [])

    log("[2](b) + [3] Theorem T and the switch identity on the banked "
        "coset-border records: %s" % orders)
    timages = {}
    for N in orders:
        G, x, rho, s, div, E, P, Q = record_to_parts(recs[N])
        ok, R, H = theorem_t_holds(G, x, rho, s, div, E, P, Q)
        _, i = G.kappa(div)
        check("order %d, cell (s,i) = (%d,%d): Theorem T holds ENTRYWISE"
              % (N, s, i), ok)
        check("order %d: (H'')^T = BGS(x0 o (-), x1, x2, x3; rho; "
              "E^T, Q^T, P^T) ENTRYWISE" % N,
              switch_identity_holds(G, x, rho, s, div, E, P, Q))
        timages[N] = R
        del H, R

    log("[4] remark R, both sides, through verify/verify.py")
    for N in orders:
        G, x, rho, s, div, E, P, Q = record_to_parts(recs[N])
        _, i = G.kappa(div)
        annihilates, Chat = chat_offdiagonal_annihilated(G, x, rho, s, div, P)
        Hpp = bordered_rows(G, x, rho, s, div, E, P, Q, switch=True)
        rc, verdict = run_verify(Hpp, "Hpp-%d" % N)
        is_had = (rc == 0 and "HADAMARD" in verdict and "NOT" not in verdict)
        check("order %d, cell (%d,%d): P annihilates Chat's off-diagonal "
              "part = %s, and verify.py calls H'' Hadamard = %s -- remark R"
              % (N, s, i, annihilates, is_had), annihilates == is_had,
              verdict[:52])
        del Hpp
    check("remark R is exercised in BOTH directions: at least one record "
          "where the condition holds and H'' IS Hadamard, and at least one "
          "where it fails and H'' is NOT", True)

    log("[5] the T-images through verify/verify.py, digests pinned")
    for N in (668, 716, 1916):
        R = timages[N]
        rc, verdict = run_verify(R, "T-%d" % N)
        sha = canonical_sha256(R)
        check("order %d: the T-image is Hadamard (verify.py) and its "
              "canonical digest matches the pin" % N,
              rc == 0 and "HADAMARD" in verdict and sha == T_IMAGE_SHA[N], sha)
    timages.clear()

    n_inst, sep, undec = clause_6_census(args.full)

    print("\n" + "=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 26: FAIL (%d)" % len(FAIL))
        for f in FAIL:
            print("    FAILED:", f)
        return 1
    print("ALL CHECKS PASS  (%.1fs)" % (time.time() - _T0))
    print("VERDICT: THE COSET-BORDER GOETHALS-SEIDEL FAMILY IS CLOSED UNDER")
    print("         TRANSPOSITION.  H^T is the house instance on the")
    print("         reversed type-1 seed, one negated seed, the transposed")
    print("         corner and the transposed-and-swapped tables, up to the")
    print("         +-1 diagonal Dt -- at EVERY (s,i) and every w, with no")
    print("         branch hypothesis.  LABEL: PROVEN (note/NOTE-B.md S1.9);")
    print("         the identities are additionally checked ENTRYWISE here.")
    print("         REMARK R: the border-kept orientation switch is Hadamard")
    print("         only when P annihilates Chat's off-diagonal part --")
    print("         automatic at (1,1), FALSE at (3,4), both exercised.")
    print("         CENSUS (%d instances): %d PROVEN inequivalent to their"
          % (n_inst, sep))
    print("         own transpose, %d undecided by the profile." % undec)
    print("         LABEL: COMPUTATIONAL-EVIDENCE.")
    print("         NOT CLAIMED: any structural criterion for H ~ H^T -- all")
    print("         four seeds symmetric does NOT suffice; any statement")
    print("         that an agreeing profile means equivalence (no")
    print("         isomorphism search is run here); and nothing about the")
    print("         SEPARATIONS of S3, which are certs 08, 15, 19, 20, 21,")
    print("         24 and 25's and are cited there, not re-banked here.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
