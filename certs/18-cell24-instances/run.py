#!/usr/bin/env python3
"""cert 18 -- the even-s branch: four (2,4) coset-border Hadamard matrices,
the (2,4) border proposition, and Theorem 3.

  THE CELL.  Theorem E's GENERAL branch (note/NOTE-B.md S1.2.1) admits
  Gbar = Z4, S = {chi, chi^3} (the two faithful characters), s = 2, and
  the Gram M = 4(chi + chi^3) = (8, 0, -8, 0).  s is EVEN, so the cell
  lies outside the house branch that Theorem C classifies.  (H2) reads
  sum_q PAF_q = 4n at 0, -8 on K\\0, 0 on the odd cosets, +8 on the coset
  kappa^{-1}(2).

  FOUR MATRICES -- note/NOTE-B.md S2.3.  H(88) twice (n = 20, w = 5,
  rho-bar 0 and 1) and H(56) twice (n = 12, w = 3, rho-bar 0 and 3), each
  satisfying (H1)-(H4) literally, each assembled here and accepted by
  verify/verify.py, each at a pinned canonical digest.  None is Hadamard-
  equivalent to H(2) (x) H' -- the only Kronecker factorisation an order
  88 or 56 admits -- and S generates the dual of Gbar, so none is a
  collapse of a smaller-index construction.  These are the first
  coset-border Hadamard matrices with even s known to this laboratory;
  the ORDERS were never open, and no novelty of existence is claimed.

  THE (2,4) BORDER PROPOSITION -- note/NOTE-B.md S2.4.  (a) (H1) holds
  iff Q[4I+c+2] = -Q[4I+c] and the eight rows Q[4I+c], c in {0,1}, form
  a Hadamard matrix Q' of order 8.  (b) if w > 2s = 4 then in every
  complete instance the row table is anti-periodic, P' in H(8), and
  E = -(1/16) P Chat^T Q is in H(8).  (c) for EVERY admissible S-part
  (112 of them), EVERY kappa(rho) (4) and EVERY Q' (480 right-orbits of
  H(8)) a kit exists -- 215 040 classes.  So the border is never the
  obstruction at (2,4), at any w.

  THEOREM 3 -- note/NOTE-B.md S2.4.  Under (H1), (H2), s >= 1, w > s
  (Theorem E'), if S contains no real character then 4 | i, s is even
  and N == 0 (mod 8).  Hence at N == 4 (mod 8) -- in particular at
  N = 2092 -- every realisable Gram has a live real character.

(a), (b) and Theorem 3 are PAPER proofs (note/NOTE-B.md S2.4).  This
certificate carries the instances, the census (c), and a small-case
exhaustion around Theorem 3.

WHAT THIS SCRIPT DOES  (standard library only, exact integers only)

  [0] data/cell24-records.json, SHA-256 pinned here.

  [1] THE FOUR INSTANCES.  Per record: the coset sums and the S-part are
      RECOMPUTED from the seeds and required to equal the banked ones;
      (H2) checked lag by lag; Parseval at the live characters
      (sum_q |sigma-hat_q(chi)|^2 = N - w Mhat(chi) = 8); the row-sum law
      sum_q r_q^2 = N; Q' and P' Hadamard of order 8; the anti-periodic
      doubling to Q (16x8) and P (8x16); (H1) Q Q^T = I4 (x) M;
      P P^T = 16 I; E RECOMPUTED here as -(1/16) P Chat^T Q from this
      certificate's own Goethals-Seidel array of the coset sums over Z4
      (the compression lemma) and required to equal the banked corner;
      E in H(8); (H3) E E^T + w P P^T = N I; (H4) E Q^T + P Chat^T = 0;
      the Sigma-bar law Chat Chat^T = Chat^T Chat = I4 (x) dev(Sigma-bar);
      then assembly by this certificate's own block-explicit assembler,
      verify/verify.py (the trust chain), and the canonical digest pinned
      three ways -- computed here, reported by verify.py, and carried by
      the data file.

  [2] STRUCTURE, with controls.  (i) NOT a Kronecker product: if
      H ~ H(2) (x) H' then the rows split into N/2 disjoint pairs whose
      pointwise products all equal +-z for one vector z, and pointwise
      products are invariant under Hadamard equivalence up to global sign
      and column permutation; so the largest sign-normalised
      pointwise-product class must have >= N/2 members.  Controls:
      Sylvester H(16) and a scrambled copy pass the test (largest class
      >= 8); Paley H(12), which admits no Kronecker factorisation, fails
      it.  (ii) NOT a collapse: the common kernel of S in Gbar = Z4 is
      {0}, so S generates the dual.

  [3] THEOREM 3, small cases.  For every abelian Gbar of order i <= 16
      and every Galois-stable S containing no real character: |S| is even
      (the theorem's first clause, checked on every case); and where
      4 does not divide i and s <= 2, an EXHAUSTIVE depth-first search
      for one block Q_0 in {+-1}^{i x 4s} with Q_0 Q_0^T = M_S finds
      nothing -- (H1) is unsatisfiable, as the theorem says.  At the
      (2,4) cell itself, where 4 | i, the same search is allowed to
      succeed and does (a sanity control: the theorem forbids nothing
      there, and this certificate's own instances exist).

  [4] THE (2,4) BORDER CENSUS, EXHAUSTIVE.  The 480 right-orbit
      representatives of H(8) are built from the 30 affine (AG(3,2))
      structures on eight labels times 16 plane-sign classes, checked
      pairwise distinct by canonical form, checked to catch 3000
      deterministically generated labelled H(8), and the count 480 is
      re-derived independently from a backtracking count of NORMALISED
      H(8) (151 200; 151 200 * 2^15 / (8! * 2^8) = 480).  The 112
      S-parts are enumerated from Parseval and parity.  For each of the
      112 x 4 x 480 = 215 040 classes the admissible anti-periodic rows
      are found and an 8-clique in their orthogonality graph is searched
      for; a kit must be found for every class.  A deterministic sample
      of the kits is then re-verified in exact integers with a FULL sigma
      (an arbitrary S^c part added) at w = 5 and at w = 130.

Usage:
  python certs/18-cell24-instances/run.py
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
import time
from math import gcd

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
VERIFY = os.path.join(ROOT, "verify", "verify.py")
DATA = os.path.join(ROOT, "data", "cell24-records.json")
OUT = os.path.join(HERE, "out")

_T0 = time.time()
FAIL = []
NCHECK = [0]

S_CELL, I_CELL = 2, 4
M_GRAM = {0: 8, 1: 0, 2: -8, 3: 0}

# ---------------------------------------------------------------- the pins
FILE_PINS = {
    "data/cell24-records.json":
        "9727b392940d416d3f25dca5d51d2db71cd499bc73c3b8dc4efd22801180f179",
}

# Canonical SHA-256 of the four matrices (the digest verify/verify.py
# reports).  Also carried by the data file, and compared against it.
SHA = {
    "H88_cell24_n20_seed0_rho0":
        "942b3f32fcd75e72a64f92d9c294b0d0cedbbd0965fe5a14213c30b8b66ffc8a",
    "H88_cell24_n20_seed1_rho1":
        "4cae47d1c5054a86ca48154c5e9cd99845294be4d275dd8560bf6b05fe5f08e7",
    "H56_cell24_n12_seed0_rho0":
        "ad67ee2c9d1f4d0343b250824dec301759e097d72d949f1b9c91f22cab026b85",
    "H56_cell24_n12_seed2_rho3":
        "fa82808c8dbc0245f3d427312975183fb947438bdcafba56107ee2280d6e4aff",
}

# The census shape, pinned so a change of engine cannot quietly turn this
# certificate into a smaller statement.
N_ORBITS = 480
N_SPARTS = 112
N_CLASSES = 215040
NORMALISED_H8 = 151200
ADM_HIST = {8: 114688, 16: 86016, 32: 14336}


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


def pm(s):
    assert set(s) <= {"+", "-"}, s
    return [1 if ch == "+" else -1 for ch in s]


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


# ======================================================================
# small exact linear algebra
# ======================================================================

def mmT(A, B):
    """A * B^T."""
    return [[sum(x * y for x, y in zip(ra, rb)) for rb in B] for ra in A]


def mm(A, B):
    Bt = list(zip(*B))
    return [[sum(x * y for x, y in zip(ra, cb)) for cb in Bt] for ra in A]


def T_(A):
    return [list(r) for r in zip(*A)]


def paf(x, t):
    n = len(x)
    return sum(x[u] * x[(u + t) % n] for u in range(n))


def is_hadamard_small(X):
    n = len(X)
    G = mmT(X, X)
    return all(G[a][b] == (n if a == b else 0)
               for a in range(n) for b in range(n))


def gs_blocks(seqs, m, rho):
    """The four block types over Z_m, standard orientation (NOTE-B S1.0):
    A = dev(x0); X R has (g,h) -> x(rho-g-h); X^T R has (g,h) -> x(g+h-rho)."""
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


def gs_array(seqs, m, rho):
    B = gs_blocks(seqs, m, rho)
    out = []
    for I in range(4):
        for g in range(m):
            out.append([B[I][J][g][h] for J in range(4) for h in range(m)])
    return out


def double_Q(Qp):
    """Q[4I+c] = (+1 if c < 2 else -1) * Q'[2I + (c mod 2)]  (16 x 8)."""
    return [[(1 if c < 2 else -1) * v for v in Qp[2 * I + (c % 2)]]
            for I in range(4) for c in range(4)]


def double_P(Pp):
    """P[r][4J+c] = (+1 if c < 2 else -1) * P'[r][2J + (c mod 2)]  (8 x 16)."""
    return [[(1 if c < 2 else -1) * row[2 * J + (c % 2)]
             for J in range(4) for c in range(4)] for row in Pp]


# ======================================================================
# clause 1 -- the four instances
# ======================================================================

def write_and_verify(rows, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name + ".txt")
    text = "\n".join("".join("+" if v == 1 else "-" for v in r)
                     for r in rows) + "\n"
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write(text)
    sha = hashlib.sha256(text.encode("ascii")).hexdigest()
    proc = subprocess.run([sys.executable, VERIFY, path],
                          capture_output=True, text=True,
                          env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    verdict = [l for l in (proc.stdout + proc.stderr).splitlines()
               if l.startswith("VERDICT")]
    verdict = verdict[-1] if verdict else "(no verdict)"
    print("      verify.py: exit %d :: %s" % (proc.returncode, verdict[:100]))
    os.unlink(path)
    return proc.returncode, verdict, sha


def clause_1_instance(rec):
    n, rho = int(rec["group"][0]), int(rec["r_shift"][0])
    i, s = I_CELL, S_CELL
    w = n // i
    N = 4 * (n + s)
    log("[1] instance %s  (n = %d, w = %d, rho = %d, N = %d)"
        % (rec["name"], n, w, rho, N))
    xs = [pm(x) for x in rec["seeds"]]
    check("four +-1 seeds of length n = %d; the record's declared order, "
          "group, s, coset divisors and w are consistent" % n,
          len(xs) == 4 and all(len(x) == n for x in xs)
          and int(rec["order"]) == N and int(rec["s"]) == s
          and rec["coset_divisors"] == [4] and int(rec["w"]) == w
          and int(rec["rho_bar"]) == rho % 4)

    # --- the coset sums and the S-part, recomputed from the seeds
    sigma = [[sum(x[g] for g in range(n) if g % 4 == c) for c in range(4)]
             for x in xs]
    tau = [[sg[0] - sg[2], sg[1] - sg[3]] for sg in sigma]
    check("coset sums and S-part RECOMPUTED from the seeds == the banked "
          "'coset_sums' and 'tau'",
          sigma == [list(map(int, r)) for r in rec["coset_sums"]]
          and tau == [list(map(int, r)) for r in rec["tau"]],
          "tau = %s" % (tau,))

    # --- (H2) and the two Parseval consequences
    agg = [sum(paf(x, t) for x in xs) for t in range(n)]
    want = [4 * n if t == 0 else -M_GRAM[t % 4] for t in range(n)]
    check("(H2): aggregate PAF = 4n at 0, -8 on K\\0, 0 on the odd cosets, "
          "+8 on the coset 2", agg == want,
          "first mismatch %s" % next(((t, agg[t], want[t]) for t in range(n)
                                      if agg[t] != want[t]), None))
    energy = sum(a * a + b * b for a, b in tau)
    check("Parseval at the live characters: sum_q |sigma-hat_q(chi)|^2 = "
          "N - w*Mhat(chi) = 8", energy == 8)
    check("row-sum law: sum_q r_q^2 = N - w*sum_c M(c) = %d" % N,
          sum(sum(x) ** 2 for x in xs) == N)

    # --- the border kit
    Qp = [pm(r) for r in rec["col_table_8"]]
    Pp = [pm(r) for r in rec["row_table_8"]]
    check("Q' and P' are 8x8 Hadamard matrices",
          is_hadamard_small(Qp) and is_hadamard_small(Pp))
    Q = double_Q(Qp)
    P = double_P(Pp)
    QQt = mmT(Q, Q)
    check("(H1): Q Q^T = I4 (x) M with M = (8, 0, -8, 0)",
          all(QQt[a][b] == (M_GRAM[(a - b) % 4] if a // 4 == b // 4 else 0)
              for a in range(16) for b in range(16)))
    PPt = mmT(P, P)
    check("P P^T = 4i I = 16 I",
          all(PPt[a][b] == (16 if a == b else 0)
              for a in range(8) for b in range(8)))
    Chat = gs_array(sigma, 4, rho % 4)
    PCt = mmT(P, Chat)                                    # 8 x 16
    E16 = mm(PCt, Q)                                      # 8 x 8
    check("E = -(1/16) P Chat^T Q has every entry in {+-16}, so E is +-1",
          all(abs(v) == 16 for row in E16 for v in row))
    E = [[-v // 16 for v in row] for row in E16]
    check("the RECOMPUTED E equals the banked corner (the data file's "
          "'corner' is checked, not trusted)",
          E == [pm(r) for r in rec["corner"]])
    check("E is an 8x8 Hadamard matrix", is_hadamard_small(E))
    EEt = mmT(E, E)
    check("(H3): E E^T + w P P^T = N I   (8 + 16w = %d)" % N,
          all(EEt[a][b] + w * PPt[a][b] == (N if a == b else 0)
              for a in range(8) for b in range(8)))
    EQt = mmT(E, Q)
    check("(H4): E Q^T + P Chat^T = 0",
          all(EQt[a][b] + PCt[a][b] == 0
              for a in range(8) for b in range(16)))
    CCt, CtC = mmT(Chat, Chat), mm(T_(Chat), Chat)
    sbar = [4 * n - 8 * (w - 1), 0, 8 * w, 0]
    check("Sigma-bar law: Chat Chat^T = Chat^T Chat = I4 (x) dev(%s)" % sbar,
          CCt == CtC
          and all(CCt[a][b] == (sbar[(a - b) % 4] if a // 4 == b // 4 else 0)
                  for a in range(16) for b in range(16)))

    # --- assemble, verify, pin
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
    check("assembled a %d x %d +-1 matrix" % (N, N),
          len(rows) == N and all(len(r) == N and set(r) <= {1, -1}
                                 for r in rows))
    rc, verdict, sha = write_and_verify(rows, rec["name"])
    check("verify.py exit 0, HADAMARD order=%d" % N,
          rc == 0 and "HADAMARD" in verdict and "order=%d" % N in verdict)
    check("canonical SHA-256 pinned three ways -- this run, verify.py, and "
          "the data file: %s..." % SHA[rec["name"]][:16],
          sha == SHA[rec["name"]] and SHA[rec["name"]] in verdict
          and rec["canonical_sha256"] == SHA[rec["name"]], sha)
    return rows, sha


# ======================================================================
# clause 2 -- structure, with controls
# ======================================================================

def kronecker_class(rows):
    """The largest multiplicity of a sign-normalised pointwise row product."""
    N = len(rows)
    packed = []
    for r in rows:
        bits = 0
        for j, v in enumerate(r):
            if v == 1:
                bits |= 1 << j
        packed.append(bits)
    full = (1 << N) - 1
    counts = {}
    for a in range(N):
        for b in range(a + 1, N):
            x = packed[a] ^ packed[b]        # 1 where the product is -1
            if x & 1:
                x ^= full                    # normalise: first coordinate +
            counts[x] = counts.get(x, 0) + 1
    return max(counts.values())


def clause_2_controls():
    log("[2] controls for the Kronecker invariant")

    def sylvester(k):
        h = [[1]]
        while len(h) < k:
            h = [r + r for r in h] + [r + [-v for v in r] for r in h]
        return h
    H16 = sylvester(16)
    m = kronecker_class(H16)
    check("Sylvester H(16) = H(2) (x) H(8): largest class %d >= N/2 = 8 -- "
          "the test PASSES a genuine Kronecker product" % m, m >= 8)
    rng = random.Random(22)
    perm_r = list(range(16))
    rng.shuffle(perm_r)
    perm_c = list(range(16))
    rng.shuffle(perm_c)
    sr = [rng.choice((1, -1)) for _ in range(16)]
    sc = [rng.choice((1, -1)) for _ in range(16)]
    Hs = [[sr[a] * sc[b] * H16[perm_r[a]][perm_c[b]] for b in range(16)]
          for a in range(16)]
    check("a scrambled (Hadamard-equivalent) copy gives the same largest "
          "class -- the invariant is an invariant",
          kronecker_class(Hs) == m)
    q = 11
    chi = {0: 0}
    for a in range(1, q):
        chi[a] = -1
    for a in range(1, q):
        chi[(a * a) % q] = 1
    Qm = [[chi[(b - a) % q] for b in range(q)] for a in range(q)]
    H12 = ([[1] * 12]
           + [[1] + [(Qm[a][b] if a != b else -1) for b in range(q)]
              for a in range(q)])
    check("Paley H(12) is Hadamard (control input)", is_hadamard_small(H12))
    check("Paley H(12): largest class %d < 6 -- order 12 admits no Kronecker "
          "factorisation, and the test says so" % kronecker_class(H12),
          kronecker_class(H12) < 6)


def clause_2_instance(rec, rows):
    N = len(rows)
    kc = kronecker_class(rows)
    check("%-28s (i) NOT H(2) (x) H': largest sign-normalised pointwise-"
          "product class has %d row pairs, far below N/2 = %d"
          % (rec["name"], kc, N // 2), kc < N // 2)
    ker = [c for c in range(4) if all((k * c) % 4 == 0 for k in (1, 3))]
    check("%-28s (ii) NOT a collapse: the common kernel of S = {chi, chi^3} "
          "in Gbar = Z4 is {0}, so S generates the dual" % rec["name"],
          ker == [0])


# ======================================================================
# clause 3 -- Theorem 3, small cases
# ======================================================================

class Ab:
    def __init__(self, f):
        self.f = tuple(f)
        self.elts = list(itertools.product(*[range(m) for m in self.f]))
        self.zero = self.elts[0]
        self.e = 1
        for m in self.f:
            self.e = self.e * m // gcd(self.e, m)

    def add(self, a, b):
        return tuple((x + y) % m for x, y, m in zip(a, b, self.f))

    def sub(self, a, b):
        return tuple((x - y) % m for x, y, m in zip(a, b, self.f))

    def order(self, a):
        k, x = 1, a
        while x != self.zero:
            x = self.add(x, a)
            k += 1
        return k

    def chi_exp(self, a, c):
        return sum(x * y * (self.e // m)
                   for x, y, m in zip(a, c, self.f)) % self.e


def abelian_groups(n):
    def partitions(k):
        if k == 0:
            yield ()
            return
        for first in range(k, 0, -1):
            for rest in partitions(k - first):
                if not rest or rest[0] <= first:
                    yield (first,) + rest
    fac, m, p = {}, n, 2
    while m > 1:
        while m % p == 0:
            fac[p] = fac.get(p, 0) + 1
            m //= p
        p += 1
    per = [[tuple(p ** x for x in part) for part in partitions(a)]
           for p, a in sorted(fac.items())]
    return [tuple(sorted(sum(c, ()), reverse=True))
            for c in itertools.product(*per)] or [(1,)]


def ramanujan(d, j):
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
    seen, orbits = set(), []
    for a in G.elts:
        if a in seen:
            continue
        d = G.order(a)
        orb = frozenset(tuple((k * x) % m for x, m in zip(a, G.f))
                        for k in range(1, d + 1) if gcd(k, d) == 1)
        seen |= orb
        orbits.append((d, orb))
    return orbits


def is_real_char(G, a):
    half = G.e // 2 if G.e % 2 == 0 else -1
    return all(G.chi_exp(a, c) in (0, half) for c in G.elts)


def M_S(G, orbs):
    out = {}
    for c in G.elts:
        tot = 0
        for d, orb in orbs:
            a = next(iter(orb))
            x = G.chi_exp(a, c)
            tot += ramanujan(d, x // (G.e // d))
        out[c] = 4 * tot
    return out


def dfs_block(G, Mf, s, limit_nodes=2000000):
    """Find Q_0 in {+-1}^{i x 4s} with Q_0 Q_0^T[c,c'] = M(c - c'), or prove
    that none exists.  Exhaustive."""
    i = len(G.elts)
    vecs = list(itertools.product((1, -1), repeat=4 * s))
    rows = []
    nodes = [0]
    els = G.elts

    def rec(k):
        nodes[0] += 1
        if nodes[0] > limit_nodes:
            raise TimeoutError
        if k == i:
            return True
        for v in vecs:
            if k == 0 and v != vecs[0]:
                continue
            if all(sum(x * y for x, y in zip(v, rows[j]))
                   == Mf[G.sub(els[k], els[j])] for j in range(k)):
                rows.append(v)
                if rec(k + 1):
                    return True
                rows.pop()
        return False
    try:
        ok = rec(0)
    except TimeoutError:
        return "TIMEOUT", nodes[0]
    return ok, nodes[0]


def clause_3():
    log("[3] Theorem 3: Galois-stable S with no real character, every "
        "abelian Gbar of order i <= 16")
    odd_sizes = []
    exhausted = 0
    unexpected = []
    sanity_ok = None
    for i in range(2, 17):
        for f in abelian_groups(i):
            G = Ab(f)
            orbs = galois_orbits(G)
            nonreal = [(d, o) for d, o in orbs
                       if not any(is_real_char(G, a) for a in o)]
            gname = "x".join("Z%d" % m for m in f)
            for r in range(1, len(nonreal) + 1):
                for combo in itertools.combinations(nonreal, r):
                    s = sum(len(o) for d, o in combo)
                    if s % 2:
                        odd_sizes.append((gname, s))
                    if s > 2 and (i, s) != (4, 2):
                        continue          # DFS budget: 4s <= 8 columns
                    Mf = M_S(G, list(combo))
                    if i % 4 != 0:
                        res, nodes = dfs_block(G, Mf, s)
                        exhausted += 1
                        if res is not False:
                            unexpected.append((i, gname, s, res))
                        print("      i=%2d %-9s |S|=%d  4|i: no   exhaustive "
                              "DFS for Q_0 (%dx%d): %s (%d nodes)"
                              % (i, gname, s, i, 4 * s, res, nodes),
                              flush=True)
                    elif (i, s) == (4, 2):
                        res, nodes = dfs_block(G, Mf, s)
                        sanity_ok = (res is True)
                        print("      i=%2d %-9s |S|=%d  4|i: yes  exhaustive "
                              "DFS for Q_0 (%dx%d): %s (%d nodes)  [the (2,4) "
                              "cell -- the theorem forbids nothing here]"
                              % (i, gname, s, i, 4 * s, res, nodes),
                              flush=True)
    check("|S| is EVEN in every Galois-stable S with no real character over "
          "every abelian Gbar of order <= 16 (the theorem's first clause)",
          not odd_sizes, str(odd_sizes[:5]))
    check("where 4 does not divide i (and s <= 2, the search budget), an "
          "EXHAUSTIVE DFS finds no block Q_0 with Q_0 Q_0^T = M_S: %d cases, "
          "all empty" % exhausted, not unexpected, str(unexpected[:3]))
    check("sanity at the (2,4) cell, where 4 | i: the same search DOES find "
          "a Q_0 (the theorem forbids nothing there)", sanity_ok is True)


# ======================================================================
# clause 4 -- the (2,4) border census
# ======================================================================

def sylvester8():
    return [[(-1) ** bin(x & y).count("1") for y in range(8)]
            for x in range(8)]


def count_normalised_h8():
    """Backtracking count of the 8x8 Hadamard matrices whose first row and
    first column are all +."""
    pool = [v for v in itertools.product((1, -1), repeat=8)
            if v[0] == 1 and sum(v) == 0]
    cnt = [0]

    def rec(depth, cands):
        if depth == 8:
            cnt[0] += 1
            return
        for v in cands:
            rec(depth + 1,
                [u for u in cands if sum(a * b for a, b in zip(u, v)) == 0])
    rec(1, pool)
    return cnt[0], len(pool)


def canon8(Q):
    """Canonical form under RIGHT signed permutations: normalise each column
    by its first entry, then sort the columns."""
    sg = Q[0]
    return tuple(sorted(tuple(Q[k][y] * sg[y] for k in range(1, 8))
                        for y in range(8)))


def orbit_reps():
    """The 480 right-orbit representatives of H(8), built deterministically:
    30 affine AG(3,2) structures on eight labels x 16 plane-sign classes."""
    H0 = sylvester8()
    planes = [frozenset(y for y in range(8)
                        if bin(a & y).count("1") % 2 == b)
              for a in range(1, 8) for b in (0, 1)]
    assert len(planes) == 14
    structs = {}
    for tail in itertools.permutations(range(1, 8)):
        pi = (0,) + tail                      # pi: label k -> point pi(k)
        inv = {pi[k]: k for k in range(8)}
        key = frozenset(frozenset(inv[y] for y in pl) for pl in planes)
        if key not in structs:
            structs[key] = pi
    reps = []
    for key, pi in structs.items():
        free = [k for k in range(8) if pi[k] in (3, 5, 6, 7)]
        for signs in itertools.product((1, -1), repeat=4):
            d = [1] * 8
            for k, sg in zip(free, signs):
                d[k] = sg
            reps.append([[d[k] * H0[pi[k]][y] for y in range(8)]
                         for k in range(8)])
    return reps, len(structs)


def sparts_24():
    """Every admissible (2,4) S-part: tau_q = (sigma_q(0)-sigma_q(2),
    sigma_q(1)-sigma_q(3)) with entries even and sum_q |tau_q|^2 = 4s = 8.
    Both signs kept."""
    vals = (0, 2, -2)
    return [t for t in itertools.product(itertools.product(vals, repeat=2),
                                         repeat=4)
            if sum(a * a + b * b for a, b in t) == 8]


def sigma0_from_tau(tau):
    """sigma0_q(c) = (1/2) Re((alpha + i beta) i^{-c}) = (a/2, b/2, -a/2,
    -b/2)."""
    return [[a // 2, b // 2, -a // 2, -b // 2] for a, b in tau]


def antiperiodic_rows():
    """All anti-periodic p in {+-1}^16 with p[0] = +1:
    p[4J+c] = t[2J + (c mod 2)] * (+1 if c < 2 else -1)."""
    out = []
    for t in itertools.product((1, -1), repeat=8):
        if t[0] != 1:
            continue
        out.append(tuple(t[2 * J + (c % 2)] * (1 if c < 2 else -1)
                         for J in range(4) for c in range(4)))
    return out


def clause_4(reverify_sample=200):
    log("[4] the (2,4) border census: every (S-part, kappa(rho), Q') class")

    # --- the 480 orbit representatives, two ways
    t0 = time.time()
    n0, npool = count_normalised_h8()
    total_labelled = n0 << 15               # 2^15 = 2^8 row signs x 2^7 more
    check("normalised H(8) count by backtracking: %d (pool of %d balanced "
          "rows); total labelled T = %d; T / (8! 2^8) = %d orbits under right "
          "signed permutations (%.1fs)"
          % (n0, npool, total_labelled,
             total_labelled // (40320 * 256), time.time() - t0),
          n0 == NORMALISED_H8 and total_labelled % (40320 * 256) == 0
          and total_labelled // (40320 * 256) == N_ORBITS)
    reps, nstruct = orbit_reps()
    check("orbit representatives built from %d affine AG(3,2) structures x "
          "16 plane-sign classes: %d of them" % (nstruct, len(reps)),
          nstruct == 30 and len(reps) == N_ORBITS)
    check("every representative is in fact Hadamard",
          all(is_hadamard_small(Q) for Q in reps))
    canons = {canon8(Q) for Q in reps}
    check("the %d representatives are pairwise distinct as right-orbits "
          "(%d canonical forms) -- so they are ALL of them, by the count "
          "above" % (N_ORBITS, len(canons)), len(canons) == N_ORBITS)
    rng = random.Random(24)
    H0 = sylvester8()
    stray = 0
    for _ in range(3000):
        pr = list(range(8))
        rng.shuffle(pr)
        pc = list(range(8))
        rng.shuffle(pc)
        rs = [rng.choice((1, -1)) for _ in range(8)]
        cs = [rng.choice((1, -1)) for _ in range(8)]
        Q = [[rs[k] * cs[y] * H0[pr[k]][pc[y]] for y in range(8)]
             for k in range(8)]
        if canon8(Q) not in canons:
            stray += 1
    check("3000 deterministically generated labelled H(8) all land in the "
          "representative set (every H(8) is equivalent to Sylvester)",
          stray == 0, "strays = %d" % stray)

    # --- the admissible-row characterisation, per Q'
    #   a row p is admissible for Q' iff a(p) := the 8-vector
    #   a[2I+j] = u[4I+j] - u[4I+j+2],  u = Chat0 p,
    #   satisfies (a Q')[y] = +-16 for every y.  Multiplying on the right by
    #   Q'^T (Q' Q'^T = 8I) that says exactly a = 2 v Q'^T for some
    #   v in {+-1}^8, and then E's row is -v.  So the admissible a for a
    #   given Q' are a set of 256 vectors, computed once.
    ADM = []
    for Qp in reps:
        ADM.append({tuple(2 * sum(v[y] * Qp[k][y] for y in range(8))
                          for k in range(8)): v
                    for v in itertools.product((1, -1), repeat=8)})
    PALL = antiperiodic_rows()
    check("anti-periodic sign representatives (p[0] = +1): %d" % len(PALL),
          len(PALL) == 128)
    ORTH = []
    for a in range(128):
        m = 0
        for b in range(128):
            if sum(x * y for x, y in zip(PALL[a], PALL[b])) == 0:
                m |= 1 << b
        ORTH.append(m)

    def clique8(mask):
        out = []

        def rec(cand, need):
            if need == 0:
                return True
            while cand:
                v = cand.bit_length() - 1
                cand &= ~(1 << v)
                if bin(cand).count("1") + 1 < need:
                    return False
                out.append(v)
                if rec(cand & ORTH[v], need - 1):
                    return True
                out.pop()
            return False
        return out if rec(mask, 8) else None

    taus = sparts_24()
    check("admissible (2,4) S-parts (both signs kept): %d; with kappa(rho): "
          "%d classes per Q'" % (len(taus), 4 * len(taus)),
          len(taus) == N_SPARTS)

    # --- the census
    t0 = time.time()
    hits = total = 0
    misses = []
    adm_hist = {}
    kits = {}
    for ti, tau in enumerate(taus):
        for rho in range(4):
            C0 = gs_array(sigma0_from_tau(tau), 4, rho)
            sparse = [[(m, C0[k][m]) for m in range(16) if C0[k][m]]
                      for k in range(16)]
            AV = []
            for p in PALL:
                u = [sum(c * p[m] for m, c in row) for row in sparse]
                AV.append(tuple(u[4 * I + j] - u[4 * I + j + 2]
                                for I in range(4) for j in range(2)))
            for qi in range(N_ORBITS):
                D = ADM[qi]
                mask = 0
                nadm = 0
                for pj in range(128):
                    if AV[pj] in D:
                        mask |= 1 << pj
                        nadm += 1
                adm_hist[nadm] = adm_hist.get(nadm, 0) + 1
                total += 1
                cl = clique8(mask) if nadm >= 8 else None
                if cl is None:
                    misses.append((ti, rho, qi))
                else:
                    hits += 1
                    kits[(ti, rho, qi)] = (tuple(cl),
                                           tuple(D[AV[a]] for a in cl))
    secs = time.time() - t0
    check("kits found for %d / %d (S-part, kappa(rho), Q') classes -- the "
          "border is never the obstruction at (2,4) (%.1fs)"
          % (hits, total, secs),
          total == N_CLASSES and hits == N_CLASSES and not misses,
          str(misses[:3]))
    check("admissible anti-periodic rows per class, histogram: %s"
          % dict(sorted(adm_hist.items())), adm_hist == ADM_HIST)

    # --- exact re-verification, with a full sigma, at two widths
    rng = random.Random(2404)
    keys = sorted(kits)
    picked = [keys[(k * 1009 + 17) % len(keys)]
              for k in range(reverify_sample)]
    okc = 0
    for (ti, rho, qi) in picked:
        cl, vs = kits[(ti, rho, qi)]
        Qp = reps[qi]
        Q = double_Q(Qp)
        P = [list(PALL[a]) for a in cl]
        E = [[-v for v in vs[r]] for r in range(8)]
        tau = taus[ti]
        s0 = sigma0_from_tau(tau)
        # add an arbitrary S^c component: a constant plus b*(-1)^c per seed
        rp = [(rng.randrange(-9, 10), rng.randrange(-9, 10)) for _ in range(4)]
        sigma = [[s0[q][c] + rp[q][0] + rp[q][1] * (-1) ** c
                  for c in range(4)] for q in range(4)]
        Ch = gs_array(sigma, 4, rho)
        QQt = mmT(Q, Q)
        EEt, PPt = mmT(E, E), mmT(P, P)
        H4 = mmT(E, Q)
        PCt = mmT(P, Ch)
        good = (all(QQt[a][b] == (M_GRAM[(a - b) % 4]
                                  if a // 4 == b // 4 else 0)
                    for a in range(16) for b in range(16))
                and all(EEt[a][b] == (8 if a == b else 0)
                        for a in range(8) for b in range(8))
                and all(PPt[a][b] == (16 if a == b else 0)
                        for a in range(8) for b in range(8))
                and all(abs(v) == 1 for row in E for v in row)
                and all(H4[a][b] + PCt[a][b] == 0
                        for a in range(8) for b in range(16)))
        for w in (5, 130):
            N = 4 * (4 * w + 2)
            good = good and all(
                EEt[a][b] + w * PPt[a][b] == (N if a == b else 0)
                for a in range(8) for b in range(8))
        okc += bool(good)
    check("%d sampled kits re-verified in exact integers with a FULL sigma "
          "(an arbitrary S^c part added) -- (H1), E in H(8), P P^T = 16 I, "
          "(H4), and (H3) at BOTH w = 5 and w = 130: %d / %d"
          % (reverify_sample, okc, reverify_sample), okc == reverify_sample)
    check("(H3) is w-free at (2,4): it reads 8 + 16w = N = 4(4w+2), an "
          "identity in w -- which is why one kit serves every order",
          all(8 + 16 * w == 4 * (4 * w + 2) for w in range(1, 500)))


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reverify", type=int, default=200,
                    help="how many census kits to re-verify exactly with a "
                         "full sigma (default 200; the census itself is "
                         "always exhaustive)")
    args = ap.parse_args(argv)

    print("=" * 78)
    print("cert 18 -- the even-s branch at (s,i) = (2,4): H(88) x2, H(56) x2,")
    print("           the border census, and Theorem 3")
    print("           trust chain: %s"
          % os.path.relpath(VERIFY, ROOT).replace("\\", "/"))
    print("=" * 78)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    try:
        rc = _body(args)
    finally:
        shutil.rmtree(OUT, ignore_errors=True)
    print("generated matrices deleted; nothing left in %s   (%.1fs total)"
          % (os.path.relpath(OUT, ROOT).replace("\\", "/"),
             time.time() - _T0))
    return rc


def _body(args):
    log("[0] the banked data file, SHA-256 pinned in this script")
    for name, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, name))
        check("%-30s" % name, got == want, got[:24] + "...")
    with open(DATA, encoding="ascii") as fh:
        blob = json.load(fh)
    check("data file declares the cell (s,i) = (2,4), Gbar = Z4, "
          "M = (8,0,-8,0)",
          blob["cell"]["s"] == S_CELL and blob["cell"]["i"] == I_CELL
          and blob["cell"]["gbar"] == "Z4"
          and blob["cell"]["M"] == [M_GRAM[c] for c in range(4)])
    recs = blob["records"]
    check("four records", len(recs) == 4
          and sorted(r["name"] for r in recs) == sorted(SHA))

    clause_2_controls()
    built = []
    for rec in recs:
        rows, sha = clause_1_instance(rec)
        built.append((rec, rows, sha))
    log("[2] the two structural claims, on each instance")
    for rec, rows, _sha in built:
        clause_2_instance(rec, rows)

    clause_3()
    clause_4(args.reverify)

    log("SUMMARY")
    for rec, rows, sha in built:
        print("    %-28s N=%3d  sha %s" % (rec["name"], len(rows), sha))

    print("\n" + "=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 18: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: the even-s branch of Theorem E is INHABITED.  At the")
    print("         cell (s,i) = (2,4), Gbar = Z4, S = {chi, chi^3},")
    print("         M = (8,0,-8,0): four Hadamard matrices -- H(88) x2 and")
    print("         H(56) x2 -- satisfy (H1)-(H4) literally, assemble here")
    print("         and pass verify/verify.py at their pinned digests; none")
    print("         is equivalent to H(2) (x) H' and none is a collapse; the")
    print("         border layer never obstructs (215 040 / 215 040 classes")
    print("         admit a kit, exhaustively); and Theorem 3's small-case")
    print("         exhaustion is consistent in every case tested.")
    print("         LABEL: PROVEN-BY-CERTIFICATE (the matrices, the census,")
    print("         the two structural claims) + PROVEN (the border")
    print("         proposition (a),(b) and Theorem 3, paper-grade,")
    print("         note/NOTE-B.md S2.4).")
    print("         NOT claimed: anything about H(2092) -- the (2,4) cell")
    print("         does not land there; any novelty of the ORDERS 56 and")
    print("         88, which were never open; any equivalence statement")
    print("         among the four, or to other known H(56)/H(88); that the")
    print("         anti-periodic row table is FORCED at w <= 4.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
