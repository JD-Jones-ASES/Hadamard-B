#!/usr/bin/env python3
"""cert 32 -- the good-matrix product theorem: the small-order controls, the
group-ring proof instantiated, the Williamson control, and the closure sign.

  THE THEOREM (note/NOTE-B.md S1.13; PROVEN there, paper-grade).  Let n be
  odd and let (A, B, C, D) be a Goethals-Seidel quadruple on Z_n with A of
  SKEW TYPE (a_0 = +1, a_{-k} = -a_k) and B, C, D SYMMETRIC.  Then

        a_{2k} = -(b_0 c_0 d_0) * a_k b_k c_k d_k      for every k != 0.

  Under the normalisation b_0 = c_0 = d_0 = +1 this reads
  a_{2k} = -a_k b_k c_k d_k; WITHOUT the prefactor the identity is FALSE.

  THIS CERTIFICATE PROVES NOTHING.  The proof is in the note.  What it
  carries is the finite evidence around the proof, and four controls that
  can fail:

  [A] THE SMALL-ORDER ENUMERATION, complete.  Every good quadruple at
      n = 7, 9, 11, 13, 15 is found by a meet-in-the-middle over ALL skew and
      ALL symmetric seeds with b_0 free -- no normalisation, no prefilter,
      nothing assumed about the shape of a solution.  Counts
      528 / 288 / 1440 / 3456 / 4224 = 9 936.  The corrected identity holds
      in EVERY one of the 9 936.  Each unprefixed variant -- a_{2k} =
      -a_k b_k c_k d_k and a_{2k} = +a_k b_k c_k d_k -- holds in exactly
      HALF, namely on the half where b_0 c_0 d_0 has the matching sign: the
      prefactor is not decoration.  A sample of the quadruples is assembled
      into a Goethals-Seidel array and handed to verify/verify.py, the trust
      chain, which must return exit 0 on a Hadamard matrix of order 4n.

  [B] THE GROUP-RING PROOF, INSTANTIATED.  The proof's central identity,
      over F_2[Z_n],

            T^[2] + U = e + eps * J,
            T = U + V + W + Z,  eps = 1 + |V| + |W| + |Z|  (mod 2),

      with U, V, W, Z the {0,1} indicators of the -1 positions of A, B, C, D
      and X^[2] the squaring map (X^[2])_{2a} = x_a, is re-checked coefficient
      by coefficient on every one of the 9 936 quadruples.  Its x^{2k}
      coefficient IS the theorem.  CONTROL: on quadruples that are not good
      -- random skew/symmetric quadruples with a nonzero aggregate PAF -- the
      identity must FAIL, and the run checks that it does.

  [C] THE WILLIAMSON CONTROL, which could fail and does, where it should.
      With FOUR SYMMETRIC seeds the same lemma gives U* = U for every seed,
      so T^[2] = e + (sigma + 1) J and the conclusion is the CONSTANT-PRODUCT
      corollary

            a_k b_k c_k d_k = -(a_0 b_0 c_0 d_0)   for every k != 0,

      and NOT the doubling relation.  Complete meet-in-the-middle at
      n = 7, 9, 11, 13, 15: 960 / 2112 / 1920 / 5184 / 4608 = 14 784
      Williamson quadruples; the constant-product corollary holds in ALL
      14 784, and the doubling relation in NONE of them.  (The unprefixed
      doubling variants hold in 48 of the 960 at n = 7 and in none at any
      larger n -- a small-n coincidence, recorded so the reader is not
      surprised by it.)  Skewness of A is load-bearing.

  [D] THE CLOSURE SIGN, COMPUTED AND NEVER ASSUMED.  Doubling permutes
      Z_n \\ {0}.  Along one cycle the theorem determines A from a single
      entry, and consistency around the cycle is exactly ONE parity
      condition.  In the multiplier setting -- M_0 <= Z_n^* of odd order with
      the seeds constant on the orbits of <M_0, -1>, and 2 generating
      Z_n^* / <M_0, -1> in L steps -- that condition reads

            prod_j pi_j = closure * (-(b_0 c_0 d_0))^L,
            pi_j = b_{2^j} c_{2^j} d_{2^j},
            closure = +1 if 2^L in M_0,  -1 if 2^L in -M_0.

      Computed here: at n = 7 with M_0 trivial, L = 3 and 2^3 = 1 in M_0, so
      closure = +1 and prod pi = -1 under the normalisation -- and that
      prediction is CHECKED against all 528 enumerated good quadruples at
      n = 7, with the general prefactor, so it is a control and not a
      restatement.  At n = 523 with M_0 the subgroup of order 3, L = 87 and
      2^87 = 463 lies in -M_0, so closure = -1 and prod pi = +1.  Every step
      of that arithmetic is redone here from the definition of Z_523^*.

  WHAT THIS CERTIFICATE DOES NOT SAY.  Nothing about the existence of any
  Hadamard matrix at any order not exhibited here.  Nothing about n = 523
  beyond the group arithmetic of [D]: no seed quadruple is known there, and
  the closure condition is a necessary condition on a hypothetical one.
  Nothing about good matrices at n > 15.

Standard library only, exact integers only, no floats anywhere, no network.
Generated matrices are written, verified and deleted; nothing is committed.

Usage:
  python certs/32-good-product-theorem/run.py
"""

import argparse
import hashlib
import itertools
import os
import random
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
VERIFY = os.path.join(ROOT, "verify", "verify.py")
OUT = os.path.join(HERE, "out")

ORDERS = (7, 9, 11, 13, 15)
EXPECTED_GOOD = {7: 528, 9: 288, 11: 1440, 13: 3456, 15: 4224}
EXPECTED_WILLIAMSON = {7: 960, 9: 2112, 11: 1920, 13: 5184, 15: 4608}

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
# seeds and PAFs
# ======================================================================

def skew_seeds(n):
    """a_0 = +1, a_{-k} = -a_k: 2^{(n-1)/2} of them, all of them."""
    h = (n - 1) // 2
    out = []
    for bits in itertools.product((1, -1), repeat=h):
        a = [0] * n
        a[0] = 1
        for k in range(1, h + 1):
            a[k] = bits[k - 1]
            a[n - k] = -bits[k - 1]
        out.append(tuple(a))
    return out


def sym_seeds(n):
    """b_{-k} = b_k with b_0 FREE: 2^{(n+1)/2} of them, all of them."""
    h = (n - 1) // 2
    out = []
    for s0 in (1, -1):
        for bits in itertools.product((1, -1), repeat=h):
            b = [0] * n
            b[0] = s0
            for k in range(1, h + 1):
                b[k] = bits[k - 1]
                b[n - k] = bits[k - 1]
            out.append(tuple(b))
    return out


def paf_half(x):
    """PAF_x(t) for t = 1 .. (n-1)/2; PAF_x(-t) = PAF_x(t) covers the rest."""
    n = len(x)
    return tuple(sum(x[j] * x[(j + t) % n] for j in range(n))
                 for t in range(1, (n - 1) // 2 + 1))


def paf_full(x):
    n = len(x)
    return tuple(sum(x[j] * x[(j + t) % n] for j in range(n))
                 for t in range(n))


def mitm(n, first_seeds, second_seeds):
    """Every quadruple (x0,x1,x2,x3) with x0,x1 from first_seeds (in that
    order) and x2,x3 from second_seeds whose aggregate PAF vanishes off
    zero, by meet in the middle on the (n-1)/2 independent lags."""
    h = (n - 1) // 2
    pf = {x: paf_half(x) for x in set(first_seeds) | set(second_seeds)}
    tail = {}
    for c in second_seeds:
        pc = pf[c]
        for d in second_seeds:
            pd = pf[d]
            key = tuple(pc[t] + pd[t] for t in range(h))
            tail.setdefault(key, []).append((c, d))
    out = []
    for a in first_seeds:
        pa = pf[a]
        for b in second_seeds:
            pb = pf[b]
            key = tuple(-(pa[t] + pb[t]) for t in range(h))
            for cd in tail.get(key, ()):
                out.append((a, b) + cd)
    return out


def aggregate_vanishes(q):
    n = len(q[0])
    agg = [0] * n
    for x in q:
        p = paf_full(x)
        for t in range(n):
            agg[t] += p[t]
    return agg[0] == 4 * n and all(v == 0 for v in agg[1:])


# ======================================================================
# [A] the enumeration, the identity, and the trust chain
# ======================================================================

def identity_holds(q, prefactor):
    a, b, c, d = q
    n = len(a)
    return all(a[(2 * k) % n] == prefactor(q) * a[k] * b[k] * c[k] * d[k]
               for k in range(1, n))


PRE_CORRECTED = lambda q: -(q[1][0] * q[2][0] * q[3][0])
PRE_MINUS = lambda q: -1
PRE_PLUS = lambda q: 1


def gs_array(sig, n, krho):
    """note/NOTE-B.md S1.0, standard orientation over Z_n."""
    def dev(f):
        return [[f[(h - g) % n] for h in range(n)] for g in range(n)]

    def xr(f):
        return [[f[(krho - g - h) % n] for h in range(n)] for g in range(n)]

    def xtr(f):
        return [[f[(g + h - krho) % n] for h in range(n)] for g in range(n)]

    def ng(X):
        return [[-v for v in r] for r in X]
    A = dev(sig[0])
    BR, CR, DR = xr(sig[1]), xr(sig[2]), xr(sig[3])
    BtR, CtR, DtR = xtr(sig[1]), xtr(sig[2]), xtr(sig[3])
    blk = [[A, BR, CR, DR],
           [ng(BR), A, DtR, ng(CtR)],
           [ng(CR), ng(DtR), A, BtR],
           [ng(DR), CtR, ng(BtR), A]]
    return [[blk[I][J][g][h] for J in range(4) for h in range(n)]
            for I in range(4) for g in range(n)]


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
    os.unlink(path)
    return proc.returncode, sha


def clause_a():
    log("[A] the complete small-order enumeration, and the identity")
    goods = {}
    for n in ORDERS:
        t0 = time.time()
        G = mitm(n, skew_seeds(n), sym_seeds(n))
        goods[n] = G
        check("n = %2d: %d good quadruples by complete meet-in-the-middle "
              "over all %d skew and all %d symmetric seeds, b_0 free (%.1fs)"
              % (n, len(G), 2 ** ((n - 1) // 2), 2 ** ((n + 1) // 2),
                 time.time() - t0),
              len(G) == EXPECTED_GOOD[n], len(G))
        nc = sum(1 for q in G if identity_holds(q, PRE_CORRECTED))
        check("n = %2d: a_{2k} = -(b_0 c_0 d_0) a_k b_k c_k d_k holds in "
              "%d / %d" % (n, nc, len(G)), nc == len(G))
        nm = sum(1 for q in G if identity_holds(q, PRE_MINUS))
        npl = sum(1 for q in G if identity_holds(q, PRE_PLUS))
        half = len(G) // 2
        check("n = %2d: the unprefixed variants hold in %d and %d of %d -- "
              "exactly half each, on the half where b_0 c_0 d_0 matches, so "
              "the prefactor is load-bearing" % (n, nm, npl, len(G)),
              nm == npl == half)
        # a normalised control: under b0 = c0 = d0 = +1, the "+" form is empty
        norm = [q for q in G if q[1][0] == q[2][0] == q[3][0] == 1]
        check("n = %2d: under b_0 = c_0 = d_0 = +1 there are %d quadruples; "
              "the '-' form holds in all of them and the '+' form in none"
              % (n, len(norm)),
              len(norm) > 0
              and all(identity_holds(q, PRE_MINUS) for q in norm)
              and not any(identity_holds(q, PRE_PLUS) for q in norm))
        check("n = %2d: every enumerated quadruple has aggregate PAF "
              "4n at 0 and 0 off it (re-checked lag by lag on a sample)"
              % n, all(aggregate_vanishes(q) for q in G[:40]))
    tot = sum(len(goods[n]) for n in ORDERS)
    check("9 936 good quadruples in all, and the corrected identity holds in "
          "every one", tot == 9936)

    log("[A'] the trust chain: assembled arrays are Hadamard")
    rng = random.Random(20260905)
    for n in ORDERS:
        q = rng.choice(goods[n])
        H = gs_array(list(q), n, 0)
        rc, sha = write_and_verify(H, "good-n%d" % n)
        check("n = %2d: GS(A,B,C,D; rho = 0) is Hadamard of order %d -- "
              "verify/verify.py exit %d, canonical sha256 %s..."
              % (n, 4 * n, rc, sha[:16]), rc == 0)
    return goods


# ======================================================================
# [B] the group-ring identity, instantiated
# ======================================================================

def indicator(x):
    """the {0,1} indicator of the -1 positions, as an F_2 coefficient list."""
    return [0 if v == 1 else 1 for v in x]


def sq2(t, n):
    """the squaring map on F_2[Z_n], n odd: (T^[2])_{2a} = t_a."""
    out = [0] * n
    for a in range(n):
        out[(2 * a) % n] = t[a]
    return out


def groupring_identity(q):
    """T^[2] + U = e + eps J over F_2[Z_n]; returns True iff it holds."""
    a, b, c, d = q
    n = len(a)
    U = indicator(a)
    V, W, Z = indicator(b), indicator(c), indicator(d)
    T = [(U[j] + V[j] + W[j] + Z[j]) % 2 for j in range(n)]
    lhs = [(sq2(T, n)[j] + U[j]) % 2 for j in range(n)]
    eps = (1 + sum(V) + sum(W) + sum(Z)) % 2
    rhs = [((1 if j == 0 else 0) + eps) % 2 for j in range(n)]
    return lhs == rhs


def clause_b(goods):
    log("[B] the group-ring identity T^[2] + U = e + eps J over F_2[Z_n], "
        "coefficient by coefficient")
    for n in ORDERS:
        G = goods[n]
        k = sum(1 for q in G if groupring_identity(q))
        check("n = %2d: holds in %d / %d good quadruples" % (n, k, len(G)),
              k == len(G))
    # the control: on quadruples that are NOT good, the identity must fail
    rng = random.Random(4242)
    for n in ORDERS:
        sk, sy = skew_seeds(n), sym_seeds(n)
        bad = []
        while len(bad) < 400:
            q = (rng.choice(sk), rng.choice(sy), rng.choice(sy),
                 rng.choice(sy))
            if not aggregate_vanishes(q):
                bad.append(q)
        held = sum(1 for q in bad if groupring_identity(q))
        check("n = %2d CONTROL: on 400 quadruples that are NOT good the "
              "identity fails %d / 400 times -- it is a consequence of "
              "goodness, not of the seed shapes" % (n, 400 - held),
              held < 400, "%d held" % held)


# ======================================================================
# [C] the Williamson control
# ======================================================================

def clause_c():
    log("[C] the Williamson control: four SYMMETRIC seeds")
    tot = 0
    for n in ORDERS:
        t0 = time.time()
        S = sym_seeds(n)
        W = mitm(n, S, S)
        tot += len(W)
        check("n = %2d: %d Williamson quadruples, complete over all %d "
              "symmetric seeds per slot, b_0 free (%.1fs)"
              % (n, len(W), len(S), time.time() - t0),
              len(W) == EXPECTED_WILLIAMSON[n], len(W))
        const = sum(1 for (a, b, c, d) in W
                    if all(a[k] * b[k] * c[k] * d[k]
                           == -(a[0] * b[0] * c[0] * d[0])
                           for k in range(1, n)))
        check("n = %2d: the CONSTANT-PRODUCT corollary a_k b_k c_k d_k = "
              "-(a_0 b_0 c_0 d_0) holds in %d / %d" % (n, const, len(W)),
              const == len(W))
        doub = sum(1 for q in W if identity_holds(q, PRE_CORRECTED))
        check("n = %2d: the DOUBLING relation holds in %d / %d -- it must "
              "not, because A is not skew here" % (n, doub, len(W)),
              doub == 0)
        unpref = sum(1 for q in W if identity_holds(q, PRE_MINUS))
        want = 48 if n == 7 else 0
        check("n = %2d: the unprefixed doubling variant holds in %d "
              "(a small-n coincidence at n = 7, nothing elsewhere)"
              % (n, unpref), unpref == want, unpref)
    check("14 784 Williamson quadruples in all; the constant-product "
          "corollary holds in every one and the doubling relation in none",
          tot == 14784, tot)


# ======================================================================
# [D] the closure sign
# ======================================================================

def multiplier_data(n, m):
    """M_0 = the subgroup of Z_n^* of order m (assumed unique, checked);
    returns (M_0, orbits of <M_0,-1> on Z_n \\ {0}, L, 2^L mod n)."""
    units = [u for u in range(1, n) if _gcd(u, n) == 1]
    M0 = sorted(u for u in units if pow(u, m, n) == 1)
    if len(M0) != m:
        return None
    full = sorted(set(M0) | set((-u) % n for u in M0))
    seen, orbits = set(), []
    for k in range(1, n):
        if k in seen:
            continue
        o = frozenset((k * u) % n for u in full)
        seen |= o
        orbits.append(o)
    L = 1
    p = 2 % n
    while p not in full:
        p = (p * 2) % n
        L += 1
        if L > n:
            return None
    return M0, orbits, L, p


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def clause_d(goods):
    log("[D] the closure sign, computed from the group and then checked "
        "against the enumeration")
    # n = 7, M_0 trivial
    dat = multiplier_data(7, 1)
    check("n = 7, M_0 = {1}: <M_0,-1> has %d elements and %d orbits on "
          "Z_7 \\ {0}; 2 reaches <M_0,-1> in L = %d steps, at 2^L = %d"
          % (2, len(dat[1]), dat[2], dat[3]),
          dat is not None and len(dat[1]) == 3 and dat[2] == 3
          and dat[3] == 1)
    M0, orbits, L, pL = dat
    closure = 1 if pL in M0 else -1
    check("2^3 = 1 lies in M_0, so closure = +1 and the cycle condition is "
          "prod_j pi_j = closure * (-(b_0 c_0 d_0))^L = -(b_0 c_0 d_0)",
          closure == 1)
    G = goods[7]
    reps = [1, 2, 4]
    okc = 0
    for q in G:
        a, b, c, d = q
        pre = -(b[0] * c[0] * d[0])
        prod = 1
        for j in reps:
            prod *= b[j] * c[j] * d[j]
        if prod == closure * pre ** L:
            okc += 1
    check("and that prediction holds in %d / %d enumerated good quadruples "
          "at n = 7 -- a control that could fail" % (okc, len(G)),
          okc == len(G))
    check("under the normalisation b_0 = c_0 = d_0 = +1 it reads "
          "prod_j pi_j = -1",
          all(b[1] * c[1] * d[1] * b[2] * c[2] * d[2] * b[4] * c[4] * d[4]
              == -1
              for (a, b, c, d) in G if b[0] == c[0] == d[0] == 1))

    # n = 523, M_0 of order 3
    n = 523
    check("523 is prime and |Z_523^*| = 522 = 2 * 3^2 * 29",
          all(523 % p for p in range(2, 24)) and 522 == 2 * 9 * 29)
    dat = multiplier_data(n, 3)
    check("Z_523^* has exactly one subgroup of order 3, M_0 = %s"
          % (dat[0] if dat else None), dat is not None and len(dat[0]) == 3)
    M0, orbits, L, pL = dat
    check("<M_0,-1> has order 6 and %d orbits on Z_523 \\ {0}, all of size 6"
          % len(orbits),
          len(orbits) == 87 and all(len(o) == 6 for o in orbits))
    check("2 has order exactly L = %d in Z_523^* / <M_0,-1>, so doubling is a "
          "single %d-cycle on those orbits" % (L, L), L == 87)
    check("2^87 = %d mod 523, and -%d = %d lies in M_0, so 2^87 lies in -M_0"
          % (pL, pL, (-pL) % n), pL == 463 and ((-pL) % n) in M0)
    closure523 = 1 if pL in M0 else -1
    check("hence closure = -1 at 523, and the cycle condition reads "
          "prod_j pi_j = closure * (-(b_0 c_0 d_0))^87 = +1 under the "
          "normalisation -- ONE parity condition, not 87", closure523 == -1
          and closure523 * (-1) ** L == 1)
    check("NOT claimed: that any such seed quadruple exists at n = 523; the "
          "closure condition is necessary, and its hypothesis is unmet",
          True)


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args(argv)

    print("=" * 78)
    print("cert 32 -- the good-matrix product theorem: controls, the "
          "group-ring identity,")
    print("           the Williamson control, and the closure sign")
    print("=" * 78)
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    try:
        goods = clause_a()
        clause_b(goods)
        clause_c()
        clause_d(goods)
    finally:
        shutil.rmtree(OUT, ignore_errors=True)

    print("\n" + "=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 32: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: the good-matrix product theorem's corrected form")
    print("         a_{2k} = -(b_0 c_0 d_0) a_k b_k c_k d_k holds in ALL")
    print("         9 936 good quadruples at n = 7, 9, 11, 13, 15, found by a")
    print("         complete meet-in-the-middle with b_0 free; each")
    print("         unprefixed variant holds in exactly half, so the")
    print("         prefactor is load-bearing.  The group-ring identity the")
    print("         proof turns on holds on every one of them and fails on")
    print("         quadruples that are not good.  With four SYMMETRIC seeds")
    print("         the constant-product corollary holds in all 14 784")
    print("         Williamson quadruples and the doubling relation in NONE:")
    print("         skewness of A is load-bearing too.  The closure sign is")
    print("         computed from the group at n = 7 and at n = 523, and the")
    print("         n = 7 prediction is checked against the enumeration.")
    print("         LABEL: COMPUTATIONAL-EVIDENCE beside the PROVEN theorem")
    print("         of note/NOTE-B.md S1.13; the enumerations themselves are")
    print("         PROVEN-BY-CERTIFICATE (complete, exact, no prefilter).")
    print("         NOT claimed: any existence statement at n = 523, where no")
    print("         seed quadruple is known and the closure condition is")
    print("         necessary only; anything at n > 15.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
