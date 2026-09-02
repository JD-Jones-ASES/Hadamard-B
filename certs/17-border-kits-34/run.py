#!/usr/bin/env python3
"""cert 17 -- the (3,4) border-kit census, and Theorem F on the banked records.

  THEOREM F (note/NOTE-B.md S1.8) fixes the structure of a border kit at
  i = s+1 under w > 2s: Q^T Q = 4i I_{4s}; P (I4 (x) P_{S^c}) = 0 and
  P^T P = Q Q^T = I4 (x) M; Chat^T Chat = Chat Chat^T = I4 (x) Sigma-bar;
  E = -(1/4i) P Chat^T Q; and (H3)+(H4) reduce to the single condition
  E in {+-1}.  (H4) sees the seeds only through the S-PART of the
  coset-sum table, so a kit is an ORDER-INDEPENDENT object attached to
  (s, i, Gbar, kappa(rho), S-part).

  THE (3,4) CENSUS (note/NOTE-B.md S1.8).  At (s,i) = (3,4) with the
  house Gram M = 16 I - 4 J and Gbar in {Z4, Z2xZ2}: (i) exactly one seed
  is coset-balanced and the other three carry |sigma-hat_q(chi)| = 2 at
  every nontrivial character -- the SILENT-SEED LEMMA; (ii) there are
  exactly 2048 admissible S-parts per group; (iii) EVERY S-part admits a
  kit at EVERY kappa(rho); (iv) at N = 2092 the coset-sum tables form a
  shell of 4192 x 64 = 2048 x 131 = 268 288 tables per group.

  COROLLARY (the (3,4) cell is one-layer).  Every house-profile seed
  quadruple on an abelian group of order 4w, w > 6, with an index-4
  subgroup extends to a Hadamard matrix of order 16w + 12; the border is
  never the obstruction.  At N = 2092 (w = 130): a seed quadruple IS
  H(2092).

Theorem F is a PAPER proof.  This certificate carries the finite parts:
the S-part census, the shell arithmetic, the exhaustive kit census, the
from-scratch controls through verify/verify.py, and Theorem F measured
literally on the seven banked coset-border records with s >= 1, w > 2s.

WHAT THIS SCRIPT DOES  (default path: standard library only, seconds)

  [A] S-PART CENSUS AND THE SILENT-SEED LEMMA.  All spectra at the three
      nontrivial characters with sum_q |sigma-hat_q(chi)|^2 = 4s = 12
      (Corollary E2) and entries in 2Z are enumerated by BRUTE FORCE, and
      the ones for which sigma = (r + T)/4 is integral with the right
      parity for some integer row-sum vector are the S-parts.  The
      silent-seed shape is then checked on every one of them.  The
      sigma-shell at N = 2092 -- even r with sum r_q^2 = 2092 -- is
      enumerated and cross-multiplied against the S-parts.

  [B] KIT CENSUS.  For every S-part and EVERY kappa(rho) -- all four per
      group, no symmetry reduction -- kitlib.py finds (E, P, Q), and each
      kit is re-verified exactly with a FULL sigma (r-part included) at
      w = 130: (H1) Q Q^T = I4 (x) M, (H3) E E^T + w P P^T = N I, (H4)
      E Q^T + P Chat^T = 0, E in H(12), P P^T = 16 I.  Default: a fixed
      deterministic sample (every 64th S-part, both groups, all four
      kappa(rho): 256 classes).  `--full`: all 2 x 4 x 2048 = 16384
      classes.  Both census digests are pinned in pins.json.

  [C] CONTROLS THROUGH verify/verify.py.  From-scratch (3,4) instances at
      n = 8 (w = 2): the house-profile seed quadruples are found by a full
      meet-in-the-middle over all 256 sequences on Z8 (two rho) and on
      Z2xZ4 (both quotients), the kits by this engine, the matrices by
      this certificate's own assembler -- four H(44), all green; Z2^3
      with |K| = 2 carries no house quadruple, and that is checked too.
      Then H(124) on Z2xZ2xZ7 (w = 7) twice: once with the order-1916
      record's kit transported verbatim (Theorem F(g)), once with this
      engine's own kit -- two different artifacts, both green.

  [D] THEOREM F ON THE BANKED RECORDS.  For every record in
      data/payload-records.json with s >= 1 and w > 2s -- 668, 716, 1676,
      1772 at (1,1); 1916 at (3,4); 1388 at (5,6); 1436 at (7,8) -- the
      identities of Theorem F (a)-(f) are re-multiplied entry by entry,
      and at i = s+1 the (H4)-kernel is exhibited by shifting sigma by
      per-seed constants and requiring (H4) to be unchanged.

Usage:
  python certs/17-border-kits-34/run.py
  python certs/17-border-kits-34/run.py --full
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

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

import kitlib as K                                            # noqa: E402

VERIFY = os.path.join(ROOT, "verify", "verify.py")
RECORDS = os.path.join(ROOT, "data", "payload-records.json")
OUT = os.path.join(HERE, "out")

_T0 = time.time()
FAIL = []
NCHECK = [0]
W2092 = 130
PINS = {}

# ---------------------------------------------------------------- the pins
# SHA-256 of the banked record file this certificate reads.
FILE_PINS = {
    "data/payload-records.json":
        "9afbf60efed82e97ffe229e0059f06dfa7849a2ac2434c52a6ada09559738afb",
}

# The H(124) parameters used by clause [C].  This is the order-1916 record's
# border kit (data/payload-records.json, order 1916) transported by
# Theorem F(g) onto G = Z2 x Z2 x Z7 -- same Gbar = Z2xZ2, same kappa(rho),
# same S-part, a different w.  It is embedded here rather than banked in
# data/ because it is not a public record: it is a control this certificate
# constructs, and the transport is exactly what clause [C] is testing.  The
# seeds are a house-profile quadruple on Z2 x Z2 x Z7 (n = 28, w = 7 > 2s).
H124 = dict(
    order=124, group=[2, 2, 7], s=3, coset_divisors=[2, 2, 1],
    r_shift=[0, 0, 0],
    seeds=["++-++----+--++--+--++--+--++",
           "+-+------+--+---+-+---+-+--+",
           "+++-------++++----+++----+++",
           "+++-+--+-+-+-++-+-+-++--+-++"],
    corner=["++--+++--++-", "----+-++-+-+", "-+---+-+++--", "-+--------++",
            "++-+--+-++-+", "---++---+++-", "-++-+-+-+---", "+-+------+--",
            "++-++--+----", "---+-++-----", "-+++++---+-+", "+---++--+--+"],
    row_table=["+-+-+-+---++-++-", "+-+--+-+--+++--+", "+-+-++--++----++",
               "+-+---++++--++--", "++--+-+-+--++--+", "++---+-++--+-++-",
               "++--++---++-++--", "++----++-++---++", "+--++--++-+-+-+-",
               "+--++--+-+-+-+-+", "+--+-++-+-+--+-+", "+--+-++--+-++-+-"],
    col_table=["+-+--+-++-+-++--", "+-+--+-+-+-+--++", "+-+-+-+-+--+-++-",
               "+-+-+-+--++-+--+", "++--+--+++---+-+", "++--+--+--+++-+-",
               "++---++-++--+-+-", "++---++---++-+-+", "+--+--+++-+---++",
               "+--+--++-+-+++--", "+--+++--+--++--+", "+--+++---++--++-"],
)


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
    return [1 if ch == "+" else -1 for ch in s]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


# ======================================================================
# clause A -- the S-part census, the silent-seed lemma, the shell
# ======================================================================

def sparts(gname):
    """Every admissible S-part, by brute force over the spectra.

    Returns a list of (T, spec) with T = 4*sigma0 as a 4 x 4 integer table
    (rows = seeds, columns = the group's elements) and spec the spectrum
    magnitudes per seed at the three nontrivial characters.
    """
    els = K.GROUPS[gname]["elts"]
    out = []
    if gname == "Z2xZ2":
        ch = K.chars(gname)
        vals = (0, 2, -2)
        for Y in itertools.product(itertools.product(vals, repeat=4), repeat=3):
            # Y[t][q]: character t, seed q; Corollary E2 per character
            if any(sum(v * v for v in Y[t]) != 12 for t in range(3)):
                continue
            T = [[sum(Y[t][q] * int(ch[t](c)) for t in range(3)) for c in els]
                 for q in range(4)]
            if all(len({v % 8 for v in T[q]}) == 1 for q in range(4)):
                spec = [[abs(Y[t][q]) for t in range(3)] for q in range(4)]
                out.append((T, spec))
    else:
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
                    z2 = [ab[q][0] ** 2 + ab[q][1] ** 2 for q in range(4)]
                    spec = [[abs(y[q]), z2[q], z2[q]] for q in range(4)]
                    out.append((T, spec))
    return out


def shell_2092():
    """Even r = (r_0,...,r_3) with sum r_q^2 = 2092, ordered, signs kept."""
    out = []
    rng = range(-44, 45, 2)
    sq = {r: r * r for r in rng}
    for a in rng:
        for b in rng:
            ab = sq[a] + sq[b]
            if ab > 2092:
                continue
            for c in rng:
                abc = ab + sq[c]
                if abc > 2092:
                    continue
                rem = 2092 - abc
                d = int(round(rem ** 0.5))
                for dd in (d - 1, d, d + 1):
                    if dd * dd == rem and dd % 2 == 0:
                        out.append((a, b, c, dd))
                        if dd != 0:
                            out.append((a, b, c, -dd))
    return sorted(set(out))


def clause_a():
    log("[A] the S-part census, the silent-seed lemma and the sigma-shell")
    shell = shell_2092()
    check("sigma-shell at N = 2092: %d ordered even r with sum r_q^2 = 2092 "
          "(= r_4(523), Jacobi)" % len(shell), len(shell) == 4192)
    check("every even r in the shell has exactly one coordinate == 0 (mod 8) "
          "and three == 2 (mod 4) (forced by 523 == 3 mod 8)",
          all(sum(1 for v in r if v % 8 == 0) == 1
              and sum(1 for v in r if v % 4 == 2) == 3 for r in shell))
    SP = {}
    for gname in ("Z2xZ2", "Z4"):
        t0 = time.time()
        sp = sparts(gname)
        SP[gname] = sp
        check("%-6s %d Parseval-admissible S-parts (brute force over the "
              "spectra, %.1fs)" % (gname + ":", len(sp), time.time() - t0),
              len(sp) == 2048)
        silent_ok = True
        for T, spec in sp:
            zero_seeds = [q for q in range(4) if all(v == 0 for v in spec[q])]
            live_ok = all(
                all(v == (2 if gname == "Z2xZ2" else (2 if t == 0 else 4))
                    for t, v in enumerate(spec[q]))
                for q in range(4) if q not in zero_seeds)
            if len(zero_seeds) != 1 or not live_ok:
                silent_ok = False
                break
        check("%-6s silent-seed lemma: exactly one seed has zero spectrum; "
              "the other three have |sigma-hat_q(chi)| = 2 at every "
              "nontrivial character" % (gname + ":"), silent_ok)
        taus = [tuple((-T[q][0]) % 8 for q in range(4)) for T, _ in sp]
        by_res = {}
        for r in shell:
            key = tuple(v % 8 for v in r)
            by_res[key] = by_res.get(key, 0) + 1
        per_spart = [by_res.get(t, 0) for t in taus]
        check("%-6s every S-part is compatible with exactly 131 shell "
              "vectors (total %d = 268,288 coset-sum tables)"
              % (gname + ":", sum(per_spart)),
              all(v == 131 for v in per_spart) and sum(per_spart) == 268288)
        tau_count = {}
        for t in taus:
            tau_count[t] = tau_count.get(t, 0) + 1
        per_r = [tau_count.get(tuple(v % 8 for v in r), 0) for r in shell]
        check("%-6s every shell vector is compatible with exactly 64 S-parts"
              % (gname + ":"), all(v == 64 for v in per_r))
        check("%-6s the silent seed is the seed with r_q == 0 (mod 8)"
              % (gname + ":"),
              all((taus[k][q] == 0) == all(v == 0 for v in sp[k][1][q])
                  for k in range(len(sp)) for q in range(4)))
    return SP, shell


# ======================================================================
# clause B -- the kit census
# ======================================================================

def full_sigma(T, shell_by_res, w):
    """A compatible full sigma at width w: the first shell r with
    r == -tau (mod 8)."""
    tau = tuple((-T[q][0]) % 8 for q in range(4))
    r = shell_by_res[tau][0]
    return [[(r[q] + T[q][c]) // 4 for c in range(4)] for q in range(4)], r


def clause_b(SP, shell, full):
    log("[B] the kit census (%s)"
        % ("--full: all 16384 classes" if full else
           "default: a fixed deterministic sample"))
    rng = random.Random(3403)
    QL = K.default_Q_list(rng, extra=300)
    shell_by_res = {}
    for r in shell:
        shell_by_res.setdefault(tuple(v % 8 for v in r), []).append(r)
    digest = hashlib.sha256()
    found = total = 0
    tries_hist = {}
    maxtries = 0
    bad = []
    for gname in ("Z2xZ2", "Z4"):
        sp = SP[gname]
        idxs = range(len(sp)) if full else range(7, len(sp), 64)
        for krho in K.GROUPS[gname]["elts"]:
            f0 = found
            t0 = time.time()
            for k in idxs:
                T, _ = sp[k]
                total += 1
                E, P, Q, tries = K.kit_for(T, gname, krho, QL)
                ok = E is not None
                if ok:
                    sigma, r = full_sigma(T, shell_by_res, W2092)
                    v = K.verify_kit(E, P, Q, sigma, gname, krho, W2092)
                    ok = all(v.values())
                    if not ok:
                        bad.append((gname, krho, k, v))
                    tries_hist[tries] = tries_hist.get(tries, 0) + 1
                    maxtries = max(maxtries, tries)
                else:
                    bad.append((gname, krho, k, "no kit"))
                found += ok
                digest.update(("%s|%s|%s|%d\n"
                               % (gname, krho, T, ok)).encode("ascii"))
            print("      %-6s krho=%-6s %d / %d classes admit a kit (%.1fs)"
                  % (gname, krho, found - f0, len(idxs), time.time() - t0),
                  flush=True)
    check("every class in the sweep admits a kit, re-verified exactly at "
          "w = 130 with a full sigma: %d / %d" % (found, total),
          found == total and not bad, str(bad[:3]))
    print("      Q-candidates tried before success, histogram: %s (max %d "
          "of the 301 in the fixed list)"
          % (dict(sorted(tries_hist.items())), maxtries))
    hx = digest.hexdigest()
    print("      census digest: %s" % hx)
    if full:
        check("full-census digest == the pin: %s..." % PINS["census_full"][:16],
              hx == PINS["census_full"], hx)
    else:
        check("sample-census digest == the pin: %s..."
              % PINS["census_sample"][:16], hx == PINS["census_sample"], hx)
    return hx


# ======================================================================
# this certificate's own assembler
# ======================================================================

class Grp:
    """Mixed-radix abelian group, tape convention (row-major)."""

    def __init__(self, factors):
        self.f = list(factors)
        self.elts = list(itertools.product(*[range(m) for m in self.f]))
        self.n = len(self.elts)
        self.idx = {e: k for k, e in enumerate(self.elts)}

    def add(self, a, b):
        return tuple((x + y) % m for x, y, m in zip(a, b, self.f))

    def neg(self, a):
        return tuple((-x) % m for x, m in zip(a, self.f))

    def sub(self, a, b):
        return self.add(a, self.neg(b))


def kappa_of(G, divisors):
    def kap(g):
        return tuple(x % d for x, d in zip(g, divisors))
    return kap, Grp(divisors)


def gs_entry(seeds, G, rho, I, J, g, h):
    """Entry (I,g),(J,h) of the standard-orientation GS array
    (note/NOTE-B.md S1.0)."""
    x = seeds
    if I == J:
        return x[0][G.idx[G.sub(h, g)]]
    table = {
        (0, 1): (1, "R", 1), (0, 2): (2, "R", 1), (0, 3): (3, "R", 1),
        (1, 0): (1, "R", -1), (1, 2): (3, "TR", 1), (1, 3): (2, "TR", -1),
        (2, 0): (2, "R", -1), (2, 1): (3, "TR", -1), (2, 3): (1, "TR", 1),
        (3, 0): (3, "R", -1), (3, 1): (2, "TR", 1), (3, 2): (1, "TR", -1),
    }
    q, form, sgn = table[(I, J)]
    if form == "R":
        return sgn * x[q][G.idx[G.sub(rho, G.add(g, h))]]
    return sgn * x[q][G.idx[G.sub(G.add(g, h), rho)]]


def assemble(G, seeds, rho, s, kap, Gq, E, P, Q):
    """The bordered matrix, as +/- strings.  Border strips constant on the
    K-cosets inside each superblock; column table Q indexed
    [i*I + kappa(g)][c]."""
    n, i = G.n, Gq.n
    N = 4 * (n + s)
    ch = {1: "+", -1: "-"}
    rows = []
    for r in range(4 * s):
        line = [ch[v] for v in E[r]]
        for J in range(4):
            line += [ch[P[r][i * J + Gq.idx[kap(h)]]] for h in G.elts]
        rows.append("".join(line))
    for I in range(4):
        for g in G.elts:
            line = [ch[Q[i * I + Gq.idx[kap(g)]][c]] for c in range(4 * s)]
            for J in range(4):
                line += [ch[gs_entry(seeds, G, rho, I, J, g, h)]
                         for h in G.elts]
            rows.append("".join(line))
    assert len(rows) == N and all(len(r) == N for r in rows)
    return rows


def write_and_verify(rows, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name + ".txt")
    text = "\n".join(rows) + "\n"
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write(text)
    sha = hashlib.sha256(text.encode("ascii")).hexdigest()
    proc = subprocess.run([sys.executable, VERIFY, path],
                          capture_output=True, text=True,
                          env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    verdict = [l for l in (proc.stdout + proc.stderr).splitlines()
               if l.startswith("VERDICT")]
    verdict = verdict[-1] if verdict else "(no verdict)"
    print("      verify.py %-26s exit %d :: %s"
          % (name, proc.returncode, verdict[:96]))
    os.unlink(path)
    return proc.returncode, verdict, sha


def general_gs(Gq, sig, rho_bar):
    """The GS array over the quotient group of the integer functions sig."""
    return [[gs_entry(sig, Gq, rho_bar, I, J, g, h)
             for J in range(4) for h in Gq.elts]
            for I in range(4) for g in Gq.elts]


# ======================================================================
# clause C -- controls through verify/verify.py
# ======================================================================

def paf(G, x, t):
    return sum(x[G.idx[u]] * x[G.idx[G.add(u, t)]] for u in G.elts)


def to_kitlib_group(Gq):
    """Identify the order-4 quotient with kitlib's Z4 or Z2xZ2 labelling."""
    orders = [next(k for k in range(1, 5)
                   if all(((k * x) % m) == 0 for x, m in zip(e, Gq.f)))
              for e in Gq.elts]
    if 4 in orders:
        gen = Gq.elts[orders.index(4)]
        lab = {Gq.elts[0]: 0, gen: 1, Gq.add(gen, gen): 2,
               Gq.add(Gq.add(gen, gen), gen): 3}
        return "Z4", lab
    nz = [e for e in Gq.elts if e != Gq.elts[0]]
    lab = {Gq.elts[0]: (0, 0), nz[0]: (0, 1), nz[1]: (1, 0), nz[2]: (1, 1)}
    return "Z2xZ2", lab


def clause_c(QL):
    log("[C] controls: from-scratch H(44) at n = 8, and H(124) two ways")
    cases = [([8], [4], (0,)), ([8], [4], (7,)),
             ([2, 4], [1, 4], (0, 0)), ([2, 4], [2, 2], (0, 0)),
             ([2, 2, 2], [1, 2, 2], (0, 0, 0))]
    for factors, divisors, rho in cases:
        G = Grp(factors)
        kap, Gq = kappa_of(G, divisors)
        n = G.n
        Kset = [g for g in G.elts if kap(g) == Gq.elts[0]]
        target = tuple(4 * n if t == G.elts[0]
                       else (-12 if t in Kset else 4) for t in G.elts)
        seqs = list(itertools.product((1, -1), repeat=n))
        PF = [tuple(paf(G, x, t) for t in G.elts) for x in seqs]
        table = {}
        for a in range(len(seqs)):
            for b in range(a, len(seqs)):
                table.setdefault(tuple(u + v for u, v in zip(PF[a], PF[b])),
                                 []).append((a, b))
        sol = None
        nsol = 0
        for a in range(len(seqs)):
            for b in range(a, len(seqs)):
                need = tuple(t - u - v
                             for t, u, v in zip(target, PF[a], PF[b]))
                lst = table.get(need)
                if lst:
                    nsol += len(lst)
                    if sol is None:
                        sol = (a, b) + lst[0]
        tag = "G=%s K=ker(%s)" % ("x".join("Z%d" % m for m in factors),
                                  divisors)
        if factors == [2, 2, 2]:
            check("%s: NO house-profile quadruple at n = 8 (full MITM over "
                  "all 256 sequences)" % tag, sol is None)
            continue
        check("%s: house-profile quadruples exist at n = 8 (full MITM, %d "
              "ordered pair-pairs)" % (tag, nsol), sol is not None)
        seeds = [list(seqs[q]) for q in sol]
        gname, lab = to_kitlib_group(Gq)
        kl = {lab[e]: e for e in Gq.elts}
        els4 = K.GROUPS[gname]["elts"]
        sigma = [[sum(seeds[q][G.idx[g]] for g in G.elts if kap(g) == kl[c])
                  for c in els4] for q in range(4)]
        r = [sum(x) for x in seeds]
        T = [[4 * sigma[q][c] - r[q] for c in range(4)] for q in range(4)]
        krho = lab[kap(rho)]
        E, P, Q, tries = K.kit_for(T, gname, krho, QL)
        check("%s rho=%s (Gbar=%s, krho=%s): this engine finds a kit (%d "
              "Q-candidates)" % (tag, rho, gname, krho, tries), E is not None)
        if E is None:
            continue
        v = K.verify_kit(E, P, Q, sigma, gname, krho, n // 4)
        check("  kit re-verified exactly at w = %d: %s" % (n // 4, v),
              all(v.values()))
        Qg = [[Q[4 * I + els4.index(lab[e])][c] for c in range(12)]
              for I in range(4) for e in Gq.elts]
        Pg = [[P[rr][4 * J + els4.index(lab[e])]
               for J in range(4) for e in Gq.elts] for rr in range(12)]
        rows = assemble(G, seeds, rho, 3, kap, Gq, E, Pg, Qg)
        rc, verdict, sha = write_and_verify(
            rows, "H44_%s_rho%s" % ("x".join(map(str, factors)),
                                    "".join(map(str, rho))))
        check("  the assembled H(44) is verify.py green",
              rc == 0 and "HADAMARD order=44" in verdict)
    # H(124), twice
    rec = H124
    G = Grp(rec["group"])
    kap, Gq = kappa_of(G, rec["coset_divisors"])
    seeds = [pm(x) for x in rec["seeds"]]
    rho = tuple(rec["r_shift"])
    E = [pm(x) for x in rec["corner"]]
    P = [pm(x) for x in rec["row_table"]]
    colT = [pm(x) for x in rec["col_table"]]
    Q = [[colT[c][k] for c in range(12)] for k in range(16)]
    Kset = [g for g in G.elts if kap(g) == Gq.elts[0]]
    prof = all(sum(paf(G, x, t) for x in seeds)
               == (4 * G.n if t == G.elts[0] else (-12 if t in Kset else 4))
               for t in G.elts)
    check("H(124) control: the house profile on Z2xZ2xZ7 (n = 28, w = 7 > "
          "2s = 6)", prof)
    rows = assemble(G, seeds, rho, 3, kap, Gq, E, P, Q)
    rc, verdict, sha = write_and_verify(rows, "H124_transported_kit")
    check("H(124) with the order-1916 kit transported verbatim (Theorem "
          "F(g)): verify.py green, digest == the pin %s..."
          % PINS["h124_lane"][:16],
          rc == 0 and "HADAMARD order=124" in verdict
          and sha == PINS["h124_lane"], sha)
    gname, lab = to_kitlib_group(Gq)
    kl = {lab[e]: e for e in Gq.elts}
    els4 = K.GROUPS[gname]["elts"]
    sigma = [[sum(seeds[q][G.idx[g]] for g in G.elts if kap(g) == kl[c])
              for c in els4] for q in range(4)]
    r = [sum(x) for x in seeds]
    T = [[4 * sigma[q][c] - r[q] for c in range(4)] for q in range(4)]
    krho = lab[kap(rho)]
    E2, P2, Q2, tries = K.kit_for(T, gname, krho, QL)
    check("H(124) S-part: this engine finds its OWN kit (%d Q-candidates)"
          % tries, E2 is not None)
    if E2 is not None:
        v = K.verify_kit(E2, P2, Q2, sigma, gname, krho, 7)
        check("  own kit re-verified exactly at w = 7: %s" % v,
              all(v.values()))
        Qg = [[Q2[4 * I + els4.index(lab[e])][c] for c in range(12)]
              for I in range(4) for e in Gq.elts]
        Pg = [[P2[rr][4 * J + els4.index(lab[e])]
               for J in range(4) for e in Gq.elts] for rr in range(12)]
        rows2 = assemble(G, seeds, rho, 3, kap, Gq, E2, Pg, Qg)
        rc2, verdict2, sha2 = write_and_verify(rows2, "H124_own_kit")
        check("  a second H(124) -- same seeds, this engine's border -- is "
              "verify.py green; digest %s..." % sha2[:16],
              rc2 == 0 and "HADAMARD order=124" in verdict2)
        check("  the two H(124) are different artifacts (different borders)",
              sha2 != sha)


# ======================================================================
# clause D -- Theorem F on the banked records
# ======================================================================

def mmT(A, B):
    return [[dot(a, b) for b in B] for a in A]


def mm(A, B):
    Bt = list(zip(*B))
    return [[dot(a, b) for b in Bt] for a in A]


def clause_d():
    log("[D] Theorem F, measured on the banked coset-border records")
    with open(RECORDS, encoding="ascii") as fh:
        recs = json.load(fh)["orders"]
    n_checked = 0
    for rec in recs:
        s = int(rec["s"])
        if s == 0:
            continue
        G = Grp(rec["group"])
        div = rec["coset_divisors"]
        kap, Gq = kappa_of(G, div)
        i = Gq.n
        w = G.n // i
        N = int(rec["order"])
        seeds = [pm(x) for x in rec["seeds"]]
        rho = tuple(rec["r_shift"])
        E = [pm(x) for x in rec["corner"]]
        P = [pm(x) for x in rec["row_table"]]
        colT = [pm(x) for x in rec["col_table"]]
        Q = [[colT[c][k] for c in range(4 * s)] for k in range(4 * i)]
        sigma = [[sum(seeds[q][G.idx[g]] for g in G.elts if kap(g) == c)
                  for c in Gq.elts] for q in range(4)]
        rb = kap(rho)
        Ch = general_gs(Gq, sigma, rb)
        M = [[(4 * s + 4 if a == b else 0) - 4 for b in range(i)]
             for a in range(i)]                                # house form
        I4M = [[M[a % i][b % i] if a // i == b // i else 0
                for b in range(4 * i)] for a in range(4 * i)]
        tag = ("N=%d (s,i,w)=(%d,%d,%d) Gbar=%s"
               % (N, s, i, w,
                  "x".join("Z%d" % d for d in div if d > 1) or "1"))
        if w <= 2 * s:
            print("      [info] %s: w <= 2s, Theorem F not applicable, "
                  "skipped" % tag)
            continue
        n_checked += 1
        QtQ = mm([list(x) for x in zip(*Q)], Q)
        ok_a = all(QtQ[a][b] == (4 * i if a == b else 0)
                   for a in range(4 * s) for b in range(4 * s))
        PtP = mm([list(x) for x in zip(*P)], P)
        QQt = mmT(Q, Q)
        ok_b = (PtP == I4M) and (QQt == I4M)
        CCt, CtC = mmT(Ch, Ch), mm([list(x) for x in zip(*Ch)], Ch)
        sbar0 = 4 * G.n - 4 * s * (w - 1)
        Sb = [[0] * (4 * i) for _ in range(4 * i)]
        for a in range(4 * i):
            for b in range(4 * i):
                if a // i == b // i:
                    Sb[a][b] = sbar0 if a == b else 4 * w
        ok_c = (CCt == Sb) and (CtC == Sb)
        PCQ = mm(mmT(P, Ch), Q)
        ok_e = all(PCQ[a][b] == -4 * i * E[a][b]
                   for a in range(4 * s) for b in range(4 * s))
        res1 = mmT(E, Q)
        PC1 = mmT(P, Ch)
        h4 = all(res1[a][b] + PC1[a][b] == 0
                 for a in range(4 * s) for b in range(4 * i))
        if i == s + 1:
            ok_bal = all(sum(P[r][i * J + c] for c in range(i)) == 0
                         for r in range(4 * s) for J in range(4))
            rng = random.Random(N)
            shifts = [rng.randrange(-7, 8) for _ in range(4)]
            sig2 = [[sigma[q][c] + shifts[q] for c in range(i)]
                    for q in range(4)]
            h4k = (PC1 == mmT(P, general_gs(Gq, sig2, rb)))
            extra = "P superblock-balanced | (H4)-kernel = per-seed constants"
        else:
            # (1,1): S = {1} is the whole dual, S^c empty, the kernel is
            # 0-dimensional and the row sums enter (H4) directly -- NOTE-B
            # Theorem D (D-e).  Theorem F(f) is vacuous there.
            ok_bal = h4k = True
            extra = "(i = s: S^c empty, kernel 0-dimensional, r direct)"
        check("%s: Q^T Q = 4i I | P^T P = Q Q^T = I4(x)M | Chat^T Chat = "
              "Chat Chat^T = I4(x)Sigma-bar | E = -(1/4i) P Chat^T Q | (H4) "
              "| %s" % (tag, extra),
              ok_a and ok_b and ok_c and ok_e and ok_bal and h4 and h4k,
              str([ok_a, ok_b, ok_c, ok_e, h4, ok_bal, h4k]))
    check("banked records with s >= 1 and w > 2s checked: %d -- the four "
          "(1,1) records, and one each at (3,4), (5,6), (7,8)" % n_checked,
          n_checked == 7)


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--full", action="store_true",
                    help="census all 2 x 4 x 2048 = 16384 (group, kappa(rho), "
                         "S-part) classes instead of the deterministic "
                         "sample; the full-census digest is pinned")
    args = ap.parse_args(argv)
    with open(os.path.join(HERE, "pins.json"), encoding="ascii") as fh:
        PINS.update(json.load(fh))

    print("=" * 78)
    print("cert 17 -- the (3,4) border-kit census, and Theorem F on the "
          "records")
    print("            trust chain: %s"
          % os.path.relpath(VERIFY, ROOT).replace("\\", "/"))
    print("=" * 78)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    try:
        rc = _body(args)
    finally:
        shutil.rmtree(OUT, ignore_errors=True)
    return rc


def _body(args):
    log("[0] the banked record file, SHA-256 pinned in this script")
    for name, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, name))
        check("%-32s" % name, got == want, got[:24] + "...")

    SP, shell = clause_a()
    clause_b(SP, shell, args.full)
    rng = random.Random(3403)
    QL = K.default_Q_list(rng, extra=300)
    clause_c(QL)
    clause_d()

    print("\n" + "=" * 78)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 17: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: at (s,i) = (3,4) with the house Gram, over both abelian")
    print("         groups of order 4, there are exactly 2048 admissible")
    print("         S-parts per group, all of the silent-seed shape; the")
    print("         sigma-shell at N = 2092 is 4192 x 64 = 2048 x 131 =")
    print("         268 288 tables per group; and every class in the sweep")
    print("         run above admits a border kit, re-verified exactly at")
    print("         w = 130.  LABEL: PROVEN (Theorem F, paper-grade,")
    print("         note/NOTE-B.md S1.8) + PROVEN-BY-CERTIFICATE (the")
    print("         censuses) + MEASURED (Theorem F on 7 / 7 banked records).")
    if args.full:
        print("         CENSUS: ALL 16384 classes, digest matched to the pin.")
    else:
        print("         CENSUS: the deterministic 256-class sample.  The")
        print("         16384-class census is `--full` (see NOTES.md).")
    print("         NOT claimed: anything about the SEED layer -- no")
    print("         house-profile quadruple on a group of order 520 is")
    print("         known; the number or structure of kits per S-part;")
    print("         anything at w <= 2s beyond the H(44) controls.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
