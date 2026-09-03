#!/usr/bin/env python3
"""cert 24 -- the three classes at order 2060 hold with the TRANSPOSE added
to the group.

  THE RELATION.  A ~ B is Hadamard equivalence: B = D_r P_r A P_c D_c.
  The TRANSPOSE-EXTENDED relation is  A ~~ B  iff  A ~ B or A ~ B^T.
  Refuting A ~~ B takes TWO refutations, one for each disjunct.  The
  first is profile(A) != profile(B).  For the second, note that ~ is
  symmetric and transposition an involution, so A ~ B^T iff A^T ~ B;
  EITHER of profile(A) != profile(B^T) and profile(A^T) != profile(B)
  therefore refutes it, and one of the two suffices.  The |T4| 4-profile
  is a Hadamard-equivalence invariant (note/NOTE-B.md S3.1, invariant
  I5), and the transpose of a Hadamard matrix is Hadamard, so each of
  those profiles is that same invariant computed on that other matrix.

  THEOREM (2060).  Let P be '2060-plain', the plain Goethals-Seidel
  array over the raw Z_515 seed; G be '2060-gist', the x104-twisted
  array that is byte-for-byte the publicly posted H(2060) (both cert
  07); and H'' be P with its TWELVE OFF-DIAGONAL 515-BLOCKS NEGATED --
  the orientation switch in its UNBORDERED form, 2060 being the
  degenerate s = 0 layer where the array is a plain 4x4 GS array of
  circulant blocks and there is no border to leave alone (cert 22).
  Then P, G and H'' are pairwise inequivalent UNDER THE
  TRANSPOSE-EXTENDED RELATION.  So ORDER 2060 CARRIES AT LEAST THREE
  HADAMARD EQUIVALENCE CLASSES with the transpose in the group -- the
  statement cert 22 explicitly withheld.  2060 and 1772 were the note's
  last two row-side separation statements, and cert 25 discharges 1772's
  the same day: with these two, NO separation statement in
  note/NOTE-B.md is row-side any longer.

  PROOF (finite, exact).  Every leg is an exact 4-profile comparison
  over all C(2060,4) = 748 155 697 135 row 4-subsets, in two
  arithmetics that agree bin for bin.

    P ~~ G   is refuted by  P !~ G   (146 of 147, cert 07)
                       and  P^T !~ G (134 of 134, new here);
    P ~~ H'' is refuted by  P !~ H'' (107 of 145, cert 22)
                       and  P !~ (H'')^T (145 of 146, new)
                       -- with P^T !~ H'' (145 of 146, new) as the
                          redundant second route;
    G ~~ H'' is refuted by  G !~ H'' (146 of 147, cert 22)
                       and  G !~ (H'')^T (134 of 134, new).

  An invariant that differs is a separation.  []

  WHICH ROUTE EACH PAIR TAKES, AND WHY IT MATTERS HERE.  Unlike cert 21
  at 1676, the two routes are NOT interchangeable at this order,
  because G^T is not banked (below).  The pair {P, G} is refuted only
  by the A^T vs B route (P^T against G); the pair {G, H''} only by the
  A vs B^T route (G against (H'')^T); the pair {P, H''} by both.  Every
  pair has at least one, which is all the relation asks, and clause [4]
  asserts exactly which entries are n/a and that each n/a pair is
  carried by its other route.

  REMARK, not the headline.  Under PLAIN Hadamard equivalence the five
  matrices profiled at this order -- P, G, H'', P^T, (H'')^T -- are
  pairwise inequivalent: all ten profile comparisons separate, the
  least separated pair by 92 bins.  FIVE classes are therefore
  EXHIBITED at 2060 by two constructions, the orientation switch and
  transposition.  The house counts THREE, because the
  transpose-extended relation is the one under which a matrix and its
  transpose are the same object, and three is the count that survives
  either convention.

  WHAT IS NOT PROFILED.  G^T -- the posted matrix's own transpose --
  was NOT enumerated by the campaign and is NOT banked here.  This
  certificate therefore says NOTHING about G vs G^T at 2060 (cert 19
  decided the analogous question at 668 only), the "A vs B^T" route is
  n/a for the two pairs whose B is G, and the plain-equivalence count
  above is "at least FIVE exhibited", not eight.  The theorem is
  unaffected: each pair's two refutations are both in hand.

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The two new C(2060,4) enumerations -- (H'')^T and P^T -- were not run
  inside this repository.  They ran in the source laboratory
  (Hadamard-2060, experiments/inequiv/exact_profile_big.py, the engine
  unchanged since the 2060 registration, numpy, 16 threads on a rented
  c2d-highcpu-16, prof42-2, us-east1-b, 2026-09-03) under the
  pre-registration experiments/pr0042/REGISTRATION.md, FLUSHED 10:17
  UTC on 2026-09-02 before any matrix it governs was built, which fixed
  the objects (S2: H_2060-orient-T and H_2060-plain-T are named there,
  and no third 2060 transpose is), the decision rules (S4: "differs in
  any bin => inequivalent"; "equal in every bin => MEASURED, nothing
  proved"; and that the two transposes give the transpose-extended form
  for all three pairs) and the kill criteria (S5: blas != bits in any
  bin is a hard stop) in advance.  The matrices were built AND verified
  at the desk (experiments/pr0042/build_matrices.py, every digest in
  its manifest.json, each through this repository's verify/verify.py,
  the plain source re-verified against cert 07's pin first); the rented
  machine enumerated and nothing else.  The DEFAULT path of this script
  AUDITS all ten banks.  Say "banked exact computation AUDITED" of a
  default run; the word "replayed" belongs to --full, which is OFFERED
  AND PRICED BELOW BUT HAS NOT BEEN RUN IN THIS REPOSITORY AT THIS
  ORDER (as in cert 22).

WHAT THIS SCRIPT DOES  (default path: standard library only, ~20 s)

  (0) Pins the SHA-256 of all ten banked files it reads.
  (1) Rebuilds P and G from data/sep2060-records.json -- the four
      normalised seeds, the normalising shifts, the twist multiplier and
      the CRT relabelling -- after re-verifying the classical
      Goethals-Seidel condition sum_q PAF_q(t) = 4v*[t = 0] on the raw
      seeds (the s = 0 layer of the master theorem, checked, not
      assumed); verifies both through verify/verify.py and pins their
      canonical digests against cert 07's.  Forms H'' by negating the
      twelve off-diagonal 515-blocks, checks that exactly 12*515^2 cells
      changed, checks the alternate-orientation identity BOTH ways --
      as a sign pattern cell by cell, and against the alternate array
      ASSEMBLED from the same raw seeds -- and pins H''.  Then
      TRANSPOSES P and H'' in-process, hands all five matrices to
      verify/verify.py, and pins all five canonical digests.  G^T is
      neither built nor pinned: no profile of it exists.
  (1b) Control C7 -- the dim V / dim W trap on the real objects.
  (2) Loads and AUDITS ten banked exact 4-profiles (five matrices x
      two implementations) -- it does not recompute them -- asserting
      in exact integer arithmetic: every populated bin = 4 (mod 8)
      (2060 = 4 mod 8); the bin counts (145 / 133 / 145 / 123 / 123);
      every bin key canonical in [0, 2060] and every count a positive
      integer; the counts total C(2060,4); the second moment equals
      n^3(n-1)(n-2)/24, recomputed here AND compared against the field
      the bank declares; the schema, folding, arithmetic and matrix
      each bank declares; and each bank's declared matrix digest
      against the in-process digest of the matrix rebuilt in THIS run.
      Then blas == bits bin for bin on each of the five matrices.  Two
      schemas are in play -- cert 07's sep2060-exact-profile/1 and the
      campaign's exact-4-profile/1 -- and each is checked on its own
      terms.
  (3) All ten pair comparisons, in both arithmetics, each asserted to
      its exact differing-bin count and union support size.
  (4) The TRANSPOSE-EXTENDED verdicts, DERIVED IN CODE from those
      counts: a pair is separated only if profile(A) != profile(B) and
      at least one of the two second-refutation routes is banked and
      nonzero.  Then the plain-equivalence remark, also derived.
  (5) Controls: five small Hadamard matrices profiled two ways, one of
      them the route --full takes; the UNBORDERED orientation switch
      exercised end to end on the GS controls H(28) and H(36), which
      are the exact structural analogue of the 2060 switch (s = 0
      there too); the transposed-profile route on matrices small
      enough for straight enumeration (symmetric Sylvester H(8)/H(16),
      where profile(M^T) = profile(M) is FORCED and is checked; Paley
      H(20), which is not symmetric and is MEASURED, never asserted);
      the support control (at 2060 the transposes DROP 23 of the
      originals' bins and ADD one, where at 1676 they swap two and at
      668/716 they drop one); the comparator in the null direction; a
      total-preserving corruption only the second moment can catch;
      and the dim-V trap on Sylvester H(16).
  (6) --full: recompute a transposed 2060 profile here from the rows
      clause [1] verified (certs/06-668-separation/full_recompute.py,
      imported by path, not copied, so certs
      06/08/11/13/14/15/19/20/21/22/23/24 cannot drift apart), after a
      smoke test against the forced profile of Sylvester H(128).
      numpy is imported only under this flag, is finder-side only and
      is never in the trust chain; BLAS threads are capped at three
      before it loads.  PRICE, and why it has not been run: exactly
      cert 22's, since the order and the module are the same -- one
      2060 leg here is of order 15 h (the source laboratory's MEASURED
      716->2060 ratio of 137 applied to the 400.3 s 716 leg this same
      module took in this repository) or 22 h on the Theta(n^5) law
      quoted elsewhere here, and the blas route materialises a
      2 120 770 x 2060 pair matrix -- 4.4 GB as int8, 17.5 GB as the
      float32 copy -- far past this desk.  The banked legs' own
      seconds fields price the campaign's engine for comparison:
      6 645.1 s blas and 17 685.3 s bits for (H'')^T, 6 650.4 s and
      17 678.4 s for P^T, at 16 rented threads.  NO --full LEG HAS
      BEEN RUN IN THIS REPOSITORY AT 2060.

Usage:
  python certs/24-transpose-extended-2060/run.py
  python certs/24-transpose-extended-2060/run.py --full --impl bits --matrix or-T
  python certs/24-transpose-extended-2060/run.py --full
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

N = 2060
V = 515
OUT = os.path.join(HERE, "out")

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the five matrices (the digest verify/verify.py
# reports).  ALL FIVE are rebuilt in clause [1] of THIS run: P and G carry
# the pins cert 07 fixed at this order (G's is additionally the digest of
# the publicly posted artifact), H'' is formed here from P (cert 22's pin),
# and the two transposes are formed here from P and H''.  G^T is
# deliberately absent: the campaign never profiled it, so this certificate
# neither builds nor pins it.
SHA = {
    "plain":  "510f89b7b423c85da0c7cada52cb0f62e0415d736d2040d6701f66c4e524cf6a",
    "gist":   "c7a145d86210740dd3f8ea21ca896a54d6916007a042638f17c8c47f097200f7",
    "orient": "4e1891b095b8aafa21176e494038f199b495c96a840bdb003e231c160870b801",
    "plain-T": "5a980e3ad69f02fdece4f0aca40f46b9afafe15703e87b82941ff661db1a2960",
    "orient-T": "8558904d9d61c7547b835c25791da97f7d3e0bc1cd852de082b8475c01b34337",
}
# The alternate-orientation array assembled directly from the seeds in
# clause [1](b).  It is not a campaign object and no profile of it is
# banked; it is pinned because it is built and verified on every run.
SHA_ALT = \
    "40e1d1c8cd40e94016c453f12e520a8518e7d29b773d3adaae3f484eca64398d"

# SHA-256 of the banked data FILES themselves, so a silently edited bank is
# a hard error rather than a different theorem.  The two *-T-* pairs are
# this certificate's own (banked 2026-09-03 under cert 24); the two orient
# banks are cert 22's and the four plain/gist banks cert 07's, all reused
# verbatim and re-pinned here at exactly the values those certificates
# carry, so the certificates cannot drift apart on shared evidence.
# data/sep2060-records.json is NOT file-pinned here: it is shared with
# certs 07 and 22, and the binding pin on it is the canonical digest of
# each matrix it produces, checked in clause [1].
FILE_PINS = {
    "data/sep2060-plain-T-exact-blas.json":
        "1f8aa2469f22d3fee20ce3cf3618dbc275c81d4f0a6493271a0532eb77e9ae57",
    "data/sep2060-plain-T-exact-bits.json":
        "8a5389e0b74c0d27c5981ac1ca5605b834517f751e83a81b14d8759aa9e557d6",
    "data/sep2060-orient-T-exact-blas.json":
        "b21f217eaa70dfe01da599cac3c260f4cef9436d398833962bc46224298a771f",
    "data/sep2060-orient-T-exact-bits.json":
        "468386dab0c5b6c7c317cdc4f3694113f1a2d8cfac449f2a19ba618209c56694",
    "data/sep2060-orient-exact-blas.json":
        "38135aef205b4428760dc0439b29196b6776f215397fb34880156c66ee283f00",
    "data/sep2060-orient-exact-bits.json":
        "18dcf4c2e5d603324182eeef45c13e89ef80d5c0c0a5add3a7f90c333e4e5e87",
    "data/sep2060-exact-blas-plain.json":
        "5428aeac7b570fff55975c2b737fae9e8d0b717ec511735b68893e609a0037d8",
    "data/sep2060-exact-bits-plain.json":
        "e6c3af94712d0ba5cf3a3047796ccd474970036fec211b41a5579b7ff892ca49",
    "data/sep2060-exact-blas-gist.json":
        "a20b9a63cd3d93046c251b5c19aabeeac412b8f7933bbafa82d0210320e3aef0",
    "data/sep2060-exact-bits-gist.json":
        "9d8cc4b55c297c7e948df3e7639613a0580fc3e54af9eb12399bc010337f8a93",
}

# (tag, impl) -> (file, schema, the name the bank must declare)
PROFILES = {
    ("plain", "blas"): ("data/sep2060-exact-blas-plain.json",
                        "sep2060-exact-profile/1", "plain"),
    ("plain", "bits"): ("data/sep2060-exact-bits-plain.json",
                        "sep2060-exact-profile/1", "plain"),
    ("gist", "blas"): ("data/sep2060-exact-blas-gist.json",
                       "sep2060-exact-profile/1", "gist"),
    ("gist", "bits"): ("data/sep2060-exact-bits-gist.json",
                       "sep2060-exact-profile/1", "gist"),
    ("orient", "blas"): ("data/sep2060-orient-exact-blas.json",
                         "exact-4-profile/1", "H_2060-orient"),
    ("orient", "bits"): ("data/sep2060-orient-exact-bits.json",
                         "exact-4-profile/1", "H_2060-orient"),
    ("plain-T", "blas"): ("data/sep2060-plain-T-exact-blas.json",
                          "exact-4-profile/1", "H_2060-plain-T"),
    ("plain-T", "bits"): ("data/sep2060-plain-T-exact-bits.json",
                          "exact-4-profile/1", "H_2060-plain-T"),
    ("orient-T", "blas"): ("data/sep2060-orient-T-exact-blas.json",
                           "exact-4-profile/1", "H_2060-orient-T"),
    ("orient-T", "bits"): ("data/sep2060-orient-T-exact-bits.json",
                           "exact-4-profile/1", "H_2060-orient-T"),
}

LABEL = {"plain": "P", "gist": "G", "orient": "H''",
         "plain-T": "P^T", "orient-T": "(H'')^T"}

TAGS = ("plain", "gist", "orient", "plain-T", "orient-T")

# The three classes the theorem separates, and the transpose of each that
# is BANKED.  "gist" has no entry: G^T was never enumerated.
CLASSES = ("plain", "gist", "orient")
TRANSPOSE_OF = {"plain": "plain-T", "orient": "orient-T"}

# The bin counts, pinned so a drifting bank cannot quietly turn this
# certificate into a different (weaker) statement.
NBINS = {"plain": 145, "gist": 133, "orient": 145,
         "plain-T": 123, "orient-T": 123}

# Every pair, with its differing-bin count and the size of the union of the
# two supports.  NEW = a comparison that needs one of the two profiles
# banked under cert 24; cert 07 / cert 22 = re-affirmed from those banks.
SEP = {
    ("plain", "gist"):        (146, 147, "cert 07"),
    ("plain", "orient"):      (107, 145, "cert 22"),
    ("gist", "orient"):       (146, 147, "cert 22"),
    ("plain", "plain-T"):     (145, 146, "NEW"),
    ("plain", "orient-T"):    (145, 146, "NEW"),
    ("gist", "plain-T"):      (134, 134, "NEW"),
    ("gist", "orient-T"):     (134, 134, "NEW"),
    ("orient", "plain-T"):    (145, 146, "NEW"),
    ("orient", "orient-T"):   (145, 146, "NEW"),
    ("plain-T", "orient-T"):  (92, 123, "NEW"),
}

# The support structure, asserted.  P and H'' share one 145-bin support
# (cert 22: the switch moves counts, not support) and P^T and (H'')^T
# share another of 123 bins.  Unlike 1676, where the transposes SWAP two
# bins, and unlike 668/716, where each transpose DROPS one, at 2060 the
# transposes drop TWENTY-THREE of the originals' bins and add ONE.
ORIG_ONLY = [844, 860, 876, 924, 940, 956, 972, 988, 1004, 1020, 1036,
             1052, 1068, 1076, 1084, 1100, 1116, 1124, 1148, 1156, 1164,
             1180, 1236]
TRANS_ONLY = [1204]
NUNION_MIXED = 146              # |support(original) U support(transpose)|
NUNION_ALL = 148                # all five supports together
# G's own support structure, as cert 22 asserts it.
GIST_MISSING = [940, 972, 988, 1004, 1020, 1036, 1052, 1068, 1084, 1100,
                1116, 1148, 1164, 1180]
GIST_ONLY = [892, 908]
# The ONE bin where an original (P or H'') agrees with a transpose, and
# its count there; G reads 1380 in that bin, which is why G against a
# transpose separates in EVERY bin of the union.
AGREE_MIXED = [900]
AGREE_MIXED_COUNT = 300
GIST_AT_AGREE = 1380
# Where P^T vs (H'')^T stops separating: the top differing bin, and how
# many bins above it agree.
TOP_DIFF_TT = 788
TAIL_TT = 24
# The originals' top bin, which is populated in no transpose -- so at
# 2060, unlike 668/716/1676/1772, the extreme tail DOES separate on every
# original-versus-transpose comparison.
TOP_BIN_ORIG = 1236

HEXDIGITS = set("0123456789abcdef")
FLIP = str.maketrans("+-", "-+")
FAIL = []
NCHECK = [0]


def is_sha256(v):
    return (isinstance(v, str) and len(v) == 64
            and all(c in HEXDIGITS for c in v))


def check(label, cond, extra=""):
    NCHECK[0] += 1
    print("  [%s] %s%s" % ("PASS" if cond else "FAIL", label,
                           ("  " + extra) if extra else ""))
    if not cond:
        FAIL.append(label)
    return cond


def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def rows_sha256(rows):
    """The canonical digest verify.py reports: '+/-' rows, newline-joined."""
    h = hashlib.sha256()
    for r in rows:
        h.update((r + "\n").encode("ascii"))
    return h.hexdigest()


if hasattr(int, "bit_count"):  # 3.10+
    def popcount(x):
        return x.bit_count()
else:  # 3.9 fallback -- kept local so this certificate imports nothing
    def popcount(x):        # from the trust chain
        return bin(x).count("1")


# ======================================================================
# exact-integer profile arithmetic
# ======================================================================

def c_n_4(n):
    return n * (n - 1) * (n - 2) * (n - 3) // 24


def second_moment_want(n):
    """sum over 4-subsets of T4^2 = n^3 (n-1)(n-2) / 24  (NOTE-B.md S3.1).

    Proof.  Summing T4^2 over ORDERED 4-tuples gives
    sum_{c,c'} (sum_i H[i][c]H[i][c'])^4 = n^5 by column orthogonality.
    The only degenerate tuples with T4 != 0 are those whose index multiset
    has all multiplicities even: n tuples with all four indices equal and
    3n(n-1) with two values twice each, each contributing n^2.  Hence
    24 * sum_{4-subsets} T4^2 = n^5 - 3n^4 + 2n^3.
    """
    return n ** 3 * (n - 1) * (n - 2) // 24


def audit(prof, n, where):
    """The two forced identities, as ASSERTS.  Exact integers throughout.

    Neither identity is tuned for by anything that produced `prof`, and the
    second one is sharp: at n = 2060 it pins a 16-digit number to the unit.
    """
    tot = sum(prof.values())
    want = c_n_4(n)
    assert tot == want, (
        "%s: profile totals %d, not C(%d,4) = %d" % (where, tot, n, want))
    m2 = sum(k * k * v for k, v in prof.items())
    m2w = second_moment_want(n)
    assert m2 == m2w, (
        "%s: second-moment identity FAILED -- sum T4^2 = %d, want "
        "n^3(n-1)(n-2)/24 = %d" % (where, m2, m2w))
    return tot, m2


def compare(pa, pb):
    """(divergent bins, delta sum, first-moment delta, union support size)."""
    ks = sorted(set(pa) | set(pb))
    diff = [(k, pa.get(k, 0), pb.get(k, 0)) for k in ks
            if pa.get(k, 0) != pb.get(k, 0)]
    return (diff,
            sum(q - p for _k, p, q in diff),
            sum(k * (q - p) for k, p, q in diff),
            len(ks))


# ======================================================================
# clause 1 -- rebuild the pair from the seed record, switch, transpose,
# verify, pin
#
# The block algebra below is cert 07's, restated here (as cert 22 restates
# it) rather than imported so that this certificate stands on its own
# file; every convention it uses is written out in
# data/sep2060-records.json's `conventions` block and is checked against
# that block on every run.
# ======================================================================

def circ(sq):
    """circ(x)[i][j] = x[(i - j) mod v]."""
    v = len(sq)
    r0 = "".join(sq[(-j) % v] for j in range(v))
    return [r0[v - i:] + r0[:v - i] for i in range(v)]


def rev_seq(sq):
    v = len(sq)
    return "".join(sq[(-k) % v] for k in range(v))


def twist(sq, k):
    """twist(x, k)[i][j] = circ(x)[i][(k*j) mod v] -- the multiplier twist."""
    v = len(sq)
    sig = [(k * j) % v for j in range(v)]
    return ["".join(row[t] for t in sig) for row in circ(sq)]


def neg(rs):
    return [r.translate(FLIP) for r in rs]


def mulJ(rs):
    return [r[::-1] for r in rs]


def conjJ(rs):
    return [r[::-1] for r in rs[::-1]]


def relabel(rs, inv):
    return ["".join(rs[i][t] for t in inv) for i in inv]


def grid_rows(grid, v):
    """Flatten a 4x4 grid of v x v blocks row-major (the record's rule)."""
    out = []
    for r in range(4):
        b0, b1, b2, b3 = grid[r]
        for p in range(v):
            out.append(b0[p] + b1[p] + b2[p] + b3[p])
    return out


def raw_seeds(recs):
    """raw[q][u] = normalised_seed[q][(u + shift_q) mod v]."""
    ns = recs["normalised_seed"]
    e = list(recs["normalising_shifts"])
    return ["".join(ns[q][(u + e[q]) % V] for u in range(V))
            for q in range(4)]


def gs_condition(seqs):
    """The classical Goethals-Seidel condition on the seed quadruple:
    sum_q PAF_q(t) = 4v * [t = 0].  This is the s = 0 (unbordered) layer
    of the master theorem -- the hypothesis that makes the plain array
    Hadamard -- and it is re-verified here rather than assumed."""
    v = len(seqs[0])
    xs = [[1 if ch == "+" else -1 for ch in s] for s in seqs]
    for t in range(v):
        tot = 0
        for x in xs:
            tot += sum(x[h] * x[(h + t) % v] for h in range(v))
        if tot != (4 * v if t == 0 else 0):
            return False, t
    return True, None


def build_pair(recs):
    """P and G from the banked seed record.  Returns (plain, gist, raw)."""
    raw = raw_seeds(recs)
    a, b, c, d = raw

    A = circ(a)
    Bj, Cj, Dj = mulJ(circ(b)), mulJ(circ(c)), mulJ(circ(d))
    Bt, Ct, Dt = (mulJ(circ(rev_seq(b))), mulJ(circ(rev_seq(c))),
                  mulJ(circ(rev_seq(d))))
    plain = [[A, Bj, Cj, Dj],
             [neg(Bj), A, Dt, neg(Ct)],
             [neg(Cj), neg(Dt), A, Bt],
             [neg(Dj), Ct, neg(Bt), A]]

    k = int(recs["gist_array"]["twist"])
    B_, C_, D_ = mulJ(twist(b, k)), mulJ(twist(c, k)), mulJ(twist(d, k))
    X12, X13, X23 = neg(conjJ(D_)), conjJ(C_), neg(conjJ(B_))
    tw = [[A, B_, C_, D_],
          [neg(B_), A, X12, X13],
          [neg(C_), neg(X12), A, X23],
          [neg(D_), neg(X13), neg(X23), A]]
    ORD = [5 * (g % 103) + (g % 5) for g in range(V)]   # the CRT relabelling
    inv = [0] * V
    for g in range(V):
        inv[ORD[g]] = g
    tw = [[relabel(blk, inv) for blk in row] for row in tw]
    return grid_rows(plain, V), grid_rows(tw, V), raw


def build_alternate(raw):
    """The ALTERNATE Goethals-Seidel orientation over the same raw seeds:
    NOTE-B.md S1.0's standard array with its SIX TRANSPOSED BLOCKS
    negated.  Assembled from the seeds, independently of H''."""
    a, b, c, d = raw
    A = circ(a)
    Bj, Cj, Dj = mulJ(circ(b)), mulJ(circ(c)), mulJ(circ(d))
    Bt, Ct, Dt = (mulJ(circ(rev_seq(b))), mulJ(circ(rev_seq(c))),
                  mulJ(circ(rev_seq(d))))
    alt = [[A, Bj, Cj, Dj],
           [neg(Bj), A, neg(Dt), Ct],
           [neg(Cj), Dt, A, neg(Bt)],
           [neg(Dj), neg(Ct), Bt, A]]
    return grid_rows(alt, V)


def transpose_rows(rows):
    return ["".join(col) for col in zip(*rows)]


def verify_rows(tag, rows, want_sha):
    """Hand a matrix to the trust chain and pin its canonical digest."""
    path = os.path.join(OUT, "H2060_%s.txt" % tag)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify", "verify.py"), path],
        capture_output=True, text=True)
    verdict = ((proc.stdout or proc.stderr).strip().splitlines()
               or ["(no output)"])[-1]
    lab = "%-8s" % tag
    check("%s verify/verify.py" % lab, proc.returncode == 0, verdict[:70])
    dig = rows_sha256(rows)
    check("%s canonical sha256 == pin" % lab, dig == want_sha,
          dig[:24] + "...")
    check("%s verify.py reports the same digest" % lab, dig in verdict)
    os.remove(path)                    # 4.3 MB apiece; six of them
    return dig


def orientation_switch(rows, n):
    """H'' = H with the twelve off-diagonal blocks negated, UNBORDERED.

    At 2060 there is no border: superblock I = k // n for every row and
    column index k, and every off-diagonal block is negated.  (The
    bordered form of certs 13/14/20/23 takes the first 4s rows and columns
    out first; here s = 0 and that step is empty.)
    """
    Nn = len(rows)
    assert Nn == 4 * n
    out = []
    for r, row in enumerate(rows):
        I = r // n
        chars = list(row)
        for c in range(Nn):
            if c // n != I:
                chars[c] = "+" if chars[c] == "-" else "-"
        out.append("".join(chars))
    return out


def alternate_orientation_identity(rows_H, rows_Hpp, n):
    """S H'' S == H_alt, where H_alt negates the six transposed blocks of
    the standard array (blocks (I,J) with I,J >= 1, I != J) and
    S = diag(1,-1,-1,-1) (x) I_n.  A sign-pattern identity, checked cell
    by cell.  This is certs 13/14/20/23's clause with the border term
    gone, s = 0 having emptied it."""
    Nn = len(rows_H)
    sg = [1 if k // n == 0 else -1 for k in range(Nn)]
    for r in range(Nn):
        hr, pr = rows_H[r], rows_Hpp[r]
        I = r // n
        for c in range(Nn):
            J = c // n
            alt = -1 if (I >= 1 and J >= 1 and I != J) else 1
            lhs = sg[r] * sg[c] * (1 if pr[c] == "+" else -1)
            rhs = alt * (1 if hr[c] == "+" else -1)
            if lhs != rhs:
                return False
    return True


def signed_conjugate(rows, n):
    """S M S with S = diag(1,-1,-1,-1) (x) I_n -- an element of the
    Hadamard equivalence group (a diagonal +-1 on each side), so it moves
    no invariant."""
    Nn = len(rows)
    sg = [1 if k // n == 0 else -1 for k in range(Nn)]
    out = []
    for r, row in enumerate(rows):
        if sg[r] == 1:
            out.append("".join(ch if sg[c] == 1 else
                               ("+" if ch == "-" else "-")
                               for c, ch in enumerate(row)))
        else:
            out.append("".join(ch if sg[c] == -1 else
                               ("+" if ch == "-" else "-")
                               for c, ch in enumerate(row)))
    return out


def cells_changed(a, b):
    """How many entries differ between two sign-string matrices."""
    tot = 0
    for ra, rb in zip(a, b):
        tot += sum(1 for x, y in zip(ra, rb) if x != y)
    return tot


# ======================================================================
# F2 difference code -- dim V (NOT invariant) and dim W (invariant)
# ======================================================================

def f2_rows(rows):
    """'+' -> 0, '-' -> 1, one Python int per row (bit c = column c)."""
    tab = str.maketrans("+-", "01")
    return [int(r.translate(tab)[::-1], 2) for r in rows]


def f2_rank(vectors):
    """Rank over F2 by highest-set-bit pivoting.  Exact, no floats."""
    piv = {}
    rk = 0
    for v in vectors:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                rk += 1
                break
    return rk


def dim_V_W(rows):
    """V = span{r_i + r_j}, W = V + <all-ones>.

    V = span{r_i + r_0 : i > 0} because r_i + r_j = (r_i+r_0) + (r_j+r_0),
    so the n-1 generators below span the same space as all C(n,2) of them.
    """
    f = f2_rows(rows)
    n = len(rows)
    gens = [f[i] ^ f[0] for i in range(1, n)]
    ones = (1 << len(rows[0])) - 1
    return f2_rank(gens), f2_rank(gens + [ones])


# ======================================================================
# control matrices and two independent full-profile routes
# ======================================================================

def sylvester(k):
    n = 1 << k
    return ["".join("+" if bin(x & y).count("1") % 2 == 0 else "-"
                    for y in range(n)) for x in range(n)]


def sylvester_profile_forced(n):
    """The 4-profile of Sylvester H(n) is forced, so it is a control with a
    PREDICTED answer rather than merely a self-consistent one.

    Rows are indexed by F2^k with H[x][y] = (-1)^<x,y>, so
    T4({a,b,c,d}) = sum_y (-1)^<a+b+c+d, y> = n if a+b+c+d = 0 and 0
    otherwise.  Choosing any three distinct a, b, c forces d = a+b+c, and
    d is automatically outside {a,b,c} (d = a would force b = c), so the
    4-subsets summing to zero are counted n(n-1)(n-2) times over ordered
    triples and each 4-subset supplies 4*3*2 = 24 of them.
    """
    hit = n * (n - 1) * (n - 2) // 24
    return {0: c_n_4(n) - hit, n: hit}


def paley1(q):
    """Paley type I Hadamard matrix of order q+1, q prime = 3 (mod 4)."""
    res = {(x * x) % q for x in range(1, q)}

    def chi(a):
        a %= q
        return 0 if a == 0 else (1 if a in res else -1)

    n = q + 1
    M = [[0] * n for _ in range(n)]
    M[0][0] = 1
    for j in range(1, n):
        M[0][j] = 1
        M[j][0] = -1
    for i in range(1, n):
        for j in range(1, n):
            M[i][j] = 1 if i == j else chi(i - j)
    return ["".join("+" if v == 1 else "-" for v in row) for row in M]


# Four +-1 sequences of length v with sum_q PAF_q(t) = 4v*[t = 0] -- the
# classical Goethals-Seidel condition, i.e. exactly the s = 0 layer the
# 2060 arrays inhabit.  The PAF condition is RE-VERIFIED below, never
# assumed, and each assembled matrix goes through the same C0 Hadamard
# check as every other control.  These two give profiles with three and
# four populated bins, which Sylvester and Paley matrices do not.
GS_CONTROLS = [
    (7, ["+------", "+++----", "++--+--", "+-+-+--"]),            # -> H(28)
    (9, ["++-------", "++-+-----", "++-+-+---", "++--+-+--"]),    # -> H(36)
]


def gs_control(v, seqs):
    """Assemble the plain (unbordered) GS array over Z_v from four seeds,
    with this file's own block algebra -- the same code path that builds
    the 2060 plain array, at a size where the Hadamard property can be
    checked by brute force.  Returns (rows, paf_ok, raw)."""
    a, b, c, d = seqs
    A = circ(a)
    Bj, Cj, Dj = mulJ(circ(b)), mulJ(circ(c)), mulJ(circ(d))
    Bt, Ct, Dt = (mulJ(circ(rev_seq(b))), mulJ(circ(rev_seq(c))),
                  mulJ(circ(rev_seq(d))))
    grid = [[A, Bj, Cj, Dj],
            [neg(Bj), A, Dt, neg(Ct)],
            [neg(Cj), neg(Dt), A, Bt],
            [neg(Dj), Ct, neg(Bt), A]]
    ok, _t = gs_condition(seqs)
    return grid_rows(grid, v), ok, list(seqs)


def gs_control_alternate(v, seqs):
    """The same four seeds in the ALTERNATE orientation."""
    a, b, c, d = seqs
    A = circ(a)
    Bj, Cj, Dj = mulJ(circ(b)), mulJ(circ(c)), mulJ(circ(d))
    Bt, Ct, Dt = (mulJ(circ(rev_seq(b))), mulJ(circ(rev_seq(c))),
                  mulJ(circ(rev_seq(d))))
    grid = [[A, Bj, Cj, Dj],
            [neg(Bj), A, neg(Dt), Ct],
            [neg(Cj), Dt, A, neg(Bt)],
            [neg(Dj), neg(Ct), Bt, A]]
    return grid_rows(grid, v)


def is_hadamard(rows):
    """C0: a control that is not what it says it is invalidates the control."""
    n = len(rows)
    S = [[1 if ch == "+" else -1 for ch in r] for r in rows]
    if any(len(r) != n for r in S):
        return False
    for i in range(n):
        for j in range(i + 1, n):
            if sum(a * b for a, b in zip(S[i], S[j])) != 0:
                return False
    return True


def profile_straight(rows):
    """Route 1: straight O(C(n,4)) enumeration over the +-1 entries.

    No packing, no popcount, no pair vectors -- T4 is summed column by
    column from the signs, exactly as it is defined.
    """
    n = len(rows)
    S = [[1 if ch == "+" else -1 for ch in r] for r in rows]
    hist = {}
    for a in range(n):
        Sa = S[a]
        for b in range(a + 1, n):
            Sb = S[b]
            ab = [Sa[c] * Sb[c] for c in range(n)]
            for c_ in range(b + 1, n):
                Sc = S[c_]
                abc = [ab[c] * Sc[c] for c in range(n)]
                for d in range(c_ + 1, n):
                    Sd = S[d]
                    t = 0
                    for c in range(n):
                        t += abc[c] * Sd[c]
                    t = abs(t)
                    hist[t] = hist.get(t, 0) + 1
    return {k: v for k, v in hist.items() if v}


def profile_pairvec(rows):
    """Route 2: the pair-vector / Gram-triangle route, in stdlib integers.

    This is the SAME bookkeeping certs/06-668-separation/full_recompute.py
    relies on -- histogram |<u_P,u_Q>| over the upper triangle of U U^T
    including the diagonal, drop the m diagonal terms at bin n, drop the
    n*C(n-1,2) index-sharing pairs at bin 0, divide the rest by 3 -- and
    running it here on matrices small enough for route 1 is what validates
    the route `--full` takes.  Bit-packed popcount arithmetic, i.e. a
    different language of arithmetic from route 1.
    """
    n = len(rows)
    pk = []
    for r in rows:
        v = 0
        for c, ch in enumerate(r):
            if ch == "-":
                v |= 1 << c
        pk.append(v)
    U = [pk[i] ^ pk[j] for i in range(n) for j in range(i + 1, n)]
    m = len(U)
    hist = {}
    for a in range(m):
        ua = U[a]
        for b in range(a, m):
            t = abs(n - 2 * popcount(ua ^ U[b]))
            hist[t] = hist.get(t, 0) + 1
    hist[n] = hist.get(n, 0) - m                       # drop P == Q
    inter = n * (n - 1) * (n - 2) // 2                 # pairs sharing an index
    hist[0] = hist.get(0, 0) - inter
    out = {}
    for k, v in hist.items():
        if v == 0:
            continue
        assert v % 3 == 0, "triangle histogram not divisible by 3 at %d" % k
        out[k] = v // 3
    return out


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="RECOMPUTE a transposed exact 2060 profile here "
                         "with numpy and compare to the banks bin for bin; "
                         "without it the banked profiles are audited, not "
                         "recomputed.  OF ORDER 15 HOURS per leg at this "
                         "order, and the blas path wants ~17.5 GB; not yet "
                         "run in this repository")
    ap.add_argument("--impl", choices=("blas", "bits", "both"),
                    default="both",
                    help="which arithmetic path --full recomputes")
    ap.add_argument("--matrix", choices=("pl-T", "or-T", "all"),
                    default="all",
                    help="which NEW profile --full recomputes (default both; "
                         "at this order one matrix one way is already many "
                         "hours)")
    args = ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 24 -- order 2060: the three classes hold with the TRANSPOSE")
    print("           added to the group")
    print("=" * 72)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    # The scratch directory holds one ~4.3 MB generated matrix at a time.
    # Wrapping the body means an exception -- in the optional numpy path or
    # anywhere else -- cannot leave one on disk.
    try:
        rc = _body(args, t_start)
    finally:
        shutil.rmtree(OUT, ignore_errors=True)
    print("generated matrices deleted; nothing left in %s   (%.1fs)"
          % (rel(OUT), time.time() - t_start))
    return rc


def _audit_bank(tag, impl, p, blob, schema, name, built):
    """The forced identities, the declared headers, and the matrix binding
    -- against the in-process digest of the matrix rebuilt in THIS run."""
    lab = "%-7s %-4s" % (LABEL[tag], impl)
    tot, m2 = audit(p, N, "%s/%s" % (tag, impl))
    check("%s  %d bins, all |T4| = 4 (mod 8)" % (lab, NBINS[tag]),
          len(p) == NBINS[tag] and all(k % 8 == 4 for k in p))
    check("%s  total == C(2060,4) == %d" % (lab, c_n_4(N)), tot == c_n_4(N))
    check("%s  second moment == n^3(n-1)(n-2)/24 == %d"
          % (lab, second_moment_want(N)), m2 == second_moment_want(N))
    # Every count is a positive integer and every key a canonical decimal
    # in [0, 2060]: a bank cannot smuggle a negative or an empty bin past
    # the two identities by cancellation.
    check("%s  every bin key is canonical in [0, %d] and every count a "
          "positive integer" % (lab, N),
          all(isinstance(v, int) and v > 0 for v in p.values())
          and all(0 <= k <= N for k in p),
          "min bin %d, max bin %d, min count %d"
          % (min(p), max(p), min(p.values())))
    # The declared headers.  The two schemas carry different header sets
    # (cert 07's banks predate the campaign's), so each is checked on its
    # own terms and neither is allowed to be silently absent.
    if schema == "exact-4-profile/1":
        hdr_ok = (int(blob.get("second_moment", -1)) == m2
                  and int(blob.get("second_moment_want", -1)) == m2
                  and int(blob.get("total", -1)) == tot
                  and int(blob.get("n", -1)) == N
                  and int(blob.get("C_n_4", -1)) == c_n_4(N))
        named = (blob.get("matrix") == name
                 and blob.get("producer_filename") == name + ".txt"
                 and bool(blob.get("engine")))
        extra = "matrix=%r engine=%r" % (blob.get("matrix"),
                                         blob.get("engine"))
    else:
        hdr_ok = (int(blob.get("second_moment", -1)) == m2
                  and int(blob.get("total", -1)) == tot
                  and int(blob.get("n", -1)) == N
                  and "C_n_4" not in blob)
        named = blob.get("tag") == name and bool(blob.get("producer"))
        extra = "tag=%r" % (blob.get("tag"),)
    check("%s  banked second_moment / total / n headers agree with the "
          "recomputation" % lab, hdr_ok)
    check("%s  bank declares schema %s, the |T4| folding, its arithmetic "
          "and its matrix" % (lab, schema),
          blob.get("schema") == schema
          and str(blob.get("folded", "")).startswith("|T4|")
          and blob.get("impl") == impl and named, extra)
    # Matrix identity.  The declared digest is compared against the digest
    # of the matrix REBUILT IN THIS RUN -- not against a static string --
    # so a bank cannot drift onto a different object.  Keyed on PRESENCE,
    # not truthiness: a declared digest that is empty, null, or not 64 hex
    # digits is a FAILURE.
    declared = (blob.get("matrix_canonical_sha256")
                or blob.get("matrix_sha256"))
    producer = (blob.get("matrix_sha256")
                or blob.get("producer_matrix_sha256"))
    check("%s  bank names the matrix rebuilt in THIS run (%s)"
          % (lab, name),
          is_sha256(declared) and declared == built[tag],
          (declared[:24] + "...") if is_sha256(declared)
          else "declared = %r" % (declared,))
    check("%s  the producer's own matrix_sha256 agrees with it" % lab,
          is_sha256(producer) and producer == declared)


def _body(args, t_start):
    # ---------------------------------------------------------- clause 0
    print("\n[0] the ten banked data files, SHA-256 pinned in this script")
    for name, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, name))
        check("%-42s" % name, got == want, got[:24] + "...")

    # ---------------------------------------------------------- clause 1
    print("\n[1] rebuild P and G from the banked seed record, re-verify the")
    print("    Goethals-Seidel condition, verify, pin; form H'' by the")
    print("    UNBORDERED orientation switch and check it two ways; then")
    print("    TRANSPOSE P and H'' in-process and put both through the trust")
    print("    chain.  G^T is NOT built: the campaign never profiled it, so")
    print("    this certificate neither pins nor uses it.")
    with open(os.path.join(ROOT, "data", "sep2060-records.json"),
              "r", encoding="utf-8") as fh:
        recs = json.load(fh)
    check("record is the order-2060 seed pair, v = %d, N = 4v = %d"
          % (V, N),
          int(recs["v"]) == V and int(recs["order"]) == N
          and recs["schema"] == "gs-2060-pair/1"
          and len(recs["normalised_seed"]) == 4
          and all(len(s) == V for s in recs["normalised_seed"]),
          "s = 0: four circulant blocks, no border")
    check("the record declares the two canonical digests this run must "
          "reproduce",
          recs["plain_array"]["canonical_sha256"] == SHA["plain"]
          and recs["gist_array"]["canonical_sha256"] == SHA["gist"]
          and int(recs["gist_array"]["twist"]) == 104)

    t0 = time.time()
    raw = raw_seeds(recs)
    paf_ok, bad_t = gs_condition(raw)
    check("the raw seeds satisfy sum_q PAF_q(t) = 4v*[t = 0] -- the s = 0 "
          "layer of the master theorem, re-verified here",
          paf_ok, "v = %d, row sums %s" % (V, list(recs["row_sums"]))
          if paf_ok else "fails at t = %s" % (bad_t,))
    check("the row sums the record declares are the seeds' own",
          [sum(1 if ch == "+" else -1 for ch in s) for s in raw]
          == list(recs["row_sums"]))
    plain, gist, raw2 = build_pair(recs)
    check("assembly is deterministic: the raw seeds used by the builder "
          "are the ones just checked", raw2 == raw,
          "assembled in %.1fs" % (time.time() - t0))
    check("layout N = 4*v = %d, four %d-blocks, no border (s = 0)" % (N, V),
          len(plain) == N and len(gist) == N
          and all(len(r) == N for r in plain))
    dig_p = verify_rows("plain", plain, SHA["plain"])
    dig_g = verify_rows("gist", gist, SHA["gist"])

    rows_o = orientation_switch(plain, V)
    check("H''      differs from P in exactly the twelve off-diagonal "
          "515-blocks", cells_changed(plain, rows_o) == 12 * V * V,
          "%d cells = 12*515^2" % (12 * V * V))
    check("H''      (a) S H'' S == the alternate-orientation sign pattern, "
          "cell by cell", alternate_orientation_identity(plain, rows_o, V),
          "S = diag(1,-1,-1,-1) (x) I_515; no border term at s = 0")
    alt = build_alternate(raw)
    conj = signed_conjugate(rows_o, V)
    check("H''      (b) S H'' S IS the alternate GS array assembled from "
          "the same seeds, cell for cell", conj == alt,
          "the six transposed blocks negated -- NOTE-B.md S1.0")
    dig_a = verify_rows("alt", alt, SHA_ALT)
    check("H''      the alternate array is a DIFFERENT matrix from P and "
          "from H''", alt != plain and alt != rows_o,
          "it is H'' up to the signed conjugation S, which is in the "
          "equivalence group")
    dig_o = verify_rows("orient", rows_o, SHA["orient"])

    built = {"plain": dig_p, "gist": dig_g, "orient": dig_o}
    mats = {"plain": plain, "gist": gist, "orient": rows_o}
    for tag in ("plain", "orient"):
        rows = mats[tag]
        rt = transpose_rows(rows)
        ttag = TRANSPOSE_OF[tag]
        check("%-8s is NOT the original -- the transpose is a different "
              "matrix here" % LABEL[ttag], rt != rows,
              "%d of %d rows differ"
              % (sum(1 for a, b in zip(rt, rows) if a != b), N))
        check("%-8s transposing twice returns the original, cell for cell"
              % LABEL[ttag], transpose_rows(rt) == rows)
        built[ttag] = verify_rows(ttag, rt, SHA[ttag])
        mats[ttag] = rt
    check("the five matrices carry five DISTINCT canonical digests",
          len(set(built.values())) == 5, "%d matrices" % len(built))
    check("the alternate array's digest is a sixth, distinct from all five",
          dig_a not in set(built.values()))

    # ---------------------------------------------------------- clause 1b
    # The dim-V trap on the real objects, while the rows are still in hand.
    print("\n[1b] control C7 -- the dim-V trap on the real 2060 objects")
    dims = {t: dim_V_W(mats[t]) for t in TAGS}
    check("C7  dim W (INVARIANT) is %d on all five 2060 objects, which "
          "clause [4] proves pairwise inequivalent" % (N - 1),
          set(w for _v, w in dims.values()) == {N - 1} and len(dims) == 5,
          "matching invariants prove nothing: this one separates none of "
          "the ten pairs")
    check("C7  dim V (NOT invariant) takes more than one value across them "
          "-- and is worthless",
          len(set(v for v, _w in dims.values())) > 1,
          ", ".join("%s %d" % (LABEL[t], dims[t][0]) for t in TAGS)
          + "  <- do NOT read this as a separation")

    # The rebuilt rows are wanted again only by --full, which recomputes a
    # profile FROM THEM.  On the default path they go now.
    if args.full:
        sel = {"pl-T": "plain-T", "or-T": "orient-T"}
        want = (("plain-T", "orient-T") if args.matrix == "all"
                else (sel[args.matrix],))
        ROWS = {t: mats[t] for t in want}
    else:
        ROWS = {}
    mats.clear()
    del plain, gist, rows_o, alt, conj

    # ---------------------------------------------------------- clause 2
    print("\n[2] the ten banked exact 4-profiles, AUDITED in exact integers")
    print("    (audited, not recomputed: see --full)")
    prof = {}
    for (tag, impl), (name, schema, mname) in sorted(PROFILES.items()):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(tag, impl)] = p
        _audit_bank(tag, impl, p, blob, schema, mname, built)
    for tag in TAGS:
        check("%-7s blas == bits, bin for bin (two independent "
              "implementations)" % LABEL[tag],
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    # ---------------------------------------------------------- clause 3
    print("\n[3] the separation table -- every pair, in both arithmetics")
    keys = sorted(SEP, key=lambda k: (SEP[k][2] == "NEW", k[0], k[1]))
    print("      %-22s %10s %9s %18s  %s"
          % ("pair", "differing", "union", "sum |T4|*delta", "source"))
    for (a, b) in keys:
        want, wantu, src = SEP[(a, b)]
        diff, dsum, m1, u = compare(prof[(a, "blas")], prof[(b, "blas")])
        dbits = compare(prof[(a, "bits")], prof[(b, "bits")])[0]
        print("      %-22s %10d %9d %18d  %s"
              % ("%s vs %s" % (LABEL[a], LABEL[b]), len(diff), u, m1, src))
        check("%-7s vs %-7s  differ in exactly %d of the %d bins of their "
              "union, deltas summing to zero, blas and bits alike"
              % (LABEL[a], LABEL[b], want, wantu),
              len(diff) == want and u == wantu and dsum == 0
              and len(dbits) == want, "%d differing" % len(diff))
        check("%-7s vs %-7s  the FIRST moment, which nothing forces, does "
              "differ" % (LABEL[a], LABEL[b]), m1 != 0,
              "sum |T4|*delta = %d" % m1)

    print("\n    the support structure, asserted")
    Sp = set(prof[("plain", "blas")])
    Sg = set(prof[("gist", "blas")])
    So = set(prof[("orient", "blas")])
    SpT = set(prof[("plain-T", "blas")])
    SoT = set(prof[("orient-T", "blas")])
    check("H'' and P populate the SAME %d bins -- the switch moves counts, "
          "not the support" % NBINS["plain"],
          So == Sp and len(So) == NBINS["plain"])
    check("P^T and (H'')^T populate the SAME %d bins -- another single "
          "support, on the transpose side" % NBINS["plain-T"],
          SpT == SoT and len(SpT) == NBINS["plain-T"])
    check("the two supports are NOT the same: the transposes DROP %d of "
          "the originals' bins" % len(ORIG_ONLY),
          sorted(Sp - SpT) == ORIG_ONLY,
          "counts in P %s"
          % [prof[("plain", "blas")][k] for k in ORIG_ONLY[:6]] + " ...")
    check("and ADD one, |T4| = %s, populated in both transposes and "
          "neither original" % TRANS_ONLY,
          sorted(SpT - Sp) == TRANS_ONLY,
          "counts %s" % [(k, prof[("plain-T", "blas")].get(k, 0),
                          prof[("orient-T", "blas")].get(k, 0))
                         for k in TRANS_ONLY])
    check("so an original-versus-transpose union support is %d bins -- "
          "which is why those comparisons read n of %d"
          % (NUNION_MIXED, NUNION_MIXED),
          len(Sp | SpT) == NUNION_MIXED and len(So | SoT) == NUNION_MIXED)
    check("the union of all five supports is %d bins" % NUNION_ALL,
          len(Sp | Sg | So | SpT | SoT) == NUNION_ALL)
    check("G populates %d bins: it lacks %d of the originals' and has %s "
          "of its own (cert 22's structure, re-asserted)"
          % (NBINS["gist"], len(GIST_MISSING), GIST_ONLY),
          sorted(Sp - Sg) == GIST_MISSING and sorted(Sg - Sp) == GIST_ONLY
          and len(Sp | Sg) == 147, "union support 147 bins")
    check("G's union with a transpose is only %d bins, and G separates in "
          "EVERY one of them" % SEP[("gist", "plain-T")][1],
          len(Sg | SpT) == SEP[("gist", "plain-T")][1]
          and len(Sg | SoT) == SEP[("gist", "orient-T")][1]
          and SEP[("gist", "plain-T")][0] == len(Sg | SpT)
          and SEP[("gist", "orient-T")][0] == len(Sg | SoT),
          "no bin of either union agrees -- the strongest form the "
          "comparison can take")

    print("\n    where the invariant does and does not separate")
    for a in ("plain", "orient"):
        for b in ("plain-T", "orient-T"):
            A, B = prof[(a, "blas")], prof[(b, "blas")]
            ks = sorted(set(A) | set(B))
            ag = [k for k in ks if A.get(k, 0) == B.get(k, 0)]
            check("%-7s vs %-7s  exactly one bin agrees: |T4| = %s, at %d "
                  "counts on each side"
                  % (LABEL[a], LABEL[b], AGREE_MIXED, AGREE_MIXED_COUNT),
                  ag == AGREE_MIXED
                  and A[AGREE_MIXED[0]] == AGREE_MIXED_COUNT
                  and B[AGREE_MIXED[0]] == AGREE_MIXED_COUNT)
    check("G reads %d in that same bin, not %d -- which is exactly why G "
          "against a transpose separates in every bin"
          % (GIST_AT_AGREE, AGREE_MIXED_COUNT),
          prof[("gist", "blas")].get(AGREE_MIXED[0]) == GIST_AT_AGREE)
    check("the EXTREME TAIL separates on every original-versus-transpose "
          "leg: the top bin |T4| = %d is populated in the originals and in "
          "no transpose" % TOP_BIN_ORIG,
          max(Sp) == TOP_BIN_ORIG and max(So) == TOP_BIN_ORIG
          and TOP_BIN_ORIG not in SpT and TOP_BIN_ORIG not in SoT,
          "as at 2060's own row-side legs (cert 22), and unlike 668, 716, "
          "1676 and 1772, where the tail always agreed")
    dTT = compare(prof[("plain-T", "blas")],
                  prof[("orient-T", "blas")])[0]
    if not dTT:                            # guard: an empty diff has no max
        check("P^T vs (H'')^T  the extreme tail does NOT separate them",
              False, "no divergent bin at all -- the separation FAILED")
    else:
        topTT = max(k for k, _p, _q in dTT)
        tailTT = [k for k in sorted(SpT) if k > topTT]
        check("P^T     vs (H'')^T  the extreme tail does NOT separate them: "
              "the top differing bin is |T4| = %d and the %d bins above it "
              "agree" % (TOP_DIFF_TT, TAIL_TT),
              topTT == TOP_DIFF_TT and len(tailTT) == TAIL_TT
              and all(prof[("plain-T", "blas")][k]
                      == prof[("orient-T", "blas")][k] for k in tailTT),
              "up to |T4| = %d" % max(SpT))

    print("\n      the three new legs the theorem rests on (first eight "
          "divergent bins each)")
    for a, b in (("gist", "plain-T"), ("plain", "orient-T"),
                 ("gist", "orient-T")):
        diff = compare(prof[(a, "blas")], prof[(b, "blas")])[0]
        if not diff:
            check("%-7s vs %-7s  has divergent bins at all"
                  % (LABEL[a], LABEL[b]), False, "the separation FAILED")
            continue
        big = max(diff, key=lambda t: abs(t[2] - t[1]))
        print("      %s vs %s: largest |delta| = %d at |T4| = %d "
              "(%.2e of that bin)"
              % (LABEL[a], LABEL[b], abs(big[2] - big[1]), big[0],
                 abs(big[2] - big[1]) / big[1] if big[1] else 0.0))
        print("      %6s %18s %18s %15s" % ("|T4|", LABEL[a], LABEL[b],
                                            "delta"))
        for k, p, q in diff[:8]:
            print("      %6d %18d %18d %+15d" % (k, p, q, q - p))

    # ---------------------------------------------------------- clause 4
    print("\n[4] the TRANSPOSE-EXTENDED verdicts, derived here from those "
          "counts")
    print("    A ~~ B requires A ~ B or A ~ B^T, so each pair needs TWO")
    print("    refutations.  The second disjunct is refuted by EITHER")
    print("    profile(A) != profile(B^T) or profile(A^T) != profile(B),")
    print("    since A ~ B^T iff A^T ~ B; a route is available only where")
    print("    that transpose is banked, and G^T is not.")

    def ndiff(x, y):
        return len(compare(prof[(x, "blas")], prof[(y, "blas")])[0])

    verdicts = {}
    routes = {}
    print("\n      %-16s %8s %10s %12s %10s"
          % ("pair", "A vs B", "A vs B^T", "A^T vs B", "verdict"))
    for i, a in enumerate(CLASSES):
        for b in CLASSES[i + 1:]:
            d0 = ndiff(a, b)
            d1 = ndiff(a, TRANSPOSE_OF[b]) if b in TRANSPOSE_OF else None
            d2 = ndiff(TRANSPOSE_OF[a], b) if a in TRANSPOSE_OF else None
            second = [d for d in (d1, d2) if d is not None and d > 0]
            ok = d0 > 0 and bool(second)
            verdicts[(a, b)] = ok
            routes[(a, b)] = (d1, d2)
            print("      %-16s %8d %10s %12s %10s"
                  % ("%s vs %s" % (LABEL[a], LABEL[b]), d0,
                     ("%d" % d1) if d1 is not None else "n/a",
                     ("%d" % d2) if d2 is not None else "n/a",
                     "SEPARATED" if ok else "OPEN"))
            check("%-7s ~~ %-7s is REFUTED: A ~ B by %d bins, and A ~ B^T "
                  "by %d bins on the '%s' route"
                  % (LABEL[a], LABEL[b], d0, (second[0] if second else 0),
                     "A vs B^T" if (d1 is not None and d1 > 0)
                     else "A^T vs B"), ok)
    check("order 2060: all %d pairs of the %d classes are separated under "
          "the transpose-extended relation"
          % (len(verdicts), len(CLASSES)), all(verdicts.values()))
    check("the ONE n/a in the A vs B^T column is the pair whose B is G, "
          "because G^T was never profiled at 2060",
          sorted(k for k in routes if routes[k][0] is None)
          == [("plain", "gist")],
          "P vs G is carried by its A^T vs B route instead: P^T vs G, "
          "%d bins" % routes[("plain", "gist")][1])
    check("the ONE n/a in the A^T vs B column is the pair whose A is G, for "
          "the same reason",
          sorted(k for k in routes if routes[k][1] is None)
          == [("gist", "orient")],
          "G vs H'' is carried by its A vs B^T route instead: G vs "
          "(H'')^T, %d bins" % routes[("gist", "orient")][0])
    check("so exactly one pair -- P vs H'' -- has BOTH routes banked, and "
          "both separate",
          sorted(k for k in routes if None not in routes[k])
          == [("plain", "orient")]
          and all(d > 0 for d in routes[("plain", "orient")]),
          "%d and %d bins" % routes[("plain", "orient")])
    check("no pair is left with only one refutation: each has at least one "
          "of the two second-refutation routes banked and nonzero",
          all(any(d is not None and d > 0 for d in routes[k])
              for k in routes))

    print("\n    the remark, stated and not headlined: under PLAIN Hadamard")
    print("    equivalence the five matrices profiled here are pairwise")
    print("    inequivalent.")
    pairs5 = [(a, b) for i, a in enumerate(TAGS) for b in TAGS[i + 1:]]
    worst = min(ndiff(a, b) for a, b in pairs5)
    check("all %d pairs of {P, G, H'', P^T, (H'')^T} separate (least "
          "separated pair: %d bins)" % (len(pairs5), worst), worst > 0,
          "at least FIVE classes exhibited by two constructions, the "
          "orientation switch and transposition; the house counts the three")
    check("the plain count is 'at least five EXHIBITED', not eight: G^T is "
          "not profiled at 2060",
          "gist-T" not in {t for (t, _i) in prof},
          "nothing is said here about G vs G^T at this order")

    # ---------------------------------------------------------- clause 5
    print("\n[5] controls")

    print("\n  C1 -- full |T4| profiles of small Hadamard matrices, two ways")
    controls = [(sylvester(3), "Sylvester H(8)"),
                (sylvester(4), "Sylvester H(16)"),
                (paley1(19), "Paley I H(20)")]
    gs_built = {}
    for v, seqs in GS_CONTROLS:
        rows, paf_ok, seeds = gs_control(v, seqs)
        check("C0  %-16s GS condition sum_q PAF_q(t) = 4v*[t=0]"
              % ("GS H(%d)" % (4 * v)), paf_ok, "v = %d" % v)
        gs_built[v] = (rows, seeds)
        controls.append((rows, "GS H(%d)" % (4 * v)))
    small = {}
    for rows, name in controls:
        nn = len(rows)
        check("C0  %-16s is in fact Hadamard" % name, is_hadamard(rows))
        t0 = time.time()
        p1 = profile_straight(rows)
        p2 = profile_pairvec(rows)
        audit(p1, nn, "%s straight" % name)
        audit(p2, nn, "%s pairvec" % name)
        small[name] = (rows, p1)
        check("C1  %-16s straight enumeration == pair-vector route" % name,
              p1 == p2, "%d bins (%.1fs)" % (len(p1), time.time() - t0))
        check("C1  %-16s |T4| = n (mod 8) on every populated bin" % name,
              all(k % 8 == nn % 8 for k in p1), "n mod 8 = %d" % (nn % 8))
        if name.startswith("Sylvester"):
            check("C1  %-16s matches the FORCED Sylvester profile" % name,
                  p1 == sylvester_profile_forced(nn),
                  "predicted %s" % sylvester_profile_forced(nn))

    print("\n  C2 -- the UNBORDERED orientation switch, exercised end to end")
    print("        on GS controls: s = 0 there too, so these are the exact")
    print("        structural analogue of the 2060 switch")
    for v, (rows, seeds) in sorted(gs_built.items()):
        nn = 4 * v
        sw = orientation_switch(rows, v)
        check("C2  GS H(%d) with its twelve off-diagonal blocks negated is "
              "Hadamard" % nn, is_hadamard(sw))
        check("C2  GS H(%d): the switch moved exactly 12*v*v = %d cells"
              % (nn, 12 * v * v), cells_changed(rows, sw) == 12 * v * v)
        check("C2  GS H(%d): (a) S H'' S == the alternate-orientation sign "
              "pattern" % nn, alternate_orientation_identity(rows, sw, v))
        alt_c = gs_control_alternate(v, seeds)
        check("C2  GS H(%d): (b) S H'' S IS the alternate GS array over the "
              "same seeds, cell for cell" % nn,
              signed_conjugate(sw, v) == alt_c)
        check("C2  GS H(%d): the alternate array is Hadamard too" % nn,
              is_hadamard(alt_c))
        check("C2  GS H(%d): the |T4| profile of H'' equals that of the "
              "alternate array -- the conjugation is in the group" % nn,
              profile_straight(sw) == profile_straight(alt_c),
              "so (b) exhibits H'' as the other orientation, not a "
              "coincidence of counts")

    print("\n  C3 -- the transposed-profile route, on matrices small enough")
    print("        for straight O(C(n,4)) enumeration")
    for name in ("Sylvester H(8)", "Sylvester H(16)", "Paley I H(20)"):
        rows, p1 = small[name]
        nn = len(rows)
        rt = transpose_rows(rows)
        p2 = profile_straight(rt)
        audit(p2, nn, name + "^T")
        check("C3  %-16s the transpose is Hadamard; profile(M^T) audited"
              % name, is_hadamard(rt))
        if rt == rows:
            check("C3  %-16s is SYMMETRIC, so profile(M^T) == profile(M) is "
                  "FORCED -- and holds" % name, p1 == p2, "%d bins" % len(p1))
        else:
            print("      [MEAS] C3  %-16s is not symmetric; profile(M^T) == "
                  "profile(M) is %s (measured, never asserted)"
                  % (name, p1 == p2))
    print("      (C3 exists because clause [4] rests on transposed "
          "profiles.")
    print("       At this order the route is emphatically not vacuous: the")
    print("       transposes populate a bin the originals do not and lose")
    print("       twenty-three the originals have -- clause [3].)")

    print("\n  C4 -- the comparator, exercised in the null direction")
    check("C4  every banked profile against itself: 0 differing bins",
          all(compare(prof[k], prof[k])[0] == [] for k in prof),
          "%d profiles" % len(prof))

    print("\n  C5 -- negative control: the second-moment assert must FIRE")
    victim = dict(prof[("plain-T", "blas")])
    kk = sorted(victim)
    victim[kk[0]] -= 1
    victim[kk[-1]] += 1                    # total preserved, m2 moved
    fired = False
    try:
        audit(victim, N, "C5-corrupted")
    except AssertionError as exc:
        fired = "second-moment identity FAILED" in str(exc)
    check("C5  a total-preserving corruption of the new P^T bank is "
          "rejected by audit()", fired,
          "moved one count from |T4|=%d to |T4|=%d" % (kk[0], kk[-1]))
    check("C5  the corrupted profile still totals C(2060,4) -- so ONLY the "
          "second moment could catch it",
          sum(victim.values()) == c_n_4(N))

    print("\n  C7b -- the dim-V trap, demonstrated on Sylvester H(16)")
    import random                                            # noqa: E402
    rng = random.Random(20260903)
    h16 = sylvester(4)
    v0, w0 = dim_V_W(h16)
    flip = [rng.random() < 0.5 for _ in range(16)]
    if all(flip) or not any(flip):          # a global flip is not a witness
        flip[0] = not flip[0]
    h16b = [(r.translate(FLIP) if f else r) for r, f in zip(h16, flip)]
    check("C7b the negated matrix is still Hadamard", is_hadamard(h16b))
    v1, w1 = dim_V_W(h16b)
    check("C7b dim V MOVES under signed row negation", v0 != v1,
          "%d -> %d  (seed 20260903, %d of 16 rows negated)"
          % (v0, v1, sum(flip)))
    check("C7b dim W does NOT move", w0 == w1, "%d -> %d" % (w0, w1))
    check("C7b the |T4| profile does not move either (it is an invariant)",
          profile_straight(h16) == profile_straight(h16b))

    # ---------------------------------------------------------- clause 6
    replayed = []
    if args.full:
        impls = (("blas", "bits") if args.impl == "both" else (args.impl,))
        print("\n[6] --full: RECOMPUTING transposed exact 2060 profiles "
              "here, from the rows")
        print("    clause [1] verified, with numpy (%s path%s, %s).  This is "
              "the step" % (" and ".join(impls),
                            "" if len(impls) == 1 else "s",
                            ", ".join(LABEL[t] for t in sorted(ROWS))))
        print("    that earns the word 'replayed'.  It is also a THIRD "
              "arithmetic route:")
        print("    the banks came from the canonical-split engine, this one "
              "enumerates the")
        print("    U U^T triangle.  Of order 15 h per leg here -- see "
              "NOTES.md.")
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS"):
            os.environ[var] = "3"           # set BEFORE numpy is imported
        # The recompute machinery is cert 06's full_recompute.py, imported by
        # path rather than copied, so the certificates cannot drift apart.
        sys.path.insert(0, os.path.join(ROOT, "certs", "06-668-separation"))
        import full_recompute as FR                          # noqa: E402
        # Smoke the ported numpy paths BEFORE spending the hours on them, on
        # a matrix whose profile is forced (so this is a positive control
        # with a predicted answer) and big enough that the packed path needs
        # more than one uint64 word per row.
        h128 = sylvester(7)
        want128 = sylvester_profile_forced(128)
        for impl in impls:
            got = FR.profile(h128, 128, impl, progress=False)
            audit(got, 128, "full/H128/%s" % impl)
            check("[full] Sylvester H(128) %-4s == the forced profile "
                  "(%d uint64 words/row)" % (impl, (128 + 63) // 64),
                  got == want128)
        for tag in sorted(ROWS):
            rows = ROWS[tag]
            # the rows recomputed from are exactly the rows clause [1]
            # verified and pinned -- re-derived here, not assumed
            check("[full] %-7s the rows about to be enumerated are the ones "
                  "verify.py accepted" % LABEL[tag],
                  rows_sha256(rows) == built[tag] and built[tag] is not None)
            for impl in impls:
                t0 = time.time()
                got = FR.profile(rows, N, impl)
                audit(got, N, "full/%s/%s" % (tag, impl))
                secs = time.time() - t0
                for bimpl in ("blas", "bits"):
                    check("[full] %-7s recomputed %-4s == banked %-4s, bin "
                          "for bin" % (LABEL[tag], impl, bimpl),
                          got == prof[(tag, bimpl)],
                          "%d bins, %.0fs" % (len(got), secs))
                replayed.append("%s/%s" % (LABEL[tag], impl))
            ROWS[tag] = None
            del rows
        ROWS.clear()
    else:
        print("\n[6] --full not requested: the banked profiles were AUDITED, "
              "not recomputed.")
        print("    Nothing above shows that a banked histogram was computed "
              "from the")
        print("    matrices clause [1] rebuilt.  `--full` re-derives them "
              "here with")
        print("    numpy; at order 2060 one leg is of order 15 h (137x the "
              "716 leg this")
        print("    same module took here, the source laboratory's MEASURED "
              "ratio; 197x")
        print("    and about 22 h on the n^5 law quoted elsewhere here) and "
              "the blas")
        print("    path materialises a 2 120 770 x 2060 pair matrix -- "
              "4.4 GB as int8,")
        print("    17.5 GB as its float32 copy -- so NO --full leg has been "
              "run in this")
        print("    repository at 2060, for cert 22's reasons and at cert "
              "22's price.")
        print("    For scale, the two new legs' own seconds fields: "
              "(H'')^T took")
        print("    6 645.1 s blas and 17 685.3 s bits, P^T 6 650.4 s and "
              "17 678.4 s, on")
        print("    16 rented threads.  numpy is finder-side only and never "
              "in the")
        print("    trust chain.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 24: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: with the TRANSPOSE added to the group --")
    print("         ORDER 2060 CARRIES AT LEAST THREE HADAMARD EQUIVALENCE")
    print("         CLASSES (P, G, H''), the statement cert 22 explicitly")
    print("         withheld.  LABEL: PROVEN.")
    print("         Every pair carries two refutations: %d of %d bins and, "
          "by" % (SEP[("plain", "gist")][0], SEP[("plain", "gist")][1]))
    print("         the A^T vs B route, %d of %d (P vs G); %d of %d and %d "
          "of %d"
          % (SEP[("gist", "plain-T")][0], SEP[("gist", "plain-T")][1],
             SEP[("plain", "orient")][0], SEP[("plain", "orient")][1],
             SEP[("plain", "orient-T")][0], SEP[("plain", "orient-T")][1]))
    print("         (P vs H''); %d of %d and %d of %d (G vs H'') -- each an"
          % (SEP[("gist", "orient")][0], SEP[("gist", "orient")][1],
             SEP[("gist", "orient-T")][0], SEP[("gist", "orient-T")][1]))
    print("         exact |T4| 4-profile comparison over all C(2060,4) =")
    print("         %d row 4-subsets in two arithmetics that" % c_n_4(N))
    print("         agree bin for bin, on profiles hitting the second")
    print("         moment %d to the unit." % second_moment_want(N))
    print("         NO SEPARATION STATEMENT IN note/NOTE-B.md IS ROW-SIDE")
    print("         ANY LONGER: 2060 and 1772 were the last two, and")
    print("         cert 25 settled 1772 the same day.")
    if replayed:
        print("         PROFILES: RECOMPUTED in this run from the rebuilt")
        print("         matrices and matched to the bank bin for bin (%s)."
              % ", ".join(replayed))
    else:
        print("         PROFILES: banked exact computation AUDITED in this")
        print("         run, not recomputed.  The two new C(2060,4)")
        print("         enumerations ran in the source laboratory")
        print("         (Hadamard-2060, experiments/inequiv/")
        print("         exact_profile_big.py) on 2026-09-03, on a rented")
        print("         16-thread machine, under a pre-registration flushed")
        print("         before the matrices were built.  `--full` would be")
        print("         the replay; it has NOT been run here at 2060.")
    print("         REMARK: under PLAIN Hadamard equivalence the five")
    print("         matrices profiled here are pairwise inequivalent, so")
    print("         at least FIVE classes are exhibited; the house counts")
    print("         the three that survive either convention.  NOT claimed:")
    print("         anything about G vs G^T at 2060 (G^T was never")
    print("         profiled); any general theorem about orientation; that")
    print("         three is the number of classes at 2060; and NO novelty")
    print("         or priority of any kind -- order 2060 was settled by the")
    print("         publicly posted matrix, which is G itself.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
