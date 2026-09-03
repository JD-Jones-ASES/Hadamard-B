#!/usr/bin/env python3
"""cert 25 -- the three classes at order 1772 hold with the TRANSPOSE added
to the group.

  THE RELATION.  A ~ B is Hadamard equivalence: B = D_r P_r A P_c D_c.
  The TRANSPOSE-EXTENDED relation is  A ~~ B  iff  A ~ B or A ~ B^T.
  Since ~ is symmetric and transposition is an involution,
  A ~ B^T iff A^T ~ B, so refuting A ~~ B takes TWO refutations:
  profile(A) != profile(B), and profile(A) != profile(B^T) -- or
  equivalently profile(A^T) != profile(B), the same statement reached
  from the other side.  The |T4| 4-profile is a Hadamard-equivalence
  invariant (note/NOTE-B.md S3.1, invariant I5), and the transpose of a
  Hadamard matrix is Hadamard, so profile(B^T) is that same invariant
  computed on that other matrix.

  THEOREM (1772).  Let H be the decoded (s,i) = (1,1) bordered
  Goethals-Seidel record at order 1772 (certs 01/23), H' its Lemma-T
  i = 2 rebuild (certs 02/23), and H'' the orientation switch of H --
  the twelve off-diagonal core blocks negated, border unchanged
  (cert 23).  Then H, H' and H'' are pairwise inequivalent UNDER THE
  TRANSPOSE-EXTENDED RELATION.  So ORDER 1772 CARRIES AT LEAST THREE
  HADAMARD EQUIVALENCE CLASSES with the transpose in the group -- the
  statement cert 23 explicitly withheld.

  PROOF (finite, exact).  Every leg is an exact 4-profile comparison
  over all C(1772,4) = 409 422 905 815 row 4-subsets, in two
  arithmetics that agree bin for bin.  H ~~ H' is refuted by H !~ H'
  (57 of 89 bins, cert 23) and H !~ (H')^T (91 of 92, new here);
  H ~~ H'' by H !~ H'' (58, cert 23) and H !~ (H'')^T (91, new);
  H' ~~ H'' by H' !~ H'' (53, cert 23) and H' !~ (H'')^T (91, new) --
  and, the other route, (H')^T !~ H'' (91, new).  An invariant that
  differs is a separation.  []

  REMARK, not the headline.  Under PLAIN Hadamard equivalence the five
  matrices profiled at this order -- H, H', H'', (H')^T, (H'')^T --
  are pairwise inequivalent: all ten profile comparisons separate, the
  least separated pair by 52 bins.  FIVE classes are therefore
  EXHIBITED at 1772 by three constructions and transposition.  The
  house counts THREE, because the transpose-extended relation is the
  one under which a matrix and its transpose are the same object, and
  three is the count that survives either convention.

  WHAT IS NOT PROFILED.  H^T -- the decoded record's own transpose at
  1772 -- was NOT enumerated by the campaign and is NOT banked here.
  This certificate therefore says NOTHING about H vs H^T at 1772
  (cert 19 decided that question at 668 only), the "A^T vs B" route is
  n/a for the two pairs whose A is H, and the plain-equivalence count
  above is "at least FIVE exhibited", not eight.  The theorem is
  unaffected: each pair's two refutations are both in hand.

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The two new C(1772,4) enumerations -- (H')^T and (H'')^T -- were not
  run inside this repository.  They ran in the source laboratory
  (Hadamard-2060, experiments/inequiv/exact_profile_big.py, the engine
  unchanged since the 2060 registration, numpy, 16 threads on a rented
  c2d-highcpu-16, prof42-1, us-east1-b, 2026-09-03) under the
  pre-registration experiments/pr0042/REGISTRATION.md, FLUSHED 10:17
  UTC on 2026-09-02 before any matrix it governs was built, which fixed
  the objects (S2: H_1772-twisted-T and H_1772-orient-T are named
  there, and no third 1772 transpose is), the decision rules (S4:
  "differs in any bin => inequivalent"; "equal in every bin =>
  MEASURED, nothing proved") and the kill criteria (S5: blas != bits in
  any bin is a hard stop) in advance.  The matrices were built AND
  verified at the desk (experiments/pr0042/build_matrices.py, every
  digest in its manifest.json, each through this repository's
  verify/verify.py); the rented machine enumerated and nothing else.
  The DEFAULT path of this script AUDITS all ten banks.  Say "banked
  exact computation AUDITED" of a default run; the word "replayed"
  belongs to --full, which is OFFERED AND PRICED BELOW BUT HAS NOT BEEN
  RUN IN THIS REPOSITORY AT THIS ORDER (as in cert 23).

WHAT THIS SCRIPT DOES  (default path: standard library only, seconds)

  (0) Pins the SHA-256 of all ten banked files it reads.
  (1) Rebuilds H from the banked record through tools/bordered_gs.py
      (every master-theorem hypothesis re-checked, not merely
      assembled); forms H'' by negating the twelve off-diagonal core
      blocks and checks the alternate-orientation identity cell by
      cell; re-derives the twisted seeds as the psi-twist of the
      decoded seeds and rebuilds H'; then TRANSPOSES H' and H''
      in-process, hands all five matrices to verify/verify.py, and
      pins all five canonical digests.  H^T is neither built nor
      pinned: no profile of it exists.
  (1b) Control C6 -- the dim V / dim W trap on the real objects.
  (2) Loads and AUDITS ten banked exact 4-profiles (five matrices x
      two implementations) -- it does not recompute them -- asserting
      in exact integer arithmetic: every populated bin = 4 (mod 8);
      89 populated bins on each original and 88 on each transpose;
      every bin key canonical in [0, 1772] and every count a positive
      integer; total = C(1772,4); second moment =
      n^3(n-1)(n-2)/24, recomputed here AND compared against the
      fields the bank declares; the schema, folding, arithmetic and
      matrix the bank declares; and each bank's declared matrix digest
      against the in-process digest of the matrix rebuilt in THIS run.
      Then blas == bits bin for bin on each of the five matrices.
  (3) All ten pair comparisons, in both arithmetics, each asserted to
      its exact differing-bin count and union support size.
  (4) The TRANSPOSE-EXTENDED verdicts, DERIVED IN CODE from those
      counts: a pair is separated only if both of its refutations are
      nonzero.  Then the plain-equivalence remark, also derived.
  (5) Controls: five small Hadamard matrices profiled two ways, one of
      them the route --full takes; the transposed-profile route on
      matrices small enough for straight enumeration (symmetric
      Sylvester H(8)/H(16), where profile(M^T) = profile(M) is FORCED
      and is checked; Paley H(20), which is not symmetric and is
      MEASURED, never asserted); the support control (at 1772 the
      transposes DROP FOUR of the originals' bins and ADD THREE, where
      at 1676 the swap is two for two and at 668/716 each transpose
      simply drops one); the comparator in the null direction; a
      total-preserving corruption only the second moment can catch;
      and the dim-V trap on Sylvester H(16).
  (6) --full: recompute a transposed 1772 profile here from the rows
      clause [1] verified (certs/06-668-separation/full_recompute.py,
      imported by path, not copied, so certs
      06/08/11/13/14/15/19/20/21/22/23/24/25 cannot drift apart),
      after a smoke test against the forced profile of Sylvester
      H(128).  numpy is imported only under this flag, is finder-side
      only and is never in the trust chain; BLAS threads are capped at
      three before it loads.  PRICE, and why it has not been run:
      exactly cert 23's, since the order and the module are the same
      -- one 1772 leg is about 68x the 716 leg the same module took
      here (400.3 s), i.e. of order 7-8 h for one blas matrix at three
      threads, from the source laboratory's MEASURED sub-n^5 scaling
      (its 716->2060 ratio of 137, PR-0042 REGISTRATION.md
      Amendment 1, exponent 4.66); on the Theta(n^5) law quoted
      elsewhere in this repository the same leg is 93x, about 10.3 h.
      Both say hours; the smaller measured one is quoted so the price
      is not inflated in this certificate's own favour.  And
      full_recompute.py materialises a C(n,2) x n pair matrix -- at
      n = 1772 that is 1 569 106 x 1 772, i.e. 2.78 GB as int8 and
      11.1 GB as the float32 copy the blas path makes, past this desk.
      The bits path (a 1 569 106 x 28 uint64 packing) is the tractable
      one here.  NO --full LEG HAS BEEN RUN IN THIS REPOSITORY AT 1772.

Usage:
  python certs/25-transpose-extended-1772/run.py
  python certs/25-transpose-extended-1772/run.py --full --impl bits --matrix or-T
  python certs/25-transpose-extended-1772/run.py --full
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
sys.path.insert(0, os.path.join(ROOT, "tools"))

import bordered_gs as BGS                                    # noqa: E402

N = 1772
OUT = os.path.join(HERE, "out")

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the five matrices (the digest verify/verify.py
# reports).  ALL FIVE are rebuilt in clause [1] of THIS run: H and H' carry
# the pins certs 01 and 02 fixed at this order, H'' is formed here from H
# (cert 23's pin), and the two transposes are formed here from H' and H''.
# H^T is deliberately absent: the campaign never profiled it, so this
# certificate neither builds nor pins it.
SHA = {
    "dec":  "1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2",
    "tw":   "82484769a28ac93201f208ca3256bfd491f8edc8bb5e0309764c4f609a113378",
    "or":   "7f1fae050def5b9b7bdc491c05b24551465cbea8d3d9482a9cd23c98ba607e53",
    "tw-T": "471f705168cc87b4a1256625ce345a0308f356c5a4e9b1807154163078773238",
    "or-T": "0dffc98fbb6e290a6592bf4e253cf7ab973add8be8b642af83434386a2568864",
}

# SHA-256 of the banked data FILES themselves, so a silently edited bank is
# a hard error rather than a different theorem.  The four *-T-* banks are
# this certificate's own (banked 2026-09-03 under cert 25); the six others
# are cert 23's, reused verbatim and re-pinned here.
# data/payload-records.json and data/twisted-i2-records.json are NOT
# file-pinned: they are shared with certs 01, 02 and 23, and the binding
# pin on each is the canonical digest of the matrix it produces, checked in
# clause [1] -- reinforced, for the twisted record, by the psi-twist
# re-derivation that binds it outright to payload-records.json.
FILE_PINS = {
    "data/sep1772-decoded-exact-blas.json":
        "5985d5f9e1a7ceb54d12dc65e7d5179412eaffb78d94ac3ecf8366db5edbc0d4",
    "data/sep1772-decoded-exact-bits.json":
        "f4b2522d2b8d8ff06ed0195f72051ef28e1f47d7aa691d0afd425a1c84cf98e6",
    "data/sep1772-twisted-exact-blas.json":
        "3b010bf5406916977f060a033990427633bebff0891895b72a0e943afde0f76e",
    "data/sep1772-twisted-exact-bits.json":
        "481ba9f78be1c1d0d2072cfad67598d30c74423ec12e12702350dbba6fb72b35",
    "data/sep1772-orient-exact-blas.json":
        "f58bb4d7db1106950e6506899ffc1329787740cc35b4e0c48574b1d6f06bbea7",
    "data/sep1772-orient-exact-bits.json":
        "57ea8694738d3a6bffe5fa56bf7bab5c6e0d354425b7678f57f6b0d691472f5d",
    "data/sep1772-twisted-T-exact-blas.json":
        "20507070c8cf28702bc9093b2f5324e0736744959cfed3c7f44ea9900df3d101",
    "data/sep1772-twisted-T-exact-bits.json":
        "673d558086ac938e75ba28c2a2b240edf3737482a07714cdaa33aa419b1300dc",
    "data/sep1772-orient-T-exact-blas.json":
        "960c9a1893accfd9e29eb44febd19226340f4db612188090daaa5a50d222e0b7",
    "data/sep1772-orient-T-exact-bits.json":
        "e2901d9976067293c12794b6cf7cf004af9f044acd2b92d7f1abbad25374e94a",
}

PROFILES = {
    ("dec", "blas"):  "data/sep1772-decoded-exact-blas.json",
    ("dec", "bits"):  "data/sep1772-decoded-exact-bits.json",
    ("tw", "blas"):   "data/sep1772-twisted-exact-blas.json",
    ("tw", "bits"):   "data/sep1772-twisted-exact-bits.json",
    ("or", "blas"):   "data/sep1772-orient-exact-blas.json",
    ("or", "bits"):   "data/sep1772-orient-exact-bits.json",
    ("tw-T", "blas"): "data/sep1772-twisted-T-exact-blas.json",
    ("tw-T", "bits"): "data/sep1772-twisted-T-exact-bits.json",
    ("or-T", "blas"): "data/sep1772-orient-T-exact-blas.json",
    ("or-T", "bits"): "data/sep1772-orient-T-exact-bits.json",
}

LABEL = {"dec": "H", "tw": "H'", "or": "H''",
         "tw-T": "(H')^T", "or-T": "(H'')^T"}

TAGS = ("dec", "tw", "or", "tw-T", "or-T")

# The producer's own name for each matrix, as the bank declares it; checked,
# so a bank cannot be silently swapped for another matrix's with the right
# digest edited in.
MATRIX_NAME = {"dec": "H_1772-decoded", "tw": "H_1772-twisted",
               "or": "H_1772-orient", "tw-T": "H_1772-twisted-T",
               "or-T": "H_1772-orient-T"}

# The three originals share one support of 89 bins and the two transposes
# share another of 88 -- and unlike 1676, where the two supports differ by
# two bins each way, at 1772 the transposes DROP FOUR of the originals'
# bins and ADD THREE.  The union is therefore 92.  (At 668 and 716 each
# transpose simply drops one bin: cert 15, C3.)
NBINS = {"dec": 89, "tw": 89, "or": 89, "tw-T": 88, "or-T": 88}
NBINS_ORIG = 89
NBINS_TRANS = 88
NUNION = 92
ORIG_ONLY = [636, 668, 708, 772]      # populated in every original, no transpose
TRANS_ONLY = [620, 724, 916]          # populated in every transpose, no original

# Every pair, with its differing-bin count and the size of the union of the
# two supports.  NEW = a comparison that needs one of the two profiles
# banked under cert 25; cert 23 = re-affirmed from that certificate's banks.
SEP = {
    ("dec", "tw"):    (57, 89, "cert 23"),
    ("dec", "or"):    (58, 89, "cert 23"),
    ("tw", "or"):     (53, 89, "cert 23"),
    ("dec", "tw-T"):  (91, 92, "NEW"),
    ("dec", "or-T"):  (91, 92, "NEW"),
    ("tw", "tw-T"):   (91, 92, "NEW"),
    ("tw", "or-T"):   (91, 92, "NEW"),
    ("or", "tw-T"):   (91, 92, "NEW"),
    ("or", "or-T"):   (91, 92, "NEW"),
    ("tw-T", "or-T"): (52, 88, "NEW"),
}

# The ONE bin that agrees on every original-versus-transpose comparison:
# the isolated top bin |T4| = 1764, at a single 4-subset in all five
# profiles.  It is also the top of every support, so the "extreme tail
# does not separate" statement at this order is exactly that one bin.
TOP_BIN = 1764
TOP_BIN_COUNT = 1
# Where each original-versus-transpose comparison stops separating, and
# where the transpose pair does.
TOP_DIFF_MIXED = 916           # a transpose-only bin, hence a difference
TOP_DIFF_TT = 428
TAIL_TT = 34

# The three classes the theorem separates, and the transpose of each that is
# BANKED.  "dec" has no entry: H^T was never enumerated (see the docstring).
CLASSES = ("dec", "tw", "or")
TRANSPOSE_OF = {"tw": "tw-T", "or": "or-T"}

HEXDIGITS = set("0123456789abcdef")
CH = {1: "+", -1: "-"}
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
    second one is sharp: at n = 1772 it pins a 15-digit number to the unit.
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
# clause 1 -- rebuild, switch, transpose, verify, pin
# ======================================================================

def load_record(path, order):
    with open(path, "r", encoding="ascii") as fh:
        blob = json.load(fh)
    if "orders" in blob:
        return [r for r in blob["orders"] if int(r["order"]) == order][0]
    return blob


def psi_twist_seeds(src):
    """Lemma T's twist, re-derived: x'_q[g] = psi(g) x_q[g] with
    psi(g) = (-1)^g, the character of Z_v whose kernel is the index-2
    subgroup.  Returns the four twisted seed strings."""
    v = int(src["group"][0])
    psi = [1 if g % 2 == 0 else -1 for g in range(v)]
    base = [BGS.signs(x) for x in src["seeds"]]
    return ["".join(CH[psi[g] * base[q][g]] for g in range(v))
            for q in range(len(base))]


def transpose_rows(rows):
    return ["".join(col) for col in zip(*rows)]


def verify_rows(tag, rows):
    """Hand a matrix to the trust chain and pin its canonical digest."""
    path = os.path.join(OUT, "H1772_%s.txt" % tag)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify", "verify.py"), path],
        capture_output=True, text=True)
    verdict = ((proc.stdout or proc.stderr).strip().splitlines()
               or ["(no output)"])[-1]
    lab = "%-7s" % LABEL[tag]
    check("%s verify/verify.py" % lab, proc.returncode == 0, verdict[:70])
    dig = rows_sha256(rows)
    check("%s canonical sha256 == pin" % lab, dig == SHA[tag],
          dig[:24] + "...")
    check("%s verify.py reports the same digest" % lab, dig in verdict)
    os.remove(path)                    # 3.1 MB apiece; five of them
    return dig


def build_and_verify(tag, rec):
    """Re-check the master-theorem hypotheses, assemble, verify, pin."""
    t0 = time.time()
    rep, rows = BGS.check_record(rec)
    lab = "%-7s" % LABEL[tag]
    if not rows:
        check("%s hypotheses" % lab, False, str(rep.get("failures")))
        return rep, None, None
    check("%s hypotheses H0-H4, D1/D3/D5, Sigma-bar, compression lemma"
          % lab, rep["hypotheses_ok"] and rep["compression_lemma_crosscheck"],
          "s=%d i=%d w=%d" % (rep["s"], rep["i"], rep["w"]))
    dig = verify_rows(tag, rows)
    print("       rebuilt and verified in %.1fs" % (time.time() - t0))
    return rep, rows, dig


def orientation_switch(rows, n, s):
    """H'' = H with the twelve off-diagonal core blocks negated.

    Layout (tools/bordered_gs.py): the first 4s rows and columns are the
    border (corner E, row strips); the core follows as four superblocks of n
    rows each, superblock I = (k - 4s) // n.  Border rows/columns untouched.
    """
    b = 4 * s
    Nn = len(rows)
    assert Nn == b + 4 * n

    def blk(k):
        return -1 if k < b else (k - b) // n

    out = []
    for r, row in enumerate(rows):
        I = blk(r)
        if I < 0:
            out.append(row)
            continue
        chars = list(row)
        for c in range(b, Nn):
            if blk(c) != I:
                chars[c] = "+" if chars[c] == "-" else "-"
        out.append("".join(chars))
    return out


def alternate_orientation_identity(rows_H, rows_Hpp, n, s):
    """S H'' S == H_alt, where H_alt negates the six transposed core blocks
    of the standard array (blocks (I,J) with I,J >= 1, I != J), and signs
    the border strips by superblock: P[a][J](-1)^[J != 0], Q[I](-1)^[I != 0].
    A sign-pattern identity, checked cell by cell."""
    b = 4 * s
    Nn = len(rows_H)

    def blk(k):
        return -1 if k < b else (k - b) // n

    def alt_sign(r, c):
        I, J = blk(r), blk(c)
        if I >= 1 and J >= 1 and I != J:
            return -1
        if I == -1 and J >= 1:
            return -1
        if J == -1 and I >= 1:
            return -1
        return 1

    sg = [1 if blk(k) <= 0 else -1 for k in range(Nn)]
    for r in range(Nn):
        hr, pr = rows_H[r], rows_Hpp[r]
        for c in range(Nn):
            lhs = sg[r] * sg[c] * (1 if pr[c] == "+" else -1)
            rhs = alt_sign(r, c) * (1 if hr[c] == "+" else -1)
            if lhs != rhs:
                return False
    return True


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


# Four +-1 sequences of length v with sum_q PAF_q(g) = 4v*[g = 0] -- the
# classical Goethals-Seidel condition, i.e. the degenerate s = 0 layer of
# this repository's master theorem.  The PAF condition is RE-VERIFIED
# below, never assumed, and the assembled matrix is put through the same C0
# Hadamard check as every other control.  These two give profiles with
# three and four populated bins, which Sylvester and Paley matrices do not.
GS_CONTROLS = [
    (7, ["+------", "+++----", "++--+--", "+-+-+--"]),            # -> H(28)
    (9, ["++-------", "++-+-----", "++-+-+---", "++--+-+--"]),    # -> H(36)
]


def gs_control(v, seqs):
    """Assemble the plain GS array over Z_v; returns (rows, paf_ok)."""
    xs = [BGS.signs(s) for s in seqs]
    paf_ok = True
    for g in range(v):
        tot = sum(x[h] * x[(h + g) % v] for x in xs for h in range(v))
        if tot != (4 * v if g == 0 else 0):
            paf_ok = False
    G = BGS.AbelianGroup([v])
    rows = BGS.assemble(G, xs, G.sub_table(), G.idx((v - 1,)),
                        0, 1, [0] * v, None, None, None)
    return rows, paf_ok


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
                    help="RECOMPUTE a transposed exact 1772 profile here "
                         "with numpy and compare to the banks bin for bin; "
                         "without it the banked profiles are audited, not "
                         "recomputed.  HOURS per leg at this order, and the "
                         "blas path wants ~11.1 GB; not yet run in this "
                         "repository")
    ap.add_argument("--impl", choices=("blas", "bits", "both"),
                    default="both",
                    help="which arithmetic path --full recomputes")
    ap.add_argument("--matrix", choices=("tw-T", "or-T", "all"), default="all",
                    help="which NEW profile --full recomputes (default both; "
                         "at this order one matrix one way is already hours)")
    args = ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 25 -- order 1772: the three classes hold with the TRANSPOSE")
    print("           added to the group")
    print("=" * 72)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    # The scratch directory holds one ~3.1 MB generated matrix at a time.
    # Wrapping the body means an exception -- in the optional numpy path or
    # anywhere else -- cannot leave one on disk.
    try:
        rc = _body(args, t_start)
    finally:
        shutil.rmtree(OUT, ignore_errors=True)
    print("generated matrices deleted; nothing left in %s   (%.1fs)"
          % (rel(OUT), time.time() - t_start))
    return rc


def _audit_bank(tag, impl, p, blob, built):
    """The forced identities, the declared headers, and the matrix binding
    -- against the in-process digest of the matrix rebuilt in THIS run."""
    lab = "%-7s %-4s" % (LABEL[tag], impl)
    tot, m2 = audit(p, N, "%s/%s" % (tag, impl))
    check("%s  %d bins, all |T4| = 4 (mod 8)" % (lab, NBINS[tag]),
          len(p) == NBINS[tag] and all(k % 8 == 4 for k in p))
    check("%s  total == C(1772,4) == %d" % (lab, c_n_4(N)), tot == c_n_4(N))
    check("%s  second moment == n^3(n-1)(n-2)/24 == %d"
          % (lab, second_moment_want(N)), m2 == second_moment_want(N))
    # Every count is a positive integer and every key a canonical decimal
    # in [0, 1772]: a bank cannot smuggle a negative or an empty bin past
    # the two identities by cancellation.
    check("%s  every bin key is canonical in [0, %d] and every count a "
          "positive integer" % (lab, N),
          all(isinstance(v, int) and v > 0 for v in p.values())
          and all(0 <= k <= N for k in p),
          "min bin %d, max bin %d, min count %d"
          % (min(p), max(p), min(p.values())))
    check("%s  banked second_moment / second_moment_want / total / n / "
          "C_n_4 headers agree with the recomputation" % lab,
          int(blob.get("second_moment", -1)) == m2
          and int(blob.get("second_moment_want", -1)) == m2
          and int(blob.get("total", -1)) == tot
          and int(blob.get("n", -1)) == N
          and int(blob.get("C_n_4", -1)) == c_n_4(N))
    # The bank must declare WHAT it is: the schema, the folding (the SIGNED
    # T4 histogram is not an invariant -- NOTE-B.md S3.1), which of the two
    # arithmetics produced it, and which matrix.
    check("%s  bank declares schema exact-4-profile/1, the |T4| folding, "
          "its arithmetic and its matrix" % lab,
          blob.get("schema") == "exact-4-profile/1"
          and str(blob.get("folded", "")).startswith("|T4|")
          and blob.get("impl") == impl
          and blob.get("matrix") == MATRIX_NAME[tag]
          and blob.get("producer_filename") == MATRIX_NAME[tag] + ".txt"
          and bool(blob.get("engine")),
          "matrix=%r engine=%r" % (blob.get("matrix"), blob.get("engine")))
    # Matrix identity.  The declared digest is compared against the digest of
    # the matrix REBUILT IN THIS RUN -- not against a static string -- so a
    # bank cannot drift onto a different object.  Keyed on PRESENCE, not
    # truthiness: a declared digest that is empty, null, or not 64 hex digits
    # is a FAILURE.
    declared = (blob.get("matrix_canonical_sha256")
                or blob.get("matrix_sha256"))
    producer = blob.get("matrix_sha256")
    check("%s  bank names the matrix rebuilt in THIS run (%s)"
          % (lab, MATRIX_NAME[tag]),
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
    print("\n[1] rebuild H, re-check the theorem's hypotheses, verify, pin; "
          "form H''; re-derive and rebuild H';")
    print("    then TRANSPOSE H' and H'' in-process and put both through the "
          "trust chain.  H^T is NOT built: the")
    print("    campaign never profiled it, so this certificate neither pins "
          "nor uses it.")
    dec_rec = load_record(os.path.join(ROOT, "data", "payload-records.json"),
                          N)
    tw_rec = load_record(os.path.join(ROOT, "data",
                                      "twisted-i2-records.json"), N)
    rep_d, rows_d, dig_d = build_and_verify("dec", dec_rec)
    if rows_d is None:
        print("\nFATAL: the decoded matrix did not rebuild.")
        return 1
    n, s = int(rep_d["w"]) * int(rep_d["i"]), int(rep_d["s"])
    check("H       layout n = |G| = %d, s = %d, N = 4(n+s) = %d"
          % (n, s, 4 * (n + s)), 4 * (n + s) == N and len(rows_d) == N)

    rows_o = orientation_switch(rows_d, n, s)
    check("H''     differs from H in exactly the twelve off-diagonal core "
          "blocks", cells_changed(rows_d, rows_o) == 12 * n * n,
          "%d cells = 12*%d^2" % (12 * n * n, n))
    check("H''     S H'' S == the alternate-orientation array with the "
          "signed border",
          alternate_orientation_identity(rows_d, rows_o, n, s),
          "S = diag(I_4, diag(1,-1,-1,-1) (x) I_%d)" % n)
    dig_o = verify_rows("or", rows_o)

    # H' is cert 02's shared bank rather than a file of this certificate's
    # own, so its seeds are RE-DERIVED here as the psi-twist of the decoded
    # seeds: the record is bound to payload-records.json by computation.
    check("the two records are genuinely different instances",
          dec_rec["seeds"] != tw_rec["seeds"]
          and int(dec_rec.get("coset_divisors", [1])[0]) == 1
          and int(tw_rec["coset_divisors"][0]) == 2,
          "decoded i=1, rebuild i=2, different seeds")
    check("the twisted seeds are the psi-twist of the decoded seeds, "
          "re-derived here",
          list(tw_rec["seeds"]) == psi_twist_seeds(dec_rec)
          and list(tw_rec["group"]) == list(dec_rec["group"])
          and list(tw_rec["r_shift"]) == list(dec_rec["r_shift"]),
          "psi(g) = (-1)^g on Z_%d; rho = %d is odd, so psi(rho) = -1"
          % (int(dec_rec["group"][0]), int(dec_rec["r_shift"][0])))
    rep_t, rows_t, dig_t = build_and_verify("tw", tw_rec)
    if rows_t is None:
        print("\nFATAL: the twisted matrix did not rebuild.")
        return 1

    built = {"dec": dig_d, "or": dig_o, "tw": dig_t}
    mats = {"dec": rows_d, "or": rows_o, "tw": rows_t}
    for tag in ("tw", "or"):
        rows = mats[tag]
        rt = transpose_rows(rows)
        ttag = TRANSPOSE_OF[tag]
        lab = "%-7s" % LABEL[ttag]
        check("%s is NOT the original -- the transpose is a different matrix "
              "here" % lab, rt != rows,
              "%d of %d rows differ"
              % (sum(1 for a, b in zip(rt, rows) if a != b), N))
        check("%s transposing twice returns the original, cell for cell"
              % lab, transpose_rows(rt) == rows)
        built[ttag] = verify_rows(ttag, rt)
        mats[ttag] = rt
    check("the five matrices carry five DISTINCT canonical digests",
          len(set(built.values())) == 5, "%d matrices" % len(built))

    # ---------------------------------------------------------- clause 1b
    # The dim-V trap on the real objects, while the rows are still in hand.
    print("\n[1b] control C6 -- the dim-V trap on the real 1772 objects")
    dims = {t: dim_V_W(mats[t]) for t in TAGS}
    check("C6  dim W (INVARIANT) is %d on all five 1772 objects, which "
          "clause [4] proves pairwise inequivalent" % (N - 1),
          set(w for _v, w in dims.values()) == {N - 1} and len(dims) == 5,
          "matching invariants prove nothing: this one separates none of the "
          "ten pairs")
    check("C6  dim V (NOT invariant) takes more than one value across them "
          "-- and is worthless",
          len(set(v for v, _w in dims.values())) > 1,
          ", ".join("%s %d" % (LABEL[t], dims[t][0]) for t in TAGS)
          + "  <- do NOT read this as a separation")

    # The rebuilt rows are wanted again only by --full, which recomputes a
    # profile FROM THEM.  On the default path they go now.
    if args.full:
        want = ("tw-T", "or-T") if args.matrix == "all" else (args.matrix,)
        ROWS = {t: mats[t] for t in want}
    else:
        ROWS = {}
    mats.clear()
    del rows_d, rows_t, rows_o

    # ---------------------------------------------------------- clause 2
    print("\n[2] the ten banked exact 4-profiles, AUDITED in exact integers")
    print("    (audited, not recomputed: see --full)")
    prof = {}
    for (tag, impl), name in sorted(PROFILES.items()):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(tag, impl)] = p
        _audit_bank(tag, impl, p, blob, built)
    for tag in TAGS:
        check("%-7s blas == bits, bin for bin (two independent "
              "implementations)" % LABEL[tag],
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    # ---------------------------------------------------------- clause 3
    print("\n[3] the separation table -- every pair, in both arithmetics")
    keys = sorted(SEP, key=lambda k: (SEP[k][2] == "NEW", k[0], k[1]))
    print("      %-22s %10s %9s %16s  %s"
          % ("pair", "differing", "union", "sum |T4|*delta", "source"))
    for (a, b) in keys:
        want, wantu, src = SEP[(a, b)]
        diff, dsum, m1, u = compare(prof[(a, "blas")], prof[(b, "blas")])
        dbits = compare(prof[(a, "bits")], prof[(b, "bits")])[0]
        print("      %-22s %10d %9d %16d  %s"
              % ("%s vs %s" % (LABEL[a], LABEL[b]), len(diff), u, m1, src))
        check("%-7s vs %-7s  differ in exactly %d of the %d bins of their "
              "union, deltas summing to zero, blas and bits alike"
              % (LABEL[a], LABEL[b], want, wantu),
              len(diff) == want and u == wantu and dsum == 0
              and len(dbits) == want, "%d differing" % len(diff))
        check("%-7s vs %-7s  the FIRST moment, which nothing forces, does "
              "differ" % (LABEL[a], LABEL[b]), m1 != 0,
              "sum |T4|*delta = %d" % m1)

    print("\n    on every original-versus-transpose leg exactly ONE bin "
          "agrees, and it")
    print("    is the isolated top bin |T4| = %d, at %d count apiece"
          % (TOP_BIN, TOP_BIN_COUNT))
    for a in CLASSES:
        for b in ("tw-T", "or-T"):
            A, B = prof[(a, "blas")], prof[(b, "blas")]
            ks = sorted(set(A) | set(B))
            ag = [k for k in ks if A.get(k, 0) == B.get(k, 0)]
            check("%-7s vs %-7s  the agreeing bins are exactly [%d]"
                  % (LABEL[a], LABEL[b], TOP_BIN),
                  ag == [TOP_BIN]
                  and A[TOP_BIN] == B[TOP_BIN] == TOP_BIN_COUNT,
                  "%d of %d bins differ" % (len(ks) - len(ag), len(ks)))

    print("\n      the divergent bins of the three new legs the theorem "
          "rests on (first eight each)")
    for a, b in (("dec", "tw-T"), ("dec", "or-T"), ("tw", "or-T")):
        diff = compare(prof[(a, "blas")], prof[(b, "blas")])[0]
        ks = sorted(set(prof[(a, "blas")]) | set(prof[(b, "blas")]))
        if not diff:                       # guard: an empty diff has no max
            check("%-7s vs %-7s  the extreme tail does NOT separate them"
                  % (LABEL[a], LABEL[b]), False,
                  "no divergent bin at all -- the separation FAILED")
            continue
        big = max(diff, key=lambda t: abs(t[2] - t[1]))
        top = max(k for k, _p, _q in diff)
        check("%-7s vs %-7s  the extreme tail does NOT separate them"
              % (LABEL[a], LABEL[b]),
              top == TOP_DIFF_MIXED
              and all(prof[(a, "blas")].get(k, 0)
                      == prof[(b, "blas")].get(k, 0)
                      for k in ks if k > top),
              "every bin above |T4| = %d agrees (%d of them, up to %d)"
              % (top, sum(1 for k in ks if k > top), ks[-1]))
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
    print("    refutations.  A ~ B^T is refuted by profile(A) != "
          "profile(B^T)")
    print("    and, equivalently, by profile(A^T) != profile(B); the second")
    print("    route is available only where A^T is banked, and H^T is not.")

    def ndiff(x, y):
        return len(compare(prof[(x, "blas")], prof[(y, "blas")])[0])

    verdicts = {}
    print("\n      %-16s %8s %10s %12s %10s"
          % ("pair", "A vs B", "A vs B^T", "A^T vs B", "verdict"))
    for i, a in enumerate(CLASSES):
        for b in CLASSES[i + 1:]:
            d0 = ndiff(a, b)
            d1 = ndiff(a, TRANSPOSE_OF[b])
            have_aT = a in TRANSPOSE_OF
            d2 = ndiff(TRANSPOSE_OF[a], b) if have_aT else None
            ok = d0 > 0 and d1 > 0
            verdicts[(a, b)] = ok
            print("      %-16s %8d %10d %12s %10s"
                  % ("%s vs %s" % (LABEL[a], LABEL[b]), d0, d1,
                     ("%d" % d2) if have_aT else "n/a",
                     "SEPARATED" if ok else "OPEN"))
            check("%-4s ~~ %-4s is REFUTED: %d bins against B and %d against "
                  "B^T" % (LABEL[a], LABEL[b], d0, d1), ok)
            if have_aT:
                check("%-4s vs %-4s  the other route, A^T vs B, separates too "
                      "(%d bins)" % (LABEL[a], LABEL[b], d2), d2 > 0)
    check("order 1772: all %d pairs of the %d classes are separated under "
          "the transpose-extended relation"
          % (len(verdicts), len(CLASSES)), all(verdicts.values()))
    check("the 'A^T vs B' column is n/a exactly where A = H, because H^T was "
          "never profiled at 1772",
          sorted(k for k in CLASSES if k not in TRANSPOSE_OF) == ["dec"],
          "each such pair is refuted by its middle column alone")

    print("\n    the remark, stated and not headlined: under PLAIN Hadamard")
    print("    equivalence the five matrices profiled here are pairwise")
    print("    inequivalent.")
    pairs5 = [(a, b) for i, a in enumerate(TAGS) for b in TAGS[i + 1:]]
    worst = min(ndiff(a, b) for a, b in pairs5)
    check("all %d pairs of {H, H', H'', (H')^T, (H'')^T} separate (least "
          "separated pair: %d bins)" % (len(pairs5), worst), worst > 0,
          "at least FIVE classes exhibited by three constructions and "
          "transposition; the house counts the three")
    check("the plain count is 'at least five EXHIBITED', not eight: H^T is "
          "not profiled at 1772",
          "dec-T" not in {t for (t, _i) in prof},
          "nothing is said here about H vs H^T at this order")

    # ---------------------------------------------------------- clause 5
    print("\n[5] controls")

    print("\n  C1 -- full |T4| profiles of small Hadamard matrices, two ways")
    controls = [(sylvester(3), "Sylvester H(8)"),
                (sylvester(4), "Sylvester H(16)"),
                (paley1(19), "Paley I H(20)")]
    for v, seqs in GS_CONTROLS:
        rows, paf_ok = gs_control(v, seqs)
        check("C0  %-16s GS condition sum_q PAF_q(g) = 4v*[g=0]"
              % ("GS H(%d)" % (4 * v)), paf_ok, "v = %d" % v)
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
    # the orientation switch on a control: still Hadamard (it is a GS array
    # in the other orientation), so the switch is exercised at a size where
    # the Hadamard property can be checked by brute force
    rows28, _ = gs_control(7, GS_CONTROLS[0][1])
    sw28 = orientation_switch(rows28, 7, 0)
    check("C0  GS H(28) with its twelve off-diagonal blocks negated is "
          "Hadamard", is_hadamard(sw28))
    check("C0  GS H(28): the switch moved exactly 12*n*n = %d cells"
          % (12 * 7 * 7), cells_changed(rows28, sw28) == 12 * 7 * 7)
    check("C0  GS H(28): S H'' S == the alternate-orientation sign pattern, "
          "cell by cell",
          alternate_orientation_identity(rows28, sw28, 7, 0))

    print("\n  C2 -- the transposed-profile route, on matrices small enough")
    print("        for straight O(C(n,4)) enumeration")
    for name in ("Sylvester H(8)", "Sylvester H(16)", "Paley I H(20)"):
        rows, p1 = small[name]
        nn = len(rows)
        rt = transpose_rows(rows)
        p2 = profile_straight(rt)
        audit(p2, nn, name + "^T")
        check("C2  %-16s the transpose is Hadamard; profile(M^T) audited"
              % name, is_hadamard(rt))
        if rt == rows:
            check("C2  %-16s is SYMMETRIC, so profile(M^T) == profile(M) is "
                  "FORCED -- and holds" % name, p1 == p2, "%d bins" % len(p1))
        else:
            print("      [MEAS] C2  %-16s is not symmetric; profile(M^T) == "
                  "profile(M) is %s (measured, never asserted)"
                  % (name, p1 == p2))
    print("      (C2 exists because clause [4] rests on transposed profiles.")
    print("       At this order the route is not vacuous: the transposes")
    print("       populate three bins the originals do not and lose four")
    print("       the originals have -- C3.)")

    print("\n  C3 -- the transposes are genuinely different objects, and at")
    print("        1772 they DROP FOUR of the originals' bins and ADD THREE")
    print("        (at 1676 the swap is two for two, cert 21; at 668 and 716")
    print("        each transpose simply drops one, cert 15)")
    orig = [set(prof[(t, "blas")]) for t in ("dec", "tw", "or")]
    trans = [set(prof[(t, "blas")]) for t in ("tw-T", "or-T")]
    check("C3  the three originals share one support of %d bins" % NBINS_ORIG,
          orig[0] == orig[1] == orig[2] and len(orig[0]) == NBINS_ORIG)
    check("C3  the two transposes share another support of %d bins"
          % NBINS_TRANS,
          trans[0] == trans[1] and len(trans[0]) == NBINS_TRANS)
    check("C3  |T4| = %s populated in every original and no transpose"
          % ORIG_ONLY, sorted(orig[0] - trans[0]) == ORIG_ONLY,
          "counts %s" % [prof[("dec", "blas")].get(k, 0) for k in ORIG_ONLY])
    check("C3  and |T4| = %s populated in every transpose and no original"
          % TRANS_ONLY, sorted(trans[0] - orig[0]) == TRANS_ONLY,
          "counts %s" % [prof[("or-T", "blas")][k] for k in TRANS_ONLY])
    check("C3  so the union support of the five profiles is %d bins" % NUNION,
          len(orig[0] | trans[0]) == NUNION)
    check("C3  |T4| = %d is the isolated TOP bin of all five supports, at %d "
          "count apiece -- the one bin every original-versus-transpose "
          "comparison agrees on" % (TOP_BIN, TOP_BIN_COUNT),
          all(prof[(t, "blas")].get(TOP_BIN) == TOP_BIN_COUNT for t in TAGS)
          and all(max(prof[(t, "blas")]) == TOP_BIN for t in TAGS),
          "the next bin down is %d in an original and %d in a transpose"
          % (max(k for k in orig[0] if k < TOP_BIN),
             max(k for k in trans[0] if k < TOP_BIN)))
    dTT = compare(prof[("tw-T", "blas")], prof[("or-T", "blas")])[0]
    if not dTT:
        check("C3  (H')^T vs (H'')^T  the extreme tail does NOT separate "
              "them", False, "no divergent bin at all -- FAILED")
    else:
        topTT = max(k for k, _p, _q in dTT)
        tailTT = [k for k in sorted(trans[0]) if k > topTT]
        check("C3  (H')^T vs (H'')^T  the bulk separates and the tail does "
              "not: the top differing bin is |T4| = %d and the %d bins above "
              "it agree" % (TOP_DIFF_TT, TAIL_TT),
              topTT == TOP_DIFF_TT and len(tailTT) == TAIL_TT
              and all(prof[("tw-T", "blas")][k]
                      == prof[("or-T", "blas")][k] for k in tailTT),
              "up to |T4| = %d" % max(trans[0]))

    print("\n  C4 -- the comparator, exercised in the null direction")
    check("C4  every banked profile against itself: 0 differing bins",
          all(compare(prof[k], prof[k])[0] == [] for k in prof),
          "%d profiles" % len(prof))

    print("\n  C5 -- negative control: the second-moment assert must FIRE")
    victim = dict(prof[("or-T", "blas")])
    kk = sorted(victim)
    victim[kk[0]] -= 1
    victim[kk[-1]] += 1                    # total preserved, m2 moved
    fired = False
    try:
        audit(victim, N, "C5-corrupted")
    except AssertionError as exc:
        fired = "second-moment identity FAILED" in str(exc)
    check("C5  a total-preserving corruption of the new (H'')^T bank is "
          "rejected by audit()", fired,
          "moved one count from |T4|=%d to |T4|=%d" % (kk[0], kk[-1]))
    check("C5  the corrupted profile still totals C(1772,4) -- so ONLY the "
          "second moment could catch it",
          sum(victim.values()) == c_n_4(N))

    print("\n  C6b -- the dim-V trap, demonstrated on Sylvester H(16)")
    import random                                            # noqa: E402
    rng = random.Random(20260903)
    h16 = sylvester(4)
    v0, w0 = dim_V_W(h16)
    flip = [rng.random() < 0.5 for _ in range(16)]
    if all(flip) or not any(flip):          # a global flip is not a witness
        flip[0] = not flip[0]
    h16b = [(r.translate(FLIP) if f else r) for r, f in zip(h16, flip)]
    check("C6b the negated matrix is still Hadamard", is_hadamard(h16b))
    v1, w1 = dim_V_W(h16b)
    check("C6b dim V MOVES under signed row negation", v0 != v1,
          "%d -> %d  (seed 20260903, %d of 16 rows negated)"
          % (v0, v1, sum(flip)))
    check("C6b dim W does NOT move", w0 == w1, "%d -> %d" % (w0, w1))
    check("C6b the |T4| profile does not move either (it is an invariant)",
          profile_straight(h16) == profile_straight(h16b))

    # ---------------------------------------------------------- clause 6
    replayed = []
    if args.full:
        impls = (("blas", "bits") if args.impl == "both" else (args.impl,))
        print("\n[6] --full: RECOMPUTING transposed exact 1772 profiles here, "
              "from the rows")
        print("    clause [1] verified, with numpy (%s path%s, %s).  This is "
              "the step" % (" and ".join(impls),
                            "" if len(impls) == 1 else "s",
                            ", ".join(LABEL[t] for t in sorted(ROWS))))
        print("    that earns the word 'replayed'.  It is also a THIRD "
              "arithmetic route:")
        print("    the banks came from the canonical-split engine, this one "
              "enumerates the")
        print("    U U^T triangle.  Hours per leg at this order -- see "
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
        print("    numpy; at order 1772 one leg is hours (about 68x the 716 "
              "leg this")
        print("    same module took here on the measured sub-n^5 scaling, "
              "93x on the")
        print("    n^5 law cert 14 quotes -- 7-10 h) and the blas path wants "
              "about")
        print("    11.1 GB for its float32 pair matrix, so NO --full leg has "
              "been run")
        print("    in this repository at 1772 -- cert 23's position, "
              "unchanged.  For")
        print("    scale, the two new legs' own seconds fields: (H'')^T took "
              "3 401.2 s")
        print("    blas and 8 511.5 s bits, (H')^T 3 405.4 s and 8 501.5 s, "
              "on 16 rented")
        print("    threads.  numpy is finder-side only and never in the "
              "trust chain.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 25: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: with the TRANSPOSE added to the group --")
    print("         ORDER 1772 CARRIES AT LEAST THREE HADAMARD EQUIVALENCE")
    print("         CLASSES (H, H', H''), the statement cert 23 explicitly")
    print("         withheld.  LABEL: PROVEN.")
    print("         Every pair carries two refutations: %d and %d bins"
          % (SEP[("dec", "tw")][0], SEP[("dec", "tw-T")][0]))
    print("         (H vs H'), %d and %d (H vs H''), %d and %d (H' vs H''),"
          % (SEP[("dec", "or")][0], SEP[("dec", "or-T")][0],
             SEP[("tw", "or")][0], SEP[("tw", "or-T")][0]))
    print("         each an exact |T4| 4-profile comparison over all")
    print("         C(1772,4) = %d row 4-subsets in two" % c_n_4(N))
    print("         arithmetics that agree bin for bin, on profiles hitting")
    print("         the second moment %d to the unit."
          % second_moment_want(N))
    print("         NO SEPARATION STATEMENT IN note/NOTE-B.md IS ROW-SIDE")
    print("         ANY LONGER: 1772 and 2060 were the last two, and cert 24")
    print("         settled 2060 the same day.")
    if replayed:
        print("         PROFILES: RECOMPUTED in this run from the rebuilt")
        print("         matrices and matched to the bank bin for bin (%s)."
              % ", ".join(replayed))
    else:
        print("         PROFILES: banked exact computation AUDITED in this")
        print("         run, not recomputed.  The two new C(1772,4)")
        print("         enumerations ran in the source laboratory")
        print("         (Hadamard-2060, experiments/inequiv/")
        print("         exact_profile_big.py) on 2026-09-03, on a rented")
        print("         16-thread machine, under a pre-registration flushed")
        print("         before the matrices were built.  `--full` would be")
        print("         the replay; it has NOT been run here at 1772.")
    print("         REMARK: under PLAIN Hadamard equivalence the five")
    print("         matrices profiled here are pairwise inequivalent, so")
    print("         at least FIVE classes are exhibited; the house counts")
    print("         the three that survive either convention.  NOT claimed:")
    print("         anything about H vs H^T at 1772 (H^T was never")
    print("         profiled); any general theorem about psi(rho) = -1 or")
    print("         about orientation; that three is the number of classes")
    print("         at 1772; any novelty or priority.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
