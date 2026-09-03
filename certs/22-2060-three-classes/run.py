#!/usr/bin/env python3
"""cert 22 -- order 2060 carries at least THREE Hadamard equivalence
classes: the orientation switch of the plain Goethals-Seidel realisation
is a third class at the founding order.

  THEOREM.  Let P be '2060-plain', the plain Goethals-Seidel array over
  the raw Z_515 seed, and G be '2060-gist', the x104-twisted array that
  is byte-for-byte the publicly posted H(2060) (both cert 07, which
  proves P !~ G by the same invariant used here).  Let H'' be P with its
  TWELVE OFF-DIAGONAL 515-BLOCKS NEGATED -- the orientation switch, in
  its UNBORDERED form.  At 2060 the array is a plain 4x4 Goethals-Seidel
  array of circulant 515-blocks, N = 4*515 with s = 0: there is no
  border, so the switch negates twelve blocks and there is nothing to
  leave alone.  Then H'' is a Hadamard matrix, and it is
  Hadamard-inequivalent to each of P and G: there is no
  H'' = D_r P_r X P_c D_c with P permutation matrices and D diagonal +-1,
  for X either of them.  ORDER 2060 THEREFORE CARRIES AT LEAST THREE
  HADAMARD EQUIVALENCE CLASSES.

  ROW-SIDE ONLY.  The transposed profiles at 2060 -- H''^T and P^T -- are
  separate legs of the same campaign, still running when this certificate
  was written and so not banked here, and NOTHING is claimed under the
  transpose-extended relation at this order.  This is exactly cert 20's
  caveat at 1676, which cert 21 later discharged there; 2060 is now the
  only order in note/NOTE-B.md whose separation statement is row-side.

  PROOF (finite, exact).  The multiset {|T4(i,j,k,l)|} over all
  C(2060,4) = 748 155 697 135 row 4-subsets, with
  T4 = sum_c H[i][c]H[j][c]H[k][c]H[l][c], is a Hadamard-equivalence
  invariant (note/NOTE-B.md S3.1, invariant I5: each row negation
  contributes one sign to T4 so |T4| is fixed, each column negation
  contributes d_c^4 = 1, and permutations relabel).  H'' and P populate
  the SAME 145 bins and 107 of those bin counts differ; H'' and G have
  union support 147 and 146 of those bins differ.  (Cert 07's own leg,
  P against G, is 146 of 147 and is re-derived here.)  An invariant that
  differs is a separation.  []

  WHAT H'' IS.  S H'' S, with S = diag(1,-1,-1,-1) (x) I_515, is the same
  four seeds assembled in the ALTERNATE Goethals-Seidel orientation --
  the six transposed blocks negated (note/NOTE-B.md S1.0).  Certs 13, 14
  and 20 check that identity at 668, 716 and 1676 as a sign pattern, cell
  by cell, and carry a border-signing term (P[a][J](-1)^[J != 0],
  Q[I](-1)^[I != 0]) because those orders are bordered.  At 2060 there is
  no border term and the conjugation is by S alone; and because the seeds
  themselves are in hand here, this certificate checks the identity TWICE:

    (a) as a sign-pattern identity, cell by cell, in the unbordered form
        of certs 13/14/20 -- S[r] S[c] H''[r][c] == alt(r,c) H[r][c] with
        alt(r,c) = -1 exactly on the six transposed blocks;
    (b) by ASSEMBLING the alternate-orientation array directly from the
        same raw seeds -- [[A,B,C,D],[-B,A,-Dt,Ct],[-C,Dt,A,-Bt],
        [-D,-Ct,Bt,A]], the standard array of NOTE-B S1.0 with its six
        transposed blocks negated -- putting it through verify/verify.py,
        and comparing it to S H'' S cell for cell.

  (b) is what (a) can only assert: it exhibits H'' as a signed
  conjugation of an independently constructed Goethals-Seidel array, so
  H'' is not merely "twelve blocks flipped" but the other orientation of
  the founding seed.

  REMARK, as at 668, 716 and 1676.  The GS orientation, a "convention" in
  NOTE-B.md S1.0, is NOT a gauge for Hadamard equivalence at the founding
  order either.  One seed quadruple; three classes; two constructions
  (the plain array and the x104 twist that is the public artifact) plus
  the orientation switch.  Four orders with the same verdict and NO
  general theorem: nothing is claimed at any other order, nothing about
  orientation in general, and nothing to the effect that three is the
  number of classes at 2060.

  PRIORITY.  None is claimed, here or anywhere.  Order 2060 was settled
  by the publicly posted matrix; '2060-gist' IS that matrix, and cert
  07's disclosure language binds -- this certificate counts classes among
  the artifacts banked here and says nothing about existence at 2060 or
  about who first exhibited a Hadamard matrix of this order.

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The H'' enumeration was not run inside this repository.  It ran in the
  source laboratory -- Hadamard-2060, experiments/inequiv/
  exact_profile_big.py (the unchanged engine), 16 threads on a rented
  c2d-highcpu-16 (prof42-2, us-east1-b), 2026-09-02/03, under the
  pre-registration experiments/pr0042/REGISTRATION.md, flushed 10:17Z
  on 2026-09-02 before any matrix it governs was built (Amendment 1
  ~11:05Z added the second instance and re-priced) -- in two arithmetics
  that agree bin for bin, and its output is banked in data/.  The four
  comparison profiles are cert 07's, reused verbatim and re-pinned here.
  The matrix itself was built AND verified at the desk
  (build_matrices.py, manifest.json, this repository's verify/verify.py,
  the plain source reproducing cert 07's pin); the rented machine
  enumerated and nothing else.  The DEFAULT path of this script AUDITS
  all six banks.  Say "banked exact computation AUDITED" of a default
  run; the word "replayed" belongs to --full, which is OFFERED AND PRICED
  BELOW BUT HAS NOT BEEN RUN IN THIS REPOSITORY AT THIS ORDER.

WHAT THIS SCRIPT DOES  (default path: standard library only, ~15 s)

  (0) Pins the SHA-256 of all six banked files it reads.
  (1) Rebuilds P and G from data/sep2060-records.json -- the four
      normalised seeds, the normalising shifts, the twist multiplier and
      the CRT relabelling -- after re-verifying the classical
      Goethals-Seidel condition sum_q PAF_q(t) = 4v*[t = 0] on the raw
      seeds (the s = 0 layer of the master theorem, checked, not
      assumed); verifies both through verify/verify.py and pins their
      canonical digests against cert 07's.  Forms H'' by negating the
      twelve off-diagonal 515-blocks, checks that exactly 12*515^2 cells
      changed, checks the alternate-orientation identity both ways
      (a) and (b) above, verifies the alternate array AND H'' through
      verify/verify.py, and pins H''.
  (1b) The dim V / dim W trap on the three real 2060 objects.
  (2) Loads and AUDITS the six banked exact 4-profiles (three matrices x
      two implementations) -- it does not recompute them -- pinning each
      file's SHA-256 in code and asserting, in exact integer arithmetic:
      every populated bin = 4 (mod 8) (2060 = 4 mod 8); the bin counts
      (145 / 145 / 133); the counts total C(2060,4); the second moment
      equals n^3(n-1)(n-2)/24, recomputed here AND compared against the
      field the bank declares; the schema, folding and arithmetic each
      bank declares; and each bank's declared matrix digest compared
      against the in-process digest of the matrix rebuilt in THIS run.
      Then blas == bits bin for bin on each matrix.  Two schemas are in
      play -- cert 07's sep2060-exact-profile/1 and the campaign's
      exact-4-profile/1 -- and each is checked on its own terms.
  (3) The three separations -- H'' vs P (107 of 145), H'' vs G (146 of
      147), P vs G (146 of 147, cert 07's) -- with the support structure
      and the divergent bins printed.
  (4) The theorem, DERIVED IN CODE from those counts, and the row-side
      caveat asserted: no transposed 2060 profile is banked.
  (5) Controls, all standard library: five small Hadamard matrices
      profiled two ways (one of them the route --full takes); the
      UNBORDERED orientation switch on the GS controls H(28) and H(36),
      which are the exact structural analogue of the 2060 switch (s = 0
      there too), still Hadamard by brute force and satisfying both forms
      of the identity; the comparator in the null direction; a
      total-preserving corruption of the new bank that only the
      second-moment identity can catch; and the dim-V trap on Sylvester
      H(16).
  (6) --full: recompute a profile here from the rows clause [1] verified
      (certs/06-668-separation/full_recompute.py, imported by path, not
      copied, so certs 06/08/11/13/14/15/19/20/21/22 cannot drift apart),
      after a smoke test against the forced profile of Sylvester H(128).
      numpy is imported only under this flag, is finder-side only and is
      never in the trust chain; BLAS threads are capped at three before
      it loads.  PRICE, and why it has not been run: see NOTES.md.  One
      leg here is of order 15 h (the source laboratory's MEASURED
      716->2060 ratio of 137 applied to the 400.3 s 716 leg this same
      module took in this repository) or 22 h on the Theta(n^5) law
      quoted elsewhere here, and the blas route materialises a
      2 120 770 x 2060 pair matrix -- 4.4 GB as int8, 17.5 GB as the
      float32 copy -- far past this desk.  The banked legs' own seconds
      fields price the campaign's engine for comparison: 6 631.0 s blas
      and 17 654.7 s bits at 16 rented threads for H'', against cert 07's
      15 691 s / 30 646 s at three desk threads for P.  NO --full LEG HAS
      BEEN RUN IN THIS REPOSITORY AT 2060.

Usage:
  python certs/22-2060-three-classes/run.py
  python certs/22-2060-three-classes/run.py --full --impl bits --matrix orient
  python certs/22-2060-three-classes/run.py --full
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

# The shape of the separations, pinned so a drifting bank cannot quietly
# turn this certificate into a different (weaker) statement.
NBINS = {"plain": 145, "gist": 133, "orient": 145}
# (a, b) -> (differing bins, union support size, source)
SEP = {
    ("orient", "plain"): (107, 145, "NEW"),
    ("orient", "gist"): (146, 147, "NEW"),
    ("plain", "gist"): (146, 147, "cert 07"),
}
# The support structure, asserted: H'' and P populate the same 145 bins;
# G is missing fourteen of them and populates two of its own.
GIST_MISSING = [940, 972, 988, 1004, 1020, 1036, 1052, 1068, 1084, 1100,
                1116, 1148, 1164, 1180]
GIST_ONLY = [892, 908]
# H'' vs P: where the profiles agree.  Unlike 668/716/1676 the extreme
# tail does NOT agree here -- the top bin |T4| = 1236 is one of the 107.
AGREE_LOW = [644, 852]          # the only bins below 868 that agree
AGREE_BAND = (868, 1180)        # every bin in this closed range agrees
TOP_BIN = 1236                  # differs: 12 (H'') against 6 (P)
COMMON_BIN = 1108               # the one bin where all three agree, at 30

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the matrices (the digest verify/verify.py reports).
# P and G carry the pins cert 07 fixed at this order -- both measured
# before this repository existed, and G's is additionally the digest of
# the publicly posted artifact.  H'' and the alternate-orientation array
# are formed here from P's rows; H''s digest is the desk's, recorded in
# the source laboratory's experiments/pr0042/manifest.json when the matrix
# was built and verified there.
SHA_PLAIN = \
    "510f89b7b423c85da0c7cada52cb0f62e0415d736d2040d6701f66c4e524cf6a"
SHA_GIST = \
    "c7a145d86210740dd3f8ea21ca896a54d6916007a042638f17c8c47f097200f7"
SHA_ORIENT = \
    "4e1891b095b8aafa21176e494038f199b495c96a840bdb003e231c160870b801"
# The alternate-orientation array assembled directly from the seeds in
# clause [1](b).  It is not a campaign object and no profile of it is
# banked; it is pinned because it is built and verified on every run.
SHA_ALT = \
    "40e1d1c8cd40e94016c453f12e520a8518e7d29b773d3adaae3f484eca64398d"

# SHA-256 of the banked data FILES themselves, so a silently edited bank
# is a hard error rather than a different theorem.  The two orient banks
# are this certificate's own (banked 2026-09-03 under cert 22); the four
# others are cert 07's, reused verbatim and re-pinned here.
# data/sep2060-records.json is NOT file-pinned here: it is shared with
# cert 07, which pins it, and the binding pin on it is the canonical
# digest of each matrix it produces, checked in clause [1].
FILE_PINS = {
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

# tag -> impl -> (file, schema, the name the bank must declare)
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
}

LABEL = {"plain": "P (plain)", "gist": "G (gist)", "orient": "H''(orient)"}

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
# clause 1 -- rebuild the pair from the seed record, switch, verify, pin
#
# The block algebra below is cert 07's, restated here rather than imported
# so that this certificate stands on its own file; every convention it
# uses is written out in data/sep2060-records.json's `conventions` block
# and is checked against that block on every run.
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
    lab = "%-6s" % tag
    check("%s verify/verify.py" % lab, proc.returncode == 0, verdict[:70])
    dig = rows_sha256(rows)
    check("%s canonical sha256 == pin" % lab, dig == want_sha,
          dig[:24] + "...")
    check("%s verify.py reports the same digest" % lab, dig in verdict)
    os.remove(path)                    # 4.3 MB apiece; four of them
    return dig


def orientation_switch(rows, n):
    """H'' = H with the twelve off-diagonal blocks negated, UNBORDERED.

    At 2060 there is no border: superblock I = k // n for every row and
    column index k, and every off-diagonal block is negated.  (The
    bordered form of certs 13/14/20 takes the first 4s rows and columns
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
    by cell.  This is certs 13/14/20's clause with the border term gone,
    s = 0 having emptied it."""
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
                    help="RECOMPUTE an exact 2060 profile here with numpy "
                         "and compare to the banks bin for bin; without it "
                         "the banked profiles are audited, not recomputed. "
                         "OF ORDER 15 HOURS per leg at this order, and the "
                         "blas path wants ~17.5 GB; not yet run here")
    ap.add_argument("--impl", choices=("blas", "bits", "both"),
                    default="both",
                    help="which arithmetic path --full recomputes")
    ap.add_argument("--matrix", choices=("plain", "gist", "orient", "all"),
                    default="all",
                    help="which matrix --full recomputes (default all three; "
                         "at this order one matrix one way is already many "
                         "hours)")
    args = ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 22 -- order 2060: the orientation switch is a THIRD class")
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
    lab = "%-11s %-4s" % (LABEL[tag], impl)
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
    # (cert 07's bank predates the campaign's), so each is checked on its
    # own terms and neither is allowed to be silently absent.
    if schema == "exact-4-profile/1":
        hdr_ok = (int(blob.get("second_moment", -1)) == m2
                  and int(blob.get("second_moment_want", -1)) == m2
                  and int(blob.get("total", -1)) == tot
                  and int(blob.get("n", -1)) == N
                  and int(blob.get("C_n_4", -1)) == c_n_4(N))
        named = blob.get("matrix") == name
        extra = "matrix=%r" % (blob.get("matrix"),)
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
    check("%s  bank names the matrix rebuilt in THIS run" % lab,
          is_sha256(declared) and declared == built[tag],
          (declared[:24] + "...") if is_sha256(declared)
          else "declared = %r" % (declared,))
    check("%s  the producer's own matrix_sha256 agrees with it" % lab,
          is_sha256(producer) and producer == declared)


def _body(args, t_start):
    # ---------------------------------------------------------- clause 0
    print("\n[0] the six banked data files, SHA-256 pinned in this script")
    for name, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, name))
        check("%-42s" % name, got == want, got[:24] + "...")

    # ---------------------------------------------------------- clause 1
    print("\n[1] rebuild P and G from the banked seed record, re-verify the")
    print("    Goethals-Seidel condition, verify, pin; form H'' by the")
    print("    UNBORDERED orientation switch and check it two ways")
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
          recs["plain_array"]["canonical_sha256"] == SHA_PLAIN
          and recs["gist_array"]["canonical_sha256"] == SHA_GIST
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
    dig_p = verify_rows("plain", plain, SHA_PLAIN)
    dig_g = verify_rows("gist", gist, SHA_GIST)

    rows_o = orientation_switch(plain, V)
    check("H''    differs from P in exactly the twelve off-diagonal "
          "515-blocks", cells_changed(plain, rows_o) == 12 * V * V,
          "%d cells = 12*515^2" % (12 * V * V))
    check("H''    (a) S H'' S == the alternate-orientation sign pattern, "
          "cell by cell", alternate_orientation_identity(plain, rows_o, V),
          "S = diag(1,-1,-1,-1) (x) I_515; no border term at s = 0")
    alt = build_alternate(raw)
    conj = signed_conjugate(rows_o, V)
    check("H''    (b) S H'' S IS the alternate GS array assembled from the "
          "same seeds, cell for cell", conj == alt,
          "the six transposed blocks negated -- NOTE-B.md S1.0")
    dig_a = verify_rows("alt", alt, SHA_ALT)
    check("H''    the alternate array is a DIFFERENT matrix from P and "
          "from H''", alt != plain and alt != rows_o,
          "it is H'' up to the signed conjugation S, which is in the "
          "equivalence group")
    dig_o = verify_rows("orient", rows_o, SHA_ORIENT)
    built = {"plain": dig_p, "gist": dig_g, "orient": dig_o}
    check("the four matrices carry four DISTINCT canonical digests",
          len({dig_p, dig_g, dig_o, dig_a}) == 4)

    # ---------------------------------------------------------- clause 1b
    # The dim-V trap on the real objects, while the rows are still in hand.
    print("\n[1b] control C7 -- the dim-V trap on the real 2060 objects")
    dims = {t: dim_V_W(m) for t, m in
            (("plain", plain), ("gist", gist), ("orient", rows_o))}
    check("C7  dim W (INVARIANT) is %d on all three 2060 objects, which "
          "clause [3] proves pairwise inequivalent" % (N - 1),
          set(w for _v, w in dims.values()) == {N - 1},
          "matching invariants prove nothing: this one separates none of "
          "the three pairs")
    check("C7  dim V (NOT invariant) takes more than one value across them "
          "-- and is worthless",
          len(set(v for v, _w in dims.values())) > 1,
          ", ".join("%s %d" % (LABEL[t], dims[t][0]) for t in
                    ("plain", "gist", "orient"))
          + "  <- do NOT read this as a separation")

    # The rebuilt rows are wanted again only by --full, which recomputes a
    # profile FROM THEM.  On the default path they go now.
    if args.full:
        want = (("plain", "gist", "orient") if args.matrix == "all"
                else (args.matrix,))
        ROWS = {k: v for k, v in
                (("plain", plain), ("gist", gist), ("orient", rows_o))
                if k in want}
    else:
        ROWS = {}
    del plain, gist, rows_o, alt, conj

    # ---------------------------------------------------------- clause 2
    print("\n[2] the six banked exact 4-profiles, AUDITED in exact integers")
    print("    (audited, not recomputed: see --full)")
    prof = {}
    for (tag, impl), (name, schema, mname) in sorted(PROFILES.items()):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(tag, impl)] = p
        _audit_bank(tag, impl, p, blob, schema, mname, built)
    for tag in ("plain", "gist", "orient"):
        check("%-11s blas == bits, bin for bin (two independent "
              "implementations)" % LABEL[tag],
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    # ---------------------------------------------------------- clause 3
    print("\n[3] the three separations, in both arithmetics")
    print("      %-26s %10s %9s %18s  %s"
          % ("pair", "differing", "union", "sum |T4|*delta", "source"))
    for (a, b) in (("orient", "plain"), ("orient", "gist"),
                   ("plain", "gist")):
        want, wantu, src = SEP[(a, b)]
        diff, dsum, m1, u = compare(prof[(a, "blas")], prof[(b, "blas")])
        dbits = compare(prof[(a, "bits")], prof[(b, "bits")])[0]
        print("      %-26s %10d %9d %18d  %s"
              % ("%s vs %s" % (LABEL[a], LABEL[b]), len(diff), u, m1, src))
        check("%-11s vs %-11s  differ in exactly %d of the %d bins of "
              "their union, deltas summing to zero, blas and bits alike"
              % (LABEL[a], LABEL[b], want, wantu),
              len(diff) == want and u == wantu and dsum == 0
              and len(dbits) == want, "%d differing" % len(diff))
        check("%-11s vs %-11s  the FIRST moment, which nothing forces, "
              "does differ" % (LABEL[a], LABEL[b]), m1 != 0,
              "sum |T4|*delta = %d" % m1)

    print("\n    the support structure, asserted")
    Sp = set(prof[("plain", "blas")])
    Sg = set(prof[("gist", "blas")])
    So = set(prof[("orient", "blas")])
    check("H'' and P populate the SAME %d bins -- the switch moves counts, "
          "not the support" % NBINS["plain"], So == Sp and len(So) == 145)
    check("G populates %d bins: it lacks %d of theirs and has %s of its own"
          % (NBINS["gist"], len(GIST_MISSING), GIST_ONLY),
          sorted(Sp - Sg) == GIST_MISSING and sorted(Sg - Sp) == GIST_ONLY
          and len(Sp | Sg) == 147,
          "union support 147 bins")
    check("|T4| = %d is the ONE bin where all three profiles agree, at 30 "
          "counts apiece" % COMMON_BIN,
          prof[("plain", "blas")][COMMON_BIN]
          == prof[("gist", "blas")][COMMON_BIN]
          == prof[("orient", "blas")][COMMON_BIN] == 30)

    print("\n    H'' vs P: where the invariant does and does not separate")
    D = prof[("orient", "blas")]
    P_ = prof[("plain", "blas")]
    ks = sorted(So)
    agree = [k for k in ks if D[k] == P_[k]]
    band = [k for k in ks if AGREE_BAND[0] <= k <= AGREE_BAND[1]]
    check("the bulk separates: %d of the %d bins below |T4| = %d differ"
          % (sum(1 for k in ks if k < AGREE_BAND[0] and D[k] != P_[k]),
             sum(1 for k in ks if k < AGREE_BAND[0]), AGREE_BAND[0]),
          sorted(k for k in agree if k < AGREE_BAND[0]) == AGREE_LOW,
          "the only agreeing bins below it are |T4| = %s" % AGREE_LOW)
    check("a %d-bin band agrees exactly: every |T4| in [%d, %d]"
          % (len(band), AGREE_BAND[0], AGREE_BAND[1]),
          all(D[k] == P_[k] for k in band) and len(band) == 36)
    check("but the EXTREME TAIL separates here, unlike 668/716/1676: the "
          "top bin |T4| = %d differs" % TOP_BIN,
          D[TOP_BIN] != P_[TOP_BIN] and max(ks) == TOP_BIN,
          "H'' %d against P %d" % (D[TOP_BIN], P_[TOP_BIN]))

    for a, b in (("orient", "plain"), ("orient", "gist")):
        A, B = prof[(a, "blas")], prof[(b, "blas")]
        diff = compare(A, B)[0]
        big = max(diff, key=lambda t: abs(t[2] - t[1]))
        print("\n      %s vs %s: largest |delta| = %d at |T4| = %d "
              "(%.2e of that bin); first eight divergent bins:"
              % (LABEL[a], LABEL[b], abs(big[2] - big[1]), big[0],
                 abs(big[2] - big[1]) / big[1]))
        print("      %6s %18s %18s %15s"
              % ("|T4|", LABEL[a], LABEL[b], "delta"))
        for k, p, q in diff[:8]:
            print("      %6d %18d %18d %+15d" % (k, p, q, q - p))

    # ---------------------------------------------------------- clause 4
    print("\n[4] the theorem, derived here from those counts")

    def ndiff(x, y):
        return len(compare(prof[(x, "blas")], prof[(y, "blas")])[0])

    classes = ("plain", "gist", "orient")
    verdicts = {}
    for i, a in enumerate(classes):
        for b in classes[i + 1:]:
            d = ndiff(a, b)
            verdicts[(a, b)] = d > 0
            check("%-11s !~ %-11s : the invariant differs in %d bins"
                  % (LABEL[a], LABEL[b], d), d > 0)
    check("order 2060: all %d pairs of the %d matrices are separated, so "
          "the order carries AT LEAST THREE classes (row-side)"
          % (len(verdicts), len(classes)), all(verdicts.values()))
    check("ROW-SIDE ONLY: no transposed 2060 profile is banked in this "
          "repository, so nothing is claimed under the transpose-extended "
          "relation here",
          not any("-T" in os.path.basename(n)
                  for (n, _s, _m) in PROFILES.values()),
          "H''^T and P^T are separate, pending legs of the same campaign")

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
    for rows, name in controls:
        nn = len(rows)
        check("C0  %-16s is in fact Hadamard" % name, is_hadamard(rows))
        t0 = time.time()
        p1 = profile_straight(rows)
        p2 = profile_pairvec(rows)
        audit(p1, nn, "%s straight" % name)
        audit(p2, nn, "%s pairvec" % name)
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

    print("\n  C3 -- the comparator, exercised in the null direction")
    check("C3  every banked profile against itself: 0 differing bins",
          all(compare(prof[k], prof[k])[0] == [] for k in prof),
          "%d profiles" % len(prof))

    print("\n  C4 -- negative control: the second-moment assert must FIRE")
    victim = dict(prof[("orient", "blas")])
    kk = sorted(victim)
    victim[kk[0]] -= 1
    victim[kk[-1]] += 1                    # total preserved, m2 moved
    fired = False
    try:
        audit(victim, N, "C4-corrupted")
    except AssertionError as exc:
        fired = "second-moment identity FAILED" in str(exc)
    check("C4  a total-preserving corruption of the new H'' bank is "
          "rejected by audit()", fired,
          "moved one count from |T4|=%d to |T4|=%d" % (kk[0], kk[-1]))
    check("C4  the corrupted profile still totals C(2060,4) -- so ONLY the "
          "second moment could catch it",
          sum(victim.values()) == c_n_4(N))

    print("\n  C5 -- the dim-V trap, demonstrated on Sylvester H(16)")
    import random                                            # noqa: E402
    rng = random.Random(20260903)
    h16 = sylvester(4)
    v0, w0 = dim_V_W(h16)
    flip = [rng.random() < 0.5 for _ in range(16)]
    if all(flip) or not any(flip):          # a global flip is not a witness
        flip[0] = not flip[0]
    h16b = [(r.translate(FLIP) if f else r) for r, f in zip(h16, flip)]
    check("C5  the negated matrix is still Hadamard", is_hadamard(h16b))
    v1, w1 = dim_V_W(h16b)
    check("C5  dim V MOVES under signed row negation", v0 != v1,
          "%d -> %d  (seed 20260903, %d of 16 rows negated)"
          % (v0, v1, sum(flip)))
    check("C5  dim W does NOT move", w0 == w1, "%d -> %d" % (w0, w1))
    check("C5  the |T4| profile does not move either (it is an invariant)",
          profile_straight(h16) == profile_straight(h16b))

    # ---------------------------------------------------------- clause 6
    replayed = []
    if args.full:
        impls = (("blas", "bits") if args.impl == "both" else (args.impl,))
        print("\n[6] --full: RECOMPUTING exact 2060 profiles here, from the "
              "rows clause [1]")
        print("    verified, with numpy (%s path%s, %s).  This is the step "
              "that earns" % (" and ".join(impls),
                              "" if len(impls) == 1 else "s",
                              ", ".join(sorted(ROWS))))
        print("    the word 'replayed'.  It is also a THIRD arithmetic "
              "route: the banks")
        print("    came from the canonical-split engine, this one "
              "enumerates the U U^T")
        print("    triangle.  Of order 15 h per leg here -- see NOTES.md.")
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
            check("[full] %-11s the rows about to be enumerated are the "
                  "ones verify.py accepted" % LABEL[tag],
                  rows_sha256(rows) == built[tag] and built[tag] is not None)
            for impl in impls:
                t0 = time.time()
                got = FR.profile(rows, N, impl)
                audit(got, N, "full/%s/%s" % (tag, impl))
                secs = time.time() - t0
                for bimpl in ("blas", "bits"):
                    check("[full] %-11s recomputed %-4s == banked %-4s, bin "
                          "for bin" % (LABEL[tag], impl, bimpl),
                          got == prof[(tag, bimpl)],
                          "%d bins, %.0fs" % (len(got), secs))
                replayed.append("%s/%s" % (LABEL[tag], impl))
            ROWS[tag] = None
            del rows
        ROWS.clear()
    else:
        print("\n[6] --full not requested: the banked profiles were "
              "AUDITED, not recomputed.")
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
        print("    repository at 2060.  For scale, the banked legs' own "
              "seconds fields:")
        print("    H'' took 6 631.0 s blas and 17 654.7 s bits on 16 rented "
              "threads;")
        print("    cert 07's P took 15 691 s and 30 646 s on three desk "
              "threads.  numpy")
        print("    is finder-side only and never in the trust chain.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 22: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: at order 2060 the plain Goethals-Seidel array P over "
          "the raw")
    print("         Z_515 seed, the x104-twisted array G that is the "
          "publicly")
    print("         posted H(2060), and H'' -- P with its twelve "
          "off-diagonal")
    print("         515-blocks negated, the UNBORDERED orientation switch "
          "-- are")
    print("         pairwise Hadamard-INEQUIVALENT: %d of %d, %d of %d and "
          "%d of %d"
          % (SEP[("orient", "plain")][0], SEP[("orient", "plain")][1],
             SEP[("orient", "gist")][0], SEP[("orient", "gist")][1],
             SEP[("plain", "gist")][0], SEP[("plain", "gist")][1]))
    print("         bins of the exact |T4| 4-profile over all C(2060,4) = "
          "%d" % c_n_4(N))
    print("         row 4-subsets differ.  Two independent implementations")
    print("         agree bin for bin on each matrix; all three hit the")
    print("         second moment %d to the unit."
          % second_moment_want(N))
    print("         ORDER 2060 CARRIES AT LEAST THREE HADAMARD EQUIVALENCE")
    print("         CLASSES.  LABEL: PROVEN + PROVEN-BY-CERTIFICATE.")
    if replayed:
        print("         PROFILES: RECOMPUTED in this run from the rebuilt")
        print("         matrices and matched to the bank bin for bin (%s)."
              % ", ".join(replayed))
    else:
        print("         PROFILES: banked exact computation AUDITED in this")
        print("         run, not recomputed.  The H'' enumeration ran in "
              "the")
        print("         source laboratory (Hadamard-2060, "
              "experiments/inequiv/")
        print("         exact_profile_big.py) on 2026-09-02/03, on a rented")
        print("         16-thread machine, under a pre-registration flushed")
        print("         before the matrix was built; the four comparison")
        print("         profiles are cert 07's.  `--full` would be the")
        print("         replay; it has NOT been run here at 2060.")
    print("         ROW-SIDE ONLY: the transposed profiles at 2060 (H''^T,")
    print("         P^T) are separate, pending legs and nothing here is")
    print("         claimed under the transpose-extended relation -- cert")
    print("         20's caveat at 1676, which cert 21 discharged there.")
    print("         REMARK: the GS orientation is not a gauge for Hadamard")
    print("         equivalence at the founding order either -- one seed,")
    print("         three classes, from two constructions plus the switch.")
    print("         NOT claimed: any general theorem about orientation;")
    print("         that three is the number of classes at 2060; anything")
    print("         under the transpose-extended relation; and NO novelty")
    print("         or priority of any kind -- order 2060 was settled by")
    print("         the publicly posted matrix, which is G itself.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
