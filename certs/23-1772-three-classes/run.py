#!/usr/bin/env python3
"""cert 23 -- order 1772 carries at least THREE Hadamard equivalence
classes: the Lemma-T twist leaves the class at a FOURTH order, and the
orientation switch leaves both classes.

  THEOREM 1 (the twist, a fourth instance).  Let H be the decoded
  (s, i) = (1, 1) bordered Goethals-Seidel record at order 1772
  (data/payload-records.json; cert 01) and let H' be its Lemma-T i = 2
  rebuild (data/twisted-i2-records.json, order 1772; cert 02).  Then H and
  H' are not Hadamard-equivalent: there is no H' = D_r P_r H P_c D_c with
  P permutation matrices and D diagonal +-1.  So the Lemma-T construction
  at psi(rho) = -1 leaves the equivalence class at a FOURTH order -- 668
  (cert 06), 716 (cert 11) and 1676 (cert 20) were the first three.  Four
  instances; NO general theorem is claimed.

  THEOREM 2 (three classes).  Let H'' be H with its TWELVE OFF-DIAGONAL
  CORE BLOCKS NEGATED and the 4-row/4-column border unchanged.  Then H''
  is a Hadamard matrix, and it is Hadamard-inequivalent to each of H and
  H'.  With Theorem 1, ORDER 1772 CARRIES AT LEAST THREE HADAMARD
  EQUIVALENCE CLASSES.  ROW-SIDE ONLY IN THIS CERTIFICATE: the transposed
  profiles at 1772 (H''^T, H'^T) were separate legs of the same campaign
  -- their blas legs in, their bits legs still running when this
  certificate was written, and this repository certifies a profile only
  once blas and bits agree bin for bin -- so THIS certificate banks
  neither and claims NOTHING under the transpose-extended relation.  That
  is exactly cert 20's caveat at 1676, which cert 21 discharged there.
  DISCHARGED HERE TOO, 2026-09-03: the two bits legs landed later the same
  day and CERT 25 carries this theorem across to the transpose-extended
  relation (every original vs every transpose, 91 of 92 bins).  Nothing
  below depends on cert 25, and cert 25 depends on the counts below; with
  it and cert 24, no separation statement in note/NOTE-B.md is row-side
  any longer.

  PROOF (finite, exact).  The multiset {|T4(i,j,k,l)|} over all
  C(1772,4) = 409 422 905 815 row 4-subsets, with
  T4 = sum_c H[i][c]H[j][c]H[k][c]H[l][c], is a Hadamard-equivalence
  invariant (note/NOTE-B.md S3.1, invariant I5: each row negation
  contributes one sign to T4 so |T4| is fixed, each column negation
  contributes d_c^4 = 1, and permutations relabel).  All three profiles
  populate the same 89 bins.  H differs from H' in 57 of them, H'' from H
  in 58, and H'' from H' in 53.  An invariant that differs is a
  separation.  []

  WHAT H'' IS.  S H'' S, with S = diag(I_4, diag(1,-1,-1,-1) (x) I_n) and
  n = 442, is the same seeds and border assembled in the ALTERNATE
  Goethals-Seidel orientation (the six transposed blocks negated) with the
  border signed by superblock (P[a][J] (-1)^[J != 0], Q[I] (-1)^[I != 0])
  -- checked in clause [1] as an identity of sign patterns, exactly as at
  668 (cert 13), 716 (cert 14) and 1676 (cert 20).  So the theorem says:
  at 1772, as at 668, 716, 1676 and 2060 (cert 22), GS orientation is not
  a gauge for Hadamard equivalence, and the orientation switch is a
  different class from the psi(rho) = -1 twist.  Four bordered orders and
  the unbordered founding order; no general statement is made.

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The three C(1772,4) enumerations were not run inside this repository.
  They ran in the source laboratory -- Hadamard-2060, experiments/inequiv/
  exact_profile_big.py (the unchanged engine), 16 threads on a rented
  c2d-highcpu-16 (prof42-1, us-east1-b), 2026-09-02 23:33Z - 2026-09-03
  11:24Z, under the pre-registration experiments/pr0042/REGISTRATION.md,
  flushed 2026-09-02 10:17Z before any matrix it governs was built -- in
  two arithmetics that agree bin for bin, and their output is banked in
  data/.  The matrices themselves were built AND verified at the desk
  (build_matrices.py, manifest.json, this repository's verify/verify.py,
  the pinned ones reproducing their cert pins); the rented machine
  enumerated and nothing else.  The DEFAULT path of this script AUDITS
  all six banks.  Say "banked exact computation AUDITED" of a default
  run; the word "replayed" belongs to --full, which is OFFERED AND PRICED
  BELOW BUT HAS NOT BEEN RUN IN THIS REPOSITORY AT THIS ORDER.

WHAT THIS SCRIPT DOES  (default path: standard library only, ~tens of s)

  (0) Pins the SHA-256 of all six banked files it reads.
  (1) Rebuilds H from the banked record through tools/bordered_gs.py
      (every master-theorem hypothesis re-checked, not merely assembled),
      verifies it with verify/verify.py, pins its canonical digest; forms
      H'' by negating the twelve off-diagonal core blocks, checks that
      exactly 12*n*n entries changed and that S H'' S is the
      alternate-orientation array, verifies H'', pins ITS digest; rebuilds
      H' from the cert-02 twisted record, whose seeds are RE-DERIVED here
      as the psi-twist of the decoded seeds (psi(g) = (-1)^g on Z_442) and
      compared character for character against the bank, verifies it, pins
      its digest.
  (1b) The dim V / dim W trap on the three real 1772 objects.
  (2) Loads and AUDITS the six banked exact 4-profiles (three matrices x
      two implementations) -- it does not recompute them -- pinning each
      file's SHA-256 in code and asserting, in exact integer arithmetic:
      every populated bin = 4 (mod 8) (1772 = 4 mod 8); the bin count
      (89 on each); every bin key canonical in [0, 1772] and every count a
      positive integer; the counts total C(1772,4); the second moment
      equals n^3(n-1)(n-2)/24 (NOTE-B.md S3.1), recomputed here AND
      compared against the two fields the bank declares; the schema, the
      folding, the arithmetic and the matrix name each bank declares; and
      each bank's declared matrix digest compared against the in-process
      digest of the matrix rebuilt in THIS run.  Then blas == bits bin for
      bin on each matrix.
  (3) The three separations -- H vs H' (57 bins), H'' vs H (58), H'' vs H'
      (53) -- computed in BOTH arithmetics, with the support structure
      (identical 89-bin support; the tail that does not separate) and the
      divergent bins printed.
  (4) The theorem, DERIVED IN CODE from those counts, and the row-side
      caveat asserted: no transposed 1772 profile is banked here.
  (5) Controls, all standard library: five small Hadamard matrices
      profiled two ways (one of them the route --full takes); the
      orientation switch applied to GS H(28) and H(36) controls; the
      comparator exercised in the null direction; a total-preserving
      corruption of a banked profile that only the second-moment identity
      can catch, required to be caught; and the dim V / dim W trap on
      Sylvester H(16).
  (6) --full: recompute a profile here from the rows clause [1] verified
      (certs/06-668-separation/full_recompute.py, imported by path, not
      copied, so certs 06/08/11/13/14/15/19/20/21/22/23 cannot drift
      apart), after a smoke test against the forced profile of Sylvester
      H(128).  numpy is imported only under this flag, is finder-side only
      and is never in the trust chain; BLAS threads are capped at three
      before it loads.  PRICE, and why it has not been run: the banked
      legs' own seconds fields are the campaign engine's price at this
      order -- 3 391.6 / 3 398.6 / 3 411.7 s blas and 8 544.6 / 8 509.1 /
      8 586.5 s bits, on SIXTEEN rented threads.  This desk has three, and
      a different engine: a single 1772 leg here is about 68x the 716 leg
      the same module took in this repository (400.3 s, cert 14), i.e. of
      order 7-8 h for one blas matrix at three threads.  The 68x is the
      source laboratory's MEASURED sub-n^5 scaling (its 716->2060 ratio of
      137, PR-0042 REGISTRATION.md Amendment 1, exponent 4.66); on the
      Theta(n^5) law quoted elsewhere in this repository the same leg is
      93x, about 10.3 h, which is the figure cert 14's notes already carry
      for 1772.  Both say hours; the smaller measured one is quoted so the
      price is not inflated in this certificate's own favour.  And
      full_recompute.py materialises a C(n,2) x n pair matrix -- at
      n = 1772 that is 1 569 106 x 1 772, i.e. 2.78 GB as int8 and 11.1 GB
      as the float32 copy the blas path makes, past this desk.  The bits
      path (a 1 569 106 x 28 uint64 packing) is the tractable one here.
      NO --full LEG HAS BEEN RUN IN THIS REPOSITORY AT 1772.

Usage:
  python certs/23-1772-three-classes/run.py
  python certs/23-1772-three-classes/run.py --full --impl bits --matrix orient
  python certs/23-1772-three-classes/run.py --full
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

# The shape of the separations, pinned so a drifting bank cannot quietly
# turn this certificate into a different (weaker) statement.
NBINS = 89
NDIFF = {"decoded": 58, "twisted": 53}     # H'' against each
NDIFF_TWIST = 57                           # H vs H', theorem 1

# Where each comparison stops separating: the top bin whose counts differ,
# and the bins BELOW it that nevertheless agree.  Asserted so that "the
# bulk separates and the tail does not" is a checked statement rather than
# a remark.  Support is identical (89 bins) on all three matrices.
TOP_DIFF = {("decoded", "twisted"): 476,
            ("decoded", "orient"): 476,
            ("twisted", "orient"): 452}
AGREE_BELOW_TOP = {("decoded", "twisted"): [444, 460, 468],
                   ("decoded", "orient"): [460, 468],
                   ("twisted", "orient"): [380, 404, 428, 436]}

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the three matrices (the digest verify/verify.py
# reports).  All three are rebuilt in clause [1] of THIS run: the decoded
# record and the Lemma-T rebuild carry the pins cert 01 and cert 02 already
# fixed at this order, and H'' is formed here from the decoded rows.
SHA_DECODED = \
    "1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2"
SHA_TWISTED = \
    "82484769a28ac93201f208ca3256bfd491f8edc8bb5e0309764c4f609a113378"
SHA_ORIENT = \
    "7f1fae050def5b9b7bdc491c05b24551465cbea8d3d9482a9cd23c98ba607e53"

# SHA-256 of the banked data FILES themselves, so a silently edited bank is
# a hard error rather than a different theorem.  All six are this
# certificate's own.  data/payload-records.json and
# data/twisted-i2-records.json are NOT file-pinned here: they are shared
# with certs 01 and 02, and the binding pin on each is the canonical digest
# of the matrix it produces, checked in clause [1] -- reinforced, for the
# twisted record, by the psi-twist re-derivation that binds it outright to
# payload-records.json.
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
}

# (tag, impl) -> (file, the matrix name the bank must declare)
PROFILES = {
    ("decoded", "blas"): ("data/sep1772-decoded-exact-blas.json",
                          "H_1772-decoded"),
    ("decoded", "bits"): ("data/sep1772-decoded-exact-bits.json",
                          "H_1772-decoded"),
    ("twisted", "blas"): ("data/sep1772-twisted-exact-blas.json",
                          "H_1772-twisted"),
    ("twisted", "bits"): ("data/sep1772-twisted-exact-bits.json",
                          "H_1772-twisted"),
    ("orient", "blas"): ("data/sep1772-orient-exact-blas.json",
                         "H_1772-orient"),
    ("orient", "bits"): ("data/sep1772-orient-exact-bits.json",
                         "H_1772-orient"),
}

LABEL = {"decoded": "H  (decoded)", "twisted": "H' (twist)",
         "orient": "H''(orient)"}

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
# clause 1 -- rebuild, switch, verify, pin
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


def verify_rows(tag, rows, want_sha):
    path = os.path.join(OUT, "H1772_%s.txt" % tag)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify", "verify.py"), path],
        capture_output=True, text=True)
    verdict = ((proc.stdout or proc.stderr).strip().splitlines()
               or ["(no output)"])[-1]
    check("%-8s verify/verify.py" % tag, proc.returncode == 0, verdict[:78])
    dig = rows_sha256(rows)
    check("%-8s canonical sha256 == pin" % tag, dig == want_sha,
          dig[:24] + "...")
    check("%-8s verify.py reports the same digest" % tag, dig in verdict)
    os.remove(path)                    # 3.1 MB apiece; three of them
    return dig


def build_and_verify(tag, rec, want_sha):
    """Re-check the master-theorem hypotheses, assemble, verify, pin."""
    t0 = time.time()
    rep, rows = BGS.check_record(rec)
    if not rows:
        check("%-8s hypotheses" % tag, False, str(rep.get("failures")))
        return rep, None, None
    check("%-8s hypotheses H0-H4, D1/D3/D5, Sigma-bar, compression lemma"
          % tag, rep["hypotheses_ok"] and rep["compression_lemma_crosscheck"],
          "s=%d i=%d w=%d" % (rep["s"], rep["i"], rep["w"]))
    dig = verify_rows(tag, rows, want_sha)
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
                    help="RECOMPUTE an exact 1772 profile here with numpy "
                         "and compare to the banks bin for bin; without it "
                         "the banked profiles are audited, not recomputed. "
                         "HOURS per leg at this order, and the blas path "
                         "wants ~11.1 GB; not yet run in this repository")
    ap.add_argument("--impl", choices=("blas", "bits", "both"),
                    default="both",
                    help="which arithmetic path --full recomputes")
    ap.add_argument("--matrix", choices=("decoded", "twisted", "orient",
                                         "all"), default="all",
                    help="which matrix --full recomputes (default all three; "
                         "at this order one matrix one way is already hours)")
    args = ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 23 -- order 1772: the twist's fourth instance, and a THIRD "
          "class")
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


def _audit_bank(tag, impl, p, blob, mname, built):
    """The forced identities, the declared headers, and the matrix binding
    -- against the in-process digest of the matrix rebuilt in THIS run."""
    lab = "%-8s %-4s" % (tag, impl)
    tot, m2 = audit(p, N, "%s/%s" % (tag, impl))
    check("%s  %d bins, all |T4| = 4 (mod 8)" % (lab, NBINS),
          len(p) == NBINS and all(k % 8 == 4 for k in p))
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
    check("%s  bank declares schema exact-4-profile/1, the |T4| folding, "
          "its arithmetic and its matrix" % lab,
          blob.get("schema") == "exact-4-profile/1"
          and str(blob.get("folded", "")).startswith("|T4|")
          and blob.get("impl") == impl
          and blob.get("matrix") == mname
          and blob.get("producer_filename") == mname + ".txt"
          and bool(blob.get("engine")),
          "matrix=%r engine=%r" % (blob.get("matrix"), blob.get("engine")))
    # Matrix identity.  The declared digest is compared against the digest
    # of the matrix REBUILT IN THIS RUN -- not against a static string --
    # so a bank cannot drift onto a different object.  Keyed on PRESENCE,
    # not truthiness: a declared digest that is empty, null, or not 64 hex
    # digits is a FAILURE.
    declared = (blob.get("matrix_canonical_sha256")
                or blob.get("matrix_sha256"))
    producer = blob.get("matrix_sha256")
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
    print("\n[1] rebuild H, re-check the theorem's hypotheses, verify, pin; "
          "form H'', verify, pin; re-derive and rebuild H'")
    dec_rec = load_record(os.path.join(ROOT, "data", "payload-records.json"),
                          N)
    tw_rec = load_record(os.path.join(ROOT, "data",
                                      "twisted-i2-records.json"), N)
    rep_d, rows_d, dig_d = build_and_verify("decoded", dec_rec, SHA_DECODED)
    if rows_d is None:
        print("\nFATAL: the decoded matrix did not rebuild.")
        return 1
    n, s = int(rep_d["w"]) * int(rep_d["i"]), int(rep_d["s"])
    check("decoded  layout n = |G| = %d, s = %d, N = 4(n+s) = %d"
          % (n, s, 4 * (n + s)), 4 * (n + s) == N and len(rows_d) == N)
    rows_o = orientation_switch(rows_d, n, s)
    check("orient   H'' differs from H in exactly the twelve off-diagonal "
          "core blocks", cells_changed(rows_d, rows_o) == 12 * n * n,
          "%d cells = 12*%d^2" % (12 * n * n, n))
    check("orient   S H'' S == the alternate-orientation array with the "
          "signed border",
          alternate_orientation_identity(rows_d, rows_o, n, s),
          "S = diag(I_4, diag(1,-1,-1,-1) (x) I_%d)" % n)
    dig_o = verify_rows("orient", rows_o, SHA_ORIENT)

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
    rep_t, rows_t, dig_t = build_and_verify("twisted", tw_rec, SHA_TWISTED)
    if rows_t is None:
        print("\nFATAL: the twisted matrix did not rebuild.")
        return 1
    built = {"decoded": dig_d, "orient": dig_o, "twisted": dig_t}
    check("the three matrices carry three DISTINCT canonical digests",
          len(set(built.values())) == 3)

    # ---------------------------------------------------------- clause 1b
    # The dim-V trap on the real objects, while the rows are still in hand.
    print("\n[1b] control C6 -- the dim-V trap on the real 1772 objects")
    dims = {t: dim_V_W(m) for t, m in (("decoded", rows_d),
                                       ("twisted", rows_t),
                                       ("orient", rows_o))}
    check("C6  dim W (INVARIANT) is %d on all three, which clause [3] "
          "proves pairwise inequivalent" % (N - 1),
          set(w for _v, w in dims.values()) == {N - 1},
          "matching invariants prove nothing: this one separates none of "
          "the three pairs")
    check("C6  dim V (NOT invariant) takes more than one value across them "
          "-- and is worthless",
          len(set(v for v, _w in dims.values())) > 1,
          ", ".join("%s %d" % (t, dims[t][0]) for t in
                    ("decoded", "twisted", "orient"))
          + "  <- do NOT read this as a separation")

    # The rebuilt rows are wanted again only by --full, which recomputes a
    # profile FROM THEM.  On the default path they go now.
    if args.full:
        want = (("decoded", "twisted", "orient") if args.matrix == "all"
                else (args.matrix,))
        ROWS = {k: v for k, v in (("decoded", rows_d), ("twisted", rows_t),
                                  ("orient", rows_o)) if k in want}
    else:
        ROWS = {}
    del rows_d, rows_t, rows_o

    # ---------------------------------------------------------- clause 2
    print("\n[2] the six banked exact 4-profiles, AUDITED in exact integers")
    print("    (audited, not recomputed: see --full)")
    prof = {}
    for (tag, impl), (name, mname) in sorted(PROFILES.items()):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(tag, impl)] = p
        _audit_bank(tag, impl, p, blob, mname, built)
    for tag in ("decoded", "twisted", "orient"):
        check("%-8s  blas == bits, bin for bin (two independent "
              "implementations)" % tag,
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    # ---------------------------------------------------------- clause 3
    print("\n[3] the three separations, in both arithmetics")

    def separate(a, b, want):
        la, lb = LABEL[a], LABEL[b]
        pair = "%s vs %s" % (la, lb)
        A, B = prof[(a, "blas")], prof[(b, "blas")]
        diff, dsum, m1, u = compare(A, B)
        dbits = compare(prof[(a, "bits")], prof[(b, "bits")])[0]
        ks = sorted(set(A) | set(B))
        check("%-26s identical support: the same %d populated bins"
              % (pair, NBINS), set(A) == set(B) and u == NBINS)
        check("%-26s bin counts differ in exactly %d of the %d bins, blas "
              "and bits alike" % (pair, want, NBINS),
              len(diff) == want and len(dbits) == want,
              "%d differing (blas), %d (bits)" % (len(diff), len(dbits)))
        check("%-26s the differences sum to zero (both totals are C(1772,4))"
              % pair, dsum == 0)
        check("%-26s the FIRST moment, which nothing forces, does differ"
              % pair, m1 != 0, "sum |T4|*delta = %d" % m1)
        if not diff:                       # guard: an empty diff has no max
            check("%-26s the extreme tail does NOT separate them" % pair,
                  False, "no divergent bin at all -- the separation FAILED")
            return diff
        top = max(k for k, _p, _q in diff)
        tail = [k for k in ks if k > top]
        agree_below = [k for k in ks if k < top and A[k] == B[k]]
        check("%-26s the extreme tail does NOT separate them" % pair,
              top == TOP_DIFF[(a, b)] and all(A[k] == B[k] for k in tail),
              "every bin above |T4| = %d agrees (%d of them, up to %d)"
              % (top, len(tail), ks[-1]))
        check("%-26s the BULK separates: below |T4| = %d only %s agree"
              % (pair, top, agree_below or "no bins"),
              agree_below == AGREE_BELOW_TOP[(a, b)],
              "%d of the %d bins up to |T4| = %d differ"
              % (len(diff), len([k for k in ks if k <= top]), top))
        big = max(diff, key=lambda t: abs(t[2] - t[1]))
        print("      largest |delta| = %d at |T4| = %d, i.e. %.2e of that "
              "bin; first eight divergent bins:"
              % (abs(big[2] - big[1]), big[0],
                 abs(big[2] - big[1]) / big[1]))
        print("      %6s %18s %18s %14s" % ("|T4|", la, lb, "delta"))
        for k, p, q in diff[:8]:
            print("      %6d %18d %18d %+14d" % (k, p, q, q - p))
        return diff

    print("\n    THEOREM 1 -- the Lemma-T twist leaves the class at 1772")
    separate("decoded", "twisted", NDIFF_TWIST)
    print("\n    THEOREM 2 -- the orientation switch leaves both classes")
    separate("decoded", "orient", NDIFF["decoded"])
    separate("twisted", "orient", NDIFF["twisted"])

    # ---------------------------------------------------------- clause 4
    print("\n[4] the theorems, derived here from those counts")

    def ndiff(x, y):
        return len(compare(prof[(x, "blas")], prof[(y, "blas")])[0])

    classes = ("decoded", "twisted", "orient")
    verdicts = {}
    for i, a in enumerate(classes):
        for b in classes[i + 1:]:
            d = ndiff(a, b)
            verdicts[(a, b)] = d > 0
            check("%-12s !~ %-12s : the invariant differs in %d bins"
                  % (LABEL[a], LABEL[b], d), d > 0)
    check("THEOREM 1: the psi(rho) = -1 Lemma-T twist leaves the class at "
          "1772 -- a FOURTH order after 668, 716 and 1676",
          verdicts[("decoded", "twisted")],
          "four instances; no general theorem is claimed")
    check("THEOREM 2: all %d pairs of the %d matrices are separated, so "
          "order 1772 carries AT LEAST THREE classes (row-side)"
          % (len(verdicts), len(classes)), all(verdicts.values()))
    check("ROW-SIDE ONLY IN THIS CERTIFICATE: it opens no transposed 1772 "
          "profile, so nothing HERE is claimed under the transpose-extended "
          "relation",
          not any("-T" in os.path.basename(nm)
                  for (nm, _m) in PROFILES.values()),
          "(H')^T and (H'')^T landed in both arithmetics later on "
          "2026-09-03; CERT 25 banks them and discharges this caveat")

    # ---------------------------------------------------------- clause 5
    print("\n[5] controls")

    print("\n  C1 -- full |T4| profiles of small Hadamard matrices, two ways")
    controls = [(sylvester(3), "Sylvester H(8)"),
                (sylvester(4), "Sylvester H(16)"),
                (paley1(19), "Paley I H(20)")]
    gs_built = {}
    for v, seqs in GS_CONTROLS:
        rows, paf_ok = gs_control(v, seqs)
        check("C0  %-16s GS condition sum_q PAF_q(g) = 4v*[g=0]"
              % ("GS H(%d)" % (4 * v)), paf_ok, "v = %d" % v)
        gs_built[v] = rows
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

    print("\n  C2 -- the orientation switch, exercised where the Hadamard")
    print("        property can be checked by brute force")
    for v, rows in sorted(gs_built.items()):
        nn = 4 * v
        sw = orientation_switch(rows, v, 0)
        check("C2  GS H(%d) with its twelve off-diagonal blocks negated is "
              "Hadamard" % nn, is_hadamard(sw))
        check("C2  GS H(%d): the switch moved exactly 12*v*v = %d cells"
              % (nn, 12 * v * v), cells_changed(rows, sw) == 12 * v * v)
        check("C2  GS H(%d): S H'' S == the alternate-orientation sign "
              "pattern, cell by cell" % nn,
              alternate_orientation_identity(rows, sw, v, 0))

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
    check("C4  a total-preserving corruption of a banked profile is "
          "rejected by audit()", fired,
          "moved one count from |T4|=%d to |T4|=%d" % (kk[0], kk[-1]))
    check("C4  the corrupted profile still totals C(1772,4) -- so ONLY the "
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
        print("\n[6] --full: RECOMPUTING exact 1772 profiles here, from the "
              "rows clause [1]")
        print("    verified, with numpy (%s path%s, %s).  This is the step "
              "that earns" % (" and ".join(impls),
                              "" if len(impls) == 1 else "s",
                              ", ".join(sorted(ROWS))))
        print("    the word 'replayed'.  It is also a THIRD arithmetic "
              "route: the banks")
        print("    came from the canonical-split engine, this one enumerates "
              "the U U^T")
        print("    triangle.  Hours per leg at this order -- see NOTES.md.")
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
        # more than one uint64 word per row -- which the n <= 36 controls
        # above cannot exercise.
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
            check("[full] %-8s the rows about to be enumerated are the ones "
                  "verify.py accepted" % tag,
                  rows_sha256(rows) == built[tag] and built[tag] is not None)
            for impl in impls:
                t0 = time.time()
                got = FR.profile(rows, N, impl)
                audit(got, N, "full/%s/%s" % (tag, impl))
                secs = time.time() - t0
                for bimpl in ("blas", "bits"):
                    check("[full] %-8s recomputed %-4s == banked %-4s, bin "
                          "for bin" % (tag, impl, bimpl),
                          got == prof[(tag, bimpl)],
                          "%d bins, %.0fs" % (len(got), secs))
                replayed.append("%s/%s" % (tag, impl))
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
        print("    numpy.  The campaign's own price at this order is in the "
              "banks'")
        print("    seconds fields: 3 391.6 / 3 398.6 / 3 411.7 s blas and "
              "8 544.6 /")
        print("    8 509.1 / 8 586.5 s bits for H / H' / H'' on SIXTEEN "
              "rented threads.")
        print("    This desk has three and a different engine: one leg is "
              "about 68x")
        print("    the 716 leg this same module took here (400.3 s) on the "
              "measured")
        print("    sub-n^5 scaling, 93x on the n^5 law cert 14 quotes -- "
              "7-10 hours --")
        print("    and the blas path wants about 11.1 GB for its float32 "
              "pair matrix,")
        print("    so NO --full leg has been run in this repository at 1772.")
        print("    numpy is finder-side only and never in the trust chain.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 23: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: at order 1772 the decoded (1,1) record H, its Lemma-T")
    print("         i = 2 rebuild H', and H'' -- the decoded record with its")
    print("         twelve off-diagonal core blocks negated -- are pairwise")
    print("         Hadamard-INEQUIVALENT: %d, %d and %d of the %d bins of "
          "the" % (NDIFF_TWIST, NDIFF["decoded"], NDIFF["twisted"], NBINS))
    print("         exact |T4| 4-profile over all C(1772,4) = %d"
          % c_n_4(N))
    print("         row 4-subsets differ.  Two independent implementations")
    print("         agree bin for bin on each matrix; all three hit the")
    print("         second moment %d to the unit."
          % second_moment_want(N))
    print("         ORDER 1772 CARRIES AT LEAST THREE HADAMARD EQUIVALENCE")
    print("         CLASSES, and the psi(rho) = -1 twist provably leaves the")
    print("         class at a FOURTH order (668, 716, 1676, 1772).")
    print("         LABEL: PROVEN + PROVEN-BY-CERTIFICATE.")
    if replayed:
        print("         PROFILES: RECOMPUTED in this run from the rebuilt")
        print("         matrices and matched to the bank bin for bin (%s)."
              % ", ".join(replayed))
    else:
        print("         PROFILES: banked exact computation AUDITED in this")
        print("         run, not recomputed.  The three C(1772,4)")
        print("         enumerations ran in the source laboratory")
        print("         (Hadamard-2060, experiments/inequiv/")
        print("         exact_profile_big.py) on 2026-09-02/03, on a rented")
        print("         16-thread machine, under a pre-registration flushed")
        print("         before the matrices were built.  `--full` would be")
        print("         the replay; it has NOT been run here at 1772.")
    print("         ROW-SIDE ONLY IN THIS CERTIFICATE: the transposed 1772")
    print("         profiles were a pending leg when it was written, so")
    print("         nothing HERE is claimed under the transpose-extended")
    print("         relation -- cert 20's caveat at 1676, which cert 21")
    print("         discharged there.  DISCHARGED HERE TOO: both legs")
    print("         landed in both arithmetics on 2026-09-03 and CERT 25")
    print("         makes the transpose-extended statement at 1772.")
    print("         REMARK: the GS orientation is not a gauge for Hadamard")
    print("         equivalence at a FOURTH order.  Four instances; no")
    print("         general theorem.  NOT claimed: any general statement")
    print("         about psi(rho) = -1 or about orientation; that three is")
    print("         the number of classes at 1772; anything under the")
    print("         transpose-extended relation; any novelty or priority.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
