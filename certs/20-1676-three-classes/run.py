#!/usr/bin/env python3
"""cert 20 -- order 1676 carries at least THREE Hadamard equivalence
classes: the Lemma-T twist leaves the class at a third order, and the
orientation switch leaves both classes.

  THEOREM 1 (the twist, a third instance).  Let H be the decoded
  (s, i) = (1, 1) bordered Goethals-Seidel record at order 1676
  (data/payload-records.json; cert 01) and let H' be its Lemma-T i = 2
  rebuild (data/twisted-i2-records.json, order 1676; cert 02).  Then H and
  H' are not Hadamard-equivalent: there is no H' = D_r P_r H P_c D_c with
  P permutation matrices and D diagonal +-1.  So the Lemma-T construction
  at psi(rho) = -1 leaves the equivalence class at a THIRD order -- 668
  (cert 06) and 716 (cert 11) were the first two.  Three instances; NO
  general theorem is claimed.

  THEOREM 2 (three classes).  Let H'' be H with its TWELVE OFF-DIAGONAL
  CORE BLOCKS NEGATED and the 4-row/4-column border unchanged.  Then H''
  is a Hadamard matrix, and it is Hadamard-inequivalent to each of H and
  H'.  With Theorem 1, ORDER 1676 CARRIES AT LEAST THREE HADAMARD
  EQUIVALENCE CLASSES.  ROW-SIDE ONLY (THIS CERTIFICATE): the transposed
  profiles at 1676 were a separate leg of the same campaign, not complete
  in both arithmetics when this certificate was written, so THIS
  certificate banks neither and claims NOTHING under the
  transpose-extended relation.  DISCHARGED by CERT 21 (2026-09-02), which
  banks both and makes the transpose-extended statement at 1676; nothing
  below depends on it.

  PROOF (finite, exact).  The multiset {|T4(i,j,k,l)|} over all
  C(1676,4) = 327 588 749 775 row 4-subsets, with
  T4 = sum_c H[i][c]H[j][c]H[k][c]H[l][c], is a Hadamard-equivalence
  invariant (note/NOTE-B.md S3.1, invariant I5: each row negation
  contributes one sign to T4 so |T4| is fixed, each column negation
  contributes d_c^4 = 1, and permutations relabel).  All three profiles
  populate the same 142 bins.  H differs from H' in 68 of them, H'' from
  H in 70, and H'' from H' in 66.  An invariant that differs is a
  separation.  []

  WHAT H'' IS.  S H'' S, with S = diag(I_4, diag(1,-1,-1,-1) (x) I_n) and
  n = 418, is the same seeds and border assembled in the ALTERNATE
  Goethals-Seidel orientation (the six transposed blocks negated) with the
  border signed by superblock (P[a][J] (-1)^[J != 0], Q[I] (-1)^[I != 0])
  -- checked in clause [1] as an identity of sign patterns, exactly as at
  668 (cert 13) and 716 (cert 14).  So the theorem says: at 1676, as at
  668 and at 716, GS orientation is not a gauge for Hadamard equivalence,
  and the orientation switch is a different class from the psi(rho) = -1
  twist.  Three orders; no general statement is made.

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The three C(1676,4) enumerations were not run inside this repository.
  They ran in the source laboratory -- Hadamard-2060, experiments/inequiv/
  exact_profile_big.py (the unchanged engine), 16 threads on a rented
  c2d-highcpu-16 (prof42-1, us-east1-b), 2026-09-02 10:36Z-19:55Z, under
  the pre-registration experiments/pr0042/REGISTRATION.md, flushed 10:17Z
  before any matrix it governs was built -- in two arithmetics that agree
  bin for bin, and their output is banked in data/.  The matrices
  themselves were built AND verified at the desk (build_matrices.py,
  manifest.json, this repository's verify/verify.py, the pinned ones
  reproducing their cert pins); the rented machine enumerated and nothing
  else.  The DEFAULT path of this script AUDITS all six banks.  Say
  "banked exact computation AUDITED" of a default run; the word "replayed"
  belongs to --full, which is OFFERED AND PRICED BELOW BUT HAS NOT BEEN
  RUN IN THIS REPOSITORY AT THIS ORDER.

WHAT THIS SCRIPT DOES  (default path: standard library only, ~tens of s)

  (1) Rebuilds H from the banked record through tools/bordered_gs.py
      (every master-theorem hypothesis re-checked, not merely assembled),
      verifies it with verify/verify.py, pins its canonical digest; forms
      H'' by negating the twelve off-diagonal core blocks, checks that
      exactly 12*n*n entries changed and that S H'' S is the
      alternate-orientation array, verifies H'', pins ITS digest; rebuilds
      H' from the cert-02 twisted record, whose seeds are RE-DERIVED here
      as the psi-twist of the decoded seeds (psi(g) = (-1)^g on Z_418) and
      compared character for character against the bank, verifies it, pins
      its digest.
  (2) Loads and AUDITS the six banked exact 4-profiles (three matrices x
      two implementations) -- it does not recompute them -- pinning each
      file's SHA-256 in code and asserting, in exact integer arithmetic:
      every populated bin = 4 (mod 8) (1676 = 4 mod 8); the counts total
      C(1676,4); the second moment equals n^3(n-1)(n-2)/24 (NOTE-B.md
      S3.1), recomputed here AND compared against the field the bank
      declares; each bank's declared matrix digest compared against the
      in-process digest of the matrix rebuilt in THIS run; then blas ==
      bits bin for bin on each matrix.
  (3) The three separations -- H vs H' (68 bins), H'' vs H (70), H'' vs H'
      (66) -- with the divergent bins printed.
  (4) Controls, all standard library: five small Hadamard matrices
      profiled two ways (one of them the route --full takes); the
      orientation switch applied to a GS H(28) control; a total-preserving
      corruption of a banked profile that only the second-moment identity
      can catch, required to be caught; and the dim V / dim W trap, on
      Sylvester H(16) and on the real 1676 objects.
  (5) --full: recompute a profile here from the rows clause [1] verified
      (certs/06-668-separation/full_recompute.py, imported by path, not
      copied, so certs 06/08/11/13/14/15/19/20 cannot drift apart), after
      a smoke test against the forced profile of Sylvester H(128).  numpy
      is imported only under this flag, is finder-side only and is never
      in the trust chain; BLAS threads are capped at three before it
      loads.  PRICE, and why it has not been run: a single 1676 leg is
      about 52x the 716 leg the same module took here (400 s), i.e. of
      order 6-7 h for one blas matrix at three threads.  The 52x is the
      source laboratory's MEASURED sub-n^5 scaling (its 716->2060 ratio
      of 137, PR-0042 REGISTRATION.md Amendment 1, exponent 4.66); on
      the Theta(n^5) law quoted elsewhere in this repository the same
      leg is 70x, about 7.8 h, which is the figure cert 14's notes
      carry.  Both say hours; the smaller measured one is quoted so the
      price is not inflated in this certificate's own favour.  And
      full_recompute.py materialises a C(n,2) x n pair matrix -- at
      n = 1676 that is 2.35 GB as int8 and 9.4 GB as the float32 copy the
      blas path makes, past this desk.  The bits path is the tractable
      one here.  NO --full LEG HAS BEEN RUN IN THIS REPOSITORY AT 1676.

Usage:
  python certs/20-1676-three-classes/run.py
  python certs/20-1676-three-classes/run.py --full --impl bits --matrix orient
  python certs/20-1676-three-classes/run.py --full
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

N = 1676
OUT = os.path.join(HERE, "out")

# The shape of the separations, pinned so a drifting bank cannot quietly
# turn this certificate into a different (weaker) statement.
NBINS = 142
NDIFF = {"decoded": 70, "twisted": 66}     # H'' against each
NDIFF_TWIST = 68                           # H vs H', theorem 1

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the three matrices (the digest verify/verify.py
# reports).  All three are rebuilt in clause [1] of THIS run: the decoded
# record and the Lemma-T rebuild carry the pins cert 01 and cert 02 already
# fixed at this order, and H'' is formed here from the decoded rows.
SHA_DECODED = \
    "8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99"
SHA_TWISTED = \
    "6a4938371ddbe4ad8bd35f21d7e61dad683b15f8f2ec1c88e88ce579c4907405"
SHA_ORIENT = \
    "16d1617cc62532b26c010f3b174c741f0b9388089516759834030d9056a84346"

# SHA-256 of the banked data FILES themselves, so a silently edited bank is
# a hard error rather than a different theorem.  All six are this
# certificate's own.  data/payload-records.json and
# data/twisted-i2-records.json are NOT file-pinned here: they are shared
# with certs 01 and 02, and the binding pin on each is the canonical digest
# of the matrix it produces, checked in clause [1] -- reinforced, for the
# twisted record, by the psi-twist re-derivation that binds it outright to
# payload-records.json.
FILE_PINS = {
    "data/sep1676-decoded-exact-blas.json":
        "57b9a43caf5246de779ad3205a45642c98f7a211be47e4ed12d718fe098781c9",
    "data/sep1676-decoded-exact-bits.json":
        "469e0b0382d479a6d917316246807cadfa1f113bb2bfcd1429ec1712622e7b94",
    "data/sep1676-twisted-exact-blas.json":
        "328ed05c9614a223d95bd35583c83433a8700249c0b50872fbfc8d846e9b5a49",
    "data/sep1676-twisted-exact-bits.json":
        "311ef88606d0543967e2b0cf46aad4f3fb3f1353cc59b018c5369e447c0c2bb1",
    "data/sep1676-orient-exact-blas.json":
        "a83b239695a3bd820de222e829e65a10a5dd66a432858af57cc950eb4ff40be2",
    "data/sep1676-orient-exact-bits.json":
        "af198e51aecd165e8a2a22ee5ece8dfa73d8ddedf314fa94684ee367db14e9d5",
}

PROFILES = {
    ("decoded", "blas"): "data/sep1676-decoded-exact-blas.json",
    ("decoded", "bits"): "data/sep1676-decoded-exact-bits.json",
    ("twisted", "blas"): "data/sep1676-twisted-exact-blas.json",
    ("twisted", "bits"): "data/sep1676-twisted-exact-bits.json",
    ("orient", "blas"): "data/sep1676-orient-exact-blas.json",
    ("orient", "bits"): "data/sep1676-orient-exact-bits.json",
}

HEXDIGITS = set("0123456789abcdef")
CH = {1: "+", -1: "-"}
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
    second one is sharp: at n = 1676 it pins a 15-digit number to the unit.
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
    path = os.path.join(OUT, "H1676_%s.txt" % tag)
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
    os.remove(path)                    # 2.8 MB apiece; three of them
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
    out = []
    for r in rows:
        v = 0
        for c, ch in enumerate(r):
            if ch == "-":
                v |= 1 << c
        out.append(v)
    return out


def f2_rank(vectors):
    basis = []
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return len(basis)


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
    return hist


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
                    help="RECOMPUTE an exact 1676 profile here with numpy "
                         "and compare to the banks bin for bin; without it "
                         "the banked profiles are audited, not recomputed. "
                         "HOURS per leg at this order, and the blas path "
                         "wants ~9.4 GB; not yet run in this repository")
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
    print("cert 20 -- order 1676: the twist's third instance, and a THIRD "
          "class")
    print("=" * 72)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    # The scratch directory holds one ~2.8 MB generated matrix at a time.
    # Wrapping the body means an exception -- in the optional numpy path or
    # anywhere else -- cannot leave one on disk.
    try:
        rc = _body(args, t_start)
    finally:
        shutil.rmtree(OUT, ignore_errors=True)
    print("generated matrices deleted; nothing left in %s   (%.1fs)"
          % (rel(OUT), time.time() - t_start))
    return rc


def _body(args, t_start):
    # ---------------------------------------------------------- clause 0
    print("\n[0] banked data files, SHA-256 pinned in this script")
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
          "%d cells" % (12 * n * n))
    check("orient   S H'' S == the alternate-orientation array with the "
          "signed border",
          alternate_orientation_identity(rows_d, rows_o, n, s))
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

    # ---------------------------------------------------------- clause 1b
    # The dim-V trap on the real objects, while the rows are still in hand.
    print("\n[1b] control C4 -- the dim-V trap on the real 1676 objects")
    vd, wd = dim_V_W(rows_d)
    vt, wt = dim_V_W(rows_t)
    vo, wo = dim_V_W(rows_o)
    check("C4  dim W (INVARIANT) agrees on all three",
          wd == wt == wo == N - 1,
          "decoded %d, twisted %d, orient %d" % (wd, wt, wo))
    check("C4  dim V (NOT invariant) differs across the pair -- and is "
          "worthless", vd != vt,
          "decoded %d, twisted %d, orient %d  <- do NOT read this as a "
          "separation" % (vd, vt, vo))

    # The rebuilt rows are wanted again only by --full, which recomputes a
    # profile FROM THEM.  On the default path they go now.
    if args.full:
        want = (("decoded", "twisted", "orient") if args.matrix == "all"
                else (args.matrix,))
        ROWS = {"decoded": rows_d, "twisted": rows_t, "orient": rows_o}
        ROWS = {k: v for k, v in ROWS.items() if k in want}
    else:
        ROWS = {}
    del rows_d, rows_t, rows_o

    # ---------------------------------------------------------- clause 2
    print("\n[2] the six banked exact 4-profiles, AUDITED in exact integers")
    print("    (audited, not recomputed: see --full)")
    prof = {}
    for (tag, impl), name in sorted(PROFILES.items()):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(tag, impl)] = p
        tot, m2 = audit(p, N, "%s/%s" % (tag, impl))
        check("%-8s %-4s  %d bins, all |T4| = 4 (mod 8)"
              % (tag, impl, len(p)),
              all(k % 8 == 4 for k in p) and len(p) == NBINS)
        check("%-8s %-4s  total == C(1676,4) == %d" % (tag, impl, c_n_4(N)),
              tot == c_n_4(N))
        check("%-8s %-4s  second moment == n^3(n-1)(n-2)/24 == %d"
              % (tag, impl, second_moment_want(N)),
              m2 == second_moment_want(N))
        # every 1676 bank declares the second moment and the header fields
        check("%-8s %-4s  banked second_moment / total / n / C(n,4) headers "
              "agree with the recomputation" % (tag, impl),
              int(blob.get("second_moment", -1)) == m2
              and int(blob.get("total", -1)) == tot
              and int(blob.get("n", -1)) == N
              and int(blob.get("C_n_4", -1)) == c_n_4(N))
        # Matrix identity.  The declared digest is compared against the
        # digest of the matrix REBUILT IN THIS RUN -- not against a static
        # string -- so a bank cannot drift onto a different object.  Keyed
        # on PRESENCE, not truthiness: a declared digest that is empty,
        # null, or not 64 hex digits is a FAILURE.
        declared = (blob.get("matrix_canonical_sha256")
                    or blob.get("matrix_sha256"))
        check("%-8s %-4s  bank names the matrix rebuilt in THIS run"
              % (tag, impl),
              is_sha256(declared) and declared == built[tag],
              (declared[:24] + "...") if is_sha256(declared)
              else "declared = %r" % (declared,))
    for tag in ("decoded", "twisted", "orient"):
        check("%-8s  blas == bits, bin for bin (two independent "
              "implementations)" % tag,
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    # ---------------------------------------------------------- clause 3
    print("\n[3] the three separations")
    D = prof[("decoded", "blas")]
    T = prof[("twisted", "blas")]
    O = prof[("orient", "blas")]

    def separate(label_a, A, label_b, B, want):
        ks = sorted(set(A) | set(B))
        diff = [(k, A.get(k, 0), B.get(k, 0)) for k in ks
                if A.get(k, 0) != B.get(k, 0)]
        check("%-22s identical support: the same %d populated bins"
              % (label_a + " vs " + label_b, NBINS),
              set(A) == set(B) and len(ks) == NBINS)
        check("%-22s bin counts differ in exactly %d of the %d bins"
              % (label_a + " vs " + label_b, want, NBINS),
              len(diff) == want, "%d differing" % len(diff))
        check("%-22s the differences sum to zero (both totals are C(1676,4))"
              % (label_a + " vs " + label_b),
              sum(q - p for _k, p, q in diff) == 0)
        m1 = sum(k * (q - p) for k, p, q in diff)
        check("%-22s the FIRST moment, which nothing forces, does differ"
              % (label_a + " vs " + label_b), m1 != 0,
              "sum |T4|*delta = %d" % m1)
        top = max(k for k, _p, _q in diff)
        tail = sum(1 for k in ks if k > top)
        check("%-22s the extreme tail does NOT separate them"
              % (label_a + " vs " + label_b),
              all(A[k] == B[k] for k in ks if k > top),
              "every bin above |T4| = %d agrees (%d of them, up to %d)"
              % (top, tail, ks[-1]))
        big = max(diff, key=lambda t: abs(t[2] - t[1]))
        print("      largest |delta| = %d at |T4| = %d, i.e. %.2e of that "
              "bin; first eight divergent bins:"
              % (abs(big[2] - big[1]), big[0],
                 abs(big[2] - big[1]) / big[1]))
        print("      %6s %18s %18s %14s" % ("|T4|", label_a, label_b, "delta"))
        for k, p, q in diff[:8]:
            print("      %6d %18d %18d %+14d" % (k, p, q, q - p))
        return diff

    print("\n    THEOREM 1 -- the Lemma-T twist leaves the class at 1676")
    separate("H  (decoded)", D, "H' (twist)", T, NDIFF_TWIST)
    print("\n    THEOREM 2 -- the orientation switch leaves both classes")
    separate("H  (decoded)", D, "H''(orient)", O, NDIFF["decoded"])
    separate("H' (twist)", T, "H''(orient)", O, NDIFF["twisted"])

    # ---------------------------------------------------------- clause 4
    print("\n[4] controls")

    print("\n  C1 -- full |T4| profiles of small Hadamard matrices, two ways")
    controls = [(sylvester(3), "Sylvester H(8)"),
                (sylvester(4), "Sylvester H(16)"),
                (paley1(19), "Paley H(20)")]
    for v, seqs in GS_CONTROLS:
        rows, paf_ok = gs_control(v, seqs)
        check("C0  %-16s GS condition sum_q PAF_q(g) = 4v*[g=0]"
              % ("GS H(%d)" % (4 * v)), paf_ok, "v = %d" % v)
        controls.append((rows, "GS H(%d)" % (4 * v)))
    for rows, name in controls:
        nn = len(rows)
        check("C0  %-16s is in fact Hadamard" % name, is_hadamard(rows))
        t0 = time.time()
        p1 = {k: v for k, v in profile_straight(rows).items() if v}
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
    # the orientation switch on a control: still Hadamard (it is a GS array
    # in the other orientation), so the switch is exercised at a size where
    # the Hadamard property can be checked by brute force
    rows28, _ = gs_control(7, GS_CONTROLS[0][1])
    sw28 = orientation_switch(rows28, 7, 0)
    check("C0  GS H(28) with its twelve off-diagonal blocks negated is "
          "Hadamard", is_hadamard(sw28))
    check("C0  GS H(28): the switch moved exactly 12*n*n = %d cells"
          % (12 * 7 * 7), cells_changed(rows28, sw28) == 12 * 7 * 7)

    print("\n  C2 -- negative control: the second-moment assert must FIRE")
    victim = dict(prof[("orient", "blas")])
    kk = sorted(victim)
    victim[kk[0]] -= 1
    victim[kk[-1]] += 1                    # total preserved, m2 moved
    fired = False
    try:
        audit(victim, N, "C2-corrupted")
    except AssertionError as exc:
        fired = "second-moment identity FAILED" in str(exc)
    check("C2  a total-preserving corruption is rejected by audit()", fired,
          "moved one count from |T4|=%d to |T4|=%d" % (kk[0], kk[-1]))
    check("C2  the corrupted profile still totals C(1676,4) -- so ONLY the "
          "second moment could catch it",
          sum(victim.values()) == c_n_4(N))

    print("\n  C3 -- the dim-V trap, demonstrated on Sylvester H(16)")
    import random                                            # noqa: E402
    rng = random.Random(20260902)
    h16 = sylvester(4)
    v0, w0 = dim_V_W(h16)
    flip = [rng.random() < 0.5 for _ in range(16)]
    if all(flip) or not any(flip):          # a global flip is not a witness
        flip[0] = not flip[0]
    FL = str.maketrans("+-", "-+")
    h16b = [(r.translate(FL) if f else r) for r, f in zip(h16, flip)]
    check("C3  the negated matrix is still Hadamard", is_hadamard(h16b))
    v1, w1 = dim_V_W(h16b)
    check("C3  dim V MOVES under signed row negation", v0 != v1,
          "%d -> %d  (seed 20260902, %d of 16 rows negated)"
          % (v0, v1, sum(flip)))
    check("C3  dim W does NOT move", w0 == w1, "%d -> %d" % (w0, w1))
    check("C3  the |T4| profile does not move either (it is an invariant)",
          profile_straight(h16) == profile_straight(h16b))

    # ---------------------------------------------------------- clause 5
    replayed = []
    if args.full:
        impls = (("blas", "bits") if args.impl == "both" else (args.impl,))
        print("\n[5] --full: RECOMPUTING exact 1676 profiles here, from the "
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
        print("\n[5] --full not requested: the banked profiles were AUDITED, "
              "not recomputed.")
        print("    Nothing above shows that a banked histogram was computed "
              "from the")
        print("    matrices clause [1] rebuilt.  `--full` re-derives them "
              "here with")
        print("    numpy; at order 1676 one leg is hours (about 52x the 716 "
              "leg this")
        print("    same module took here on the measured sub-n^5 scaling, "
              "70x on the")
        print("    n^5 law cert 14 quotes) and the blas path wants about "
              "9.4 GB for its")
        print("    float32 pair matrix, so NO --full leg has been run in "
              "this repository")
        print("    at 1676.  numpy is finder-side only and never in the "
              "trust chain.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 20: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: at order 1676 the decoded (1,1) record H, its Lemma-T")
    print("         i = 2 rebuild H', and H'' -- the decoded record with its")
    print("         twelve off-diagonal core blocks negated -- are pairwise")
    print("         Hadamard-INEQUIVALENT: %d, %d and %d of the %d bins of "
          "the" % (NDIFF_TWIST, NDIFF["decoded"], NDIFF["twisted"], NBINS))
    print("         exact |T4| 4-profile over all C(1676,4) = %d"
          % c_n_4(N))
    print("         row 4-subsets differ.  Two independent implementations")
    print("         agree bin for bin on each matrix; all three hit the")
    print("         second moment %d to the unit."
          % second_moment_want(N))
    print("         ORDER 1676 CARRIES AT LEAST THREE HADAMARD EQUIVALENCE")
    print("         CLASSES, and the psi(rho) = -1 twist provably leaves the")
    print("         class at a THIRD order (668, 716, 1676).  LABEL: PROVEN.")
    if replayed:
        print("         PROFILES: RECOMPUTED in this run from the rebuilt")
        print("         matrices and matched to the bank bin for bin (%s)."
              % ", ".join(replayed))
    else:
        print("         PROFILES: banked exact computation AUDITED in this")
        print("         run, not recomputed.  The three C(1676,4)")
        print("         enumerations ran in the source laboratory")
        print("         (Hadamard-2060, experiments/inequiv/")
        print("         exact_profile_big.py) on 2026-09-02, on a rented")
        print("         16-thread machine, under a pre-registration flushed")
        print("         before the matrices were built.  `--full` would be")
        print("         the replay; it has NOT been run here at 1676.")
    print("         ROW-SIDE ONLY IN THIS CERTIFICATE: the transposed 1676")
    print("         profiles were a pending computation when it was written,")
    print("         so nothing HERE is claimed under the transpose-extended")
    print("         relation.  DISCHARGED by CERT 21.  NOT")
    print("         claimed: anything at 1772 or 2060; any general theorem")
    print("         about psi(rho) = -1 or about orientation; that three is")
    print("         the number of classes at 1676; any novelty or priority.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
