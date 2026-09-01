#!/usr/bin/env python3
"""cert 11 -- the two Hadamard matrices of order 716 banked in this
repository are Hadamard-INEQUIVALENT.

  THEOREM.  Let H1 be the decoded (s, i) = (1, 1) bordered Goethals-Seidel
  record at order 716 (data/payload-records.json) and let H2 be its Lemma-T
  i = 2 rebuild (data/twisted-i2-records.json, order 716).  Then H1 and H2
  are not Hadamard-equivalent: there is no H2 = D_r P_r H1 P_c D_c with P
  permutation matrices and D diagonal +-1.

  PROOF (finite, exact).  The multiset
  {|T4(i,j,k,l)|} over all C(716,4) = 10 859 143 295 row 4-subsets, with
  T4 = sum_c H[i][c]H[j][c]H[k][c]H[l][c], is a Hadamard-equivalence
  invariant (note/NOTE-B.md S3.1, invariant I5: each row negation
  contributes one sign to T4 so |T4| is fixed, each column negation
  contributes d_c^4 = 1, and permutations relabel).  The two profiles
  populate the same 87 bins but 27 of the 87 bin counts DIFFER.  An
  invariant that differs is a separation.  []

  CONSEQUENCE.  Order 716 carries at least two Hadamard equivalence
  classes, and the Lemma-T construction at psi(rho) = -1 leaves the
  equivalence class at a SECOND order (668, cert 06, was the first).

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The C(716,4) enumeration was not run inside this repository.  It was run
  in the source laboratory -- Hadamard-2060, experiments/inequiv/
  exact_profile_big.py, numpy, three threads, 2026-09-01, under the
  pre-registration experiments/inequiv/REGISTRATION-716-exact.md that was
  flushed before the matrices were built -- and its output is banked in
  data/.  So the DEFAULT path of this script AUDITS a banked exact
  computation.  It does NOT recompute the histograms, and no check on the
  default path shows that a banked histogram was computed from the matrix
  this script rebuilt.  Each of the four banks DECLARES the canonical
  digest of its matrix -- recorded by the producer against the file it
  enumerated, not added afterwards -- and the default path binds that
  declaration to the matrix rebuilt in THIS run; but a self-declared
  digest is metadata, not a computation.  `--full` is the recomputation --
  it re-derives both profiles here, from the rebuilt rows, and compares
  them to the bank bin for bin.  Say "banked exact computation AUDITED" of
  a default run; say "replayed" only of a --full run.

WHAT THIS SCRIPT DOES  (default path: standard library only, ~seconds)

  (1) Rebuilds both matrices from the banked records through
      tools/bordered_gs.py -- re-checking every hypothesis of the master
      theorem on the way, not merely assembling -- hands each to
      verify/verify.py (the trust chain), compares the canonical SHA-256
      against the digest pinned below, RETAINS that in-process digest for
      the bank-identity check in step (2), and DELETES the generated
      matrices.  The twisted record is cert 02's shared bank rather than a
      file of this certificate's own, so its seeds are RE-DERIVED here as
      the psi-twist of the decoded seeds and compared character for
      character: the record is bound to data/payload-records.json by
      computation, not by a file pin.
  (2) Loads and AUDITS the four banked exact 4-profiles (two matrices x
      two independent implementations) -- it does not recompute them --
      pinning each file's SHA-256 in code, reporting for each whether the
      bank names the matrix rebuilt in step (1), and asserting on each, in
      exact integer arithmetic:
        - every populated bin is  = 4 (mod 8)  (716 = 4 mod 8);
        - the counts total C(716,4);
        - the second moment equals n^3(n-1)(n-2)/24  (NOTE-B.md S3.1);
      then asserts blas == bits bin for bin on each matrix.
  (3) Compares decoded against twisted: identical 87-bin support, 27
      differing bins, differences summing to zero.  Prints the divergent
      bins and the verdict.
  (4) Runs four controls, all in the standard library:
        C1  the full |T4| profile of five small Hadamard matrices --
            Sylvester H(8) and H(16), Paley H(20), and two plain
            Goethals-Seidel arrays H(28) and H(36) built by this
            repository's own assembler at the degenerate s = 0 layer,
            whose profiles have three and four populated bins so the
            agreement is not vacuous -- each computed TWICE, by straight
            O(C(n,4)) enumeration over the +-1 rows and by the
            pair-vector / Gram-triangle route --full uses, and required
            to agree bin for bin;
        C2  a NEGATIVE control: a banked profile is corrupted in a
            total-preserving way and the second-moment assert must fire;
        C3  the dim-V trap, on Sylvester H(16) under a deterministically
            seeded random signed row negation: dim V moves, dim W does not;
        C4  the same trap on the real objects: dim V DIFFERS across the
            716 pair (714 vs 715) and is worthless, while the invariant
            dim W is 715 on both.

  --full  RECOMPUTES both exact profiles here, from the matrices rebuilt
          in step (1), with numpy (certs/06-668-separation/
          full_recompute.py, imported by path rather than copied, so this
          certificate and certs 06/08 cannot drift apart) and compares
          them to the banked JSONs bin for bin.  This is the only path on
          which the word "replayed" is earned, and the only path that
          binds a bank to a matrix.  Note that it is a THIRD arithmetic
          route as well as a replay: the banks were produced by the
          canonical-split engine, while full_recompute.py enumerates the
          upper triangle of U U^T and divides out the threefold
          overcount.  numpy is imported only under this flag; it is
          finder-side only and is never in the trust chain.  BLAS threads
          are capped at 3.  Cost: about a quarter of an hour for
          `--impl blas` alone, about an hour and a half for both
          arithmetic paths.

Usage:
  python certs/11-716-separation/run.py
  python certs/11-716-separation/run.py --full
  python certs/11-716-separation/run.py --full --impl blas
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

N = 716
OUT = os.path.join(HERE, "out")

# The shape of the separation, pinned so a drifting bank cannot quietly
# turn this certificate into a different (weaker) statement.
NBINS = 87
NDIFF = 27
TAIL_FROM = 228              # every bin at or above this one AGREES

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the two matrices, pinned by the repository's own
# earlier work (cert 01 for the decoded record; cert 02 for the Lemma-T
# i = 2 rebuild) and re-derived here from the banked parameters alone.
SHA_DECODED = \
    "3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6"
SHA_TWISTED = \
    "6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7"

# SHA-256 of the banked data FILES themselves, so a silently edited bank is
# a hard error rather than a different theorem.  Only this certificate's own
# four profile banks are pinned here: data/payload-records.json and
# data/twisted-i2-records.json are shared with certs 01 and 02, and the
# binding pin on each is the canonical digest of the matrix it produces,
# which is checked in clause [1].
FILE_PINS = {
    "data/sep716-exact-blas-decoded.json":
        "80ee1e151ec1f759d7213d500623603716b9afa6fc382a385ce6970efac35a6b",
    "data/sep716-exact-bits-decoded.json":
        "a0d5b3a65b83c39c905ec2a1d3b25ca1c58e0106b76aaa6eb3b2feee3748aeed",
    "data/sep716-exact-blas-twisted.json":
        "e2076eb890557775edaccda3c9dcbab7e585f5181e0da1bb5914913bf0749b46",
    "data/sep716-exact-bits-twisted.json":
        "c385773a7a2bf4506b94752406be787bac3d85b8a9ee18bb23668080c6afe7bc",
}

PROFILES = {
    ("decoded", "blas"): "data/sep716-exact-blas-decoded.json",
    ("decoded", "bits"): "data/sep716-exact-bits-decoded.json",
    ("twisted", "blas"): "data/sep716-exact-blas-twisted.json",
    ("twisted", "bits"): "data/sep716-exact-bits-twisted.json",
}

HEXDIGITS = set("0123456789abcdef")
CH = {1: "+", -1: "-"}


def is_sha256(v):
    return (isinstance(v, str) and len(v) == 64
            and all(c in HEXDIGITS for c in v))


def backfilled_from_pin(blob):
    """True when a bank's own banked_note says its matrix digest was
    backfilled from the digest this repository already pinned, rather than
    recorded against the producer's own file.  None of the four banks read
    here is of that kind -- all were producer-banked on 2026-09-01 -- but
    the label is computed rather than assumed, so a bank that were ever
    re-issued as a backfill would say so on this certificate's own output."""
    return "backfilled" in str(blob.get("banked_note", ""))

FAIL = []


def check(label, cond, extra=""):
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
else:  # 3.9 fallback -- same shape as verify/verify.py's, kept local so
    # this certificate imports nothing from the trust chain
    def popcount(x):
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
    second one is sharp: it pins a 13-digit number to the unit.
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
# clause 1 -- rebuild, verify, pin
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
    path = os.path.join(OUT, "H716_%s.txt" % tag)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify", "verify.py"), path],
        capture_output=True, text=True)
    verdict = ((proc.stdout or proc.stderr).strip().splitlines()
               or ["(no output)"])[-1]
    check("%-8s verify/verify.py" % tag, proc.returncode == 0, verdict[:78])
    dig = rows_sha256(rows)
    check("%-8s canonical sha256 == pin" % tag, dig == want_sha, dig[:24] + "...")
    check("%-8s verify.py reports the same digest" % tag, dig in verdict)
    print("       rebuilt and verified in %.1fs" % (time.time() - t0))
    return rep, rows, dig


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
# this repository's master theorem.  Found by exhaustive search over
# {+,-}^v offline; the PAF condition is RE-VERIFIED below, never assumed,
# and the assembled matrix is put through the same C0 Hadamard check as
# every other control.  These two give profiles with three and four
# populated bins, which Sylvester and Paley matrices do not.
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
                    help="RECOMPUTE both 716 exact profiles here with numpy "
                         "and compare to the bank bin for bin; without it "
                         "the banked profiles are audited, not recomputed")
    ap.add_argument("--impl", choices=("blas", "bits", "both"),
                    default="both",
                    help="which arithmetic path --full recomputes "
                         "(default both, about an hour and a half; blas "
                         "alone is about a quarter of an hour)")
    args = ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 11 -- order 716: two Hadamard matrices, one invariant, "
          "two classes")
    print("=" * 72)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    # The scratch directory holds two ~515 KB generated matrices.  Wrapping
    # the body means an exception -- in the optional numpy path or anywhere
    # else -- cannot leave them on disk.
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
        check("%-38s" % name, got == want, got[:24] + "...")

    # ---------------------------------------------------------- clause 1
    print("\n[1] rebuild both matrices, re-check the theorem's hypotheses, "
          "verify, pin")
    dec_rec = load_record(os.path.join(ROOT, "data", "payload-records.json"),
                          N)
    tw_rec = load_record(os.path.join(ROOT, "data",
                                      "twisted-i2-records.json"), N)
    check("the two records are genuinely different instances",
          dec_rec["seeds"] != tw_rec["seeds"]
          and int(dec_rec.get("coset_divisors", [1])[0]) == 1
          and int(tw_rec["coset_divisors"][0]) == 2,
          "decoded i=1, rebuild i=2, different seeds")
    # The twisted record is cert 02's shared bank, so it is NOT file-pinned
    # here.  Instead its seeds are re-derived from the decoded record, which
    # binds the two records to each other by computation.
    check("the twisted seeds are the psi-twist of the decoded seeds, "
          "re-derived here",
          list(tw_rec["seeds"]) == psi_twist_seeds(dec_rec)
          and list(tw_rec["group"]) == list(dec_rec["group"])
          and list(tw_rec["r_shift"]) == list(dec_rec["r_shift"]),
          "psi(g) = (-1)^g on Z_%d; rho = %d is odd, so psi(rho) = -1"
          % (int(dec_rec["group"][0]), int(dec_rec["r_shift"][0])))
    rep_d, rows_d, dig_d = build_and_verify("decoded", dec_rec, SHA_DECODED)
    rep_t, rows_t, dig_t = build_and_verify("twisted", tw_rec, SHA_TWISTED)
    built = {"decoded": dig_d, "twisted": dig_t}
    if rows_d is None or rows_t is None:
        print("\nFATAL: a matrix did not rebuild; nothing further is checked.")
        return 1

    # ---------------------------------------------------------- clause 1b
    # The dim-V trap on the real objects, while the rows are still in hand.
    print("\n[1b] control C4 -- the dim-V trap on the 716 pair itself")
    vd, wd = dim_V_W(rows_d)
    vt, wt = dim_V_W(rows_t)
    check("dim W (INVARIANT) agrees on both", wd == wt == 715,
          "decoded %d, twisted %d" % (wd, wt))
    check("dim V (NOT invariant) differs -- and is worthless", vd != vt,
          "decoded %d, twisted %d  <- do NOT read this as a separation"
          % (vd, vt))

    # The rebuilt rows are wanted again only by --full, which recomputes the
    # profiles FROM THEM.  On the default path they go now.
    ROWS = {"decoded": rows_d, "twisted": rows_t} if args.full else {}
    del rows_d, rows_t

    # ---------------------------------------------------------- clause 2
    print("\n[2] the four banked exact 4-profiles, AUDITED in exact integers")
    print("    (audited, not recomputed: see --full)")
    prof = {}
    raw = {}
    unbound = []
    for (tag, impl), name in sorted(PROFILES.items()):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        raw[(tag, impl)] = blob
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(tag, impl)] = p
        ok_mod = all(k % 8 == 4 for k in p)
        tot, m2 = audit(p, N, "%s/%s" % (tag, impl))
        check("%-8s %-4s  %d bins, all |T4| = 4 (mod 8)" % (tag, impl, len(p)),
              ok_mod)
        check("%-8s %-4s  total == C(716,4) == %d" % (tag, impl, c_n_4(N)),
              tot == c_n_4(N))
        check("%-8s %-4s  second moment == n^3(n-1)(n-2)/24 == %d"
              % (tag, impl, second_moment_want(N)), m2 == second_moment_want(N))
        if "second_moment" in blob:      # every 716 bank declares it
            check("%-8s %-4s  banked second_moment agrees" % (tag, impl),
                  int(blob["second_moment"]) == m2)
        if blob.get("n") != N or int(blob.get("C_n_4", c_n_4(N))) != c_n_4(N):
            check("%-8s %-4s  banked n / C(n,4) header" % (tag, impl), False)
        # Matrix identity.  Where a bank declares a matrix digest it is
        # compared against the digest of the matrix REBUILT IN THIS RUN --
        # not against a static string -- so a bank cannot drift onto a
        # different object.  Where it declares none, say so out loud.
        # Keyed on PRESENCE, not truthiness: a declared digest that is
        # empty, null, or not 64 hex digits is a FAILURE.  Only a bank that
        # declares no digest at all falls to the [NOTE] below.
        declared = (blob.get("matrix_canonical_sha256")
                    or blob.get("matrix_sha256"))
        if "matrix_canonical_sha256" in blob or "matrix_sha256" in blob:
            check("%-8s %-4s  bank names the matrix rebuilt in THIS run%s"
                  % (tag, impl,
                     " [digest backfilled from this repository's own pin]"
                     if backfilled_from_pin(blob) else ""),
                  is_sha256(declared) and declared == built[tag],
                  (declared[:24] + "...") if is_sha256(declared)
                  else "declared = %r" % (declared,))
        else:
            unbound.append("%s/%s" % (tag, impl))

    if unbound:
        print("      [NOTE] no matrix digest is declared by: %s."
              % ", ".join(unbound))
        print("             Those banks are bound to the rebuilt matrices by")
        print("             --full alone; the default path cannot show that")
        print("             their histograms came from these matrices.")

    for tag in ("decoded", "twisted"):
        check("%-8s  blas == bits, bin for bin (two independent "
              "implementations)" % tag,
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    # ---------------------------------------------------------- clause 3
    print("\n[3] decoded vs twisted -- the separation")
    D, T = prof[("decoded", "blas")], prof[("twisted", "blas")]
    ks = sorted(set(D) | set(T))
    check("identical support: the same %d populated bins" % len(D),
          set(D) == set(T) and len(D) == NBINS)
    diff = [(k, D[k], T[k]) for k in ks if D.get(k, 0) != T.get(k, 0)]
    check("bin counts differ in exactly %d of the %d bins" % (NDIFF, NBINS),
          len(diff) == NDIFF, "%d differing" % len(diff))
    check("the differences sum to zero (both totals are C(716,4))",
          sum(q - p for _k, p, q in diff) == 0)
    check("the extreme tail does NOT separate them",
          all(D[k] == T[k] for k in ks if k >= TAIL_FROM),
          "every one of the %d bins from |T4| = %d to %d agrees exactly"
          % (sum(1 for k in ks if k >= TAIL_FROM), TAIL_FROM, ks[-1]))
    m1 = sum(k * (q - p) for k, p, q in diff)
    check("the FIRST moment, which nothing forces, does differ", m1 != 0,
          "sum |T4|*delta = %d" % m1)

    print("\n      the %d divergent bins" % len(diff))
    print("      %6s %16s %16s %12s %10s"
          % ("|T4|", "decoded", "twisted", "delta", "rel"))
    for k, p, q in diff:
        print("      %6d %16d %16d %+12d %10.2e"
              % (k, p, q, q - p, abs(q - p) / p))
    big = max(diff, key=lambda t: abs(t[2] - t[1]))
    print("      largest |delta| = %d at |T4| = %d, i.e. %.2e of that bin"
          % (abs(big[2] - big[1]), big[0], abs(big[2] - big[1]) / big[1]))

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
        n = len(rows)
        check("C0  %-16s is in fact Hadamard" % name, is_hadamard(rows))
        t0 = time.time()
        p1 = profile_straight(rows)
        p2 = profile_pairvec(rows)
        p1 = {k: v for k, v in p1.items() if v}
        audit(p1, n, "%s straight" % name)
        audit(p2, n, "%s pairvec" % name)
        check("C1  %-16s straight enumeration == pair-vector route"
              % name, p1 == p2,
              "%d bins: %s  (%.1fs)"
              % (len(p1), " ".join("%d:%d" % (k, p1[k])
                                   for k in sorted(p1)), time.time() - t0))
        check("C1  %-16s both routes hit C(n,4) and the second moment"
              % name, True, "C(%d,4)=%d  m2=%d" % (n, c_n_4(n),
                                                   second_moment_want(n)))
        check("C1  %-16s |T4| = n (mod 8) on every populated bin" % name,
              all(k % 8 == n % 8 for k in p1),
              "n mod 8 = %d" % (n % 8))
        if name.startswith("Sylvester"):
            check("C1  %-16s matches the FORCED Sylvester profile" % name,
                  p1 == sylvester_profile_forced(n),
                  "predicted %s" % sylvester_profile_forced(n))

    print("\n  C2 -- negative control: the second-moment assert must FIRE")
    victim = dict(prof[("decoded", "blas")])
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
    check("C2  the corrupted profile still totals C(716,4) -- so ONLY the "
          "second moment could catch it",
          sum(victim.values()) == c_n_4(N))

    print("\n  C3 -- the dim-V trap, demonstrated on Sylvester H(16)")
    import random                                            # noqa: E402
    rng = random.Random(20260901)
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
          "%d -> %d  (seed 20260901, %d of 16 rows negated)"
          % (v0, v1, sum(flip)))
    check("C3  dim W does NOT move", w0 == w1, "%d -> %d" % (w0, w1))
    check("C3  the |T4| profile does not move either (it is an invariant)",
          profile_straight(h16) == profile_straight(h16b))

    # ---------------------------------------------------------- clause 5
    replayed = []
    if args.full:
        impls = (("blas", "bits") if args.impl == "both" else (args.impl,))
        print("\n[5] --full: RECOMPUTING both exact profiles here, from the "
              "matrices")
        print("    rebuilt in clause [1], with numpy (%s path%s).  This is "
              "the" % (" and ".join(impls), "" if len(impls) == 1 else "s"))
        print("    step that earns the word 'replayed'.  It is also a THIRD")
        print("    arithmetic route: the banks came from the canonical-split")
        print("    engine, this one enumerates the U U^T triangle.")
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS"):
            os.environ[var] = "3"           # set BEFORE numpy is imported
        # The recompute machinery is cert 06's full_recompute.py, imported by
        # path rather than copied, so the certificates cannot drift apart.
        sys.path.insert(0, os.path.join(ROOT, "certs", "06-668-separation"))
        import full_recompute as FR                          # noqa: E402
        # Smoke the ported numpy paths BEFORE spending the hour on them, on a
        # matrix whose profile is forced (so this is a positive control with
        # a predicted answer) and big enough that the packed path needs more
        # than one uint64 word per row -- which the n <= 36 controls above
        # cannot exercise.
        h128 = sylvester(7)
        want128 = sylvester_profile_forced(128)
        for impl in impls:
            got = FR.profile(h128, 128, impl, progress=False)
            audit(got, 128, "full/H128/%s" % impl)
            check("[full] Sylvester H(128) %-4s == the forced profile "
                  "(%d uint64 words/row)" % (impl, (128 + 63) // 64),
                  got == want128)
        for tag in ("decoded", "twisted"):
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
        print("    matrices clause [1] rebuilt.  `--full` re-derives both "
              "here with")
        print("    numpy (~1.5 h both paths, ~15 min with `--impl blas`); "
              "numpy is")
        print("    finder-side only and never in the trust chain.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    if FAIL:
        print("cert 11: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: the decoded (1,1) record and the Lemma-T i=2 rebuild at")
    print("         order 716 are Hadamard-INEQUIVALENT, so order 716")
    print("         carries at least two Hadamard equivalence classes.")
    print("         Separating invariant: the exact |T4| 4-profile over all")
    print("         C(716,4) = %d row 4-subsets." % c_n_4(N))
    print("         %d of %d bins differ; two independent implementations"
          % (NDIFF, NBINS))
    print("         agree bin for bin; both hit the second moment %d"
          % second_moment_want(N))
    print("         to the unit.  LABEL: PROVEN.")
    if replayed:
        print("         PROFILES: RECOMPUTED in this run from the rebuilt")
        print("         matrices and matched to the bank bin for bin (%s)."
              % ", ".join(replayed))
    else:
        print("         PROFILES: banked exact computation AUDITED in this")
        print("         run, not recomputed.  The C(716,4) enumeration ran")
        print("         in the source laboratory (Hadamard-2060,")
        print("         experiments/inequiv/exact_profile_big.py) on")
        print("         2026-09-01, under a pre-registration flushed before")
        print("         the build.  `--full` recomputes it here and is the")
        print("         replay.")
    print("         The twist at psi(rho) = -1 now provably leaves the")
    print("         equivalence class at a SECOND order (668 was the first,")
    print("         cert 06).  NOT claimed here: 1676, 1772, and nothing")
    print("         about the transpose-extended relation at 716.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
