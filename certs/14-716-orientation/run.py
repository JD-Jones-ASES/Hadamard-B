#!/usr/bin/env python3
"""cert 14 -- order 716 carries at least THREE Hadamard equivalence
classes: the orientation switch of the decoded record leaves both known
classes.

  THEOREM.  Let H be the decoded (s, i) = (1, 1) bordered Goethals-Seidel
  record at order 716 (data/payload-records.json; cert 01/11), and let H''
  be H with its TWELVE OFF-DIAGONAL CORE BLOCKS NEGATED and the 4-row/
  4-column border unchanged.  Then H'' is a Hadamard matrix, and it is
  Hadamard-inequivalent to each of H (decoded) and H' (the Lemma-T i = 2
  rebuild, cert 02/11).  Since cert 11 proves H and H' inequivalent,
  order 716 carries at least three Hadamard equivalence classes.

  PROOF (finite, exact).  The multiset {|T4(i,j,k,l)|} over all
  C(716,4) = 10 859 143 295 row 4-subsets is a Hadamard-equivalence
  invariant (note/NOTE-B.md S3.1, invariant I5).  H'' differs from H in
  27 of the 87 populated bins, and from H' in 25.  []

  WHAT H'' IS.  S H'' S, with S = diag(I_4, diag(1,-1,-1,-1) (x) I_n), is
  the same seeds and border assembled in the ALTERNATE Goethals-Seidel
  orientation (the six transposed blocks negated) with the border signed
  by superblock (P[a][J] (-1)^[J != 0], Q[I] (-1)^[I != 0]) -- checked in
  clause [1] as an identity of sign patterns.  So the theorem says: at
  716, as at 668 (cert 13), GS orientation is not a gauge for Hadamard
  equivalence, and the orientation switch is a different class from the
  psi(rho) = -1 twist.  Two orders; no general statement is made.

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The C(716,4) enumeration of H'' was not run inside this repository.  It
  ran in the source laboratory -- Hadamard-2060, experiments/inequiv/
  exact_profile_big.py, numpy, three threads, 2026-09-02, under the
  pre-registration experiments/inequiv/REGISTRATION-716-orientation.md,
  flushed before the matrix was built -- in two arithmetics that agree bin
  for bin, and its output is banked in data/.  The two comparison profiles
  are the banks cert 11 already audits.  The DEFAULT path of this script
  AUDITS all six banks; `--full` RECOMPUTES the H'' profile here from the
  rows clause [1] verified (numpy, finder-side).

WHAT THIS SCRIPT DOES  (default path: standard library only, ~seconds)

  (1) Rebuilds H from the banked record through tools/bordered_gs.py
      (every master-theorem hypothesis re-checked), verifies it, pins its
      canonical digest; forms H'' by negating the twelve off-diagonal core
      blocks, verifies H'' with verify/verify.py, pins ITS canonical digest,
      and checks the alternate-orientation identity; rebuilds H' from the
      cert-02 twisted record, whose seeds are RE-DERIVED here as the
      psi-twist of the decoded seeds (psi(g) = (-1)^g on Z_178), verifies
      it, pins its digest.
  (2) Loads and AUDITS the six banked exact 4-profiles (three matrices x
      two implementations): every bin = 4 (mod 8), total = C(716,4),
      second moment = n^3(n-1)(n-2)/24 to the unit and equal to the field
      the bank declares; each bank's declared matrix digest is compared
      against the in-process digest of the matrix rebuilt in THIS run;
      blas == bits on each matrix.
  (3) The separations: H'' against each of H and H', with the divergent
      bins printed; and cert 11's H vs H' re-affirmed from the same banks.
  (4) Controls: five small Hadamard matrices profiled two ways (one of
      them the route --full takes); the orientation switch applied to a
      GS H(28) control; a total-preserving corruption that only the
      second-moment identity can catch, required to be caught; and the
      dim V / dim W trap, on Sylvester H(16) and on the real 716 objects.
  (5) --full: recompute the H'' profile here (certs/06-668-separation/
      full_recompute.py, imported, not copied) and compare to both banks.

Usage:
  python certs/14-716-orientation/run.py
  python certs/14-716-orientation/run.py --full --impl blas
  python certs/14-716-orientation/run.py --full
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

# The shape of the separations, pinned so a drifting bank cannot quietly
# turn this certificate into a different statement.
NBINS = 87
NDIFF = {"decoded": 27, "twisted": 25}
NDIFF_PRIOR = 27                     # cert 11's H vs H'

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the matrices (the digest verify/verify.py reports).
# All three are rebuilt in clause [1] of THIS run.
SHA_DECODED = \
    "3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6"
SHA_ORIENT = \
    "a6b4f56ec98004e736f0ad74af52826aece4b4ab92750e4706e44486c1885fcd"
SHA_TWISTED = \
    "6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7"

# SHA-256 of the banked data FILES themselves.  The four decoded/twisted
# banks are cert 11's, reused verbatim; the two orient banks are this
# certificate's own.
FILE_PINS = {
    "data/sep716-orient-exact-blas.json":
        "b1c6b0adf393288303f780efebef7ca40bf5d21611660d70928115ff16951cc4",
    "data/sep716-orient-exact-bits.json":
        "bd1c3c23fad9b29a169de1815a4ca44d111b89bb173ef6ac29d55905e077e69c",
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
    ("orient", "blas"): "data/sep716-orient-exact-blas.json",
    ("orient", "bits"): "data/sep716-orient-exact-bits.json",
    ("decoded", "blas"): "data/sep716-exact-blas-decoded.json",
    ("decoded", "bits"): "data/sep716-exact-bits-decoded.json",
    ("twisted", "blas"): "data/sep716-exact-blas-twisted.json",
    ("twisted", "bits"): "data/sep716-exact-bits-twisted.json",
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
else:
    def popcount(x):
        return bin(x).count("1")


# ======================================================================
# exact-integer profile arithmetic
# ======================================================================

def c_n_4(n):
    return n * (n - 1) * (n - 2) * (n - 3) // 24


def second_moment_want(n):
    """sum over 4-subsets of T4^2 = n^3 (n-1)(n-2) / 24  (NOTE-B.md S3.1)."""
    return n ** 3 * (n - 1) * (n - 2) // 24


def audit(prof, n, where):
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
    return dig


def build_and_verify(tag, rec, want_sha):
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
    """V = span{r_i + r_j}, W = V + <all-ones>."""
    f = f2_rows(rows)
    n = len(rows)
    gens = [f[i] ^ f[0] for i in range(1, n)]
    ones = (1 << len(rows[0])) - 1
    return f2_rank(gens), f2_rank(gens + [ones])


# ======================================================================
# controls (as in certs 11 and 13)
# ======================================================================

def sylvester(k):
    n = 1 << k
    return ["".join("+" if bin(x & y).count("1") % 2 == 0 else "-"
                    for y in range(n)) for x in range(n)]


def sylvester_profile_forced(n):
    hit = n * (n - 1) * (n - 2) // 24
    return {0: c_n_4(n) - hit, n: hit}


def paley1(q):
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


GS_CONTROLS = [
    (7, ["+------", "+++----", "++--+--", "+-+-+--"]),            # -> H(28)
    (9, ["++-------", "++-+-----", "++-+-+---", "++--+-+--"]),    # -> H(36)
]


def gs_control(v, seqs):
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
    hist[n] = hist.get(n, 0) - m
    inter = n * (n - 1) * (n - 2) // 2
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
                    help="RECOMPUTE the H'' exact profile here with numpy "
                         "and compare to both banks bin for bin; without it "
                         "the banked profiles are audited, not recomputed")
    ap.add_argument("--impl", choices=("blas", "bits", "both"),
                    default="both",
                    help="which arithmetic path --full recomputes")
    args = ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 14 -- order 716: the orientation switch is a THIRD class")
    print("=" * 72)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
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
    print("\n[1] rebuild H, re-check the theorem's hypotheses, verify, pin; "
          "form H'', verify, pin; rebuild H'")
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
          "core blocks",
          sum(1 for r in range(N) for c in range(N) if rows_o[r][c] != rows_d[r][c])
          == 12 * n * n, "%d cells" % (12 * n * n))
    check("orient   S H'' S == the alternate-orientation array with the "
          "signed border", alternate_orientation_identity(rows_d, rows_o, n, s))
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
    print("\n[1b] control C4 -- the dim-V trap on the real 716 objects")
    vd, wd = dim_V_W(rows_d)
    vt, wt = dim_V_W(rows_t)
    vo, wo = dim_V_W(rows_o)
    check("C4  dim W (INVARIANT) agrees on all three", wd == wt == wo == 715,
          "decoded %d, twisted %d, orient %d" % (wd, wt, wo))
    check("C4  dim V (NOT invariant) differs across the pair -- and is "
          "worthless", vd != vt,
          "decoded %d, twisted %d, orient %d  <- do NOT read this as a "
          "separation" % (vd, vt, vo))

    ROWS = {"orient": rows_o} if args.full else {}
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
        check("%-8s %-4s  %d bins, all |T4| = 4 (mod 8)" % (tag, impl, len(p)),
              all(k % 8 == 4 for k in p))
        check("%-8s %-4s  total == C(716,4) == %d" % (tag, impl, c_n_4(N)),
              tot == c_n_4(N))
        check("%-8s %-4s  second moment == n^3(n-1)(n-2)/24 == %d"
              % (tag, impl, second_moment_want(N)), m2 == second_moment_want(N))
        if "second_moment" in blob:      # every 716 bank declares it
            check("%-8s %-4s  banked second_moment agrees" % (tag, impl),
                  int(blob["second_moment"]) == m2)
        if blob.get("n") != N or int(blob.get("C_n_4", c_n_4(N))) != c_n_4(N):
            check("%-8s %-4s  banked n / C(n,4) header" % (tag, impl), False)
        declared = (blob.get("matrix_canonical_sha256")
                    or blob.get("matrix_sha256"))
        check("%-8s %-4s  bank names the matrix rebuilt in THIS run"
              % (tag, impl),
              is_sha256(declared) and declared == built[tag],
              (declared[:24] + "...") if is_sha256(declared)
              else "declared = %r" % (declared,))
    for tag in ("orient", "decoded", "twisted"):
        check("%-8s  blas == bits, bin for bin (two independent "
              "implementations)" % tag,
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    # ---------------------------------------------------------- clause 3
    print("\n[3] H'' against H and H' -- the separations")
    O = prof[("orient", "blas")]
    check("H'' populates %d bins" % NBINS, len(O) == NBINS)
    for tag, label in (("decoded", "H  (decoded)"), ("twisted", "H' (twist)")):
        X = prof[(tag, "blas")]
        ks = sorted(set(O) | set(X))
        diff = [(k, X.get(k, 0), O.get(k, 0)) for k in ks
                if X.get(k, 0) != O.get(k, 0)]
        check("H'' vs %-12s support: %s" % (label,
              "identical" if set(O) == set(X) else "DIFFERENT"), True,
              "%d union bins" % len(ks))
        check("H'' vs %-12s bin counts differ in exactly %d bins"
              % (label, NDIFF[tag]), len(diff) == NDIFF[tag],
              "%d differing" % len(diff))
        check("H'' vs %-12s the differences sum to zero" % label,
              sum(q - p for _k, p, q in diff) == 0)
        big = max(diff, key=lambda t: abs(t[2] - t[1]))
        print("      largest |delta| = %d at |T4| = %d (%.2e of that bin); "
              "first eight divergent bins:" % (abs(big[2] - big[1]), big[0],
                                               abs(big[2] - big[1]) / big[1]))
        for k, p, q in diff[:8]:
            print("        |T4|=%4d  %s %16d   H'' %16d   delta %+d"
                  % (k, label[:2], p, q, q - p))
    print("\n      the prior separation, re-affirmed from the same banks")
    D, T = prof[("decoded", "blas")], prof[("twisted", "blas")]
    kd = sorted(set(D) | set(T))
    nprior = sum(1 for k in kd if D.get(k, 0) != T.get(k, 0))
    check("H vs H'  differ in exactly %d of the %d bins (cert 11)"
          % (NDIFF_PRIOR, NBINS), nprior == NDIFF_PRIOR and len(kd) == NBINS,
          "%d differing" % nprior)

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
        check("C1  %-16s straight enumeration == pair-vector route"
              % name, p1 == p2, "%d bins (%.1fs)" % (len(p1), time.time() - t0))
        check("C1  %-16s |T4| = n (mod 8) on every populated bin" % name,
              all(k % 8 == nn % 8 for k in p1), "n mod 8 = %d" % (nn % 8))
        if name.startswith("Sylvester"):
            check("C1  %-16s matches the FORCED Sylvester profile" % name,
                  p1 == sylvester_profile_forced(nn))
    # the orientation switch on a control: still Hadamard (it is a GS array
    # in the other orientation), and the profile is an invariant of the
    # matrix, not of its labelling
    rows28, _ = gs_control(7, GS_CONTROLS[0][1])
    sw28 = orientation_switch(rows28, 7, 0)
    check("C0  GS H(28) with its twelve off-diagonal blocks negated is Hadamard",
          is_hadamard(sw28))

    print("\n  C2 -- negative control: the second-moment assert must FIRE")
    victim = dict(prof[("orient", "blas")])
    kk = sorted(victim)
    victim[kk[0]] -= 1
    victim[kk[-1]] += 1
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
        print("\n[5] --full: RECOMPUTING the H'' exact profile here, from the "
              "rows clause [1] verified, with numpy (%s)." % " and ".join(impls))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS"):
            os.environ[var] = "3"           # set BEFORE numpy is imported
        sys.path.insert(0, os.path.join(ROOT, "certs", "06-668-separation"))
        import full_recompute as FR                          # noqa: E402
        h128 = sylvester(7)
        want128 = sylvester_profile_forced(128)
        for impl in impls:
            got = FR.profile(h128, 128, impl, progress=False)
            audit(got, 128, "full/H128/%s" % impl)
            check("[full] Sylvester H(128) %-4s == the forced profile" % impl,
                  got == want128)
        rows = ROWS["orient"]
        check("[full] the rows about to be enumerated are the ones verify.py "
              "accepted", rows_sha256(rows) == SHA_ORIENT)
        for impl in impls:
            t0 = time.time()
            got = FR.profile(rows, N, impl)
            audit(got, N, "full/orient/%s" % impl)
            secs = time.time() - t0
            for bimpl in ("blas", "bits"):
                check("[full] orient recomputed %-4s == banked %-4s, bin for bin"
                      % (impl, bimpl), got == prof[("orient", bimpl)],
                      "%d bins, %.0fs" % (len(got), secs))
            replayed.append("orient/%s" % impl)
        ROWS.clear()
    else:
        print("\n[5] --full not requested: the banked profiles were AUDITED, "
              "not recomputed.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 14: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: H'' -- the decoded order-716 record with its twelve")
    print("         off-diagonal core blocks negated -- is Hadamard and is")
    print("         INEQUIVALENT to each of H (decoded) and H' (Lemma-T")
    print("         twist): %d and %d of %d bins differ."
          % (NDIFF["decoded"], NDIFF["twisted"], NBINS))
    print("         With cert 11, ORDER 716 CARRIES AT LEAST THREE")
    print("         HADAMARD EQUIVALENCE CLASSES.  LABEL: PROVEN.")
    if replayed:
        print("         PROFILE: RECOMPUTED here and matched to the bank (%s)."
              % ", ".join(replayed))
    else:
        print("         PROFILE: banked exact computation AUDITED, not")
        print("         recomputed (Hadamard-2060, exact_profile_big.py,")
        print("         2026-09-02, REGISTRATION-716-orientation.md).")
    print("         Row-side only.  NOT claimed: 1676/1772 orientation")
    print("         switches; anything transpose-extended at 716; that")
    print("         these are the only three classes.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
