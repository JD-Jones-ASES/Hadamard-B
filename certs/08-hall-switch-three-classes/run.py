#!/usr/bin/env python3
"""cert 08 -- order 668 carries at least THREE Hadamard equivalence
classes, and it still does when the transpose is added to the group.

  THEOREM.  Let H be the decoded (s, i) = (1, 1) bordered Goethals-Seidel
  record at order 668 (data/payload-records.json), let H' be its Lemma-T
  i = 2 rebuild (data/sep668-twisted-record.json), and let H* be the
  paired-Hall-switch matrix of the anonymous preprint at
  hadamard-668.vercel.app (data/sep668-hall-switch.json, rebuilt here from
  H by negating the 1 328 entries the banked mask names).  Then H, H' and
  H* are pairwise Hadamard-INEQUIVALENT.  Adding the transpose to the
  group does not merge any of the three pairs either.

  PROOF (finite, exact, replayed by this script).  The multiset
  {|T4(i,j,k,l)|} over all C(668,4) = 8 222 179 035 row 4-subsets, with
  T4 = sum_c H[i][c]H[j][c]H[k][c]H[l][c], is a Hadamard-equivalence
  invariant (note/NOTE-B.md S3.1, invariant I5).  The three profiles
  populate the same 80 bins and differ in 26 (H vs H'), 27 (H vs H*) and
  27 (H' vs H*) of the bin counts.  For the transpose-extended relation,
  A is equivalent to B only if A ~ B or A ~ B^T, so each pair needs a
  second refutation; the profiles of (H')^T and (H*)^T supply them --
  50 bins for H vs (H')^T, 49 for H vs (H*)^T, 50 for H' vs (H*)^T.  An
  invariant that differs is a separation.  []

  PRIORITY.  The preprint was FIRST to publish that order 668 carries at
  least two Hadamard equivalence classes, and its first matrix is
  byte-identical to the record banked here.  See NOTES.md.

WHAT THIS SCRIPT DOES  (standard library only, seconds)

  (0) Pins the SHA-256 of every banked file it reads.
  (1) Rebuilds all five matrices -- H, H', H*, (H*)^T, (H')^T -- from the
      banked records and the banked switch mask, re-checking every
      hypothesis of the master theorem on the way for the two records,
      hands each to verify/verify.py (the trust chain), compares the
      canonical SHA-256 against the digest pinned below, and DELETES the
      generated matrices.  It also re-derives, from data/payload-
      records.json alone, the preprint's four subsets X1..X4 (as the
      negative supports of the four decoded seeds) and its blocks K, T, S
      (as the record's corner, row_table and col_table^T), and reproduces
      both SHA-256 digests the preprint publishes for its own matrices.
  (2) Loads eight banked exact 4-profiles, auditing each in exact integer
      arithmetic: every populated bin = 4 (mod 8), the counts total
      C(668,4), and the second moment equals n^3(n-1)(n-2)/24.  Then
      asserts blas == bits bin for bin on all five matrices -- H, H',
      H*, (H*)^T and (H')^T -- so every leg of the theorem, the
      transpose-extended legs included, rests on two independent
      implementations.
  (3) Compares the three matrices pairwise under the standard relation.
  (4) Compares them pairwise under the TRANSPOSE-EXTENDED relation.
  (5) Six controls, all in the standard library -- see NOTES.md.

Usage:
  python certs/08-hall-switch-three-classes/run.py
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

N = 668
OUT = os.path.join(HERE, "out")

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the five matrices (the digest verify/verify.py
# reports), each re-derived here from banked parameters alone.
SHA = {
    "decoded":
        "bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0",
    "twisted":
        "600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3",
    "hall":
        "7f6af1d9e9c80fae52c6f178e298144d209c01f6c103db119e6c12425aa1a722",
    "hall-T":
        "565a9ca5a9db739f74215474364202a4d35fd691542b18e8ace8bcbb3c190c65",
    "twisted-T":
        "32afde351e1f44aa4236cb3c406fbbcd59c5ab2cc3e02c8380e08680db2d19d7",
}

# SHA-256 of the banked data FILES themselves, so a silently edited bank is
# a hard error rather than a different theorem.  The first five are shared
# read-only with cert 06 and carry cert 06's values unchanged.
FILE_PINS = {
    "data/sep668-twisted-record.json":
        "fe8154179ba2ebfe097c82e468368cdc8a070548555bb10140949af0560611fb",
    "data/sep668-exact-blas-decoded.json":
        "22df5ce9fcd6eb307f56981c507bb46b2a18b79861d903349dc13458a6dffcbf",
    "data/sep668-exact-bits-decoded.json":
        "0bafbf8219d33b9c74786700106aeba3086bbf577ee02bcda43768f35978fdd8",
    "data/sep668-exact-blas-twisted.json":
        "c4d8db3ba40cf8c5a244607032dab6b66d878b8fe6b98784351f7b8ae70e5a17",
    "data/sep668-exact-bits-twisted.json":
        "91d154d05ccea87a6fa98a02b4fcbf275dc6b4025650116941647216a69faf5a",
    "data/sep668-hall-switch.json":
        "13efd2402b8394c62c901af4f7cfbec7b2e474832dd3055c6b9e9e220b351c85",
    "data/sep668-hall-exact-blas.json":
        "35e716ecb43bb6190d5dd6f4160e0bc2bed4f61a3aacf07a36ff9d190810c154",
    "data/sep668-hall-exact-bits.json":
        "a6f703b499d98995f6446a1aed671284c47e99cfe869f3ce8dc8b5fd9394accb",
    "data/sep668-hall-T-exact-blas.json":
        "151fb5d6e70cf56d6a1c2aa124a597a837bca0ecf5d64958b43a34c05383e0db",
    "data/sep668-hall-T-exact-bits.json":
        "48fdb26f8b1ee5135ed278ec866e204c1ab47df168c043fabff8699c0f4fd8bb",
    "data/sep668-twisted-T-exact.json":
        "38355274ec61d33fcd96e24255e4a7b02874150cd914fdfb928d28cee751fc4a",
}

# The eight banked profiles, keyed (matrix tag, implementation).
PROFILES = {
    ("decoded", "blas"): "data/sep668-exact-blas-decoded.json",
    ("decoded", "bits"): "data/sep668-exact-bits-decoded.json",
    ("twisted", "blas"): "data/sep668-exact-blas-twisted.json",
    ("twisted", "bits"): "data/sep668-exact-bits-twisted.json",
    ("hall", "blas"): "data/sep668-hall-exact-blas.json",
    ("hall", "bits"): "data/sep668-hall-exact-bits.json",
    ("hall-T", "blas"): "data/sep668-hall-T-exact-blas.json",
    ("hall-T", "bits"): "data/sep668-hall-T-exact-bits.json",
}
# The twisted-transpose bank carries BOTH implementations in one file.
TWISTED_T = "data/sep668-twisted-T-exact.json"

FLIP = str.maketrans("+-", "-+")
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


# ======================================================================
# exact-integer profile arithmetic  (same contract as cert 06)
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
    """The two forced identities, as ASSERTS.  Exact integers throughout."""
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
    """(differing bins, delta sum, first-moment delta, same support)."""
    ks = sorted(set(pa) | set(pb))
    diff = [(k, pa.get(k, 0), pb.get(k, 0)) for k in ks
            if pa.get(k, 0) != pb.get(k, 0)]
    return (diff,
            sum(q - p for _k, p, q in diff),
            sum(k * (q - p) for k, p, q in diff),
            set(pa) == set(pb))


# ======================================================================
# clause 1 -- rebuild, switch, verify, pin
# ======================================================================

def load_record(path, order):
    with open(path, "r", encoding="ascii") as fh:
        blob = json.load(fh)
    if "orders" in blob:
        return [r for r in blob["orders"] if int(r["order"]) == order][0]
    return blob


def verify_and_pin(tag, rows, want_sha):
    """Hand a matrix to the trust chain and pin its canonical digest."""
    path = os.path.join(OUT, "H668_%s.txt" % tag)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify", "verify.py"), path],
        capture_output=True, text=True)
    verdict = ((proc.stdout or proc.stderr).strip().splitlines()
               or ["(no output)"])[-1]
    check("%-9s verify/verify.py" % tag, proc.returncode == 0, verdict[:76])
    dig = rows_sha256(rows)
    check("%-9s canonical sha256 == pin" % tag, dig == want_sha,
          dig[:24] + "...")
    check("%-9s verify.py reports the same digest" % tag, dig in verdict)
    return dig


def build_mask(spec, s, n):
    """Re-derive the switch mask from the banked rectangles AND from the
    closed forms in (s, n); return the set of (row, col) positions."""
    rects = spec["switch"]["mask_rectangles"]
    want = [(0, 4 * s, 4 * s, 4 * s + n), (4 * s, 4 * s + n, 0, 4 * s)]
    got = [(r["row_start"], r["row_stop"], r["col_start"], r["col_stop"])
           for r in rects]
    mask = set()
    for r0, r1, c0, c1 in got:
        for r in range(r0, r1):
            for c in range(c0, c1):
                mask.add((r, c))
    return mask, got, want


def paper_bytes_sha(rows):
    """The preprint's format: one byte per entry, 0x01 for +1, 0x00 for -1,
    row-major.  446 224 bytes at order 668."""
    return hashlib.sha256(
        bytes(1 if ch == "+" else 0 for r in rows for ch in r)).hexdigest()


def transpose_rows(rows):
    return ["".join(col) for col in zip(*rows)]


# ======================================================================
# small-matrix controls
# ======================================================================

def sylvester(k):
    n = 1 << k
    return ["".join("+" if bin(x & y).count("1") % 2 == 0 else "-"
                    for y in range(n)) for x in range(n)]


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


def is_hadamard(rows):
    n = len(rows)
    pk = [int(r.translate(str.maketrans("+-", "01"))[::-1], 2) for r in rows]
    if any(len(r) != n for r in rows):
        return False
    for i in range(n):
        for j in range(i + 1, n):
            if (pk[i] ^ pk[j]).bit_count() != n // 2:
                return False
    return True


def profile_straight(rows):
    """Straight O(C(n,4)) enumeration over the +-1 entries."""
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
                    t = abs(sum(abc[c] * Sd[c] for c in range(n)))
                    hist[t] = hist.get(t, 0) + 1
    return {k: v for k, v in hist.items() if v}


def xor_weight_quadruple(rows, quad):
    """min(#+, #-) of the entrywise product of four rows -- 'type-1' in the
    preprint's sense is minority == 4, equivalently |T4| = n - 8."""
    n = len(rows)
    p = [1] * n
    for r in quad:
        row = rows[r]
        p = [p[c] * (1 if row[c] == "+" else -1) for c in range(n)]
    plus = sum(1 for v in p if v > 0)
    return min(plus, n - plus), abs(2 * plus - n)


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 08 -- order 668: three matrices, one invariant, three classes")
    print("=" * 72)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)

    # ---------------------------------------------------------- clause 0
    print("\n[0] banked data files, SHA-256 pinned in this script")
    for name, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, name))
        check("%-38s" % name, got == want, got[:24] + "...")

    with open(os.path.join(ROOT, "data", "sep668-hall-switch.json"),
              "r", encoding="ascii") as fh:
        spec = json.load(fh)

    # ---------------------------------------------------------- clause 1
    print("\n[1] rebuild all five matrices, verify, pin")
    dec_rec = load_record(os.path.join(ROOT, "data", "payload-records.json"),
                          N)
    tw_rec = load_record(os.path.join(ROOT, "data",
                                      "sep668-twisted-record.json"), N)
    rep_d, rows_d = BGS.check_record(dec_rec)
    check("decoded   hypotheses H0-H4, D1/D3/D5, Sigma-bar, compression",
          rep_d["hypotheses_ok"] and rep_d["compression_lemma_crosscheck"],
          "s=%d i=%d w=%d" % (rep_d["s"], rep_d["i"], rep_d["w"]))
    rep_t, rows_t = BGS.check_record(tw_rec)
    check("twisted   hypotheses H0-H4, D1/D3/D5, Sigma-bar, compression",
          rep_t["hypotheses_ok"] and rep_t["compression_lemma_crosscheck"],
          "s=%d i=%d w=%d" % (rep_t["s"], rep_t["i"], rep_t["w"]))
    if not (rows_d and rows_t):
        print("\nFATAL: a record did not rebuild; nothing further is checked.")
        return 1

    # --- the preprint's construction data IS this record's -------------
    print("\n  the preprint's data, re-derived from data/payload-records.json")
    seeds = [BGS.signs(x) for x in dec_rec["seeds"]]
    Xs = [[g for g, v in enumerate(x) if v < 0] for x in seeds]
    ser = "\n".join(",".join(str(v) for v in sorted(S)) for S in Xs) + "\n"
    dig = hashlib.sha256(ser.encode("ascii")).hexdigest()
    ident = spec["identification"]
    check("X1..X4 (the negative supports of the four seeds) cardinalities",
          [len(x) for x in Xs] == list(ident["X_cardinalities"]),
          str([len(x) for x in Xs]))
    check("X1..X4 == the preprint's four subsets (pinned digest)",
          dig == ident["X_canonical_sha256"], dig[:24] + "...")
    E = [BGS.signs(r) for r in dec_rec["corner"]]
    P = [BGS.signs(r) for r in dec_rec["row_table"]]
    colT = [BGS.signs(r) for r in dec_rec["col_table"]]
    St = [list(c) for c in zip(*colT)]            # S = col_table^T
    check("K = corner, T = row_table, S = col_table^T satisfy the "
          "preprint's eq. (6): K K^T = T T^T = S S^T = 4I and K S^T = -2T",
          BGS.mat_mul_t(E, E) == BGS.eye(4, 4)
          and BGS.mat_mul_t(P, P) == BGS.eye(4, 4)
          and BGS.mat_mul_t(St, St) == BGS.eye(4, 4)
          and BGS.mat_mul_t(E, St)
          == [[-2 * P[a][b] for b in range(4)] for a in range(4)],
          "re-derived from the record, not transcribed")

    # --- the paired Hall switch ----------------------------------------
    print("\n  the paired Hall switch, mask rebuilt from the banked data")
    s = int(dec_rec["s"])
    n_g = BGS.AbelianGroup(dec_rec["group"]).n
    mask, got_rects, want_rects = build_mask(spec, s, n_g)
    check("the banked rectangles are the closed forms in (s, n) = (%d, %d)"
          % (s, n_g), got_rects == want_rects, str(got_rects))
    check("mask size == 2*(4s)*n == %d entries" % (2 * 4 * s * n_g),
          len(mask) == 2 * 4 * s * n_g == spec["switch"]["entries_changed"]
          == 1328, "%d" % len(mask))
    check("the two rectangles are disjoint",
          len(mask) == sum(r["entries"] for r in
                           spec["switch"]["mask_rectangles"]))
    grid = [list(r) for r in rows_d]
    for r, c in mask:
        grid[r][c] = "-" if grid[r][c] == "+" else "+"
    rows_h = ["".join(r) for r in grid]

    rows_hT = transpose_rows(rows_h)
    rows_tT = transpose_rows(rows_t)

    digs = {}
    for tag, rows in (("decoded", rows_d), ("twisted", rows_t),
                      ("hall", rows_h), ("hall-T", rows_hT),
                      ("twisted-T", rows_tT)):
        digs[tag] = verify_and_pin(tag, rows, SHA[tag])
    check("H* is the matrix the banked spec names",
          digs["hall"] == spec["result"]["canonical_sha256"])
    check("(H*)^T is the matrix the banked spec names",
          digs["hall-T"] == spec["result"]["transpose_canonical_sha256"])

    # --- the preprint's own published digests, reproduced ---------------
    print("\n  the preprint's published SHA-256 digests, reproduced here")
    pub = spec["published_digests"]
    check("H   one-byte-per-entry digest == the preprint's",
          paper_bytes_sha(rows_d) == pub["H"], pub["H"][:16] + "...")
    check("H*  one-byte-per-entry digest == the preprint's",
          paper_bytes_sha(rows_h) == pub["H_star"], pub["H_star"][:16] + "...")
    check("the pinned byte length is N^2 == %d" % (N * N),
          pub["byte_length"] == N * N)

    # ------------------------------------------------------- clause 1b
    print("\n[1b] controls on the switch itself")
    diffpos = {(r, c) for r in range(N) for c in range(N)
               if rows_d[r][c] != rows_h[r][c]}
    check("C1  the positions where H and H* differ ARE the mask, exactly",
          diffpos == mask, "%d positions" % len(diffpos))
    grid2 = [list(r) for r in rows_h]
    for r, c in mask:
        grid2[r][c] = "-" if grid2[r][c] == "+" else "+"
    check("C2  the switch is an involution: applying it twice gives H back",
          rows_sha256(["".join(r) for r in grid2]) == SHA["decoded"])
    for tag, rows in (("H", rows_d), ("H*", rows_h), ("H'", rows_t)):
        mino, t4 = xor_weight_quadruple(rows, (0, 1, 2, 3))
        check("C3  %-2s the border quadruple P = {1,2,3,4} is type-1 "
              "(minority 4, |T4| = 660)" % tag, mino == 4 and t4 == 660,
              "minority %d, |T4| %d" % (mino, t4))

    del rows_d, rows_t, rows_h, rows_hT, rows_tT, grid, grid2, diffpos

    # ---------------------------------------------------------- clause 2
    print("\n[2] the banked exact 4-profiles, audited in exact integers")
    prof = {}
    for (tag, impl), name in sorted(PROFILES.items()):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(tag, impl)] = p
        tot, m2 = audit(p, N, "%s/%s" % (tag, impl))
        check("%-9s %-4s  %d bins, all |T4| = 4 (mod 8), total == C(668,4), "
              "second moment == %d" % (tag, impl, len(p),
                                       second_moment_want(N)),
              all(k % 8 == 4 for k in p) and tot == c_n_4(N)
              and m2 == second_moment_want(N))
        if "second_moment" in blob:
            check("%-9s %-4s  the banked second_moment agrees" % (tag, impl),
                  int(blob["second_moment"]) == m2)
        if "matrix_canonical_sha256" in blob:
            check("%-9s %-4s  banked matrix_canonical_sha256 == the matrix "
                  "this script rebuilt" % (tag, impl),
                  blob["matrix_canonical_sha256"] == SHA[tag]
                  and blob.get("matrix_sha256",
                               SHA[tag]) == SHA[tag])

    with open(os.path.join(ROOT, TWISTED_T), "r", encoding="ascii") as fh:
        tt = json.load(fh)
    check("twisted-T bank names the matrix this script rebuilt",
          tt["canonical_sha256"] == SHA["twisted-T"])
    for impl, blk in sorted(tt["implementations"].items()):
        p = {int(k): int(v) for k, v in blk["profile"].items()}
        prof[("twisted-T", impl)] = p
        tot, m2 = audit(p, N, "twisted-T/%s" % impl)
        check("twisted-T %-4s  %d bins, all |T4| = 4 (mod 8), total == "
              "C(668,4), second moment == %d"
              % (impl, len(p), second_moment_want(N)),
              all(k % 8 == 4 for k in p) and tot == c_n_4(N)
              and m2 == second_moment_want(N))
        check("twisted-T %-4s  banked matrix_sha256 == the rebuilt matrix"
              % impl, blk["matrix_sha256"] == SHA["twisted-T"])

    for tag in ("decoded", "twisted", "hall", "hall-T", "twisted-T"):
        check("%-9s  blas == bits, bin for bin (two independent "
              "implementations)" % tag,
              prof[(tag, "blas")] == prof[(tag, "bits")],
              "%d bins" % len(prof[(tag, "blas")]))

    D = prof[("decoded", "blas")]
    W = prof[("twisted", "blas")]
    S = prof[("hall", "blas")]
    ST = prof[("hall-T", "blas")]
    WT = prof[("twisted-T", "bits")]

    # ---------------------------------------------------------- clause 3
    print("\n[3] the three pairs under the STANDARD relation "
          "(transpose not in the group)")
    print("      %-26s %8s %10s %14s" % ("pair", "differing", "same supp",
                                         "sum |T4|*delta"))
    want3 = {("H", "H'"): 26, ("H", "H*"): 27, ("H'", "H*"): 27}
    for (a, b), pa, pb, k in ((("H", "H'"), D, W, 26),
                              (("H", "H*"), D, S, 27),
                              (("H'", "H*"), W, S, 27)):
        diff, dsum, m1, same = compare(pa, pb)
        print("      %-26s %8d %10s %14d"
              % ("%s vs %s" % (a, b), len(diff), same, m1))
        check("%-2s vs %-2s  differ in exactly %d of the 80 bins"
              % (a, b, k), len(diff) == k == want3[(a, b)],
              "%d differing" % len(diff))
        check("%-2s vs %-2s  identical 80-bin support" % (a, b),
              same and len(pa) == len(pb) == 80)
        check("%-2s vs %-2s  the deltas sum to zero" % (a, b), dsum == 0)
        check("%-2s vs %-2s  the FIRST moment, which nothing forces, differs"
              % (a, b), m1 != 0, "sum |T4|*delta = %d" % m1)
    check("the extreme tail separates NOTHING (top bins agree on all three)",
          all(D[k] == W[k] == S[k] for k in sorted(D) if k >= 604),
          ", ".join("%d:%d" % (k, D[k]) for k in sorted(D) if k >= 604))
    check("bin 660 == 1 on all three -- the unique type-1 row quadruple "
          "the preprint's Lemma 3 asserts",
          D.get(660) == W.get(660) == S.get(660) == 1)

    print("\n      the 27 divergent bins, H vs H*")
    diff, _d, _m, _s = compare(D, S)
    print("      %6s %16s %16s %12s %10s"
          % ("|T4|", "H (decoded)", "H* (Hall)", "delta", "rel"))
    for k, p, q in diff:
        print("      %6d %16d %16d %+12d %10.2e"
              % (k, p, q, q - p, abs(q - p) / p))
    big = max(diff, key=lambda t: abs(t[2] - t[1]))
    print("      largest |delta| = %d at |T4| = %d, i.e. %.2e of that bin"
          % (abs(big[2] - big[1]), big[0], abs(big[2] - big[1]) / big[1]))

    # ---------------------------------------------------------- clause 4
    print("\n[4] the same three pairs under the TRANSPOSE-EXTENDED relation")
    print("    A ~* B  requires  A ~ B  or  A ~ B^T, so each pair needs TWO")
    print("    refutations; both are exact-profile comparisons.")
    print("      %-26s %10s %10s %8s" % ("pair", "vs B", "vs B^T", "verdict"))
    ext = []
    for (a, b), pa, pb, pbT in ((("H", "H'"), D, W, WT),
                                (("H", "H*"), D, S, ST),
                                (("H'", "H*"), W, S, ST)):
        d1 = compare(pa, pb)[0]
        d2 = compare(pa, pbT)[0]
        ok = bool(d1) and bool(d2)
        ext.append(ok)
        print("      %-26s %10d %10d %8s"
              % ("%s vs %s" % (a, b), len(d1), len(d2),
                 "SEPARATED" if ok else "OPEN"))
        check("%-2s vs %-2s  separated under the transpose-extended relation"
              % (a, b), ok,
              "%d bins vs %s, %d bins vs %s^T" % (len(d1), b, len(d2), b))
    check("the cross-check that is not needed but must hold: "
          "H' vs (H*)^T and (H')^T vs H* both separate",
          bool(compare(W, ST)[0]) and bool(compare(WT, S)[0]),
          "%d and %d bins" % (len(compare(W, ST)[0]),
                              len(compare(WT, S)[0])))
    check("the transposes are genuinely different objects: (H*)^T populates "
          "%d bins, not 80" % len(ST), len(ST) == 79,
          "bin(s) present in H* and absent in (H*)^T: %s"
          % sorted(set(S) - set(ST)))

    # ---------------------------------------------------------- clause 5
    print("\n[5] controls")

    print("\n  C4 -- negative control: the second-moment assert must FIRE")
    victim = dict(S)
    kk = sorted(victim)
    victim[kk[0]] -= 1
    victim[kk[-1]] += 1                    # total preserved, m2 moved
    fired = False
    try:
        audit(victim, N, "C4-corrupted")
    except AssertionError as exc:
        fired = "second-moment identity FAILED" in str(exc)
    check("C4  a total-preserving corruption is rejected by audit()", fired,
          "moved one count from |T4|=%d to |T4|=%d" % (kk[0], kk[-1]))
    check("C4  the corrupted profile still totals C(668,4) -- so ONLY the "
          "second moment could catch it", sum(victim.values()) == c_n_4(N))

    print("\n  C5 -- the comparator, exercised in both directions")
    for tag in ("decoded", "twisted", "hall", "hall-T", "twisted-T"):
        impl = "bits" if tag == "twisted-T" else "blas"
        check("C5  %-9s compared with itself: 0 differing bins" % tag,
              compare(prof[(tag, impl)], prof[(tag, impl)])[0] == [])

    print("\n  C6 -- the |T4| profile really is invariant: a seeded signed")
    print("        row+column permutation must not move it")
    import random                                            # noqa: E402
    for rows, name in ((sylvester(4), "Sylvester H(16)"),
                       (paley1(19), "Paley I H(20)")):
        n = len(rows)
        rng = random.Random(20260831)
        pr = list(range(n))
        pc = list(range(n))
        rng.shuffle(pr)
        rng.shuffle(pc)
        er = [rng.choice((1, -1)) for _ in range(n)]
        ec = [rng.choice((1, -1)) for _ in range(n)]
        moved = ["".join(("+" if (er[r] * ec[c] *
                                  (1 if rows[pr[r]][pc[c]] == "+" else -1))
                          == 1 else "-") for c in range(n)) for r in range(n)]
        check("C6  %-16s the image is Hadamard and DIFFERENT from the "
              "original" % name, is_hadamard(moved) and moved != rows,
              "%d rows and %d columns negated, both indices permuted"
              % (sum(1 for x in er if x < 0), sum(1 for x in ec if x < 0)))
        p0 = profile_straight(rows)
        p1 = profile_straight(moved)
        audit(p0, n, name)
        audit(p1, n, name + " moved")
        check("C6  %-16s the profile is UNCHANGED under the group element"
              % name, p0 == p1, "%d bins" % len(p0))

    print("\n  C7 -- the transposed-profile route, on matrices small enough")
    print("        for straight enumeration")
    for rows, name in ((sylvester(3), "Sylvester H(8)"),
                       (sylvester(4), "Sylvester H(16)"),
                       (paley1(19), "Paley I H(20)")):
        n = len(rows)
        rt = transpose_rows(rows)
        p1 = profile_straight(rows)
        p2 = profile_straight(rt)
        audit(p1, n, name)
        audit(p2, n, name + "^T")
        sym = (rows == rt)
        check("C7  %-16s transpose is Hadamard; profile(M^T) audited"
              % name, is_hadamard(rt))
        if sym:
            check("C7  %-16s is SYMMETRIC, so profile(M^T) == profile(M) is "
                  "forced -- and holds" % name, p1 == p2)
        else:
            print("      [MEAS] C7  %-16s is not symmetric; profile(M^T) == "
                  "profile(M) is %s (measured, not asserted)"
                  % (name, p1 == p2))
    print("      (C7 exists because clause 4 rests on transposed profiles.")
    print("       At 668 the route is not vacuous: (H*)^T populates 79 bins")
    print("       where H* populates 80.)")

    # ---------------------------------------------------------- close out
    shutil.rmtree(OUT, ignore_errors=True)
    print("\n" + "=" * 72)
    if FAIL:
        print("cert 08: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: order 668 carries at least THREE Hadamard equivalence")
    print("         classes -- the decoded (1,1) record H, its Lemma-T i=2")
    print("         rebuild H', and the preprint's paired-Hall-switch H*.")
    print("         Separating invariant: the exact |T4| 4-profile over all")
    print("         C(668,4) = %d row 4-subsets." % c_n_4(N))
    print("         Standard relation:  26 / 27 / 27 of 80 bins differ.")
    print("         Transpose-extended: all three pairs stay separated")
    print("         (%d, %d and %d bins against the transposed partner)."
          % (len(compare(D, WT)[0]), len(compare(D, ST)[0]),
             len(compare(W, ST)[0])))
    print("         Two independent implementations agree bin for bin on")
    print("         each recomputed profile; every profile hits the second")
    print("         moment %d to the unit." % second_moment_want(N))
    print("         PRIORITY: the >= 2 statement is the preprint's; the")
    print("         third class and these separations are this repo's.")
    print("         NOT claimed here: 716, 1676, 1772.")
    print("=" * 72)
    print("generated matrices deleted; nothing left in %s   (%.1fs)"
          % (rel(OUT), time.time() - t_start))
    return 0


if __name__ == "__main__":
    sys.exit(main())
