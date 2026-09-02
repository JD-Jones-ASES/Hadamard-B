#!/usr/bin/env python3
"""cert 15 -- the class counts at 668 and 716 hold with the TRANSPOSE
added to the group.

  THE RELATION.  A ~ B is Hadamard equivalence: B = D_r P_r A P_c D_c.
  The TRANSPOSE-EXTENDED relation is  A ~~ B  iff  A ~ B or A ~ B^T.
  Since A ~ B^T iff A^T ~ B, refuting A ~~ B takes TWO refutations:
  profile(A) != profile(B), and profile(A) != profile(B^T) -- or
  equivalently profile(A^T) != profile(B), which is the same statement
  reached from the other side.  The |T4| 4-profile is a Hadamard-
  equivalence invariant (note/NOTE-B.md S3.1, invariant I5), and the
  transpose of a Hadamard matrix is Hadamard, so profile(B^T) is that
  same invariant computed on that other matrix.

  THEOREM (716).  Let H be the decoded (s,i) = (1,1) bordered
  Goethals-Seidel record at order 716 (certs 01/11), H' its Lemma-T
  i = 2 rebuild (certs 02/11) and H'' the orientation switch of H
  (cert 14).  Then H, H', H'' are pairwise inequivalent UNDER THE
  TRANSPOSE-EXTENDED RELATION.  So ORDER 716 CARRIES AT LEAST THREE
  HADAMARD EQUIVALENCE CLASSES with the transpose in the group -- the
  statement certs 11 and 14 explicitly withheld.

  THEOREM (668).  Let H, H', H* be the decoded record, the Lemma-T
  rebuild and the Hall switch (certs 06/08), and H'' the orientation
  switch (cert 13).  Then all six pairs are separated under the
  transpose-extended relation, so ORDER 668 CARRIES AT LEAST FOUR
  HADAMARD EQUIVALENCE CLASSES with the transpose in the group.  Cert 08
  had the first three pairs; the fourth class is what this certificate
  adds, from the single new profile of (H'')^T.

  PROOF (finite, exact).  Every leg is an exact 4-profile comparison over
  all C(n,4) row 4-subsets -- C(716,4) = 10 859 143 295,
  C(668,4) = 8 222 179 035 -- in two arithmetics that agree bin for bin.
  At 716 all three transposes were computed, so each of the three pairs
  gets BOTH refutation routes: A vs B^T and A^T vs B.  At 668 the one
  new profile is (H'')^T: it refutes X ~ (H'')^T for each X in
  {H, H', H*}, which together with cert 13's X !~ H'' is exactly the
  pair of refutations the relation needs; cert 08 already holds the
  other three pairs.  []

  REMARK, not the headline.  At 716 the six matrices H, H', H'', H^T,
  (H')^T, (H'')^T are pairwise inequivalent under PLAIN Hadamard
  equivalence: all fifteen profile comparisons separate.  Six classes are
  therefore exhibited at 716 by three constructions.  The house counts
  three, because the transpose-extended relation is the one under which
  a matrix and its transpose are the same object, and it is the stricter
  count that survives either convention.

WHERE THE PROFILES COME FROM, AND WHAT A DEFAULT RUN ESTABLISHES

  The four new C(n,4) enumerations -- (H)^T, (H')^T, (H'')^T at 716 and
  (H'')^T at 668 -- were not run inside this repository.  They ran in the
  source laboratory (Hadamard-2060, experiments/inequiv/
  exact_profile_big.py, unchanged engine, 16 threads on a rented
  c2d-highcpu-16, 2026-09-02 ~10:25-10:36 UTC) under the pre-registration
  experiments/pr0042/REGISTRATION.md, FLUSHED 10:17 UTC before any matrix
  it governs was built.  The matrices themselves were built AND verified
  at the desk (experiments/pr0042/build_matrices.py, digests in its
  manifest.json) and only then uploaded; the rented machine enumerated
  and nothing else.  The output is banked in data/.  The DEFAULT path of
  this script AUDITS twenty-five banks; `--full` RECOMPUTES one profile
  here -- (H'')^T at 668, the leg the fourth class rests on.

WHAT THIS SCRIPT DOES  (default path: standard library only, seconds)

  (0) Pins the SHA-256 of all twenty-five banked files it reads.
  (1) At 716: rebuilds H (every master-theorem hypothesis re-checked),
      forms H'' by negating the twelve off-diagonal core blocks and
      checks the alternate-orientation identity, re-derives the twisted
      seeds as the psi-twist of the decoded seeds and rebuilds H'; then
      TRANSPOSES each of the three in-process, hands all six matrices to
      verify/verify.py, and pins all six canonical digests.
  (2) At 668: the same for H, H'' and (H'')^T (H', H*, (H')^T and (H*)^T
      are certs 06/08's and are not rebuilt here; their banks are bound
      by the digests those certificates pin, and the check labels say so).
  (3) Audits twenty-six banked exact 4-profiles: every bin = 4 (mod 8),
      total = C(n,4), second moment = n^3(n-1)(n-2)/24 recomputed AND
      compared against the field the bank declares, declared matrix
      digest against the in-process digest where the matrix is rebuilt
      here and against the earlier certificate's pin otherwise, and
      blas == bits bin for bin on each of the thirteen matrices.
  (4) The separation table: all fifteen pairs at 716 (twelve of them
      transpose-related and new here, three re-affirmed from certs 11 and
      14) and all fifteen at 668 (six with (H'')^T, nine re-affirmed from
      certs 06, 08 and 13).
  (5) The transpose-extended verdicts, DERIVED IN CODE from those counts:
      a pair is separated only if both of its refutations are nonzero.
  (6) Controls: the transposed-profile route on matrices small enough for
      straight enumeration (symmetric Sylvester H(8)/H(16), where
      profile(M^T) = profile(M) is FORCED and is checked; Paley H(20),
      which is not symmetric and is MEASURED, never asserted); the
      support control; the comparator in the null direction; a
      total-preserving corruption only the second moment can catch; and
      the dim V / dim W trap on all six 716 objects.
  (7) --full: recompute the (H'')^T profile at 668 here
      (certs/06-668-separation/full_recompute.py, imported, not copied)
      and compare to both banks bin for bin.

Usage:
  python certs/15-transpose-extended-668-716/run.py
  python certs/15-transpose-extended-668-716/run.py --full --impl blas
  python certs/15-transpose-extended-668-716/run.py --full
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

OUT = os.path.join(HERE, "out")

# ---------------------------------------------------------------- the pins
# Canonical SHA-256 of the matrices (the digest verify/verify.py reports).
# "here" marks the ones rebuilt in clause [1] or [2] of THIS run; the rest
# are bound by the digests certs 06 and 08 pin, where they are rebuilt.
SHA = {
    (716, "dec"):                                               # here
        "3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6",
    (716, "dec-T"):                                             # here
        "e1c4a6fa1f3cc24f0506eddac5bbb94bcbcc5eeab8ef1881c27ad9b8a60be278",
    (716, "tw"):                                                # here
        "6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7",
    (716, "tw-T"):                                              # here
        "41fe458af7fe215e59cd98985d4c6835f2364ad13896a0777f260bbccc21ea72",
    (716, "or"):                                                # here
        "a6b4f56ec98004e736f0ad74af52826aece4b4ab92750e4706e44486c1885fcd",
    (716, "or-T"):                                              # here
        "7445c760ccaa45d1012828f845acafec6b753596798cba93529bf5f6119de3ef",
    (668, "dec"):                                               # here
        "bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0",
    (668, "or"):                                                # here
        "af1c285cbe2def88427381ab3002a267321b282a9fa78ca37e72830b602953c7",
    (668, "or-T"):                                              # here
        "49f97ecfb6bdc05c16df3f46aa360202ed29a88e0618b7e6bbbe690f958538d9",
    (668, "tw"):                                                # cert 06
        "600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3",
    (668, "tw-T"):                                              # cert 08
        "32afde351e1f44aa4236cb3c406fbbcd59c5ab2cc3e02c8380e08680db2d19d7",
    (668, "hall"):                                              # cert 08
        "7f6af1d9e9c80fae52c6f178e298144d209c01f6c103db119e6c12425aa1a722",
    (668, "hall-T"):                                            # cert 08
        "565a9ca5a9db739f74215474364202a4d35fd691542b18e8ace8bcbb3c190c65",
}
REBUILT_HERE = {(716, "dec"), (716, "dec-T"), (716, "tw"), (716, "tw-T"),
                (716, "or"), (716, "or-T"),
                (668, "dec"), (668, "or"), (668, "or-T")}
PINNED_BY = {(668, "tw"): "cert 06", (668, "tw-T"): "cert 08",
             (668, "hall"): "cert 08", (668, "hall-T"): "cert 08"}

# SHA-256 of the banked data FILES themselves, so a silently edited bank is
# a hard error rather than a different theorem.  The eight *-T-* banks
# produced 2026-09-02 are this certificate's own; the rest are certs 06,
# 08, 11, 13 and 14's, reused verbatim.
FILE_PINS = {
    "data/sep716-decoded-T-exact-blas.json":
        "b2bd98a8c5a3403273408e009f88030b66a3994d6822be0e9f008db04bc512e5",
    "data/sep716-decoded-T-exact-bits.json":
        "854e73fc748cb4b28113f267fc9266ede5a468105a3c32138978cafc4baaf040",
    "data/sep716-twisted-T-exact-blas.json":
        "135189d8fe4dd619d72ff8ea0ff3cc4b94c07334c1cde4e227f4e08c4cb91fcb",
    "data/sep716-twisted-T-exact-bits.json":
        "2c13fb2e980782b8298e0bf9a37d7bbde569f9a832c5886aa2ec7e4198c1a87b",
    "data/sep716-orient-T-exact-blas.json":
        "b287dfcb5a26e5cf47f8ac8ec445725253cc1a4a51da67eaa219431a8d2da062",
    "data/sep716-orient-T-exact-bits.json":
        "be270c974ce237f537cb39a550819169849805f691120976e71ecfc673b30704",
    "data/sep668-orient-T-exact-blas.json":
        "536e0c136d16b8271a3c7916529fd57c86f6738a8695c1d691f4fbf98455ec2d",
    "data/sep668-orient-T-exact-bits.json":
        "81fae9e0a4f23116067769bffc80451f2c1d834101f9b65f90a2cccf4aab70cf",
    "data/sep716-exact-blas-decoded.json":
        "80ee1e151ec1f759d7213d500623603716b9afa6fc382a385ce6970efac35a6b",
    "data/sep716-exact-bits-decoded.json":
        "a0d5b3a65b83c39c905ec2a1d3b25ca1c58e0106b76aaa6eb3b2feee3748aeed",
    "data/sep716-exact-blas-twisted.json":
        "e2076eb890557775edaccda3c9dcbab7e585f5181e0da1bb5914913bf0749b46",
    "data/sep716-exact-bits-twisted.json":
        "c385773a7a2bf4506b94752406be787bac3d85b8a9ee18bb23668080c6afe7bc",
    "data/sep716-orient-exact-blas.json":
        "b1c6b0adf393288303f780efebef7ca40bf5d21611660d70928115ff16951cc4",
    "data/sep716-orient-exact-bits.json":
        "bd1c3c23fad9b29a169de1815a4ca44d111b89bb173ef6ac29d55905e077e69c",
    "data/sep668-exact-blas-decoded.json":
        "370fffe6c2f5dc53c09d3b74f8c09dd2bc2a39a1ac2b27fb5167ab4d3559387b",
    "data/sep668-exact-bits-decoded.json":
        "7bace61441f17b5e95fff433bdc5939da212e2b8735e8738d7ed3078fae456b7",
    "data/sep668-exact-blas-twisted.json":
        "8526b3cfa7938a9af334e23f722b1c215ffd1e318c0c713ecc3da1b91f5b3afe",
    "data/sep668-exact-bits-twisted.json":
        "f40bbb8c3906d6fc7374e3e04c2b68eaf29393e50b5662f69eee2426ed3f1e9a",
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
    "data/sep668-orient-exact-blas.json":
        "5ad283be1baea5c191d39ff3d9219e744ae209ff1df10fea5d152cebbb06c6fd",
    "data/sep668-orient-exact-bits.json":
        "f0eeeb6c89415d09fd8c3f08a13e72666ed0ee36c0066a4d9529e4a88745e726",
}

# (order, tag, impl) -> bank file.  data/sep668-twisted-T-exact.json is the
# one bank that carries BOTH implementations in a single file and is loaded
# separately, exactly as cert 08 loads it.
PROFILES = {
    (716, "dec", "blas"): "data/sep716-exact-blas-decoded.json",
    (716, "dec", "bits"): "data/sep716-exact-bits-decoded.json",
    (716, "dec-T", "blas"): "data/sep716-decoded-T-exact-blas.json",
    (716, "dec-T", "bits"): "data/sep716-decoded-T-exact-bits.json",
    (716, "tw", "blas"): "data/sep716-exact-blas-twisted.json",
    (716, "tw", "bits"): "data/sep716-exact-bits-twisted.json",
    (716, "tw-T", "blas"): "data/sep716-twisted-T-exact-blas.json",
    (716, "tw-T", "bits"): "data/sep716-twisted-T-exact-bits.json",
    (716, "or", "blas"): "data/sep716-orient-exact-blas.json",
    (716, "or", "bits"): "data/sep716-orient-exact-bits.json",
    (716, "or-T", "blas"): "data/sep716-orient-T-exact-blas.json",
    (716, "or-T", "bits"): "data/sep716-orient-T-exact-bits.json",
    (668, "dec", "blas"): "data/sep668-exact-blas-decoded.json",
    (668, "dec", "bits"): "data/sep668-exact-bits-decoded.json",
    (668, "tw", "blas"): "data/sep668-exact-blas-twisted.json",
    (668, "tw", "bits"): "data/sep668-exact-bits-twisted.json",
    (668, "hall", "blas"): "data/sep668-hall-exact-blas.json",
    (668, "hall", "bits"): "data/sep668-hall-exact-bits.json",
    (668, "hall-T", "blas"): "data/sep668-hall-T-exact-blas.json",
    (668, "hall-T", "bits"): "data/sep668-hall-T-exact-bits.json",
    (668, "or", "blas"): "data/sep668-orient-exact-blas.json",
    (668, "or", "bits"): "data/sep668-orient-exact-bits.json",
    (668, "or-T", "blas"): "data/sep668-orient-T-exact-blas.json",
    (668, "or-T", "bits"): "data/sep668-orient-T-exact-bits.json",
}
TWISTED_T_668 = "data/sep668-twisted-T-exact.json"

LABEL = {
    (716, "dec"): "H", (716, "dec-T"): "H^T",
    (716, "tw"): "H'", (716, "tw-T"): "(H')^T",
    (716, "or"): "H''", (716, "or-T"): "(H'')^T",
    (668, "dec"): "H", (668, "tw"): "H'", (668, "tw-T"): "(H')^T",
    (668, "hall"): "H*", (668, "hall-T"): "(H*)^T",
    (668, "or"): "H''", (668, "or-T"): "(H'')^T",
}

# Populated-bin counts, pinned so a drifting bank cannot quietly turn this
# certificate into a different statement.
NBINS = {(716, False): 87, (716, True): 86, (668, False): 80, (668, True): 79}

# Every pair, with its differing-bin count and the size of the union of the
# two supports.  NEW here = the comparisons that need one of the four
# profiles banked 2026-09-02; PRIOR = re-affirmed from the earlier certs.
SEP = {
    # ---- order 716: three matrices and their three transposes ----------
    (716, "dec", "tw"): (27, 87, "cert 11"),
    (716, "dec", "or"): (27, 87, "cert 14"),
    (716, "tw", "or"): (25, 87, "cert 14"),
    (716, "dec", "dec-T"): (57, 87, "NEW"),
    (716, "dec", "tw-T"): (57, 87, "NEW"),
    (716, "dec", "or-T"): (57, 87, "NEW"),
    (716, "dec-T", "tw"): (56, 87, "NEW"),
    (716, "dec-T", "or"): (56, 87, "NEW"),
    (716, "dec-T", "tw-T"): (28, 86, "NEW"),
    (716, "dec-T", "or-T"): (28, 86, "NEW"),
    (716, "tw", "tw-T"): (57, 87, "NEW"),
    (716, "tw", "or-T"): (57, 87, "NEW"),
    (716, "tw-T", "or"): (57, 87, "NEW"),
    (716, "tw-T", "or-T"): (26, 86, "NEW"),
    (716, "or", "or-T"): (57, 87, "NEW"),
    # ---- order 668: the four matrices and three transposes -------------
    (668, "dec", "tw"): (26, 80, "cert 06"),
    (668, "dec", "hall"): (27, 80, "cert 08"),
    (668, "tw", "hall"): (27, 80, "cert 08"),
    (668, "dec", "tw-T"): (50, 80, "cert 08"),
    (668, "dec", "hall-T"): (49, 80, "cert 08"),
    (668, "tw", "hall-T"): (50, 80, "cert 08"),
    (668, "dec", "or"): (27, 80, "cert 13"),
    (668, "tw", "or"): (27, 80, "cert 13"),
    (668, "hall", "or"): (26, 80, "cert 13"),
    (668, "dec", "or-T"): (50, 80, "NEW"),
    (668, "tw", "or-T"): (50, 80, "NEW"),
    (668, "hall", "or-T"): (49, 80, "NEW"),
    (668, "or", "or-T"): (50, 80, "NEW"),
    (668, "hall-T", "or-T"): (25, 79, "NEW"),
    (668, "tw-T", "or-T"): (24, 79, "NEW"),
}

# The classes each order claims, and the transpose of each.
CLASSES = {716: ("dec", "tw", "or"), 668: ("dec", "tw", "hall", "or")}
TRANSPOSE_OF = {"dec": "dec-T", "tw": "tw-T", "hall": "hall-T", "or": "or-T"}

# The one bin each order's transposes leave empty (a support fact, checked).
DROPPED_BIN = {716: 684, 668: 644}

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
    """(differing bins, delta sum, first-moment delta, union support size)."""
    ks = sorted(set(pa) | set(pb))
    diff = [(k, pa.get(k, 0), pb.get(k, 0)) for k in ks
            if pa.get(k, 0) != pb.get(k, 0)]
    return (diff,
            sum(q - p for _k, p, q in diff),
            sum(k * (q - p) for k, p, q in diff),
            len(ks))


# ======================================================================
# clauses 1 and 2 -- rebuild, switch, transpose, verify, pin
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


def verify_rows(order, tag, rows):
    """Hand a matrix to the trust chain and pin its canonical digest."""
    path = os.path.join(OUT, "H%d_%s.txt" % (order, tag))
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify", "verify.py"), path],
        capture_output=True, text=True)
    verdict = ((proc.stdout or proc.stderr).strip().splitlines()
               or ["(no output)"])[-1]
    lab = "%d %-7s" % (order, LABEL[(order, tag)])
    check("%s verify/verify.py" % lab, proc.returncode == 0, verdict[:70])
    dig = rows_sha256(rows)
    check("%s canonical sha256 == pin" % lab, dig == SHA[(order, tag)],
          dig[:24] + "...")
    check("%s verify.py reports the same digest" % lab, dig in verdict)
    os.remove(path)
    return dig


def build_and_verify(order, tag, rec):
    t0 = time.time()
    rep, rows = BGS.check_record(rec)
    lab = "%d %-7s" % (order, LABEL[(order, tag)])
    if not rows:
        check("%s hypotheses" % lab, False, str(rep.get("failures")))
        return rep, None, None
    check("%s hypotheses H0-H4, D1/D3/D5, Sigma-bar, compression lemma"
          % lab, rep["hypotheses_ok"] and rep["compression_lemma_crosscheck"],
          "s=%d i=%d w=%d" % (rep["s"], rep["i"], rep["w"]))
    dig = verify_rows(order, tag, rows)
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
# small-matrix controls
# ======================================================================

def sylvester(k):
    n = 1 << k
    return ["".join("+" if bin(x & y).count("1") % 2 == 0 else "-"
                    for y in range(n)) for x in range(n)]


def sylvester_profile_forced(n):
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
                    t = 0
                    for c in range(n):
                        t += abc[c] * Sd[c]
                    t = abs(t)
                    hist[t] = hist.get(t, 0) + 1
    return {k: v for k, v in hist.items() if v}


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="RECOMPUTE the (H'')^T exact profile at 668 here "
                         "with numpy and compare to both banks bin for bin; "
                         "without it the banked profiles are audited, not "
                         "recomputed")
    ap.add_argument("--impl", choices=("blas", "bits", "both"),
                    default="both",
                    help="which arithmetic path --full recomputes")
    args = ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 15 -- 716: three classes, and 668: four, with the TRANSPOSE")
    print("           added to the group")
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


def _rebuild_order(order, want_twisted, keep):
    """Rebuild the decoded record, the orientation switch and (at 716) the
    Lemma-T rebuild; transpose each; verify and pin all of them.  Returns
    {tag: digest} and, in `keep`, the rows --full will need."""
    dec_rec = load_record(os.path.join(ROOT, "data", "payload-records.json"),
                          order)
    rep_d, rows_d, dig_d = build_and_verify(order, "dec", dec_rec)
    if rows_d is None:
        return None
    n, s = int(rep_d["w"]) * int(rep_d["i"]), int(rep_d["s"])
    check("%d layout n = |G| = %d, s = %d, N = 4(n+s) = %d"
          % (order, n, s, 4 * (n + s)),
          4 * (n + s) == order and len(rows_d) == order)

    rows_o = orientation_switch(rows_d, n, s)
    check("%d H''     differs from H in exactly the twelve off-diagonal "
          "core blocks" % order,
          sum(1 for r in range(order) for c in range(order)
              if rows_o[r][c] != rows_d[r][c]) == 12 * n * n,
          "%d cells" % (12 * n * n))
    check("%d H''     S H'' S == the alternate-orientation array with the "
          "signed border" % order,
          alternate_orientation_identity(rows_d, rows_o, n, s))
    dig_o = verify_rows(order, "or", rows_o)

    built = {"dec": dig_d, "or": dig_o}
    mats = {"dec": rows_d, "or": rows_o}

    if want_twisted:
        tw_rec = load_record(os.path.join(ROOT, "data",
                                          "twisted-i2-records.json"), order)
        check("%d the two records are genuinely different instances" % order,
              dec_rec["seeds"] != tw_rec["seeds"]
              and int(dec_rec.get("coset_divisors", [1])[0]) == 1
              and int(tw_rec["coset_divisors"][0]) == 2,
              "decoded i=1, rebuild i=2, different seeds")
        check("%d the twisted seeds are the psi-twist of the decoded seeds, "
              "re-derived here" % order,
              list(tw_rec["seeds"]) == psi_twist_seeds(dec_rec)
              and list(tw_rec["group"]) == list(dec_rec["group"])
              and list(tw_rec["r_shift"]) == list(dec_rec["r_shift"]),
              "psi(g) = (-1)^g on Z_%d; rho = %d is odd, so psi(rho) = -1"
              % (int(dec_rec["group"][0]), int(dec_rec["r_shift"][0])))
        rep_t, rows_t, dig_t = build_and_verify(order, "tw", tw_rec)
        if rows_t is None:
            return None
        built["tw"] = dig_t
        mats["tw"] = rows_t

    # ---- the transposes, formed in-process and put through the chain ----
    for tag in [t for t in ("dec", "tw", "or") if t in mats
                and (order, TRANSPOSE_OF[t]) in SHA]:
        rows = mats[tag]
        rt = transpose_rows(rows)
        ttag = TRANSPOSE_OF[tag]
        lab = "%d %-7s" % (order, LABEL[(order, ttag)])
        check("%s is NOT the original -- the transpose is a different "
              "matrix here" % lab, rt != rows,
              "%d of %d rows differ"
              % (sum(1 for a, b in zip(rt, rows) if a != b), order))
        check("%s transposing twice returns the original, cell for cell"
              % lab, transpose_rows(rt) == rows)
        dig_t = verify_rows(order, ttag, rt)
        built[ttag] = dig_t
        mats[ttag] = rt

    keep.update({(order, t): r for t, r in mats.items()})
    return built


def _body(args, t_start):
    # ---------------------------------------------------------- clause 0
    print("\n[0] the twenty-five banked data files, SHA-256 pinned in this "
          "script")
    for name, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, name))
        check("%-42s" % name, got == want, got[:24] + "...")

    # ------------------------------------------------------- clauses 1, 2
    keep = {}
    print("\n[1] order 716: rebuild H, form H'', re-derive and rebuild H'; "
          "transpose all three; verify and pin all six")
    b716 = _rebuild_order(716, True, keep)
    if b716 is None:
        print("\nFATAL: a 716 matrix did not rebuild.")
        return 1

    print("\n[2] order 668: rebuild H, form H'', transpose H''; verify and "
          "pin.  H', H*, (H')^T and (H*)^T are certs 06/08's and are NOT")
    print("    rebuilt here -- their banks are bound by the digests those "
          "certificates pin.")
    b668 = _rebuild_order(668, False, keep)
    if b668 is None:
        print("\nFATAL: a 668 matrix did not rebuild.")
        return 1
    for tag, who in (("tw", "cert 06"), ("tw-T", "cert 08"),
                     ("hall", "cert 08"), ("hall-T", "cert 08")):
        b668[tag] = SHA[(668, tag)]
    built = {716: b716, 668: b668}

    # --------------------------------------------------------- clause 2b
    print("\n[2b] control C1 -- the dim-V trap, on the real objects")
    dims = {}
    for order in (716, 668):
        for tag in sorted(built[order]):
            if (order, tag) in keep:
                dims[(order, tag)] = dim_V_W(keep[(order, tag)])
    w716 = {k: v[1] for k, v in dims.items() if k[0] == 716}
    v716 = {k: v[0] for k, v in dims.items() if k[0] == 716}
    check("C1  dim W (INVARIANT) is 715 on all six 716 objects, which "
          "clause [5] proves pairwise inequivalent",
          set(w716.values()) == {715} and len(w716) == 6,
          "matching invariants prove nothing: this one separates none of "
          "the fifteen pairs")
    check("C1  dim V (NOT invariant) takes two values across them -- and is "
          "worthless", len(set(v716.values())) > 1,
          ", ".join("%s %d" % (LABEL[k], v) for k, v in sorted(v716.items()))
          + "  <- do NOT read this as a separation")
    w668 = {k: v[1] for k, v in dims.items() if k[0] == 668}
    check("C1  dim W is 667 on H, H'' and (H'')^T at 668", set(w668.values())
          == {667} and len(w668) == 3,
          ", ".join("%s %d" % (LABEL[k], v) for k, v in sorted(w668.items())))

    ROWS = {}
    if args.full:
        ROWS[(668, "or-T")] = keep[(668, "or-T")]
    keep.clear()

    # ---------------------------------------------------------- clause 3
    print("\n[3] the twenty-six banked exact 4-profiles, AUDITED in exact "
          "integers")
    print("    (audited, not recomputed: see --full)")
    prof = {}
    for (order, tag, impl), name in sorted(PROFILES.items(),
                                           key=lambda kv: (kv[0][0], kv[0][1],
                                                           kv[0][2])):
        with open(os.path.join(ROOT, name), "r", encoding="ascii") as fh:
            blob = json.load(fh)
        p = {int(k): int(v) for k, v in blob["profile"].items()}
        prof[(order, tag, impl)] = p
        _audit_bank(order, tag, impl, p, blob, built)

    # The 668 twisted-transpose bank carries BOTH implementations in one file.
    with open(os.path.join(ROOT, TWISTED_T_668), "r", encoding="ascii") as fh:
        tt = json.load(fh)
    check("668 (H')^T  bank's canonical_sha256 == the digest cert 08 pins",
          tt["canonical_sha256"] == SHA[(668, "tw-T")])
    check("668 (H')^T  bank names the SOURCE matrix cert 06 pins",
          tt.get("source_matrix_canonical_sha256",
                 SHA[(668, "tw")]) == SHA[(668, "tw")])
    for impl, blk in sorted(tt["implementations"].items()):
        p = {int(k): int(v) for k, v in blk["profile"].items()}
        prof[(668, "tw-T", impl)] = p
        _audit_bank(668, "tw-T", impl, p, blk, built,
                    name_field="matrix_sha256")

    for order in (716, 668):
        for tag in sorted(set(t for (o, t, _i) in prof if o == order)):
            a, b = prof[(order, tag, "blas")], prof[(order, tag, "bits")]
            check("%d %-7s blas == bits, bin for bin (two independent "
                  "implementations)" % (order, LABEL[(order, tag)]), a == b,
                  "%d bins" % len(a))
            istrans = tag.endswith("-T")
            check("%d %-7s populates %d bins, all |T4| = %d (mod 8)"
                  % (order, LABEL[(order, tag)], NBINS[(order, istrans)],
                     order % 8),
                  len(a) == NBINS[(order, istrans)]
                  and all(k % 8 == order % 8 for k in a))

    # ---------------------------------------------------------- clause 4
    print("\n[4] the separation table -- every pair, in both arithmetics")
    for order in (716, 668):
        keys = sorted((k for k in SEP if k[0] == order),
                      key=lambda k: (SEP[k][2] == "NEW", k[1], k[2]))
        print("\n    order %d  (%d pairs: the ones the theorems and their "
              "re-affirmations use)" % (order, len(keys)))
        print("      %-22s %10s %9s %14s  %s"
              % ("pair", "differing", "union", "sum |T4|*delta", "source"))
        for (_o, a, b) in keys:
            want, wantu, src = SEP[(order, a, b)]
            diff, dsum, m1, u = compare(prof[(order, a, "blas")],
                                        prof[(order, b, "blas")])
            dbits = compare(prof[(order, a, "bits")],
                            prof[(order, b, "bits")])[0]
            print("      %-22s %10d %9d %14d  %s"
                  % ("%s vs %s" % (LABEL[(order, a)], LABEL[(order, b)]),
                     len(diff), u, m1, src))
            check("%d %-7s vs %-7s  differ in exactly %d of %d bins, "
                  "deltas summing to zero, blas and bits alike"
                  % (order, LABEL[(order, a)], LABEL[(order, b)],
                     want, wantu),
                  len(diff) == want and u == wantu and dsum == 0
                  and len(dbits) == want,
                  "%d differing" % len(diff))

    print("\n      the divergent bins of the two new legs the theorems rest "
          "on (first six each)")
    for order, a, b in ((716, "dec", "or-T"), (668, "dec", "or-T")):
        diff = compare(prof[(order, a, "blas")],
                       prof[(order, b, "blas")])[0]
        big = max(diff, key=lambda t: abs(t[2] - t[1]))
        print("      %d  %s vs %s: largest |delta| = %d at |T4| = %d "
              "(%.2e of that bin)"
              % (order, LABEL[(order, a)], LABEL[(order, b)],
                 abs(big[2] - big[1]), big[0],
                 abs(big[2] - big[1]) / big[1]))
        for k, p, q in diff[:6]:
            print("        |T4|=%4d  %-8s %16d   %-8s %16d   delta %+d"
                  % (k, LABEL[(order, a)], p, LABEL[(order, b)], q, q - p))

    # ---------------------------------------------------------- clause 5
    print("\n[5] the TRANSPOSE-EXTENDED verdicts, derived here from those "
          "counts")
    print("    A ~~ B requires A ~ B or A ~ B^T, so each pair needs TWO")
    print("    refutations.  A ~ B^T is refuted by profile(A) != profile(B^T)")
    print("    and, equivalently, by profile(A^T) != profile(B); where both")
    print("    transposes are banked, both routes are shown and both must")
    print("    hold.")

    def ndiff(order, x, y):
        return len(compare(prof[(order, x, "blas")],
                           prof[(order, y, "blas")])[0])

    verdicts = {}
    for order in (716, 668):
        cls = CLASSES[order]
        print("\n    order %d" % order)
        print("      %-16s %8s %10s %12s %10s"
              % ("pair", "A vs B", "A vs B^T", "A^T vs B", "verdict"))
        for i, a in enumerate(cls):
            for b in cls[i + 1:]:
                d0 = ndiff(order, a, b)
                d1 = ndiff(order, a, TRANSPOSE_OF[b])
                have_aT = (order, TRANSPOSE_OF[a]) in SHA
                d2 = ndiff(order, TRANSPOSE_OF[a], b) if have_aT else None
                ok = d0 > 0 and d1 > 0
                verdicts[(order, a, b)] = ok
                print("      %-16s %8d %10d %12s %10s"
                      % ("%s vs %s" % (LABEL[(order, a)], LABEL[(order, b)]),
                         d0, d1, ("%d" % d2) if have_aT else "n/a",
                         "SEPARATED" if ok else "OPEN"))
                check("%d %-4s ~~ %-4s is REFUTED: %d bins against B and %d "
                      "against B^T"
                      % (order, LABEL[(order, a)], LABEL[(order, b)], d0, d1),
                      ok)
                if have_aT:
                    check("%d %-4s vs %-4s  the other route, A^T vs B, "
                          "separates too (%d bins)"
                          % (order, LABEL[(order, a)], LABEL[(order, b)], d2),
                          d2 > 0)
        n_cls = len(cls)
        check("order %d: all %d pairs of the %d classes are separated under "
              "the transpose-extended relation"
              % (order, n_cls * (n_cls - 1) // 2, n_cls),
              all(v for (o, _a, _b), v in verdicts.items() if o == order))

    print("\n    the remark, stated and not headlined: under PLAIN Hadamard")
    print("    equivalence the six 716 objects are pairwise inequivalent.")
    tags716 = ["dec", "dec-T", "tw", "tw-T", "or", "or-T"]
    pairs716 = [(a, b) for i, a in enumerate(tags716) for b in tags716[i + 1:]]
    worst = min(ndiff(716, a, b) for a, b in pairs716)
    check("716 all %d pairs of {H, H^T, H', (H')^T, H'', (H'')^T} separate "
          "(least separated pair: %d bins)" % (len(pairs716), worst),
          worst > 0,
          "six classes exhibited by three constructions; the house counts "
          "the three")

    # ---------------------------------------------------------- clause 6
    print("\n[6] controls")

    print("\n  C2 -- the transposed-profile route, on matrices small enough")
    print("        for straight O(C(n,4)) enumeration")
    for rows, name in ((sylvester(3), "Sylvester H(8)"),
                       (sylvester(4), "Sylvester H(16)"),
                       (paley1(19), "Paley I H(20)")):
        n = len(rows)
        rt = transpose_rows(rows)
        p1 = profile_straight(rows)
        p2 = profile_straight(rt)
        audit(p1, n, name)
        audit(p2, n, name + "^T")
        check("C2  %-16s the transpose is Hadamard; profile(M^T) audited"
              % name, is_hadamard(rt))
        if rows == rt:
            check("C2  %-16s is SYMMETRIC, so profile(M^T) == profile(M) is "
                  "FORCED -- and holds" % name, p1 == p2, "%d bins" % len(p1))
            check("C2  %-16s matches the forced Sylvester profile" % name,
                  p1 == sylvester_profile_forced(n))
        else:
            print("      [MEAS] C2  %-16s is not symmetric; profile(M^T) == "
                  "profile(M) is %s (measured, never asserted)"
                  % (name, p1 == p2))
    print("      (C2 exists because clause [5] rests on transposed profiles.")
    print("       At these orders the route is not vacuous: every transpose")
    print("       banked here populates one bin FEWER than its original.)")

    print("\n  C3 -- the transposes are genuinely different objects")
    for order in (716, 668):
        tags = sorted(set(t for (o, t, _i) in prof if o == order))
        originals = [t for t in tags if not t.endswith("-T")]
        common = set(prof[(order, originals[0], "blas")])
        for t in originals[1:]:
            common &= set(prof[(order, t, "blas")])
        check("C3  %d |T4| = %d is populated in every original at this order"
              % (order, DROPPED_BIN[order]), DROPPED_BIN[order] in common,
              "%d originals" % len(originals))
        for tag in [t for t in tags if t.endswith("-T")]:
            gone = sorted(common - set(prof[(order, tag, "blas")]))
            check("C3  %d %-7s populates %d bins, not %d -- and the bin it "
                  "loses is |T4| = %d"
                  % (order, LABEL[(order, tag)], NBINS[(order, True)],
                     NBINS[(order, False)], DROPPED_BIN[order]),
                  gone == [DROPPED_BIN[order]]
                  and len(prof[(order, tag, "blas")]) == NBINS[(order, True)],
                  "missing: %s" % gone)

    print("\n  C4 -- the comparator, exercised in the null direction")
    for order in (716, 668):
        tags = sorted(set(t for (o, t, _i) in prof if o == order))
        check("C4  %d every banked profile against itself: 0 differing bins"
              % order,
              all(compare(prof[(order, t, i)], prof[(order, t, i)])[0] == []
                  for t in tags for i in ("blas", "bits")),
              "%d profiles" % (2 * len(tags)))

    print("\n  C5 -- negative control: the second-moment assert must FIRE")
    victim = dict(prof[(668, "or-T", "blas")])
    kk = sorted(victim)
    victim[kk[0]] -= 1
    victim[kk[-1]] += 1                     # total preserved, m2 moved
    fired = False
    try:
        audit(victim, 668, "C5-corrupted")
    except AssertionError as exc:
        fired = "second-moment identity FAILED" in str(exc)
    check("C5  a total-preserving corruption of the new (H'')^T bank is "
          "rejected by audit()", fired,
          "moved one count from |T4|=%d to |T4|=%d" % (kk[0], kk[-1]))
    check("C5  the corrupted profile still totals C(668,4) -- so ONLY the "
          "second moment could catch it", sum(victim.values()) == c_n_4(668))

    # ---------------------------------------------------------- clause 7
    replayed = []
    if args.full:
        impls = (("blas", "bits") if args.impl == "both" else (args.impl,))
        print("\n[7] --full: RECOMPUTING the (H'')^T exact profile at 668 "
              "here, from")
        print("    the rows clause [2] verified, with numpy (%s).  This is "
              "the leg" % " and ".join(impls))
        print("    the 668 fourth class's transpose-robustness rests on.")
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
        rows = ROWS[(668, "or-T")]
        check("[full] the rows about to be enumerated are the ones verify.py "
              "accepted", rows_sha256(rows) == SHA[(668, "or-T")])
        for impl in impls:
            t0 = time.time()
            got = FR.profile(rows, 668, impl)
            audit(got, 668, "full/668-or-T/%s" % impl)
            secs = time.time() - t0
            for bimpl in ("blas", "bits"):
                check("[full] (H'')^T at 668 recomputed %-4s == banked %-4s, "
                      "bin for bin" % (impl, bimpl),
                      got == prof[(668, "or-T", bimpl)],
                      "%d bins, %.0fs" % (len(got), secs))
            replayed.append("668 (H'')^T/%s" % impl)
        ROWS.clear()
    else:
        print("\n[7] --full not requested: the banked profiles were AUDITED, "
              "not recomputed.")

    # ---------------------------------------------------------- close out
    print("\n" + "=" * 72)
    print("checks run: %d" % NCHECK[0])
    if FAIL:
        print("cert 15: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: with the TRANSPOSE added to the group --")
    print("         ORDER 716 CARRIES AT LEAST THREE HADAMARD EQUIVALENCE")
    print("         CLASSES (H, H', H''), the statement certs 11 and 14")
    print("         withheld; and ORDER 668 CARRIES AT LEAST FOUR")
    print("         (H, H', H*, H''), extending cert 08's form to cert 13's")
    print("         fourth class.  LABEL: PROVEN.")
    print("         Every pair carries two refutations, each an exact")
    print("         |T4| 4-profile comparison in two arithmetics that agree")
    print("         bin for bin, on profiles hitting %d (716)"
          % second_moment_want(716))
    print("         and %d (668) to the unit." % second_moment_want(668))
    if replayed:
        print("         PROFILE: the new 668 leg RECOMPUTED here and matched")
        print("         to the bank (%s)." % ", ".join(replayed))
    else:
        print("         PROFILES: banked exact computation AUDITED, not")
        print("         recomputed (Hadamard-2060, exact_profile_big.py,")
        print("         2026-09-02, experiments/pr0042/REGISTRATION.md).")
    print("         NOT claimed: anything at 2060, 1676 or 1772 (those legs")
    print("         are still running); that three and four are the counts;")
    print("         any novelty or priority at either order.")
    print("=" * 72)
    return 0


def _audit_bank(order, tag, impl, p, blob, built, name_field=None):
    """The forced identities, the declared second moment, and the matrix
    binding -- against the in-process digest where this run rebuilt the
    matrix, and against the earlier certificate's pin where it did not."""
    lab = "%d %-7s %-4s" % (order, LABEL[(order, tag)], impl)
    tot, m2 = audit(p, order, "%d/%s/%s" % (order, tag, impl))
    check("%s  total == C(%d,4) == %d and second moment == "
          "n^3(n-1)(n-2)/24 == %d"
          % (lab, order, c_n_4(order), second_moment_want(order)),
          tot == c_n_4(order) and m2 == second_moment_want(order))
    if "second_moment" in blob:
        check("%s  the banked second_moment agrees" % lab,
              int(blob["second_moment"]) == m2)
    if blob.get("n", order) != order or \
            int(blob.get("C_n_4", c_n_4(order))) != c_n_4(order):
        check("%s  banked n / C(n,4) header" % lab, False)
    declared = (blob.get(name_field) if name_field else None) or \
        blob.get("matrix_canonical_sha256") or blob.get("matrix_sha256")
    how = ("rebuilt in THIS run" if (order, tag) in REBUILT_HERE
           else "pinned by %s" % PINNED_BY[(order, tag)])
    check("%s  bank names the matrix %s" % (lab, how),
          is_sha256(declared) and declared == built[order][tag],
          (declared[:24] + "...") if is_sha256(declared)
          else "declared = %r" % (declared,))


if __name__ == "__main__":
    sys.exit(main())
