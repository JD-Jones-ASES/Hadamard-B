#!/usr/bin/env python3
"""cert 03 -- the Theorem-D gate: H(52) on G = Z2 x Z2 x Z3.

Run from the repository root:

    python certs/03-h52-gate/run.py

`data/h52-gate.json` banks a complete bordered-GS parameter record for a
from-scratch `s = 1, i = 2` instance on the NON-CYCLIC group
`G = Z2 x Z2 x Z3` (n = 12, N = 4(n+1) = 52), with `K = ker(t -> t_0)`
of index 2 (w = 6, so w > 2s) and reflection `rho = (0,0,0)`, hence
`kappa(rho) = 0` and `eps = +1`.  This is the branch of Theorem D that
the four decoded `s = 1` instances do not exercise: they are all cyclic
with odd `rho`, i.e. `eps = -1`.

This cert

  1. runs `tools/bordered_gs.py::check_record` -- every hypothesis of
     Theorem A/B and every derived law of Theorem C -- in exact stdlib
     integer arithmetic;
  2. asserts, here, the gate's declared structure (non-cyclic G, i = 2,
     K = ker of the first coordinate, w = 6 > 2s, eps = +1) and the
     clauses of Theorem D it was built to exercise: (D-a) the Gram is
     the house form, (D-b) Q is the +-doubling of a 4x4 Hadamard U,
     (D-c) P is antisymmetric across each coset pair with p 4x4
     Hadamard and E Hadamard, (D-d) the closed form
     E = -(1/4) p Lambda(d)^T U with d = (delta_0, eps*delta_1,
     eps*delta_2, eps*delta_3), and the forced Parseval value
     sum_q delta_q^2 = 4;
  3. assembles, hands the file to `verify/verify.py`, compares the
     canonical SHA-256 against the digest pinned in PINNED_SHA (coded
     comparison, hard failure on mismatch), and deletes the matrix.

Stdlib only.  Exact integers only.  No network.  Nothing is written
inside the repository.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

DATA = os.path.join("data", "h52-gate.json")
VERIFY = os.path.join("verify", "verify.py")
TOOLS = "tools"

PINNED_SHA = \
    "e2c3e48b0fc65f5283e833096824b4fec651d8c57694ae45b3842c23c87ad7ca"


def die(msg):
    print("FAIL: %s" % msg)
    sys.exit(1)


def lam(y):
    """Lambda(y) of Theorem D; Lambda(y) Lambda(y)^T = (sum y_q^2) I_4."""
    return [[y[0], y[1], y[2], y[3]],
            [-y[1], y[0], y[3], -y[2]],
            [-y[2], -y[3], y[0], y[1]],
            [-y[3], y[2], -y[1], y[0]]]


def main():
    for p in (DATA, VERIFY, os.path.join(TOOLS, "bordered_gs.py")):
        if not os.path.isfile(p):
            die("missing %s -- run this from the repository root: "
                "python certs/03-h52-gate/run.py" % p)
    sys.path.insert(0, TOOLS)
    import bordered_gs as T

    with open(DATA, "r", encoding="ascii") as fh:
        rec = json.load(fh)["record"]
    bad = []
    t0 = time.time()

    # --- declared structure ---
    if list(rec["group"]) != [2, 2, 3]:
        bad.append("group is %s, expected [2, 2, 3]" % (rec["group"],))
    if len(rec["group"]) < 2:
        bad.append("G is cyclic; the gate requires a non-cyclic G")
    if list(rec["coset_divisors"]) != [2, 1, 1]:
        bad.append("coset_divisors %s: K is not ker(t -> t_0)"
                   % (rec["coset_divisors"],))
    if list(rec["r_shift"]) != [0, 0, 0]:
        bad.append("r_shift %s: the eps = +1 branch needs kappa(rho) = 0"
                   % (rec["r_shift"],))
    if rec.get("gs_variant", "standard") != "standard":
        bad.append("gs_variant is not standard")

    # --- Theorem A/B/C, by the shared checker ---
    rep, rows = T.check_record(rec)
    if not rows or not rep["hypotheses_ok"]:
        bad.append("hypotheses failed: %s" % (rep.get("failures"),))
        for b in bad:
            print("FAIL: %s" % b)
        return 1
    if (rep["n"], rep["s"], rep["i"], rep["w"]) != (12, 1, 2, 6):
        bad.append("expected (n,s,i,w) = (12,1,2,6), got (%d,%d,%d,%d)"
                   % (rep["n"], rep["s"], rep["i"], rep["w"]))
    if not rep["D3_applicable_w_gt_2s"]:
        bad.append("w > 2s does not hold")
    if not (rep["E_is_Hadamard"] and rep["PPt_eq_4i_I"]):
        bad.append("D3 forcing failed under w > 2s")
    if not rep["compression_lemma_crosscheck"]:
        bad.append("compression-lemma cross-check failed")

    # --- Theorem D, clause by clause ---
    G = T.AbelianGroup(rec["group"])
    kappa, i = T.coset_map(G, list(rec["coset_divisors"]))
    rho = G.idx(tuple(rec["r_shift"]))
    eps = 1 if kappa[rho] == 0 else -1
    if eps != 1:
        bad.append("eps = %d; the gate is the eps = +1 branch" % eps)

    E = [T.signs(r) for r in rec["corner"]]
    P = [T.signs(r) for r in rec["row_table"]]
    colT = [T.signs(r) for r in rec["col_table"]]
    Q = [[colT[r][k] for r in range(4)] for k in range(8)]

    # (D-a) the Gram is the house form M = 8 I_2 - 4 J_2 (the genuine i=2
    # branch), not the degenerate M = 4 J_2.
    QQ = T.mat_mul_t(Q, Q)
    house = [[0] * 8 for _ in range(8)]
    for I in range(4):
        for c1 in range(2):
            for c2 in range(2):
                house[2 * I + c1][2 * I + c2] = (8 if c1 == c2 else 0) - 4
    if QQ != house:
        bad.append("(D-a) Q Q^T is not I_4 (x) (8 I_2 - 4 J_2)")

    # (D-b) Q is the +- doubling of a 4x4 Hadamard U
    if any(Q[2 * I + 1][c] != -Q[2 * I][c]
           for I in range(4) for c in range(4)):
        bad.append("(D-b) Q[2I+1] = -Q[2I] fails")
    U = [Q[2 * I] for I in range(4)]
    if T.mat_mul_t(U, U) != T.eye(4, 4):
        bad.append("(D-b) U is not a 4x4 Hadamard matrix")

    # (D-c) P[r][2J+1] = -P[r][2J]; p and E are 4x4 Hadamard
    if any(P[r][2 * J + 1] != -P[r][2 * J]
           for r in range(4) for J in range(4)):
        bad.append("(D-c) P[r][2J+1] = -P[r][2J] fails")
    p = [[P[r][2 * J] for J in range(4)] for r in range(4)]
    if T.mat_mul_t(p, p) != T.eye(4, 4):
        bad.append("(D-c) p is not a 4x4 Hadamard matrix")
    if T.mat_mul_t(E, E) != T.eye(4, 4):
        bad.append("(D-c) E is not a 4x4 Hadamard matrix")

    # (D-d) the closed form and the forced Parseval value
    sigma = rep["sigma"]                       # sigma[q][c], c in {0,1}
    delta = [sigma[q][0] - sigma[q][1] for q in range(4)]
    if sum(v * v for v in delta) != 4:
        bad.append("(D-d) sum_q delta_q^2 = %d, must be 4"
                   % sum(v * v for v in delta))
    d = [delta[0], eps * delta[1], eps * delta[2], eps * delta[3]]
    L = lam(d)
    pl = T.mat_mul_t(p, L)                     # p * Lambda(d)^T
    Eclosed = [[-sum(pl[r][k] * U[k][c] for k in range(4)) for c in range(4)]
               for r in range(4)]
    if any(v % 4 for row in Eclosed for v in row):
        bad.append("(D-d) -p Lambda(d)^T U is not divisible by 4")
    else:
        Eclosed = [[v // 4 for v in row] for row in Eclosed]
        if Eclosed != E:
            bad.append("(D-d) E != -(1/4) p Lambda(d)^T U; closed form gives "
                       "%s" % (Eclosed,))

    # --- trust chain + pinned digest ---
    tmp = tempfile.mkdtemp(prefix="cert03-")
    try:
        path = os.path.join(tmp, "H52_gate_i2.txt")
        with open(path, "w", encoding="ascii", newline="\n") as fh:
            fh.write("\n".join(rows) + "\n")
        proc = subprocess.run([sys.executable, VERIFY, path],
                              capture_output=True, text=True)
        line = ((proc.stdout or proc.stderr).strip().splitlines()
                or ["(no output)"])[-1]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if proc.returncode != 0 or "canonical_sha256=" not in line:
        bad.append("verify.py rc=%d: %s" % (proc.returncode, line))
    else:
        got = line.rsplit("canonical_sha256=", 1)[1].strip()
        if got != PINNED_SHA:
            bad.append("DIGEST MISMATCH got=%s pinned=%s" % (got, PINNED_SHA))
        if got != rec.get("pinned_sha256"):
            bad.append("banked pinned_sha256 %s disagrees with the observed "
                       "digest %s" % (rec.get("pinned_sha256"), got))

    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 03: FAIL (%d problems)" % len(bad))
        return 1
    print("  H(52)  G = Z2 x Z2 x Z3 (non-cyclic)  n=12 s=1 i=2 w=6  "
          "eps=+1  delta=%s  d=%s" % (delta, d))
    print("  Theorem D: (D-a) house Gram OK  (D-b) Q = +-doubling of a 4x4 "
          "Hadamard OK  (D-c) P antisymmetric, p and E Hadamard OK  "
          "(D-d) E = -(1/4) p Lambda(d)^T U OK  sum delta^2 = 4 OK")
    print("  %s" % line)
    print("CERT 03: PASS -- the Theorem-D gate verifies and matches its "
          "pinned digest (%.1fs)" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
