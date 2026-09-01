#!/usr/bin/env python3
"""cert 10 -- the Theorem-D census: the 768^2 border census at (s,i) = (1,2).

Run from the repository root:

    python certs/10-theorem-d-census/run.py

NOTE-B.md Sec.1.5 (D-d) reduces the whole `s = 1, i = 2` border to one 4x4
equation,

    E = -(1/4) * p * Lambda(d)^T * U,

with `p, U` 4x4 Hadamard and `d` an integer vector with sum_q d_q^2 = 4.
(D-e) then reports a census over all 768^2 ordered pairs `(p, U)` and all
24 such arguments.  This cert replays that census in full, from scratch,
and replays the two compressed-block identities that (D-a') and (D-d)
stand on.  It banks no data and reads no data: every number below is
enumerated here.

What it checks:

  1. Lambda(y) Lambda(y)^T = (sum_q y_q^2) I_4.  Each entry of the
     difference is a polynomial of degree <= 2 in y with degree <= 2 in
     each variable, so vanishing on the 5^4 grid {-2..2}^4 (5 > 2 values
     per variable) proves the identity over Z.
  2. The 24 admissible arguments: sum_q y_q^2 = 4 forces |y_q| <= 2, so
     the box {-2..2}^4 is exhaustive.  Exactly 8 even (+-2 e_j) and 16
     odd ((+-1,+-1,+-1,+-1)).
  3. The compressed-block identities over Ghat = Z_2, for both
     eps = +1 (kappa(rho) = 0) and eps = -1, built by the repository's
     own `tools/bordered_gs.py::gs_array`:
         Chat[(I,1),(J,c')] - Chat[(I,0),(J,c')] = -(-1)^c' Lambda(d)[I][J]
         Chat[(I,0),(J,0)] - Chat[(I,0),(J,1)] =            Lambda(d)[I][J]
         Chat[(I,c),(J,0)] + Chat[(I,c),(J,1)] =            Lambda(r)[I][J]
     with r_q = sigma_q(0) + sigma_q(1), delta_q = sigma_q(0) - sigma_q(1),
     d = (delta_0, eps delta_1, eps delta_2, eps delta_3).  Both sides are
     linear in the eight coset sums, so vanishing on {-1,0,1}^8 (3 > 1
     values per variable) proves the identities over Z.
  4. Exactly 768 of the 2^16 sign matrices of order 4 are Hadamard; the
     rows of each lie in one weight-parity class, so do the columns, and
     both classes split 384/384.
  5. The census itself, over all 24 x 768^2 = 14 155 776 triples
     (d, p, U), by THREE evaluators that must agree pair by pair.  They
     are one arithmetic route -- W = p Lambda(d)^T, computed once --
     read under two independent bookkeeping schemes, plus a structural
     predictor over the same W.  They are NOT three independent
     arithmetic routes; the arithmetic itself is cross-checked by the
     interpolation identities of 1 and 3 and by the exact aggregate
     counts, not by disagreement among these three:
         naive     -- test <W_r, u> in {+4, -4} for every row r of W and
                      every column u of U, pair by pair.  This is the
                      entrywise +-1 test on E = -(1/4) p Lambda(d)^T U
                      without forming E or dividing by 4; E is built
                      explicitly only in the 1 152-sample corner
                      spot-check;
         table     -- precompute <W_r, v> for all 16 sign vectors v,
                      then test the four columns of U by bitmask;
         mechanism -- the structural prediction of (D-e): rows of W/2 are
                      either +-2 e_k spikes (which accept every U) or sign
                      vectors of one weight-parity class (which accept
                      exactly the U whose columns are the other class).
     Expected counts: 768^2/2 for each even argument, 3*768^2/4 for each
     odd one.
  6. Negative controls: two single-cell corruptions of the Lambda table
     (one sign, one variable index) must move the census.  With
     `--negative-control[=sign|index]` the corruption is installed into
     the main census instead, and this cert must then FAIL.

Stdlib only.  Exact integers only.  No floating point.  No network.
Nothing is written anywhere.
"""

import itertools
import os
import sys
import time

TOOLS = "tools"

# Lambda as a (sign, variable-index) table -- NOTE-B.md Sec.1.5:
#
#           [  y0   y1   y2   y3 ]
#   Lambda= [ -y1   y0   y3  -y2 ]
#           [ -y2  -y3   y0   y1 ]
#           [ -y3   y2  -y1   y0 ]
LAM_TABLE = (
    ((1, 0), (1, 1), (1, 2), (1, 3)),
    ((-1, 1), (1, 0), (1, 3), (-1, 2)),
    ((-1, 2), (-1, 3), (1, 0), (1, 1)),
    ((-1, 3), (1, 2), (-1, 1), (1, 0)),
)


def corrupt(table, i, j, mode):
    """One deliberately wrong cell of the Lambda table."""
    rows = [list(r) for r in table]
    sg, k = rows[i][j]
    rows[i][j] = (-sg, k) if mode == "sign" else (sg, (k + 1) % 4)
    return tuple(tuple(r) for r in rows)


# the two negative controls, as (name, table, an argument that must move)
NEG_SIGN = ("sign flip at Lambda[1][3] (-y2 -> +y2)",
            corrupt(LAM_TABLE, 1, 3, "sign"), (1, 1, 1, 1))
NEG_INDEX = ("index substitution at Lambda[0][0] (y0 -> y1)",
             corrupt(LAM_TABLE, 0, 0, "idx"), (2, 0, 0, 0))

VEC = [tuple(1 if (m >> k) & 1 == 0 else -1 for k in range(4))
       for m in range(16)]
IDX = dict((v, j) for j, v in enumerate(VEC))


def die(msg):
    print("FAIL: %s" % msg)
    sys.exit(1)


def lam(y, table=LAM_TABLE):
    return [[sg * y[k] for (sg, k) in row] for row in table]


def matmul(A, B):
    return [[sum(A[a][k] * B[k][b] for k in range(len(B)))
             for b in range(len(B[0]))] for a in range(len(A))]


def transpose(A):
    return [list(col) for col in zip(*A)]


def wclass(v):
    """Weight-parity class of a +-1 vector: 0 if an even number of -1s."""
    return sum(1 for x in v if x < 0) & 1


# --------------------------------------------------------------- pieces

def check_lambda_identity(bad):
    """Lambda(y) Lambda(y)^T = (sum y^2) I_4, proved on the {-2..2}^4 grid."""
    seen = 0
    for y in itertools.product(range(-2, 3), repeat=4):
        L = lam(y)
        G = matmul(L, transpose(L))
        nrm = sum(v * v for v in y)
        for a in range(4):
            for b in range(4):
                want = nrm if a == b else 0
                if G[a][b] != want:
                    bad.append("Lambda identity fails at y=%s, (%d,%d): "
                               "%d != %d" % (y, a, b, G[a][b], want))
                    return seen
        seen += 1
    return seen


def admissible_arguments(bad):
    """The integer y with sum y_q^2 = 4.  |y_q| <= 2, so the box is total."""
    args = [y for y in itertools.product(range(-2, 3), repeat=4)
            if sum(v * v for v in y) == 4]
    even = [y for y in args if all(v % 2 == 0 for v in y)]
    odd = [y for y in args if all(v % 2 for v in y)]
    if len(args) != 24 or len(even) != 8 or len(odd) != 16:
        bad.append("arguments: %d total (%d even, %d odd), expected 24 (8, 16)"
                   % (len(args), len(even), len(odd)))
    want_even = set()
    for j in range(4):
        for sg in (1, -1):
            want_even.add(tuple(sg * 2 if k == j else 0 for k in range(4)))
    if set(even) != want_even:
        bad.append("the even arguments are not exactly +-2 e_j")
    if set(odd) != set(itertools.product((1, -1), repeat=4)):
        bad.append("the odd arguments are not exactly (+-1,+-1,+-1,+-1)")
    return args, even, odd


def check_compressed_blocks(bad):
    """The (D-a')/(D-d) compressed-block identities over Ghat = Z_2.

    Both sides are linear in the eight coset sums sigma_q(c); vanishing on
    {-1,0,1}^8 (3 > 1 values per variable) proves them over Z.
    """
    sys.path.insert(0, TOOLS)
    import bordered_gs as T
    Gq = T.QuotientGroup([2])
    seen = 0
    for rho_bar in (0, 1):
        eps = 1 if rho_bar == 0 else -1
        for sig in itertools.product((-1, 0, 1), repeat=8):
            seq = [[sig[2 * q], sig[2 * q + 1]] for q in range(4)]
            Chat = T.gs_array(Gq, seq, rho_bar)
            r = [seq[q][0] + seq[q][1] for q in range(4)]
            dl = [seq[q][0] - seq[q][1] for q in range(4)]
            d = [dl[0], eps * dl[1], eps * dl[2], eps * dl[3]]
            Ld, Lr = lam(d), lam(r)
            for I in range(4):
                for J in range(4):
                    for cp in range(2):
                        got = (Chat[2 * I + 1][2 * J + cp]
                               - Chat[2 * I][2 * J + cp])
                        want = -((-1) ** cp) * Ld[I][J]
                        if got != want:
                            bad.append("row-difference identity fails: "
                                       "eps=%d sigma=%s (I,J,c')=(%d,%d,%d) "
                                       "%d != %d"
                                       % (eps, seq, I, J, cp, got, want))
                            return seen
                    if Chat[2 * I][2 * J] - Chat[2 * I][2 * J + 1] != Ld[I][J]:
                        bad.append("(D-d) column-difference identity fails: "
                                   "eps=%d sigma=%s (I,J)=(%d,%d)"
                                   % (eps, seq, I, J))
                        return seen
                    for c in range(2):
                        got = (Chat[2 * I + c][2 * J]
                               + Chat[2 * I + c][2 * J + 1])
                        if got != Lr[I][J]:
                            bad.append("row-sum identity fails: eps=%d "
                                       "sigma=%s (I,J,c)=(%d,%d,%d) %d != %d"
                                       % (eps, seq, I, J, c, got, Lr[I][J]))
                            return seen
            seen += 1
    return seen


def hadamard4(bad):
    """All 4x4 Hadamard sign matrices, by enumeration of the 2^16 matrices."""
    had = []
    for m in range(1 << 16):
        rows = [VEC[(m >> (4 * r)) & 15] for r in range(4)]
        ok = True
        for a in range(4):
            for b in range(a + 1, 4):
                if sum(x * y for x, y in zip(rows[a], rows[b])):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            had.append(rows)
    if len(had) != 768:
        bad.append("enumerated %d Hadamard matrices of order 4, expected 768"
                   % len(had))
    rowsplit = [0, 0]
    colsplit = [0, 0]
    for H in had:
        rc = set(wclass(r) for r in H)
        cc = set(wclass(tuple(H[r][c] for r in range(4))) for c in range(4))
        if len(rc) != 1 or len(cc) != 1:
            bad.append("a 4x4 Hadamard matrix has mixed row/column classes")
            break
        rowsplit[rc.pop()] += 1
        colsplit[cc.pop()] += 1
    if rowsplit != [384, 384] or colsplit != [384, 384]:
        bad.append("class split is %s / %s, expected [384,384] / [384,384]"
                   % (rowsplit, colsplit))
    return had, rowsplit, colsplit


def classify(W):
    """The (D-e) mechanism, read off W = p Lambda(d)^T.

    Returns ("spike", None)  -- every row of W/2 is +-2 e_k: accepts every U;
            ("sign", c)      -- every row of W/2 is a sign vector of class c:
                                accepts exactly the U whose columns are the
                                other class;
            (None, None)     -- neither: the mechanism does not apply.
    """
    kinds = set()
    classes = set()
    for row in W:
        if any(v % 2 for v in row):
            return None, None
        v2 = [v // 2 for v in row]
        nz = [k for k in range(4) if v2[k]]
        if len(nz) == 1 and abs(v2[nz[0]]) == 2:
            kinds.add("spike")
        elif all(abs(v) == 1 for v in v2):
            kinds.add("sign")
            classes.add(wclass(v2))
        else:
            return None, None
    if len(kinds) != 1:
        return None, None
    kind = kinds.pop()
    if kind == "spike":
        return "spike", None
    if len(classes) != 1:
        return None, None
    return "sign", classes.pop()


def census(args, had, table, want_naive, bad):
    """The full 24 x 768^2 census.

    Returns (counts, stats).  With want_naive False only the table
    evaluator runs (used by the negative controls).
    """
    cols = [[tuple(H[r][c] for r in range(4)) for c in range(4)] for H in had]
    colmask = []
    for cs in cols:
        m = 0
        for u in cs:
            m |= 1 << IDX[u]
        colmask.append(m)
    colcls = [wclass(cs[0]) for cs in cols]
    counts = {}
    stats = {"naive_vs_table": 0, "mech_vs_table": 0, "mech_inapplicable": 0,
             "spike_p": 0, "sign_p": 0}
    for y in args:
        L = lam(y, table)
        Lt = transpose(L)
        total = 0
        for H in had:
            W = matmul(H, Lt)
            w0, w1, w2, w3 = W
            allowed = 0xFFFF
            for w in W:
                m = 0
                for jj, v in enumerate(VEC):
                    e = w[0] * v[0] + w[1] * v[1] + w[2] * v[2] + w[3] * v[3]
                    if e == 4 or e == -4:
                        m |= 1 << jj
                allowed &= m
            kind, cl = classify(W)
            if kind == "spike":
                stats["spike_p"] += 1
            elif kind == "sign":
                stats["sign_p"] += 1
            for iu in range(768):
                t_ok = (colmask[iu] & ~allowed) == 0
                if t_ok:
                    total += 1
                if want_naive:
                    n_ok = True
                    for u in cols[iu]:
                        u0, u1, u2, u3 = u
                        for w in (w0, w1, w2, w3):
                            e = (w[0] * u0 + w[1] * u1
                                 + w[2] * u2 + w[3] * u3)
                            if e != 4 and e != -4:
                                n_ok = False
                                break
                        if not n_ok:
                            break
                    if n_ok != t_ok:
                        stats["naive_vs_table"] += 1
                    if kind is None:
                        stats["mech_inapplicable"] += 1
                    else:
                        pred = True if kind == "spike" else (cl != colcls[iu])
                        if pred != t_ok:
                            stats["mech_vs_table"] += 1
        counts[y] = total
    if want_naive:
        if stats["naive_vs_table"]:
            bad.append("naive and table evaluators disagree on %d pairs"
                       % stats["naive_vs_table"])
        if stats["mech_inapplicable"]:
            bad.append("the (D-e) mechanism does not apply to %d pairs"
                       % stats["mech_inapplicable"])
        if stats["mech_vs_table"]:
            bad.append("the (D-e) mechanism mispredicts %d pairs"
                       % stats["mech_vs_table"])
    return counts, stats


def spot_check_corner(args, had, bad):
    """On a deterministic accepted sample, build E and check it is Hadamard.

    (D-d) says E E^T = 4 I_4 is automatic once E is +-1; this exhibits it.
    """
    seen = 0
    for y in args:
        L = lam(y)
        Lt = transpose(L)
        for ip in range(0, 768, 97):
            W = matmul(had[ip], Lt)
            for iu in range(0, 768, 89):
                U = had[iu]
                M = matmul(W, U)
                if any(v % 4 for row in M for v in row):
                    continue
                E = [[-v // 4 for v in row] for row in M]
                if any(abs(v) != 1 for row in E for v in row):
                    continue
                if matmul(E, transpose(E)) != [[4 if a == b else 0
                                                for b in range(4)]
                                               for a in range(4)]:
                    bad.append("an accepted corner at d=%s is not Hadamard"
                               % (y,))
                    return seen
                seen += 1
    if seen < 100:
        bad.append("corner spot-check found only %d accepted samples" % seen)
    return seen


# ----------------------------------------------------------------- main

def main(argv):
    table = LAM_TABLE
    installed = None
    for a in argv[1:]:
        if a == "--negative-control" or a == "--negative-control=sign":
            installed, table = NEG_SIGN[0], NEG_SIGN[1]
        elif a == "--negative-control=index":
            installed, table = NEG_INDEX[0], NEG_INDEX[1]
        else:
            die("unknown argument %r (expected --negative-control"
                "[=sign|index])" % a)

    if not os.path.isfile(os.path.join(TOOLS, "bordered_gs.py")):
        die("missing %s -- run this from the repository root: "
            "python certs/10-theorem-d-census/run.py"
            % os.path.join(TOOLS, "bordered_gs.py"))

    bad = []
    t0 = time.time()
    if installed:
        print("  NEGATIVE CONTROL INSTALLED: %s" % installed)
        print("  this run is expected to FAIL; a PASS here would mean the "
              "census does not see the Lambda table")

    grid = check_lambda_identity(bad)
    print("  Lambda(y) Lambda(y)^T = (sum y_q^2) I_4 on all %d points of "
          "{-2..2}^4 -- an identity of Z[y]" % grid)

    args, even, odd = admissible_arguments(bad)
    print("  arguments with sum_q y_q^2 = 4: %d = %d even (+-2 e_j) + %d odd "
          "((+-1,+-1,+-1,+-1))" % (len(args), len(even), len(odd)))

    seen = check_compressed_blocks(bad)
    print("  compressed-block identities over Ghat = Z_2 on all %d "
          "(eps, sigma) points of {-1,0,1}^8 x {+1,-1} -- identities of "
          "Z[sigma]" % seen)

    had, rowsplit, colsplit = hadamard4(bad)
    print("  4x4 Hadamard matrices: %d of the 2^16 sign matrices; row "
          "classes %s, column classes %s" % (len(had), rowsplit, colsplit))
    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 10: FAIL (%d problems)" % len(bad))
        return 1

    t1 = time.time()
    counts, stats = census(args, had, table, True, bad)
    pairs = len(had) * len(had)
    half, three_q = pairs // 2, 3 * pairs // 4
    print("  census over %d x %d^2 = %d triples (d, p, U), three evaluators "
          "(%.1fs)" % (len(args), len(had), len(args) * pairs,
                       time.time() - t1))
    for y in even:
        if counts[y] != half:
            bad.append("d=%s: %d admitting pairs, expected 768^2/2 = %d"
                       % (y, counts[y], half))
    for y in odd:
        if counts[y] != three_q:
            bad.append("d=%s: %d admitting pairs, expected 3*768^2/4 = %d"
                       % (y, counts[y], three_q))
    ev = sorted(set(counts[y] for y in even))
    od = sorted(set(counts[y] for y in odd))
    print("    even arguments (+-2 e_j), 8 of them: %s   768^2/2   = %d"
          % (ev, half))
    print("    odd  arguments, 16 of them:          %s   3*768^2/4 = %d"
          % (od, three_q))
    print("    naive vs table: %d disagreements; (D-e) mechanism vs table: "
          "%d mispredictions, %d inapplicable"
          % (stats["naive_vs_table"], stats["mech_vs_table"],
             stats["mech_inapplicable"]))
    print("    W = p Lambda(d)^T: %d (d, p) with W/2 all +-2 e_k spikes, "
          "%d with W/2 a sign matrix" % (stats["spike_p"], stats["sign_p"]))

    seen = spot_check_corner(args, had, bad)
    print("  corner spot-check: %d accepted (d, p, U) samples built E "
          "explicitly; every E is +-1 with E E^T = 4 I_4" % seen)

    if not installed:
        for name, tab, mover in (NEG_SIGN, NEG_INDEX):
            cc, _ = census(args, had, tab, False, bad)
            moved = [y for y in args if cc[y] != counts[y]]
            print("  negative control -- %s: census moves on %d of the 24 "
                  "arguments (d=%s: %d, was %d)"
                  % (name, len(moved), mover, cc[mover], counts[mover]))
            if not moved:
                bad.append("negative control %r left the census unchanged"
                           % name)
            elif mover not in moved:
                bad.append("negative control %r did not move d=%s"
                           % (name, (mover,)))

    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 10: FAIL (%d problems)" % len(bad))
        return 1
    print("CERT 10: PASS -- the 768^2 border census of NOTE-B.md Sec.1.5 "
          "(D-e) replays exactly, and the compressed-block identities of "
          "(D-a')/(D-d) hold identically (%.1fs)" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
