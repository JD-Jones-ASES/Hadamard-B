#!/usr/bin/env python3
"""bordered_gs.py -- this repository's shared checker for the bordered
Goethals-Seidel master theorem, used by certs 01, 02, 04, 06, 08 and 09.

Stdlib only.  Exact integer arithmetic only.  No numpy, no floats.

PROVENANCE.  This file is ported from the source laboratory's
deliberately INDEPENDENT second implementation of the theorem, and the
two-implementation pedigree is the point of it: there it did not import,
read, or copy that laboratory's own replay builder.  It builds the four
group-developed blocks as explicit matrices and assembles the displayed
4x4 block array by literal matrix operations (column reversal along the
reflection permutation, transposition, negation), the way the theorem
statement displays it.  The only thing shared with the replay builder is
the *data convention* by which a seed string is indexed by a group
element (mixed-radix, row-major), which is a property of the decoded
tape, not of an implementation.

LABEL MAPPING (read before quoting report keys)
===============================================
This module's H-numbering predates note/NOTE-B.md and is PERMUTED
relative to Theorem A there:

  this module          NOTE-B.md Theorem A
  H1 (PAF profile)  =  (H2)
  H2 (corner/row)   =  (H3)
  H3 (column Gram)  =  (H1)
  H4 (coupling)     =  (H4)

The four conditions are identical; only the labels differ. Report
keys carry THIS module's numbering.

WHAT IS TESTED
==============
For each parameter record (G, s, rho, K, seeds, E, P, Q) it checks the
theorem's hypotheses

  H0  order = 4*(|G| + s); seeds are +-1 of length |G|; tables are +-1
      of the declared shapes; K = ker(kappa) is a subgroup of index i.
  H1  (two-tier PAF profile)
        sum_q PAF_q(g) = 4n*[g=0] - 4s*[g in K\\{0}] + 4*[g not in K]
  H2  (corner/row-table)   E E^T + w P P^T = N I_{4s}
  H3  (column-table Gram)  Q Q^T = I_4 (x) ((4s+4) I_i - 4 J_i)
  H4  (coupling)           E Q^T + P Chat^T = 0,
      Chat = the Goethals-Seidel array over G/K built from the coset
      sums sigma_q and the reflection kappa(rho).

and the derived statements

  D3  w > 2s  =>  E is a 4s x 4s Hadamard matrix and P P^T = 4i I_{4s}
  D5  sum_q (sum_g x_q(g))^2 = 8n - 4w(s+1) + 4s
  D1  i <= s+1;   (w > 2s) => s <= i

then ASSEMBLES the matrix from the parameters exactly as the theorem
displays it and hands the result to `verify/verify.py`, the lab's sole
trust chain.  A cross-check compares the compressed core against the
G/K Goethals-Seidel array predicted by the compression lemma.

USAGE
=====
  python tools/bordered_gs.py --params data/payload-records.json
                              --verify verify/verify.py
                              --out <scratch dir>
                              [--only 668,716] [--no-verify]

Exit code 0 iff every record passed every hypothesis, the compression
lemma cross-check, and (unless --no-verify) verify.py returned HADAMARD.
"""

import argparse
import json
import os
import subprocess
import sys
import time


# ----------------------------------------------------------------- group

def prod(xs):
    p = 1
    for x in xs:
        p *= x
    return p


class AbelianGroup:
    """G = Z_{n_0} x ... x Z_{n_k}, elements indexed row-major mixed-radix.

    idx(t) = ((t_0 * n_1) + t_1) * n_2 + ...  -- the tape's convention.
    """

    def __init__(self, factors):
        self.factors = list(factors)
        for n in self.factors:
            if n < 1:
                raise ValueError("bad group factor %r" % (n,))
        self.n = prod(self.factors)
        self.elts = []
        t = [0] * len(self.factors)
        for _ in range(self.n):
            self.elts.append(tuple(t))
            for j in range(len(self.factors) - 1, -1, -1):
                t[j] += 1
                if t[j] < self.factors[j]:
                    break
                t[j] = 0

    def idx(self, t):
        k = 0
        for c, n in zip(t, self.factors):
            k = k * n + (c % n)
        return k

    def add(self, p, q):
        return self.idx(tuple(a + b for a, b in
                              zip(self.elts[p], self.elts[q])))

    def sub(self, p, q):
        return self.idx(tuple(a - b for a, b in
                              zip(self.elts[p], self.elts[q])))

    def neg(self, p):
        return self.idx(tuple(-a for a in self.elts[p]))

    def sub_table(self):
        """S[g][h] = idx(g_h - g_g) -- used to develop a matrix."""
        return [[self.sub(h, g) for h in range(self.n)]
                for g in range(self.n)]


def coset_map(G, divisors):
    """kappa: G -> Z_{d_0} x ... x Z_{d_k}, index mixed-radix over d.

    K = ker(kappa) has index i = prod(d).  Requires d_j | n_j so kappa is
    a homomorphism.
    """
    if len(divisors) != len(G.factors):
        raise ValueError("divisor arity mismatch")
    for d, n in zip(divisors, G.factors):
        if d < 1 or n % d:
            raise ValueError("d=%r does not divide n=%r" % (d, n))
    i = prod(divisors)
    kap = []
    for t in G.elts:
        c = 0
        for a, d in zip(t, divisors):
            c = c * d + (a % d)
        kap.append(c)
    return kap, i


class QuotientGroup:
    """Ghat = Z_{d_0} x ... x Z_{d_k}, same mixed-radix indexing as kappa."""

    def __init__(self, divisors):
        self.g = AbelianGroup(divisors)

    def add(self, p, q):
        return self.g.add(p, q)

    def sub(self, p, q):
        return self.g.sub(p, q)

    @property
    def n(self):
        return self.g.n


# ------------------------------------------------------- matrix helpers
# Explicit n x n integer matrices as lists of lists.  These are the
# literal operations that appear in the displayed block array.

def developed(seed, subtab):
    """X[g][h] = x(h - g)  (type-1 group development)."""
    return [[seed[k] for k in row] for row in subtab]


def transpose(X):
    return [list(col) for col in zip(*X)]


def reverse_columns(X, rperm):
    """(X R)[g][h] = X[g][ rperm[h] ], rperm[h] = index(rho - h)."""
    return [[row[rperm[h]] for h in range(len(rperm))] for row in X]


def negate(X):
    return [[-v for v in row] for row in X]


def mat_mul_t(A, B):
    """A * B^T for small matrices."""
    return [[sum(x * y for x, y in zip(ra, rb)) for rb in B] for ra in A]


def eye(k, val=1):
    return [[val if a == b else 0 for b in range(k)] for a in range(k)]


# --------------------------------------------------------- the theorem

def signs(s):
    if set(s) - {"+", "-"}:
        raise ValueError("string outside {+,-}: %r" % (s,))
    return [1 if ch == "+" else -1 for ch in s]


def check_record(rec, verbose=True):
    """Check hypotheses and assemble.  Returns (report, rows_as_strings)."""
    rep = {"order": int(rec["order"])}
    fail = []

    # The assembly below is the house ("standard") orientation ONLY.  A
    # record declaring any other gs_variant must be refused, not silently
    # assembled wrong (turn-42 skeptic item 6; D-049.1 scope clause).
    variant = rec.get("gs_variant", "standard")
    rep["gs_variant"] = variant
    if variant != "standard":
        fail.append("gs_variant %r unsupported (this checker assembles "
                    "the standard orientation only)" % (variant,))
        rep["hypotheses_ok"] = False
        rep["failures"] = fail
        return rep, []

    G = AbelianGroup(rec["group"])
    n = G.n
    s = int(rec["s"])
    N = int(rec["order"])
    rep.update({"n": n, "s": s, "group": list(rec["group"])})

    # ---- H0: shape ----
    if N != 4 * (n + s):
        fail.append("H0 order != 4(|G|+s)")
    seeds = [signs(x) for x in rec["seeds"]]
    if len(seeds) != 4 or any(len(x) != n for x in seeds):
        fail.append("H0 seed shape")

    rho = G.idx(tuple(rec["r_shift"]))

    if s > 0:
        div = list(rec["coset_divisors"])
    else:
        # s = 0: the theorem degenerates with K = G, i = 1.
        div = [1] * len(G.factors)
    kappa, i = coset_map(G, div)
    w = n // i
    rep.update({"i": i, "w": w, "K_index": i, "coset_divisors": div})

    # ---- H1: the two-tier PAF profile ----
    # sum_q PAF_q(g) computed directly from the definition.
    sub = G.sub_table()          # sub[g][h] = h - g
    prof_ok = True
    prof_bad = []
    for g in range(n):
        # PAF_q(g) = sum_h x_q(h) x_q(h+g)
        shift = [G.add(h, g) for h in range(n)]
        tot = 0
        for x in seeds:
            tot += sum(x[h] * x[shift[h]] for h in range(n))
        if g == 0:
            want = 4 * n
        elif kappa[g] == 0:
            want = -4 * s
        else:
            want = 4
        if tot != want:
            prof_ok = False
            if len(prof_bad) < 5:
                prof_bad.append([g, tot, want])
    rep["H1_two_tier_PAF"] = prof_ok
    if not prof_ok:
        fail.append("H1 PAF profile")
        rep["H1_first_failures"] = prof_bad

    # ---- D5: the row-sum law ----
    rowsums = [sum(x) for x in seeds]
    rep["seed_row_sums"] = rowsums
    d5_want = 8 * n - 4 * w * (s + 1) + 4 * s
    rep["D5_rowsum_law"] = (sum(r * r for r in rowsums) == d5_want)
    rep["D5_value"] = [sum(r * r for r in rowsums), d5_want]
    if not rep["D5_rowsum_law"]:
        fail.append("D5 row-sum law")

    # ---- border data ----
    if s > 0:
        E = [signs(r) for r in rec["corner"]]
        P = [signs(r) for r in rec["row_table"]]
        colT = [signs(r) for r in rec["col_table"]]
        if len(E) != 4 * s or any(len(r) != 4 * s for r in E):
            fail.append("H0 corner shape")
        if len(P) != 4 * s or any(len(r) != 4 * i for r in P):
            fail.append("H0 row_table shape")
        if len(colT) != 4 * s or any(len(r) != 4 * i for r in colT):
            fail.append("H0 col_table shape")
        # Q as a 4i x 4s matrix: Q[iI + c][r] = colT[r][iI + c]
        Q = [[colT[r][k] for r in range(4 * s)] for k in range(4 * i)]

        # ---- H2 ----
        EEt = mat_mul_t(E, E)
        PPt = mat_mul_t(P, P)
        h2 = all(EEt[a][b] + w * PPt[a][b] == (N if a == b else 0)
                 for a in range(4 * s) for b in range(4 * s))
        rep["H2_corner_rowtable"] = h2
        if not h2:
            fail.append("H2")

        # ---- H3 ----
        QQt = mat_mul_t(Q, Q)
        tgt = [[0] * (4 * i) for _ in range(4 * i)]
        for I in range(4):
            for c1 in range(i):
                for c2 in range(i):
                    tgt[I * i + c1][I * i + c2] = \
                        (4 * s + 4 if c1 == c2 else 0) - 4
        h3 = (QQt == tgt)
        rep["H3_coltable_Gram"] = h3
        if not h3:
            fail.append("H3")

        # ---- Chat by the compression lemma: GS array over G/K ----
        Gq = QuotientGroup(div)
        sigma = [[0] * i for _ in range(4)]
        for q in range(4):
            for g in range(n):
                sigma[q][kappa[g]] += seeds[q][g]
        rep["sigma"] = sigma
        rho_bar = kappa[rho]
        Chat = gs_array(Gq, sigma, rho_bar)   # 4i x 4i ints

        # ---- H4 ----
        h4 = all(sum(E[r][j] * Q[k][j] for j in range(4 * s)) +
                 sum(P[r][cc] * Chat[k][cc] for cc in range(4 * i)) == 0
                 for r in range(4 * s) for k in range(4 * i))
        rep["H4_coupling"] = h4
        if not h4:
            fail.append("H4")

        # ---- derived D3 / D1 ----
        rep["D3_applicable_w_gt_2s"] = (w > 2 * s)
        rep["E_is_Hadamard"] = (EEt == eye(4 * s, 4 * s))
        rep["PPt_eq_4i_I"] = (PPt == eye(4 * s, 4 * i))
        rep["D1_i_le_s_plus_1"] = (i <= s + 1)
        rep["s_le_i"] = (s <= i)
        # D4 (coset-sum Parseval bound, meaningful only for i >= 2)
        rep["D4_w_times_s1_minus_i_le_s"] = (
            (w * (s + 1 - i) <= s) if i >= 2 else "vacuous (i=1)")
        # Sigma-bar: Chat Chat^T = I_4 (x) (4(w(i-s-1)+s) I_i + 4w J_i)
        CCt = mat_mul_t(Chat, Chat)
        alpha = 4 * (w * (i - s - 1) + s)
        sbar = [[0] * (4 * i) for _ in range(4 * i)]
        for I in range(4):
            for c1 in range(i):
                for c2 in range(i):
                    sbar[I * i + c1][I * i + c2] = \
                        (alpha if c1 == c2 else 0) + 4 * w
        rep["SigmaBar_law"] = (CCt == sbar)
        if not rep["SigmaBar_law"]:
            fail.append("Sigma-bar law")
        if i >= 2 and not rep["D4_w_times_s1_minus_i_le_s"]:
            fail.append("D4")
        if w > 2 * s and not (rep["E_is_Hadamard"] and rep["PPt_eq_4i_I"]):
            fail.append("D3 (w>2s but corner/table not orthogonal)")
        if not rep["D1_i_le_s_plus_1"]:
            fail.append("D1")
        if w > 2 * s and not rep["s_le_i"]:
            fail.append("s<=i under w>2s")
        rep["i_eq_s_plus_1"] = (i == s + 1)
    else:
        E = P = Q = None
        Chat = None
        rep["H2_corner_rowtable"] = "vacuous (s=0)"
        rep["H3_coltable_Gram"] = "vacuous (s=0)"
        rep["H4_coupling"] = "vacuous (s=0)"

    # ---- assemble ----
    rows = assemble(G, seeds, sub, rho, s, i, kappa, E, P, Q)
    rep["assembled"] = "%dx%d" % (len(rows), len(rows[0]))

    # ---- cross-check: compression lemma against the assembled core ----
    if s > 0:
        ok, why = check_compression(rows, n, s, i, kappa, Chat)
        rep["compression_lemma_crosscheck"] = ok
        if not ok:
            fail.append("compression lemma: " + why)

    rep["hypotheses_ok"] = (not fail)
    rep["failures"] = fail
    return rep, rows


def gs_array(Gq, seq, rho_bar):
    """The standard Goethals-Seidel array over the (small) group Gq with
    the four integer sequences seq and reflection element rho_bar.
    Returned as a 4m x 4m integer matrix, m = |Gq|.  Built by the same
    explicit block operations used for the core."""
    m = Gq.n
    subt = [[Gq.sub(h, g) for h in range(m)] for g in range(m)]
    rperm = [Gq.sub(rho_bar, h) for h in range(m)]
    X = [developed(seq[q], subt) for q in range(4)]
    rev = [reverse_columns(X[q], rperm) for q in range(4)]
    revT = [reverse_columns(transpose(X[q]), rperm) for q in range(4)]
    B = [
        [X[0], rev[1], rev[2], rev[3]],
        [negate(rev[1]), X[0], revT[3], negate(revT[2])],
        [negate(rev[2]), negate(revT[3]), X[0], revT[1]],
        [negate(rev[3]), revT[2], negate(revT[1]), X[0]],
    ]
    out = [[0] * (4 * m) for _ in range(4 * m)]
    for I in range(4):
        for J in range(4):
            blk = B[I][J]
            for g in range(m):
                dst = out[I * m + g]
                src = blk[g]
                for h in range(m):
                    dst[J * m + h] = src[h]
    return out


def assemble(G, seeds, subtab, rho, s, i, kappa, E, P, Q):
    """Build H exactly as the theorem displays it; return +/- row strings."""
    n = G.n
    N = 4 * (n + s)
    rperm = [G.sub(rho, h) for h in range(n)]

    X = [developed(seeds[q], subtab) for q in range(4)]
    rev = [reverse_columns(X[q], rperm) for q in range(4)]
    revT = [reverse_columns(transpose(X[q]), rperm) for q in range(4)]
    blocks = [
        [X[0], rev[1], rev[2], rev[3]],
        [negate(rev[1]), X[0], revT[3], negate(revT[2])],
        [negate(rev[2]), negate(revT[3]), X[0], revT[1]],
        [negate(rev[3]), revT[2], negate(revT[1]), X[0]],
    ]
    del X, rev, revT

    ch = {1: "+", -1: "-"}
    rows = []
    if s > 0:
        for r in range(4 * s):
            parts = ["".join(ch[v] for v in E[r])]
            for J in range(4):
                parts.append("".join(ch[P[r][i * J + kappa[g]]]
                                     for g in range(n)))
            rows.append("".join(parts))
    for I in range(4):
        for g in range(n):
            parts = []
            if s > 0:
                k = i * I + kappa[g]
                parts.append("".join(ch[Q[k][c]] for c in range(4 * s)))
            for J in range(4):
                parts.append("".join(ch[v] for v in blocks[I][J][g]))
            rows.append("".join(parts))
    assert len(rows) == N and all(len(r) == N for r in rows), "bad shape"
    return rows


def check_compression(rows, n, s, i, kappa, Chat):
    """Every core row's K-coset column sums depend only on (I, kappa(g)),
    and equal the G/K Goethals-Seidel array Chat."""
    S = 4 * s
    colcls = [(j // n) * i + kappa[j % n] for j in range(4 * n)]
    seen = [None] * (4 * i)
    for rr in range(4 * n):
        line = rows[S + rr]
        vec = [0] * (4 * i)
        for j in range(4 * n):
            vec[colcls[j]] += 1 if line[S + j] == "+" else -1
        cl = (rr // n) * i + kappa[rr % n]
        if seen[cl] is None:
            seen[cl] = vec
        elif seen[cl] != vec:
            return False, "core compression not constant on class %d" % cl
    if seen != Chat:
        return False, "compressed core != G/K Goethals-Seidel array"
    return True, ""


# ------------------------------------------------------------------- CLI

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--params", required=True)
    ap.add_argument("--verify", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default="")
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--report", default="")
    ap.add_argument("--replay-report", default="",
                    help="banked replay report (pr0023); every verified "
                         "canonical_sha256 must match its banked value")
    args = ap.parse_args(argv)

    banked = {}
    if args.replay_report:
        with open(args.replay_report, "r", encoding="ascii") as fh:
            rr = json.load(fh)
        for e in rr["per_order"]:
            v = e["verify_verdict"]
            if "canonical_sha256=" not in v:
                raise ValueError("no banked sha for order %r" % e["order"])
            banked[int(e["order"])] = \
                v.rsplit("canonical_sha256=", 1)[1].strip()

    with open(args.params, "r", encoding="ascii") as fh:
        bank = json.load(fh)
    recs = bank["orders"] if "orders" in bank else [bank]
    if args.only:
        want = {int(x) for x in args.only.split(",")}
        recs = [r for r in recs if int(r["order"]) in want]

    os.makedirs(args.out, exist_ok=True)
    results = []
    all_ok = True
    for rec in recs:
        t0 = time.time()
        rep, rows = check_record(rec)
        if not rows:
            all_ok = False
            results.append(rep)
            print("order %-5d REFUSED: %s" % (rep["order"], rep["failures"]))
            continue
        path = os.path.join(args.out, "H%d_theorem.txt" % rep["order"])
        with open(path, "w", encoding="ascii", newline="\n") as fh:
            fh.write("\n".join(rows) + "\n")
        rep["artifact"] = path
        if not args.no_verify:
            out = subprocess.run(
                [sys.executable, args.verify, path],
                capture_output=True, text=True)
            verdict = (out.stdout or out.stderr).strip().splitlines()
            rep["verify"] = verdict[-1] if verdict else "(no output)"
            rep["verify_rc"] = out.returncode
            if out.returncode != 0:
                all_ok = False
            if banked:
                got = None
                if "canonical_sha256=" in rep["verify"]:
                    got = rep["verify"].rsplit(
                        "canonical_sha256=", 1)[1].strip()
                want = banked.get(rep["order"])
                rep["sha_banked"] = want
                rep["sha_match"] = (got is not None and got == want)
                if not rep["sha_match"]:
                    all_ok = False
        rep["seconds"] = round(time.time() - t0, 1)
        if not rep["hypotheses_ok"]:
            all_ok = False
        results.append(rep)
        print("order %-5d hyp_ok=%-5s %s  (%.1fs)"
              % (rep["order"], rep["hypotheses_ok"],
                 rep.get("verify", "(not verified)"), rep["seconds"]))
        sys.stdout.flush()

    summary = {"records": len(results), "all_ok": all_ok,
               "results": results}
    if args.report:
        with open(args.report, "w", encoding="ascii", newline="\n") as fh:
            json.dump(summary, fh, indent=1)
    print("ALL OK:", all_ok)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
