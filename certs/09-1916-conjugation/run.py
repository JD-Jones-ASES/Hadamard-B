#!/usr/bin/env python3
"""cert 09 -- order 1916: the character twist is a NON-HOUSE Gram and a
DIAGONAL CONJUGATION at the same time.

  THEOREM.  Let H be the decoded (s, i) = (3, 4) bordered Goethals-Seidel
  record at order 1916 (data/payload-records.json), with G = Z2 x Z2 x
  Z119, K of index 4, Gbar = G/K = Z2 x Z2, and rho = 0.  Let psi =
  psibar . kappa for a character psibar of Gbar, and let

      x'_q = psi . x_q ,   P' = P Dbar ,   Q' = Dbar Q ,   E' = E

  be the character twist (data/n1916-twist.json).  Then

  (a)  the twisted instance satisfies EVERY hypothesis of the master
       theorem in its general (Theorem A) form, with Gram
       M'(ebar) = psibar(ebar) M(ebar), and for psibar nontrivial that M'
       is NOT the house form (4s+4)I_i - 4J_i;  so "the house Gram is
       forced" is FALSE as a literal statement;  and yet

  (b)  H' = S H S exactly, with S = diag(I_{4s}, I_4 (x) diag(psi(g))) --
       952 rows and the same 952 columns negated, no permutation at all.
       So the twist manufactures NO new matrix.  This is the psi(rho) = 1
       case of note/NOTE-B.md S1.4's proposition, doing exactly what the
       proposition promises.

  Both halves are checked here in exact integer arithmetic, (b) over all
  1916^2 = 3 671 056 cells.  All three nontrivial characters of Gbar are
  swept, not just one:  all three give a non-house Gram, and all three are
  diagonal conjugations.

WHAT THIS SCRIPT DOES  (standard library only, seconds)

  (0) Pins the SHA-256 of the banked twist file.
  (1) Rebuilds the house record through tools/bordered_gs.py, which
      re-checks every hypothesis of the master theorem in its HOUSE form,
      hands it to verify/verify.py, and pins the canonical digest.
  (2) Derives the character group of Gbar from the record's own coset
      divisors, checks psi^2 = 1, psi|_K = 1 and psi(rho) = 1 for each,
      and compares against the banked psibar values.
  (3) For each of the four characters: applies the substitution,
      re-checks H1-H4 in the GENERAL form (arbitrary Gram M), measures M,
      cross-checks the compression lemma against the twisted G/K
      Goethals-Seidel array, assembles, verifies through
      verify/verify.py, and pins the canonical digest.
  (4) Asserts the Gram of each nontrivial twist is NOT house, that it is
      exactly psibar . M, and that its Fourier spectrum is the SAME
      multiset as the house one (so it sits inside NOTE-B S1.2's window
      and is not excluded by any eigenvalue argument).
  (5) Verifies H' = S H S cell by cell, all 1916^2 cells, for each twist.
  (6) Four controls -- see NOTES.md.

  tools/bordered_gs.py's check_record hard-codes the HOUSE two-tier PAF
  profile, so it REJECTS the twisted instance.  That rejection is a
  control here (C1), not a problem: it is the checker correctly saying
  "this is not a house instance".

Usage:
  python certs/09-1916-conjugation/run.py
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

N = 1916
OUT = os.path.join(HERE, "out")

SHA_HOUSE = \
    "be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9"
SHA_TWISTED = \
    "05d411faed301863d1e068651976f2e0f8e200b495af265ec20bc1bd6597175c"

FILE_PINS = {
    "data/n1916-twist.json":
        "1a3f92228074f69a7ead11d66371d18dfb39aeb2d17155f3b7fc9782b7b8d51b",
}

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
    h = hashlib.sha256()
    for r in rows:
        h.update((r + "\n").encode("ascii"))
    return h.hexdigest()


def verify_and_pin(tag, rows, want_sha=None):
    path = os.path.join(OUT, "H1916_%s.txt" % tag)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "verify", "verify.py"), path],
        capture_output=True, text=True)
    verdict = ((proc.stdout or proc.stderr).strip().splitlines()
               or ["(no output)"])[-1]
    dig = rows_sha256(rows)
    check("%-14s verify/verify.py: HADAMARD, all 1 834 570 row pairs"
          % tag, proc.returncode == 0 and dig in verdict, verdict[:74])
    if want_sha is not None:
        check("%-14s canonical sha256 == pin" % tag, dig == want_sha,
              dig[:24] + "...")
    os.remove(path)
    return dig


# ======================================================================
# the general (Theorem A) hypotheses -- arbitrary Gram M
# ======================================================================

def general_hypotheses(G, seeds, rho, s, kappa, i, Gq, E, P, Q):
    """H1-H4 with NO house assumption on M.  Exact integers.

    Labelled as tools/bordered_gs.py and note/NOTE-B.md label them, not as
    the upstream skeptic script does:

    H1  sum_q PAF_q(t) = 4n on t = 0 and -M(kappa(t)) elsewhere.  (In the
        house form M(0) = 4s and M(ebar) = -4, which is exactly the
        two-tier profile -4s on K\\{0} and +4 off K.)
    H2  E E^T + w P P^T = N I_{4s}.
    H3  Q Q^T = I_4 (x) M with M a function on Gbar -- the (a,b) entry of
        each diagonal superblock depends only on a - b, and the
        off-diagonal superblocks vanish.  This is where the Gram lives,
        and the ONLY place the house form is usually imposed.
    H4  E Q^T + P Chat^T = 0, Chat the G/K Goethals-Seidel array of the
        coset sums.
    """
    n = G.n
    w = n // i
    out = {}

    QQ = BGS.mat_mul_t(Q, Q)
    M, h3 = {}, True
    for I in range(4):
        for a in range(i):
            for b in range(i):
                e = Gq.sub(a, b)
                v = QQ[i * I + a][i * I + b]
                if e in M and M[e] != v:
                    h3 = False
                M.setdefault(e, v)
    for I in range(4):
        for J in range(4):
            if I != J and any(QQ[i * I + a][i * J + b]
                              for a in range(i) for b in range(i)):
                h3 = False
    out["H3"] = h3
    out["M"] = dict(M)

    h1 = True
    for g in range(n):
        shift = [G.add(h, g) for h in range(n)]
        tot = sum(sum(x[h] * x[shift[h]] for h in range(n)) for x in seeds)
        want = 4 * n if g == 0 else -M[kappa[g]]
        if tot != want:
            h1 = False
    out["H1"] = h1

    EE = BGS.mat_mul_t(E, E)
    PP = BGS.mat_mul_t(P, P)
    out["H2"] = all(EE[a][b] + w * PP[a][b] == (N if a == b else 0)
                    for a in range(4 * s) for b in range(4 * s))
    out["E_is_Hadamard"] = (EE == BGS.eye(4 * s, 4 * s))
    out["PPt_eq_4i_I"] = (PP == BGS.eye(4 * s, 4 * i))

    sigma = [[0] * i for _ in range(4)]
    for q in range(4):
        for g in range(n):
            sigma[q][kappa[g]] += seeds[q][g]
    Chat = BGS.gs_array(Gq, sigma, kappa[rho])
    EQ = BGS.mat_mul_t(E, Q)
    PC = BGS.mat_mul_t(P, Chat)
    out["H4"] = all(EQ[r][k] + PC[r][k] == 0
                    for r in range(4 * s) for k in range(4 * i))
    out["Chat"] = Chat
    out["sigma"] = sigma
    out["all"] = out["H1"] and out["H2"] and out["H3"] and out["H4"]
    return out


def fourier(M, Gq_group, i):
    """Eigenvalues of the Gbar-circulant M, over the real characters of
    Gbar = Z2 x Z2 (every character of an elementary abelian 2-group is
    real, so the spectrum is a multiset of integers)."""
    out = []
    for t in range(i):
        tt = Gq_group.elts[t]
        tot = 0
        for c in range(i):
            cc = Gq_group.elts[c]
            sgn = 1
            for a, b in zip(cc, tt):
                if a % 2 and b % 2:
                    sgn = -sgn
            tot += M[c] * sgn
        out.append(tot)
    return sorted(out)


# ======================================================================
# main
# ======================================================================

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.parse_args(argv)
    t_start = time.time()

    print("=" * 72)
    print("cert 09 -- order 1916: a non-house Gram that is only a diagonal")
    print("           conjugation")
    print("=" * 72)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)

    # ---------------------------------------------------------- clause 0
    print("\n[0] banked data file, SHA-256 pinned in this script")
    for name, want in sorted(FILE_PINS.items()):
        got = file_sha256(os.path.join(ROOT, name))
        check("%-30s" % name, got == want, got[:24] + "...")
    with open(os.path.join(ROOT, "data", "n1916-twist.json"),
              "r", encoding="ascii") as fh:
        spec = json.load(fh)

    # ---------------------------------------------------------- clause 1
    print("\n[1] the house record: rebuild, re-check the theorem, verify, pin")
    with open(os.path.join(ROOT, "data", "payload-records.json"),
              "r", encoding="ascii") as fh:
        rec = [r for r in json.load(fh)["orders"]
               if int(r["order"]) == N][0]
    t0 = time.time()
    rep, rows_h = BGS.check_record(rec)
    check("house hypotheses H0-H4, D1/D3/D5, Sigma-bar, compression lemma",
          rep["hypotheses_ok"] and rep["compression_lemma_crosscheck"],
          "s=%d i=%d w=%d  G=%s" % (rep["s"], rep["i"], rep["w"],
                                    rep["group"]))
    check("the record IS a house instance: its Gram is (4s+4)I_i - 4J_i",
          rep["H3_coltable_Gram"] and rep["H1_two_tier_PAF"])
    dig_h = verify_and_pin("house", rows_h, SHA_HOUSE)
    check("the banked twist file names this base matrix",
          spec["base"]["canonical_sha256"] == dig_h)
    print("       rebuilt and verified in %.1fs" % (time.time() - t0))

    # ---------------------------------------------------------- clause 2
    print("\n[2] the character group of Gbar, derived from the record")
    G = BGS.AbelianGroup(rec["group"])
    n = G.n
    s = int(rec["s"])
    div = list(rec["coset_divisors"])
    kappa, i = BGS.coset_map(G, div)
    w = n // i
    rho = G.idx(tuple(rec["r_shift"]))
    Gq = BGS.QuotientGroup(div)
    Gqg = BGS.AbelianGroup(div)
    check("the frame is (s, i) = (%d, %d), w = %d, |G| = %d, rho index %d"
          % (s, i, w, n, rho), (s, i, w, n) == (3, 4, 119, 476),
          "Gbar = G/K = %s" % (div,))
    check("Gbar is elementary abelian of order %d, so every character is "
          "real and psi^2 = 1 automatically" % i,
          all(d in (1, 2) for d in div) and Gqg.n == i,
          "Gbar elements %s" % (Gqg.elts,))
    check("kappa(rho) == 0, hence psi(rho) = 1 for EVERY character",
          kappa[rho] == 0 == spec["character"]["index_of_rho"])

    chars = []
    for t in range(i):
        tt = Gqg.elts[t]
        vals = []
        for c in range(i):
            cc = Gqg.elts[c]
            sgn = 1
            for a, b in zip(cc, tt):
                if a % 2 and b % 2:
                    sgn = -sgn
            vals.append(sgn)
        chars.append(vals)
    check("the four characters of Gbar, in the record's own coset indexing",
          chars == spec["character"]["all_psibar_values"], str(chars))
    check("the banked psibar (the first-coordinate character) is chars[%d]"
          % spec["character"]["banked_character_index"],
          chars[spec["character"]["banked_character_index"]]
          == spec["character"]["psibar_values"],
          str(spec["character"]["psibar_values"]))
    for t, vals in enumerate(chars):
        check("chi_%d: psibar^2 = 1, psibar|_K = 1 (it is a character of "
              "Gbar), psibar(kappa(rho)) = +1" % t,
              all(v in (1, -1) for v in vals) and vals[kappa[rho]] == 1,
              str(vals))

    # ---------------------------------------------------------- clause 3
    print("\n[3] each twist: hypotheses in the GENERAL form, Gram, assembly")
    seeds = [BGS.signs(x) for x in rec["seeds"]]
    E = [BGS.signs(r) for r in rec["corner"]]
    P = [BGS.signs(r) for r in rec["row_table"]]
    colT = [BGS.signs(r) for r in rec["col_table"]]
    Q = [[colT[r][k] for r in range(4 * s)] for k in range(4 * i)]
    sub = G.sub_table()

    house_M = {e: (4 * s + 4 if e == 0 else 0) - 4 for e in range(i)}
    base = general_hypotheses(G, seeds, rho, s, kappa, i, Gq, E, P, Q)
    check("the house instance satisfies the GENERAL H1-H4 too, with the "
          "house Gram", base["all"] and base["M"] == house_M,
          "M = %s" % dict(sorted(base["M"].items())))
    house_spec = fourier(base["M"], Gqg, i)
    check("the house Gram's spectrum is {0, 4s+4 (x %d)}, rank %d = s, "
          "trace %d = 4si" % (i - 1, s, 4 * s * i),
          house_spec == sorted([4 * s + 4 - 4 * i] + [4 * s + 4] * (i - 1))
          and sum(1 for v in house_spec if v) == s
          and sum(house_spec) == 4 * s * i, str(house_spec))

    results = {}
    for t, psibar in enumerate(chars):
        tag = "chi%d" % t
        xs2 = [[psibar[kappa[g]] * x[g] for g in range(n)] for x in seeds]
        Q2 = [[psibar[k % i] * v for v in Q[k]] for k in range(4 * i)]
        P2 = [[psibar[c % i] * P[r][c] for c in range(4 * i)]
              for r in range(4 * s)]
        h = general_hypotheses(G, xs2, rho, s, kappa, i, Gq, E, P2, Q2)
        rows_t = BGS.assemble(G, xs2, sub, rho, s, i, kappa, E, P2, Q2)
        ok_comp, why = BGS.check_compression(rows_t, n, s, i, kappa,
                                             h["Chat"])
        nonhouse = (h["M"] != house_M)
        want_M = {e: psibar[e] * base["M"][e] for e in range(i)}
        spec_t = fourier(h["M"], Gqg, i)
        print("\n  %s  psibar = %s" % (tag, psibar))
        check("%s  satisfies the GENERAL H1-H4 (arbitrary Gram)" % tag,
              h["all"], "H1=%s H2=%s H3=%s H4=%s"
              % (h["H1"], h["H2"], h["H3"], h["H4"]))
        check("%s  compression lemma: the assembled core compresses to the "
              "twisted G/K Goethals-Seidel array" % tag, ok_comp, why)
        check("%s  M = %s  -- %s" % (tag, dict(sorted(h["M"].items())),
                                     "HOUSE" if not nonhouse
                                     else "NOT the house form"),
              (nonhouse == (t != 0)),
              "house form is %s" % dict(sorted(house_M.items())))
        check("%s  M'(ebar) == psibar(ebar) * M(ebar), exactly" % tag,
              h["M"] == want_M)
        check("%s  the spectrum is the SAME multiset as the house Gram's "
              "-- no eigenvalue argument can exclude it" % tag,
              spec_t == house_spec, str(spec_t))
        check("%s  E is unchanged and still Hadamard; P' P'^T = 4i I" % tag,
              h["E_is_Hadamard"] and h["PPt_eq_4i_I"])
        dig = verify_and_pin(tag, rows_t,
                             SHA_TWISTED if t == spec["character"]
                             ["banked_character_index"] else None)
        if t == 0:
            check("%s  the trivial character returns the house matrix "
                  "unchanged" % tag, dig == SHA_HOUSE, dig[:24] + "...")
        check("%s  canonical sha256 == the banked value" % tag,
              dig == spec["instances"][tag]["canonical_sha256"],
              dig[:24] + "...")
        results[tag] = (psibar, rows_t, h, dig)

    # ---------------------------------------------------------- clause 4
    print("\n[4] the conjugation identity, cell by cell, all %d cells each"
          % (N * N))
    for t, psibar in enumerate(chars):
        tag = "chi%d" % t
        _pb, rows_t, _h, _d = results[tag]
        sgn = [1] * (4 * s) + [psibar[kappa[g]]
                               for _I in range(4) for g in range(n)]
        check("%s  S = diag(I_{4s}, I_4 (x) diag(psi(g))) has length N = %d "
              "and entries in {+-1}" % (tag, N),
              len(sgn) == N and set(sgn) <= {1, -1},
              "%d rows and the same %d columns negated"
              % (sum(1 for x in sgn if x < 0), sum(1 for x in sgn if x < 0)))
        bad = []
        ncells = 0
        ndiff = 0
        for r in range(N):
            src, tgt, er = rows_h[r], rows_t[r], sgn[r]
            for c in range(N):
                v = 1 if src[c] == "+" else -1
                want = "+" if er * sgn[c] * v == 1 else "-"
                if want != tgt[c]:
                    if len(bad) < 5:
                        bad.append((r, c))
                ncells += 1
                if src[c] != tgt[c]:
                    ndiff += 1
        neg = sum(1 for x in sgn if x < 0)
        check("%s  H' == S H S at every one of the %d cells" % (tag, N * N),
              not bad and ncells == N * N,
              "%d cells checked" % ncells if not bad else str(bad))
        check("%s  the map is diagonal but NOT trivial: H' differs from H "
              "in exactly 2*%d*%d = %d cells"
              % (tag, neg, N - neg, 2 * neg * (N - neg)),
              ndiff == 2 * neg * (N - neg) and (ndiff == 0) == (t == 0),
              "%d differing cells" % ndiff)

    # ---------------------------------------------------------- clause 5
    print("\n[5] controls")

    print("\n  C1 -- the house checker must REJECT the twisted instance")
    tag = "chi%d" % spec["character"]["banked_character_index"]
    psibar = results[tag][0]
    xs2 = [[psibar[kappa[g]] * x[g] for g in range(n)] for x in seeds]
    Q2 = [[psibar[k % i] * v for v in Q[k]] for k in range(4 * i)]
    P2 = [[psibar[c % i] * P[r][c] for c in range(4 * i)]
          for r in range(4 * s)]
    FL = {1: "+", -1: "-"}
    rec2 = dict(rec)
    rec2["seeds"] = ["".join(FL[v] for v in x) for x in xs2]
    rec2["row_table"] = ["".join(FL[v] for v in r) for r in P2]
    rec2["col_table"] = ["".join(FL[Q2[k][r]] for k in range(4 * i))
                         for r in range(4 * s)]
    rep2, _rows2 = BGS.check_record(rec2)
    check("C1  tools/bordered_gs.py rejects it, and for the right reasons",
          (not rep2["hypotheses_ok"])
          and not rep2["H1_two_tier_PAF"] and not rep2["H3_coltable_Gram"],
          "failures: %s" % rep2["failures"])
    check("C1  the rebuilt rows are nevertheless the twisted matrix "
          "(the assembler is orientation-agnostic)",
          rows_sha256(_rows2) == SHA_TWISTED)

    print("\n  C2 -- negative control: the cell-by-cell checker must FIRE")
    _pb, rows_t, _h, _d = results[tag]
    sgn = [1] * (4 * s) + [psibar[kappa[g]]
                           for _I in range(4) for g in range(n)]
    victim = [list(r) for r in rows_t]
    victim[7][11] = "-" if victim[7][11] == "+" else "+"
    victim = ["".join(r) for r in victim]
    caught = []
    for r in range(N):
        src, tgt, er = rows_h[r], victim[r], sgn[r]
        for c in range(N):
            v = 1 if src[c] == "+" else -1
            if ("+" if er * sgn[c] * v == 1 else "-") != tgt[c]:
                caught.append((r, c))
    check("C2  a single flipped cell is caught, and ONLY that cell",
          caught == [(7, 11)], str(caught[:4]))

    print("\n  C3 -- the conjugation is an involution")
    twice = ["".join(("+" if (sgn[r] * sgn[c] *
                              (1 if rows_t[r][c] == "+" else -1)) == 1
                      else "-") for c in range(N)) for r in range(N)]
    check("C3  S (S H S) S == H", rows_sha256(twice) == SHA_HOUSE)

    print("\n  C4 -- the twist is not a relabelling in disguise")
    neg = sum(1 for x in sgn if x < 0)
    check("C4  exactly %d of the %d rows carry psi(g) = -1, and the same "
          "%d columns" % (neg, N, neg), neg == 4 * (n // 2) == 952)
    check("C4  the row permutation and the column permutation are BOTH the "
          "identity -- the whole equivalence is diagonal", True,
          "S H S, no P_r, no P_c")

    # ---------------------------------------------------------- close out
    shutil.rmtree(OUT, ignore_errors=True)
    print("\n" + "=" * 72)
    if FAIL:
        print("cert 09: FAIL (%d): %s" % (len(FAIL), FAIL))
        return 1
    print("VERDICT: at order 1916 the character twist produces admissible")
    print("         instances whose Gram M is NOT the house form")
    print("         (4s+4)I_i - 4J_i -- so Gram-forcing is FALSE as a")
    print("         literal statement -- and every one of them is the")
    print("         diagonal conjugation S H S of the house record, so the")
    print("         twist manufactures NO new matrix.  All three")
    print("         nontrivial characters of Gbar = Z2 x Z2 behave the same")
    print("         way, each checked over all %d cells." % (N * N))
    print("         This is the psi(rho) = 1 case of NOTE-B S1.4.")
    print("         NOT claimed here: anything about psi(rho) = -1 in")
    print("         general (at 668 that case DOES leave the class --")
    print("         cert 06/08), and nothing about s >= 2 forcedness,")
    print("         which stays a CONJECTURE up to the twist.")
    print("=" * 72)
    print("generated matrices deleted; nothing left in %s   (%.1fs)"
          % (rel(OUT), time.time() - t_start))
    return 0


if __name__ == "__main__":
    sys.exit(main())
