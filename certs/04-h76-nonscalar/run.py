#!/usr/bin/env python3
"""cert 04 -- H(76) on the non-cyclic group Z2 x Z3 x Z3, seeds invariant
under a NON-SCALAR multiplier subgroup of Aut(G).

Run from the repository root:

    python certs/04-h76-nonscalar/run.py

Two things are checked, and they are independent of each other.

(1) THEOREM A, IN ITS GENERAL FORM.  This cert does NOT route the record
    through the house-form checker's H3 (`Q Q^T = I_4 (x) ((4s+4)I_i -
    4J_i)`).  It checks the four hypotheses of Theorem A as stated in
    NOTE-B.md 1.1, with the Gram `M` read off the record rather than
    assumed:

      (H1)  Q Q^T = I_4 (x) M for a symmetric, Gbar-invariant integer
            i x i matrix M (M[c,c'] = M(c-c')); all cross-superblock
            blocks of Q Q^T vanish.  M is MEASURED, not assumed.
      (H2)  sum_q PAF_q(t) = 4n*delta_{t,0} - M(kappa(t)) for all t != 0,
            with the M measured in (H1).
      (H3)  E E^T + w * P P^T = N * I_{4s}.
      (H4)  E Q^T + P Chat^T = 0, Chat = GS(sigma_0..sigma_3; kappa(rho))
            over the quotient Gbar.

    (NOTE-B.md numbers the hypotheses this way.  `tools/bordered_gs.py`
    numbers the same four conditions H1 = profile, H2 = corner/row,
    H3 = column Gram, H4 = coupling.  This cert uses the NOTE-B.md
    numbering throughout.)

    The building blocks (AbelianGroup, coset_map, QuotientGroup,
    gs_array, mat_mul_t, signs) are imported from `tools/bordered_gs.py`
    unmodified; the four hypotheses themselves are implemented here.

    As a control the cert ALSO runs the house-form checker
    `check_record` and asserts the two agree.

(2) THE MULTIPLIER.  The claim the artifact was built to gate is that
    `G` is non-cyclic and the seed quadruple is constant on the orbits
    of a NON-SCALAR multiplier subgroup `M <= Aut(G)`.  This cert
    enumerates Aut(G) from scratch (all matrices A over the invariant
    factors that give a well-defined bijective endomorphism), computes
    the exact stabiliser of the seed quadruple inside it, and asserts
    that the stabiliser is the banked group: order 2, non-scalar, 12
    orbits on G.

    "Non-scalar" here means: not of the form t -> c*t for any integer c.
    It is a property of an AUTOMORPHISM, and has nothing to do with the
    Gram matrix M of (H1) -- see NOTES.md.

Then the matrix is assembled, handed to `verify/verify.py`, digest-
compared against PINNED_SHA, and deleted.

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

DATA = os.path.join("data", "h76-nonscalar.json")
VERIFY = os.path.join("verify", "verify.py")
TOOLS = "tools"

PINNED_SHA = \
    "05ad9852f64ad2c48634105e26cf35ccb9dd2ea708caa1ce6b8ab6ca1bdb5e70"


def die(msg):
    print("FAIL: %s" % msg)
    sys.exit(1)


# ---------------------------------------------------------------- Aut(G)

def apply_aut(factors, A, t):
    return tuple(sum(A[a][b] * t[b] for b in range(len(factors))) % factors[a]
                 for a in range(len(factors)))


def all_automorphisms(G):
    """Every automorphism of G = Z_{f0} x ... x Z_{fk}, as a matrix A with
    alpha(t)_a = sum_b A[a][b] t_b mod f_a.  Brute force; G is small."""
    f = G.factors
    k = len(f)
    out = []
    entries = []

    def gen(a, b, cur):
        if a == k:
            entries.append([row[:] for row in cur])
            return
        if b == k:
            gen(a + 1, 0, cur)
            return
        for v in range(f[a]):
            # well-definedness: alpha must kill f_b * e_b, i.e.
            # A[a][b] * f_b == 0 (mod f_a)
            if (v * f[b]) % f[a]:
                continue
            cur[a][b] = v
            gen(a, b + 1, cur)
        cur[a][b] = 0

    gen(0, 0, [[0] * k for _ in range(k)])
    for A in entries:
        img = set()
        for t in G.elts:
            img.add(apply_aut(f, A, t))
        if len(img) == G.n:            # bijective => automorphism
            out.append(A)
    return out


def is_scalar(G, A):
    """alpha(t) = c*t for some integer c?"""
    f = G.factors
    for c in range(max(f)):
        if all(apply_aut(f, A, t) == tuple((c * x) % m for x, m in zip(t, f))
               for t in G.elts):
            return True
    return False


def aut_order(G, A):
    f = G.factors
    k = len(f)
    ident = [[1 if a == b else 0 for b in range(k)] for a in range(k)]
    X = [row[:] for row in A]
    o = 1
    while X != ident:
        X = [[sum(X[a][m] * A[m][b] for m in range(k)) % f[a]
              for b in range(k)] for a in range(k)]
        o += 1
        if o > 4 * G.n:
            return None
    return o


def orbit_count(G, gens):
    f = G.factors
    seen = set()
    cnt = 0
    for t in G.elts:
        if t in seen:
            continue
        cnt += 1
        stack = [t]
        seen.add(t)
        while stack:
            u = stack.pop()
            for A in gens:
                v = apply_aut(f, A, u)
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
    return cnt


# ----------------------------------------------- Theorem A, general form

def theorem_a(T, rec):
    """Check (H1)-(H4) of NOTE-B.md 1.1 with M measured, not assumed.
    Returns (ok, details, rows)."""
    G = T.AbelianGroup(rec["group"])
    n = G.n
    s = int(rec["s"])
    N = int(rec["order"])
    det = {}
    fail = []
    if N != 4 * (n + s):
        fail.append("H0 order != 4(|G|+s)")
    seeds = [T.signs(x) for x in rec["seeds"]]
    kappa, i = T.coset_map(G, list(rec["coset_divisors"]))
    w = n // i
    rho = G.idx(tuple(rec["r_shift"]))
    det.update({"n": n, "s": s, "i": i, "w": w, "N": N})

    E = [T.signs(r) for r in rec["corner"]]
    P = [T.signs(r) for r in rec["row_table"]]
    colT = [T.signs(r) for r in rec["col_table"]]
    Q = [[colT[r][k] for r in range(4 * s)] for k in range(4 * i)]

    Gq = T.QuotientGroup(list(rec["coset_divisors"]))

    # ---- (H1): Q Q^T = I_4 (x) M, M symmetric and Gbar-invariant ----
    QQ = T.mat_mul_t(Q, Q)
    h1 = True
    for I in range(4):
        for J in range(4):
            if I == J:
                continue
            if any(QQ[I * i + a][J * i + b]
                   for a in range(i) for b in range(i)):
                h1 = False
    Mfun = {}
    for I in range(4):
        for a in range(i):
            for b in range(i):
                e = Gq.sub(a, b)
                v = QQ[I * i + a][I * i + b]
                if e in Mfun and Mfun[e] != v:
                    h1 = False
                Mfun.setdefault(e, v)
    for e in range(i):                       # symmetry: M(e) = M(-e)
        if Mfun.get(e) != Mfun.get(Gq.sub(0, e)):
            h1 = False
    M = [[Mfun[Gq.sub(a, b)] for b in range(i)] for a in range(i)]
    det["M"] = M
    det["H1"] = h1
    if not h1:
        fail.append("H1 (Q Q^T = I_4 (x) M, M Gbar-invariant symmetric)")

    # ---- (H2): sum_q PAF(t) = 4n delta - M(kappa(t)) ----
    h2 = True
    bad2 = []
    for t in range(n):
        shift = [G.add(h, t) for h in range(n)]
        tot = sum(sum(x[h] * x[shift[h]] for h in range(n)) for x in seeds)
        want = 4 * n if t == 0 else -Mfun[kappa[t]]
        if tot != want:
            h2 = False
            if len(bad2) < 5:
                bad2.append([t, tot, want])
    det["H2"] = h2
    det["H2_first_failures"] = bad2
    if not h2:
        fail.append("H2 (profile against the measured M)")

    # ---- (H3): E E^T + w P P^T = N I ----
    EE = T.mat_mul_t(E, E)
    PP = T.mat_mul_t(P, P)
    h3 = all(EE[a][b] + w * PP[a][b] == (N if a == b else 0)
             for a in range(4 * s) for b in range(4 * s))
    det["H3"] = h3
    det["EEt"] = EE
    det["PPt"] = PP
    if not h3:
        fail.append("H3 (E E^T + w P P^T = N I)")

    # ---- (H4): E Q^T + P Chat^T = 0 ----
    sigma = [[0] * i for _ in range(4)]
    for q in range(4):
        for g in range(n):
            sigma[q][kappa[g]] += seeds[q][g]
    Chat = T.gs_array(Gq, sigma, kappa[rho])
    h4 = all(sum(E[r][j] * Q[k][j] for j in range(4 * s)) +
             sum(P[r][c] * Chat[k][c] for c in range(4 * i)) == 0
             for r in range(4 * s) for k in range(4 * i))
    det["H4"] = h4
    det["sigma"] = sigma
    if not h4:
        fail.append("H4 (coupling)")

    rows = T.assemble(G, seeds, G.sub_table(), rho, s, i, kappa, E, P, Q)
    det["failures"] = fail
    return (not fail), det, rows


def main():
    for p in (DATA, VERIFY, os.path.join(TOOLS, "bordered_gs.py")):
        if not os.path.isfile(p):
            die("missing %s -- run this from the repository root: "
                "python certs/04-h76-nonscalar/run.py" % p)
    sys.path.insert(0, TOOLS)
    import bordered_gs as T

    with open(DATA, "r", encoding="ascii") as fh:
        bank = json.load(fh)
    rec = bank["record"]
    bad = []
    t0 = time.time()

    # ---- (1) Theorem A, general form ----
    ok, det, rows = theorem_a(T, rec)
    if not ok:
        bad.extend(det["failures"])
    if len(rec["group"]) < 2 or min(rec["group"]) < 2:
        bad.append("G is not presented as a non-cyclic direct product")
    G = T.AbelianGroup(rec["group"])
    # non-cyclic <=> gcd of some pair of factors > 1
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a
    fs = list(rec["group"])
    noncyclic = any(gcd(fs[a], fs[b]) > 1
                    for a in range(len(fs)) for b in range(a + 1, len(fs)))
    if not noncyclic:
        bad.append("G = %s is cyclic; the gate requires a non-cyclic G" % fs)
    if det["M"] != bank.get("gram_M_measured"):
        bad.append("measured Gram M = %s disagrees with the banked "
                   "gram_M_measured = %s" % (det["M"], bank.get(
                       "gram_M_measured")))

    # control: the house-form checker must agree
    rep, rows2 = T.check_record(rec)
    if not rep["hypotheses_ok"]:
        bad.append("house-form control check_record failed: %s"
                   % (rep["failures"],))
    if rows2 != rows:
        bad.append("the two assemblies disagree")

    # ---- (2) the multiplier ----
    mult = bank["multiplier"]
    gens = mult["generators"]
    auts = all_automorphisms(G)
    seeds = [T.signs(x) for x in rec["seeds"]]
    idx = {t: k for k, t in enumerate(G.elts)}
    stab = [A for A in auts
            if all(seeds[q][idx[apply_aut(G.factors, A, t)]] ==
                   seeds[q][idx[t]] for q in range(4) for t in G.elts)]
    for A in gens:
        if A not in auts:
            bad.append("banked generator %s is not an automorphism of G" % A)
        if A not in stab:
            bad.append("the seeds are NOT invariant under the banked "
                       "generator %s" % A)
        if is_scalar(G, A):
            bad.append("the banked generator %s IS scalar" % A)
        o = aut_order(G, A)
        if o != mult["order"]:
            bad.append("banked generator has order %s, bank says %s"
                       % (o, mult["order"]))
    orb = orbit_count(G, gens)
    if orb != mult["orbits_on_G"]:
        bad.append("multiplier has %d orbits on G, bank says %d"
                   % (orb, mult["orbits_on_G"]))
    if len(stab) != mult["order"]:
        bad.append("the full stabiliser of the seed quadruple in Aut(G) has "
                   "order %d, bank claims the multiplier has order %d"
                   % (len(stab), mult["order"]))
    if not any(not is_scalar(G, A) for A in stab):
        bad.append("every element of the stabiliser is scalar -- the "
                   "non-scalar claim fails")

    # ---- trust chain + pinned digest ----
    tmp = tempfile.mkdtemp(prefix="cert04-")
    try:
        path = os.path.join(tmp, "H76_gate_2x3x3.txt")
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
        if got != bank.get("pinned_sha256"):
            bad.append("banked pinned_sha256 %s disagrees with the observed "
                       "digest %s" % (bank.get("pinned_sha256"), got))

    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 04: FAIL (%d problems)" % len(bad))
        return 1
    print("  H(76)  G = Z%s (non-cyclic)  n=%d s=%d i=%d w=%d"
          % (" x Z".join(str(x) for x in rec["group"]),
             det["n"], det["s"], det["i"], det["w"]))
    print("  Theorem A (general form, NOTE-B 1.1 numbering):")
    print("    (H1) Q Q^T = I_4 (x) M   OK   M measured = %s" % (det["M"],))
    print("         house form (4s+4)I_i - 4J_i at these parameters = [[4]] "
          "-- the measured M IS the house form")
    print("    (H2) sum_q PAF(t) = 4n.delta - M(kappa(t))   OK")
    print("    (H3) E E^T + w P P^T = N I_4                 OK")
    print("    (H4) E Q^T + P Chat^T = 0                    OK")
    print("  multiplier: |Aut(G)| = %d; stabiliser of the seed quadruple has "
          "order %d, generated by %s" % (len(auts), len(stab), gens))
    print("              non-scalar = True, orbits on G = %d" % orb)
    print("  %s" % line)
    print("CERT 04: PASS (%.1fs)" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
