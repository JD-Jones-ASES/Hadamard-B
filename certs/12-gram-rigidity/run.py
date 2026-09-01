#!/usr/bin/env python3
"""cert 12 -- Gram rigidity: the complete (3,4) classification, and the
general-branch witnesses at orders 1676 and 1772.

Run from the repository root:

    python certs/12-gram-rigidity/run.py

note/NOTE-B.md Sec.1.2.1 (Theorem E) says that under (H1)+(H2) with
s >= 1 and w > 2s, every admissible Gram is M = 4i P_S for a
conjugation- and Galois-stable S of size s, so the spectrum is
{0^(i-s), (4i)^s}; its Corollary E1 says that at i = s+1 the admissible
Grams are exactly the real-quotient-character twists of the house form
M_house = (4s+4) I_i - 4 J_i.  This cert measures both edges of that
statement.

PART A -- the complete classification at (s, i) = (3, 4).
  Every candidate Gram allowed by the conditions provable WITHOUT
  Theorem E is enumerated from first principles, over both abelian
  groups of order 4:

      A1  M(0) = 4s = 12                (diagonal of a Gram of +-1 rows
                                         of length 4s)
      A2  M(-c) = M(c)                  (M symmetric)
      A3  M(c) in Z, |M(c)| <= 12, M(c) even
                                        (<u,v> = 12 - 2d for +-1 rows)
      A4  M positive semidefinite       (it is a Gram matrix)
      A5  rank M <= s = 3               (four mutually orthogonal
                                         superblock row spans in R^{4s})
      A6  every eigenvalue < 4i + 2 = 18
                                        (the Parseval window of Theorem
                                         E's step (L1): lambda <= N/w =
                                         4i + 4s/w, and w > 2s)

  The survivors must be exactly 2 for Z_4 and 4 for Z_2 x Z_2, must all
  have spectrum {16, 16, 16, 0}, and must be exactly the real-character
  twist orbit of the house form.  The six banked witnesses of
  data/gram34-witnesses.json are then re-checked: each is a +-1 table
  Q in {+-1}^(16x12) with Q Q^T = I_4 (x) M, so every admissible Gram is
  also (H1)-realizable.  Finally the four Grams banked in
  data/n1916-twist.json -- cert 09's character-group sweep at order 1916
  -- are compared against the Z_2 x Z_2 list: they are the whole of it,
  so that sweep is an exhaustive classification at that cell.

  Eigenvalues are computed from the character formulas and then VERIFIED
  per matrix against the integer power sums tr(M^k), k = 1..4, through
  Newton's identities: a monic integer polynomial determines its root
  multiset, so nothing here rests on a transcribed formula.

  The census also reports what each hypothesis does.  Dropping A6 leaves
  48 / 1154 integer candidates (24 / 290 even); dropping A3's parity but
  keeping A6 changes nothing.  The Parseval window is doing the cutting.

PART B -- the general branch: (H1)-(H4) instances with i > s+1.
  Theorem C's D1 (i <= s+1) is proved inside the house branch.  It is
  not a consequence of (H1)-(H4) alone, and this is the certificate of
  that.  The decoded (s, i) = (1, 1) records at orders 1676 (G = Z_418,
  418 = 2*11*19) and 1772 (G = Z_442, 442 = 2*13*17) are re-read on
  subgroups K <= G of index i in {11, 19, 22, 38} and {13, 17, 26, 34}
  respectively -- every index i > 1 with w = n/i > 2s = 2 -- with

      Q'[(I,c)] := Q[I]  for every c in Gbar     (row repetition)
      P'[r][iJ + c] := P[r][J]                   (column repetition)
      M := 4 on all of Gbar,  i.e. M = 4 J_i

  and all four hypotheses re-checked in exact integer arithmetic, (H4)
  against the compressed Goethals-Seidel array over Gbar = Z_i built by
  tools/bordered_gs.py in the standard orientation.  Each is an
  (H1)-(H4) instance with i > s+1 = 2, and its Gram is 4 J_i =
  4i P_{trivial} -- exactly Theorem E's prediction, S = {trivial} being
  the unique Galois-stable singleton of the character group when i is
  prime.

  What the witnesses do NOT do is manufacture anything: the re-reading
  leaves P~ and Q~, hence the assembled matrix, unchanged.  The cert
  checks that too, so the record is what it is -- i <= s+1 fails in the
  general branch, and the failure is inert.

NEGATIVE CONTROLS.
  C1  a corrupted Parseval window (< 4i + 6 instead of < 4i + 2) must
      move the Part A census.
  C2  a perturbed Gram in Part B (M(1) := -4, the house value, with
      everything else unchanged) must fail (H1).
  Both must fire.  `--negative-control[=window|gram]` installs the
  corruption into the main path instead of the control slot, and that
  run must FAIL with a non-zero exit code.

Stdlib only.  Exact integers only.  No floating point.  No network.
Nothing is written anywhere.
"""

import hashlib
import json
import os
import sys
import time

TOOLS = "tools"
DATA_WITNESS = os.path.join("data", "gram34-witnesses.json")
DATA_TWIST = os.path.join("data", "n1916-twist.json")
DATA_RECORDS = os.path.join("data", "payload-records.json")

# data/gram34-witnesses.json is file-pinned here.  The other two are
# shared with certs 01/06/08/09 and are bound there by the canonical
# digests of the matrices they produce, so they are not file-pinned.
WITNESS_SHA = \
    "abaf4728e8ba5cd737024b9ab319640c8c634e497884a444b683fb5ee4b93307"

S = 3                       # the cell of Part A
I = 4
WINDOW = 4 * I + 2          # the true Parseval window: lambda < 18
BAD_WINDOW = 4 * I + 6      # C1's corruption
NO_WINDOW = 10 ** 9

# order -> the indices i = [G:K] swept in Part B (every divisor of n
# with i > 1 and w = n/i > 2s = 2)
GENERAL = ((1676, (11, 19, 22, 38)), (1772, (13, 17, 26, 34)))

# the relaxed counts of the Part A control table, per group:
#   (no window, integer entries) and (no window, even entries)
RELAXED = {"Z4": (48, 24), "Z2xZ2": (1154, 290)}


def die(msg):
    print("FAIL: %s" % msg)
    sys.exit(1)


def brief(ms, k=4):
    """A set of Grams, printed short: a failing control can produce
    dozens and the diagnosis is the count, not the list."""
    ms = sorted(ms)
    head = ", ".join(str(list(m)) for m in ms[:k])
    return ("%d = {%s}" % (len(ms), head) if len(ms) <= k
            else "%d = {%s, ...}" % (len(ms), head))


# ------------------------------------------------------- linear algebra

def matmul(A, B):
    return [[sum(A[a][k] * B[k][b] for k in range(len(B)))
             for b in range(len(B[0]))] for a in range(len(A))]


def mat_mul_t(A, B):
    return [[sum(x * y for x, y in zip(ra, rb)) for rb in B] for ra in A]


def eye(k, val=1):
    return [[val if a == b else 0 for b in range(k)] for a in range(k)]


def power_sums(A, upto):
    """[tr(A), tr(A^2), ..., tr(A^upto)], exact integers."""
    out = []
    Ak = eye(len(A))
    for _k in range(upto):
        Ak = matmul(Ak, A)
        out.append(sum(Ak[a][a] for a in range(len(A))))
    return out


def charpoly(A):
    """Monic integer characteristic polynomial, highest power first.

    Newton's identities from p_k = tr(A^k); every division is exact
    over Z for an integer matrix.
    """
    m = len(A)
    p = [0] + power_sums(A, m)
    e = [1] + [0] * m
    for k in range(1, m + 1):
        acc = 0
        for j in range(1, k + 1):
            acc += (-1) ** (j - 1) * e[k - j] * p[j]
        if acc % k:
            raise ValueError("Newton's identities left a remainder")
        e[k] = acc // k
    return [(-1) ** k * e[k] for k in range(m + 1)]


def poly_from_roots(roots):
    """Monic integer polynomial with exactly this integer root multiset."""
    c = [1]
    for r in roots:
        nxt = [0] * (len(c) + 1)
        for j, v in enumerate(c):
            nxt[j] += v
            nxt[j + 1] -= r * v
        c = nxt
    return c


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


# ----------------------------------------------------------- the groups

class Bar(object):
    """One of the two abelian groups of order 4, with its characters."""

    def __init__(self, name):
        self.name = name
        if name == "Z4":
            self.sub = lambda a, b: (a - b) % 4
            # chi_k(c) = i^{kc}; the real characters are k = 0 and k = 2
            self.real_chars = [[1, 1, 1, 1], [1, -1, 1, -1]]
            self.real_names = ["chi0 (real)", "chi2 (real)"]
        elif name == "Z2xZ2":
            self.sub = lambda a, b: a ^ b
            # chi_t(c) = prod_j (-1)^{t_j c_j}: all four are real
            self.real_chars = [[1 if bin(t & c).count("1") % 2 == 0 else -1
                                for c in range(4)] for t in range(4)]
            self.real_names = ["chi0=[+,+,+,+]", "chi1=[+,-,+,-]",
                               "chi2=[+,+,-,-]", "chi3=[+,-,-,+]"]
        else:
            raise ValueError(name)

    def gram(self, M):
        return [[M[self.sub(c, cp)] for cp in range(4)] for c in range(4)]

    def spectrum(self, M):
        """The four eigenvalues, from the character formulas.

        Z_2 x Z_2:  lambda_t = sum_c (-1)^{popcount(t & c)} M(c).
        Z_4:        M(1) = M(3) = a by A2, b = M(2), and
                    lambda_k = M(0) + a (i^k + i^{3k}) + b i^{2k}, i.e.
                    (M0 + 2a + b, M0 - b, M0 - 2a + b, M0 - b).
        """
        if self.name == "Z2xZ2":
            return sorted((sum(ch[c] * M[c] for c in range(4))
                           for ch in self.real_chars), reverse=True)
        if M[1] != M[3]:
            raise ValueError("Z_4 Gram is not symmetric: %r" % (M,))
        a, b = M[1], M[2]
        return sorted([M[0] + 2 * a + b, M[0] - b, M[0] - 2 * a + b,
                       M[0] - b], reverse=True)

    def symmetric(self, M):
        return all(M[self.sub(0, c)] == M[c] for c in range(4))

    def house_orbit(self):
        """The real-character twists of M_house = (4s+4) I_4 - 4 J_4."""
        house = [4 * S if c == 0 else -4 for c in range(4)]
        return [(tuple(ch[c] * house[c] for c in range(4)), nm)
                for ch, nm in zip(self.real_chars, self.real_names)]


# ------------------------------------------------- Part A: the census

def census(bar, window, bad):
    """One exhaustive pass over the candidate box.

    A1 fixes M(0) = 12 and A3 bounds |M(c)| <= 12, so {-12..12} on the
    off-origin entries is exhaustive: no Gram of +-1 rows of length 12
    lies outside it.  A2 is imposed by construction for Z_4 (M(1) =
    M(3)) and is automatic for Z_2 x Z_2 (every element is its own
    inverse); it is re-checked anyway.

    Returns (survivors, relaxed_int, relaxed_even, window_no_parity)
    where the three counts are the control table of Sec.2 of the source
    classification: A1-A2 + PSD + rank <= 3 with no window, the same
    with even entries, and the full filter with A3's parity dropped.
    """
    rng = range(-12, 13)
    if bar.name == "Z4":
        cand = ((12, a, b, a) for a in rng for b in rng)
    else:
        cand = ((12, a, b, d) for a in rng for b in rng for d in rng)
    survivors, relaxed_int, relaxed_even, no_parity = [], 0, 0, 0
    for M in cand:
        if not bar.symmetric(M):
            bad.append("%s: candidate %s is not symmetric" % (bar.name, M))
            return [], 0, 0, 0
        lam = bar.spectrum(M)
        # the character formulas, verified per matrix against tr(M^k)
        if charpoly(bar.gram(M)) != poly_from_roots(lam):
            bad.append("%s: the character formula disagrees with the "
                       "characteristic polynomial at M = %s"
                       % (bar.name, list(M)))
            return [], 0, 0, 0
        if min(lam) < 0:                        # A4
            continue
        if lam[3] != 0:                         # A5: rank <= s = 3
            continue
        even = not any(v % 2 for v in M)        # A3's parity clause
        relaxed_int += 1
        relaxed_even += 1 if even else 0
        if max(lam) >= window:                  # A6
            continue
        no_parity += 1
        if even:
            survivors.append(tuple(M))
    return sorted(survivors), relaxed_int, relaxed_even, no_parity


def check_witnesses(bank, admissible, bad):
    """Every banked Q is +-1 with Q Q^T = I_4 (x) M, and the six of them
    cover the six admissible Grams."""
    seen = {}
    for w in bank["witnesses"]:
        tag = "%s/%s" % (w["group"], w["twist"])
        bar = Bar(w["group"])
        M = list(w["M"])
        Q = []
        for row in w["Q"]:
            if len(row) != 4 * S or set(row) - set("+-"):
                bad.append("witness %s: bad Q row %r" % (tag, row))
                return seen
            Q.append([1 if ch == "+" else -1 for ch in row])
        if len(Q) != 4 * I:
            bad.append("witness %s: Q has %d rows, expected %d"
                       % (tag, len(Q), 4 * I))
            return seen
        tgt = [[0] * (4 * I) for _ in range(4 * I)]
        for a in range(4):
            for c1 in range(I):
                for c2 in range(I):
                    tgt[a * I + c1][a * I + c2] = M[bar.sub(c1, c2)]
        if mat_mul_t(Q, Q) != tgt:
            bad.append("witness %s: Q Q^T != I_4 (x) M" % tag)
        if bar.spectrum(M) != list(w["spectrum"]):
            bad.append("witness %s: banked spectrum %s, measured %s"
                       % (tag, w["spectrum"], bar.spectrum(M)))
        seen.setdefault(w["group"], set()).add(tuple(M))
    for gname, lst in sorted(admissible.items()):
        if seen.get(gname, set()) != set(lst):
            bad.append("%s: the banked witnesses cover %s; the admissible "
                       "list is %s" % (gname, brief(seen.get(gname, set())),
                                       brief(lst)))
    return seen


def check_1916(admissible, bad):
    """cert 09's four Grams at order 1916 are the whole Z_2 x Z_2 list."""
    with open(DATA_TWIST, "r", encoding="ascii") as fh:
        tw = json.load(fh)
    got = set()
    for _name, inst in sorted(tw["instances"].items()):
        gm = inst["gram_M"]
        got.add(tuple(gm[str(c)] for c in range(4)))
    if got != set(admissible["Z2xZ2"]):
        bad.append("the banked 1916 Grams %s are not the Z_2 x Z_2 "
                   "admissible list %s"
                   % (brief(got), brief(admissible["Z2xZ2"])))
    return sorted(got)


# ------------------------------------- Part B: the general-branch reads

def galois_stable_singletons(i):
    """The a in Z_i with {k a mod i : gcd(k, i) = 1} = {a}.

    Characters of Z_i are chi_a(c) = zeta_i^{ac}, and sigma_k sends
    chi_a to chi_{ka}; a singleton {chi_a} is Galois-stable iff a is
    fixed by every unit.  For i prime the only one is a = 0, the trivial
    character -- which is exactly the S that Theorem E predicts here.
    """
    units = [k for k in range(1, i) if gcd(k, i) == 1]
    return [a for a in range(i) if all((k * a) % i == a for k in units)]


def aggregate_paf(G, seeds):
    """sum_q PAF_q(t) for every t in G, straight from the definition."""
    n = G.n
    out = []
    for t in range(n):
        shift = [G.add(h, t) for h in range(n)]
        out.append(sum(sum(x[h] * x[shift[h]] for h in range(n))
                       for x in seeds))
    return out


def general_read(T, rec, idx, paf, gram_delta, bad):
    """Re-read an (s,i) = (1,1) record on K <= G of index idx."""
    G = T.AbelianGroup(rec["group"])
    n, s, N = G.n, int(rec["s"]), int(rec["order"])
    if s != 1:
        bad.append("order %d: expected an s = 1 record" % N)
        return None
    seeds = [T.signs(x) for x in rec["seeds"]]
    rho = G.idx(tuple(rec["r_shift"]))
    E = [T.signs(r) for r in rec["corner"]]
    P1 = [T.signs(r) for r in rec["row_table"]]
    C1 = [T.signs(r) for r in rec["col_table"]]
    Q1 = [[C1[r][k] for r in range(4 * s)] for k in range(4)]

    if n % idx:
        bad.append("order %d: %d does not divide n = %d" % (N, idx, n))
        return None
    kappa, i = T.coset_map(G, [idx])
    w = n // i
    rep = {"order": N, "n": n, "s": s, "i": i, "w": w,
           "rho_bar": kappa[rho]}

    # the i = 1 column table is a 4x4 Hadamard matrix
    rep["Q1_hadamard"] = (mat_mul_t(Q1, Q1) == eye(4, 4))

    # the Gram of the re-reading: M = 4 on all of Gbar, i.e. M = 4 J_i
    M = [4 + (gram_delta if c == 1 else 0) for c in range(i)]

    # (H1): Q Q^T = I_4 (x) M with Q the row repetition of Q1
    Q = [Q1[a][:] for a in range(4) for _c in range(i)]
    tgt = [[0] * (4 * i) for _ in range(4 * i)]
    for a in range(4):
        for c1 in range(i):
            for c2 in range(i):
                tgt[a * i + c1][a * i + c2] = M[(c1 - c2) % i]
    rep["H1"] = (mat_mul_t(Q, Q) == tgt)

    # (H2): sum PAF(t) = -M(kappa t) off the origin, 4n at the origin
    rep["H2"] = (paf[0] == 4 * n and
                 all(paf[t] == -M[kappa[t]] for t in range(1, n)))

    # (H3): unchanged from i = 1, since P' P'^T = i P P^T and w i = n
    P = [[P1[r][J] for J in range(4) for _c in range(i)]
         for r in range(4 * s)]
    EEt = mat_mul_t(E, E)
    PPt = mat_mul_t(P, P)
    rep["H3"] = all(EEt[a][b] + w * PPt[a][b] == (N if a == b else 0)
                    for a in range(4 * s) for b in range(4 * s))

    # (H4): against the compressed GS array over Gbar = Z_i
    Gq = T.QuotientGroup([idx])
    sigma = [[0] * i for _ in range(4)]
    for q in range(4):
        for g in range(n):
            sigma[q][kappa[g]] += seeds[q][g]
    Chat = T.gs_array(Gq, sigma, kappa[rho])
    rep["H4"] = all(sum(E[r][j] * Q[k][j] for j in range(4 * s)) +
                    sum(P[r][cc] * Chat[k][cc] for cc in range(4 * i)) == 0
                    for r in range(4 * s) for k in range(4 * i))
    rep["sigma_row_sums"] = [sum(sig) for sig in sigma]

    # Theorem E's prediction: the Gram is 4i P_{trivial} = 4 J_i, of
    # rank 1 = s and spectrum {4i, 0^(i-1)}.
    Mmat = [[M[(c - cp) % i] for cp in range(i)] for c in range(i)]
    rep["gram_is_4J"] = (Mmat == [[4] * i for _ in range(i)])
    rows = set(tuple(r) for r in Mmat)
    rep["rank_is_s"] = (len(rows) == 1 and any(next(iter(rows))))
    rep["power_sums_ok"] = (power_sums(Mmat, 4) ==
                            [(4 * i) ** k for k in range(1, 5)])
    stable = galois_stable_singletons(i)
    rep["galois_stable_singletons"] = stable
    rep["S_trivial_is_stable"] = (0 in stable)
    rep["i_gt_s_plus_1"] = (i > s + 1)
    rep["w_gt_2s"] = (w > 2 * s)

    # the re-reading is inert: P~ and Q~ do not depend on the index
    probe = sorted(set([0, 1, n // 3, n // 2, n - 1]))
    rep["strips_unchanged"] = (
        all(P[r][i * J + kappa[g]] == P1[r][J]
            for r in range(4 * s) for J in range(4) for g in probe) and
        all(Q[i * a + kappa[g]] == Q1[a]
            for a in range(4) for g in probe))

    keys = ("Q1_hadamard", "H1", "H2", "H3", "H4", "gram_is_4J",
            "rank_is_s", "power_sums_ok", "S_trivial_is_stable",
            "i_gt_s_plus_1", "w_gt_2s", "strips_unchanged")
    rep["ok"] = all(rep[k] for k in keys)
    if not rep["ok"]:
        bad.append("order %d at index %d: %s"
                   % (N, i, ", ".join(k for k in keys if not rep[k])))
    return rep


# ----------------------------------------------------------------- main

def main(argv):
    window = WINDOW
    delta = 0
    installed = None
    for a in argv[1:]:
        if a in ("--negative-control", "--negative-control=window"):
            window = BAD_WINDOW
            installed = "corrupted Parseval window (< %d)" % BAD_WINDOW
        elif a == "--negative-control=gram":
            delta = -8
            installed = "perturbed Gram in Part B (M(1) := -4)"
        else:
            die("unknown argument %r (expected --negative-control"
                "[=window|gram])" % a)

    for path in (os.path.join(TOOLS, "bordered_gs.py"), DATA_WITNESS,
                 DATA_TWIST, DATA_RECORDS):
        if not os.path.isfile(path):
            die("missing %s -- run this from the repository root: "
                "python certs/12-gram-rigidity/run.py" % path)
    sys.path.insert(0, TOOLS)
    import bordered_gs as T

    bad = []
    t0 = time.time()
    if installed:
        print("  NEGATIVE CONTROL INSTALLED: %s" % installed)
        print("  this run is expected to FAIL; a PASS here would mean the "
              "cert does not see what it claims to check")

    # ---- [0] the witness bank is pinned ----
    with open(DATA_WITNESS, "rb") as fh:
        blob = fh.read()
    dig = hashlib.sha256(blob).hexdigest()
    if dig != WITNESS_SHA:
        bad.append("%s sha256 %s != pinned %s"
                   % (DATA_WITNESS, dig, WITNESS_SHA))
    bank = json.loads(blob.decode("ascii"))
    print("  [0] %s pinned: sha256 %s..%s, %d witnesses"
          % (DATA_WITNESS, dig[:8], dig[-6:], len(bank["witnesses"])))

    # ---- [1] Part A: the census ----
    print("  [1] the complete admissible-Gram classification at "
          "(s,i) = (%d,%d), both abelian groups of order 4" % (S, I))
    admissible = {}
    for name, want in (("Z4", 2), ("Z2xZ2", 4)):
        bar = Bar(name)
        full, rint, reven, nopar = census(bar, window, bad)
        admissible[name] = full
        orbit = bar.house_orbit()
        orbit_set = set(m for m, _nm in orbit)
        if len(full) != want:
            bad.append("%s: %d admissible Grams, expected %d"
                       % (name, len(full), want))
        if set(full) != orbit_set:
            bad.append("%s: the admissible list %s is not the house twist "
                       "orbit %s" % (name, brief(full), brief(orbit_set)))
        if len(orbit_set) != len(bar.real_chars):
            bad.append("%s: the twist orbit has %d members, but Gbar has %d "
                       "real characters"
                       % (name, len(orbit_set), len(bar.real_chars)))
        off = [M for M in full if bar.spectrum(list(M)) != [16, 16, 16, 0]]
        if off:
            bad.append("%s: %s admissible Grams do not have spectrum "
                       "{16,16,16,0}, e.g. M = %s with %s"
                       % (name, len(off), list(off[0]),
                          bar.spectrum(list(off[0]))))
        if nopar != len(full):
            bad.append("%s: dropping A3's parity changes the census, "
                       "%d vs %d" % (name, nopar, len(full)))
        if (rint, reven) != RELAXED[name] and not installed:
            bad.append("%s: the no-window control table is (%d, %d), "
                       "expected %s" % (name, rint, reven, RELAXED[name]))
        print("      %-6s admissible %d   (drop A6: %d integer / %d even; "
              "keep A6, drop A3's parity: %d)   real characters of Gbar: %d"
              % (name, len(full), rint, reven, nopar, len(bar.real_chars)))
        for M, nm in orbit:
            print("        M = %-19s twist %-15s spectrum %-16s %s"
                  % (str(list(M)), nm, str(bar.spectrum(list(M))),
                     "in the census" if M in set(full) else "MISSING"))

    # ---- [2] the six witnesses ----
    seen = check_witnesses(bank, admissible, bad)
    print("  [2] the %d banked witnesses: Q in {+-1}^(%dx%d), all %d "
          "inner products of Q Q^T re-checked against I_4 (x) M; coverage "
          "%s -- every admissible Gram is (H1)-realizable"
          % (len(bank["witnesses"]), 4 * I, 4 * S, (4 * I) ** 2,
             ", ".join("%s %d/%d" % (g, len(seen.get(g, ())),
                                     len(admissible[g]))
                       for g in ("Z4", "Z2xZ2"))))

    # ---- [3] the 1916 cross-check ----
    got1916 = check_1916(admissible, bad)
    print("  [3] cert 09's four banked Grams at order 1916 vs the "
          "Z_2 x Z_2 admissible list: %s"
          % ("identical -- that sweep is the complete classification at "
             "this cell" if set(got1916) == set(admissible["Z2xZ2"])
             else "DIFFERENT"))

    # ---- [4] Part B: the general branch ----
    print("  [4] the general branch: (H1)-(H4) instances with i > s+1")
    with open(DATA_RECORDS, "r", encoding="ascii") as fh:
        records = dict((int(r["order"]), r)
                       for r in json.load(fh)["orders"])
    pafs = {}
    reports = []
    for order, indices in GENERAL:
        rec = records[order]
        G = T.AbelianGroup(rec["group"])
        pafs[order] = aggregate_paf(G, [T.signs(x) for x in rec["seeds"]])
        for idx in indices:
            rep = general_read(T, rec, idx, pafs[order], delta, bad)
            if rep is None:
                continue
            reports.append(rep)
            print("      order %d  i=%-3d w=%-3d rho_bar=%-3d  H1 %s H2 %s "
                  "H3 %s H4 %s  Gram = 4 J_%-3d = 4i P_triv %s  i > s+1 %s"
                  "  Galois-stable singletons of Z^_%d: %s"
                  % (rep["order"], rep["i"], rep["w"], rep["rho_bar"],
                     rep["H1"], rep["H2"], rep["H3"], rep["H4"], rep["i"],
                     rep["gram_is_4J"] and rep["power_sums_ok"],
                     rep["i_gt_s_plus_1"], rep["i"],
                     rep["galois_stable_singletons"]))
    if reports and all(r["strips_unchanged"] for r in reports):
        print("      every re-reading leaves P~ and Q~ unchanged: the "
              "general branch here carries no border and no matrix the "
              "i = 1 construction does not already carry")

    # ---- [5] the negative controls ----
    if not installed:
        moved = []
        for name in ("Z4", "Z2xZ2"):
            got = census(Bar(name), BAD_WINDOW, [])[0]
            if set(got) != set(admissible[name]):
                moved.append("%s %d -> %d"
                             % (name, len(admissible[name]), len(got)))
        if not moved:
            bad.append("C1: the corrupted window left the census unchanged")
        print("  [5] C1 -- corrupted Parseval window (< %d instead of "
              "< %d): the census moves (%s)"
              % (BAD_WINDOW, WINDOW,
                 "; ".join(moved) if moved else "IT DOES NOT"))

        c2 = []
        rep = general_read(T, records[1676], 11, pafs[1676], -8, c2)
        fired = rep is not None and not rep["H1"]
        if not fired:
            bad.append("C2: the perturbed Gram M(1) = -4 did not fail (H1)")
        print("      C2 -- perturbed Gram at 1676, index 11 (M(1) := -4, "
              "the house value): (H1) %s"
              % ("fails, as it must" if fired else "STILL HOLDS"))

    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 12: FAIL (%d problems)" % len(bad))
        return 1
    print("CERT 12: PASS -- at (%d,%d) the admissible Grams are exactly the "
          "real-character twists of the house form (%d for Z_4, %d for "
          "Z_2 x Z_2), all (H1)-realizable, all of spectrum {16,16,16,0}; "
          "and %d general-branch instances with i > s+1 carry the Gram "
          "4i P_trivial (%.1fs)"
          % (S, I, len(admissible["Z4"]), len(admissible["Z2xZ2"]),
             len(reports), time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
