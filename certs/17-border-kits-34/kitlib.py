"""kitlib.py -- this certificate's own (3,4) border-kit engine.  Standard library,
exact integers.  Written for this certificate; it shares no code with the
kit engine of the source laboratory, and it imports nothing from outside
this repository.

Model (note/NOTE-B.md S1.8, Theorem F (a),(b),(d),(e)):
at (s,i) = (3,4), house Gram M = 16 I - 4 J, w > 2s,
  * Q (16 x 12) has 12 mutually orthogonal block-balanced +-1 columns
    (block = the four rows iI + c, c in Gbar);  <=> (H1);
  * P (12 x 16) has 12 mutually orthogonal block-balanced +-1 rows;
  * E = -(1/16) P Chat^T Q, and E E^T = 12 I is then automatic, so the only
    remaining condition is E in {+-1}; Chat may be replaced by Chat0 (S-part)
    since P kills the per-seed constants.
Search strategy (deliberately different from the source laboratory's
biclique-over-P backtracking): FIX Q from a list of candidate 12-cliques
(first the Kronecker table H4 (x) [b1 b2 b3], then random cliques from a
fixed seed), compute V = (4 Chat0)^T Q once,
filter the 648 sign-representatives p by p.V in {+-64}^12, and look for a
12-clique in the orthogonality graph of the admissible p (bitmask DFS).
"""
import itertools
import random

# ---------------------------------------------------------------- groups of order 4
GROUPS = {
    "Z4": dict(elts=[0, 1, 2, 3],
               add=lambda a, b: (a + b) % 4,
               neg=lambda a: (-a) % 4),
    "Z2xZ2": dict(elts=[(0, 0), (0, 1), (1, 0), (1, 1)],
                  add=lambda a, b: ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2),
                  neg=lambda a: a),
}


def chars(gname):
    """Nontrivial characters as functions c -> complex, in a fixed order."""
    if gname == "Z2xZ2":
        ts = [(0, 1), (1, 0), (1, 1)]
        return [(lambda t: (lambda c: (-1) ** ((t[0] * c[0] + t[1] * c[1]) % 2)))(t) for t in ts]
    # Z4: chi_1 (complex), chi_2 (real), chi_3 = conj(chi_1)
    return [lambda c: 1j ** (c % 4), lambda c: (-1) ** (c % 4), lambda c: 1j ** ((3 * c) % 4)]


def gs_array(sig, gname, krho):
    """Standard-orientation GS array over the order-4 group of four integer
    functions sig[q] (indexed like GROUPS[gname]['elts']).  16 x 16."""
    G = GROUPS[gname]
    els, add, neg = G["elts"], G["add"], G["neg"]
    idx = {e: k for k, e in enumerate(els)}

    def dev(f):
        return [[f[idx[add(h, neg(g))]] for h in els] for g in els]

    def xr(f):    # (X R)[g,h] = x(krho - g - h)
        return [[f[idx[add(krho, neg(add(g, h)))]] for h in els] for g in els]

    def xtr(f):   # (X^T R)[g,h] = x(g + h - krho)
        return [[f[idx[add(add(g, h), neg(krho))]] for h in els] for g in els]

    def ng(X):
        return [[-v for v in r] for r in X]
    A = dev(sig[0])
    BR, CR, DR = xr(sig[1]), xr(sig[2]), xr(sig[3])
    BtR, CtR, DtR = xtr(sig[1]), xtr(sig[2]), xtr(sig[3])
    B = [[A, BR, CR, DR],
         [ng(BR), A, DtR, ng(CtR)],
         [ng(CR), ng(DtR), A, BtR],
         [ng(DR), CtR, ng(BtR), A]]
    return [[B[I][J][g][h] for J in range(4) for h in range(4)] for I in range(4) for g in range(4)]


# ---------------------------------------------------------------- the pool
BAL4 = [v for v in itertools.product((1, -1), repeat=4) if sum(v) == 0]          # 6
POOL = [sum(t, ()) for t in itertools.product(BAL4, repeat=4)]                      # 1296
REPS = [v for v in POOL if v[0] == 1]                                               # 648


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


ORTH = []
for a in range(len(REPS)):
    m = 0
    for b in range(len(REPS)):
        if dot(REPS[a], REPS[b]) == 0:
            m |= 1 << b
    ORTH.append(m)


def popcount(x):
    return bin(x).count("1")


def find_clique(mask, k):
    out = []

    def rec(cand, need):
        if need == 0:
            return True
        while cand:
            v = cand.bit_length() - 1
            cand &= ~(1 << v)
            if popcount(cand) + 1 < need:
                return False
            out.append(v)
            if rec(cand & ORTH[v], need - 1):
                return True
            out.pop()
        return False
    return list(out) if rec(mask, k) else None


# ---------------------------------------------------------------- candidate Q's
H4 = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]
B3 = [(1, 1, -1, -1), (1, -1, -1, 1), (1, -1, 1, -1)]


def kronecker_Q():
    """Q[4I + c][3b + a] = H4[I][b] * B3[a][c]: 12 orthogonal block-balanced columns."""
    return [[H4[I][b] * B3[a][c] for b in range(4) for a in range(3)] for I in range(4) for c in range(4)]


def random_Q(rng):
    """A random 12-clique of pool vectors, as a 16 x 12 matrix."""
    while True:
        order = list(range(len(REPS)))
        rng.shuffle(order)
        chosen = []
        cand = (1 << len(REPS)) - 1
        for v in order:
            if (cand >> v) & 1:
                chosen.append(v)
                cand &= ORTH[v]
                if len(chosen) == 12:
                    break
        if len(chosen) == 12:
            cols = [REPS[v] for v in chosen]
            return [[cols[j][k] for j in range(12)] for k in range(16)]


def is_valid_Q(Q):
    """(H1) at the house Gram: Q Q^T = I4 (x) (16 I - 4 J)."""
    for a in range(16):
        for b in range(16):
            want = (12 if a == b else -4) if a // 4 == b // 4 else 0
            if dot(Q[a], Q[b]) != want:
                return False
    return True


# ---------------------------------------------------------------- the search

def kit_for(T, gname, krho, Q_list):
    """T = 4*sigma0 (4x4 ints, rows q, columns c as GROUPS[gname]['elts']).
    Returns (E, P, Q, tries) or (None, None, None, tries)."""
    C4 = gs_array(T, gname, krho)                    # 4 Chat0
    C4T = [list(x) for x in zip(*C4)]
    for tries, Q in enumerate(Q_list, 1):
        # V = (4 Chat0)^T Q : 16 x 12
        V = [[sum(C4T[a][k] * Q[k][j] for k in range(16)) for j in range(12)] for a in range(16)]
        Vt = [list(x) for x in zip(*V)]
        mask = 0
        for pi, p in enumerate(REPS):
            ok = True
            for j in range(12):
                d = dot(p, Vt[j])
                if d != 64 and d != -64:
                    ok = False
                    break
            if ok:
                mask |= 1 << pi
        if popcount(mask) < 12:
            continue
        cl = find_clique(mask, 12)
        if cl is None:
            continue
        P = [list(REPS[a]) for a in cl]
        E = [[-dot(P[r], Vt[j]) // 64 for j in range(12)] for r in range(12)]
        return E, P, Q, tries
    return None, None, None, len(Q_list)


def verify_kit(E, P, Q, sigma, gname, krho, w):
    """Exact (H1),(H3),(H4) with the FULL sigma (r-part included) at width w.
    Returns dict of booleans."""
    N = 4 * (4 * w + 3)
    out = {}
    out["H1"] = is_valid_Q(Q)
    EEt = [[dot(E[a], E[b]) for b in range(12)] for a in range(12)]
    PPt = [[dot(P[a], P[b]) for b in range(12)] for a in range(12)]
    out["E_hadamard"] = all(EEt[a][b] == (12 if a == b else 0) for a in range(12) for b in range(12))
    out["PPt_16I"] = all(PPt[a][b] == (16 if a == b else 0) for a in range(12) for b in range(12))
    out["H3"] = all(EEt[a][b] + w * PPt[a][b] == (N if a == b else 0) for a in range(12) for b in range(12))
    Ch = gs_array(sigma, gname, krho)
    out["H4"] = all(dot(E[r], Q[k]) + dot(P[r], Ch[k]) == 0 for r in range(12) for k in range(16))
    out["E_pm1"] = all(abs(v) == 1 for row in E for v in row)
    return out


def default_Q_list(rng, extra=24):
    L = [kronecker_Q()]
    for _ in range(extra):
        L.append(random_Q(rng))
    for Q in L:
        assert is_valid_Q(Q)
    return L
