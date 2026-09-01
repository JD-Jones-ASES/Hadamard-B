#!/usr/bin/env python3
"""full_recompute.py -- the exact |T4| 4-profile, recomputed with numpy.

NOT IN THE TRUST CHAIN, AND NOT ON THE DEFAULT PATH.  cert 06's `run.py`
imports this module only under `--full`; nothing else in the repository
touches it, and no claim rests on it that is not already established by
the banked profiles the default path audits.  numpy is finder-side only.
The trust chain for the matrices themselves is `verify/verify.py`, which
is standard library and exact-integer throughout.

THE COMPUTATION.  For rows i != j write u_ij[c] = H[i][c]*H[j][c], and let
U be the C(n,2) x n sign matrix of those pair vectors.  For a 4-subset
{i,j,k,l},  T4 = <u_ij, u_kl>  for each of the three ways of splitting the
subset into two disjoint pairs, so the histogram of |<u_P,u_Q>| over
UNORDERED pairs {P,Q} with P != Q equals

    3 * (the 4-profile)  +  n*C(n-1,2) counts at the value 0,

the second term being the pairs of row-pairs that share an index, whose
inner product is sum_c H[j][c]H[k][c] = 0 by orthogonality.  Take the
upper triangle of U U^T, remove the diagonal and those intersecting pairs,
divide by three.

TWO INDEPENDENT ARITHMETIC PATHS, which must agree bin for bin:

  impl='blas'  float32 U U^T in row blocks, np.bincount.  This is EXACT:
               every entry is an integer of absolute value <= n < 2^24,
               and every partial sum of +-1 terms is likewise an exact
               integer below 2^24, which float32 represents exactly, so
               no rounding can occur at these sizes.
  impl='bits'  pure integer path: rows packed into uint64 words,
               |T4| = |n - 2*popcount(u_P xor u_Q)|, accumulated word by
               word.  Shares no arithmetic with the BLAS path.

Ported from the lab's experiments/inequiv/exact_profile.py (Hadamard-2060),
which is where the banked profiles in data/ were produced; `run.py` compares
what this module returns against those banks bin for bin.

Cost at n = 668 on a desktop with three BLAS threads: about 280 s for
'blas' and about 1 480 s for 'bits', per matrix.
"""

import sys

sys.dont_write_bytecode = True


def _pair_matrix(rows, n, np):
    H = np.array([[1 if c == "+" else -1 for c in r] for r in rows],
                 dtype=np.int8)
    i, j = np.triu_indices(n, 1)
    return (H[i] * H[j]).astype(np.int8)


def _finish(hist, n, m):
    """hist = |<u_P,u_Q>| over the upper triangle INCLUDING the diagonal."""
    h = dict(hist)
    h[n] = h.get(n, 0) - m                        # drop P == Q
    inter = n * (n - 1) * (n - 2) // 2            # pairs sharing one index
    h[0] = h.get(0, 0) - inter
    out = {}
    for k, v in h.items():
        if v == 0:
            continue
        assert v % 3 == 0, "histogram not divisible by 3 at %d: %d" % (k, v)
        out[k] = v // 3
    return out


def _popcount_words(np, x):
    """popcount of a uint64 array, with a fallback for numpy < 2.0."""
    if hasattr(np, "bitwise_count"):
        return np.bitwise_count(x)
    y = x
    m1 = np.uint64(0x5555555555555555)
    m2 = np.uint64(0x3333333333333333)
    m4 = np.uint64(0x0F0F0F0F0F0F0F0F)
    y = y - ((y >> np.uint64(1)) & m1)
    y = (y & m2) + ((y >> np.uint64(2)) & m2)
    y = (y + (y >> np.uint64(4))) & m4
    return ((y * np.uint64(0x0101010101010101)) >> np.uint64(56))


def _profile_blas(rows, n, block=128, progress=True):
    import numpy as np
    import time
    U = _pair_matrix(rows, n, np).astype(np.float32)
    m = U.shape[0]
    hist = np.zeros(n + 1, dtype=np.int64)
    t0 = time.time()
    for a in range(0, m, block):
        b = min(block, m - a)
        G = U[a:a + b] @ U[a:].T
        Gi = np.abs(G).astype(np.int32)
        hist += np.bincount(Gi.ravel(), minlength=n + 1)
        # the leading b x b square counts each in-block pair twice; drop the
        # strictly lower triangle so every unordered pair is counted once
        lo = Gi[:, :b][np.tril_indices(b, -1)]
        if lo.size:
            hist -= np.bincount(lo, minlength=n + 1)
        if progress and (a // block) % 200 == 0:
            print("        blas %d/%d  (%.0fs)" % (a, m, time.time() - t0),
                  flush=True)
    return _finish({int(k): int(v) for k, v in enumerate(hist) if v}, n, m)


def _profile_bits(rows, n, block=16, progress=True):
    import numpy as np
    import time
    U = _pair_matrix(rows, n, np)
    m = U.shape[0]
    nw = (n + 63) // 64
    P = np.zeros((m, nw), dtype=np.uint64)
    bits = (U < 0).astype(np.uint8)               # 1 where the entry is -1
    for w in range(nw):
        lo, hi = 64 * w, min(64 * (w + 1), n)
        chunk = bits[:, lo:hi].astype(np.uint64)
        acc = np.zeros(m, dtype=np.uint64)
        for t in range(hi - lo):
            acc |= chunk[:, t] << np.uint64(t)
        P[:, w] = acc
    hist = np.zeros(n + 1, dtype=np.int64)
    t0 = time.time()
    for a in range(0, m, block):
        b = min(block, m - a)
        d = np.zeros((b, m - a), dtype=np.int16)
        for w in range(nw):
            d += _popcount_words(
                np, P[a:a + b, w][:, None] ^ P[a:, w][None, :]
            ).astype(np.int16)
        t4 = np.abs(np.int32(n) - 2 * d.astype(np.int32))
        hist += np.bincount(t4.ravel(), minlength=n + 1)
        lo = t4[:, :b][np.tril_indices(b, -1)]
        if lo.size:
            hist -= np.bincount(lo, minlength=n + 1)
        if progress and (a // block) % 2000 == 0:
            print("        bits %d/%d  (%.0fs)" % (a, m, time.time() - t0),
                  flush=True)
    return _finish({int(k): int(v) for k, v in enumerate(hist) if v}, n, m)


def profile(rows, n, impl, progress=True):
    """The full |T4| 4-profile as {bin: count}.  impl in {'blas', 'bits'}."""
    assert len(rows) == n and all(len(r) == n for r in rows)
    if impl == "blas":
        return _profile_blas(rows, n, progress=progress)
    if impl == "bits":
        return _profile_bits(rows, n, progress=progress)
    raise ValueError("impl must be 'blas' or 'bits', not %r" % (impl,))
