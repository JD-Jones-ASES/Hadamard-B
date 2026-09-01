#!/usr/bin/env python3
"""verify.py -- the trust chain of this repository.

Checks that a matrix file encodes a Hadamard matrix: square, entries in
{+1, -1}, and H * H^T = n * I exactly.

Standard library only. Exact integer arithmetic only: rows are packed into
Python integers; the dot product of two +-1 rows equals n - 2*popcount(xor),
so orthogonality of a pair is exactly popcount(xor) == n/2. No floats
anywhere.

File format (auto-detected per line; blank lines and '#' comments ignored):
  compact:   one row per line of '+' and '-' characters
  numeric:   one row per line of whitespace/comma-separated +1/-1 entries

Exit codes: 0 = verified Hadamard (or selftest passed), 1 = verification
failed, 2 = usage/parse error.

Runs on bare python3 >= 3.9 on Windows and macOS.
"""

import argparse
import hashlib
import sys

if hasattr(int, "bit_count"):  # 3.10+
    def _popcount(x):
        return x.bit_count()
else:  # 3.9 fallback

    def _popcount(x):
        return bin(x).count("1")


def parse_matrix(text):
    """Parse matrix text -> list of rows of +-1 ints. Raises ValueError."""
    rows = []
    for ln, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if set(s) <= {"+", "-"}:
            row = [1 if c == "+" else -1 for c in s]
        else:
            row = []
            for tok in s.replace(",", " ").split():
                if tok in ("1", "+1"):
                    row.append(1)
                elif tok == "-1":
                    row.append(-1)
                else:
                    raise ValueError(
                        "line %d: entry %r not in {+1,-1}" % (ln, tok)
                    )
        rows.append(row)
    if not rows:
        raise ValueError("no matrix rows found")
    n = len(rows)
    for i, r in enumerate(rows):
        if len(r) != n:
            raise ValueError(
                "row %d has length %d, expected %d (square matrix required)"
                % (i, len(r), n)
            )
    return rows


def canonical_sha256(rows):
    """SHA-256 of the canonical '+/-' serialization (newline-joined)."""
    text = "\n".join(
        "".join("+" if v == 1 else "-" for v in row) for row in rows
    ) + "\n"
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def check_hadamard(rows, progress=False):
    """Return (ok, message, pairs_checked). Exact arithmetic throughout."""
    n = len(rows)
    if n > 2 and n % 4 != 0:
        return (False, "order %d is not 1, 2, or divisible by 4" % n, 0)
    if n == 1:
        return (True, "order 1 (trivial)", 0)
    half, rem = divmod(n, 2)
    if rem:
        return (False, "odd order %d cannot be Hadamard" % n, 0)
    packed = []
    for row in rows:
        bits = 0
        for j, v in enumerate(row):
            if v == 1:
                bits |= 1 << j
        packed.append(bits)
    total = n * (n - 1) // 2
    done = 0
    bad = []
    for i in range(n):
        ri = packed[i]
        for j in range(i + 1, n):
            if _popcount(ri ^ packed[j]) != half:
                bad.append((i, j))
                if len(bad) >= 5:
                    return (
                        False,
                        "non-orthogonal row pairs (first %d): %s"
                        % (len(bad), bad),
                        done + j - i,
                    )
        done += n - i - 1
        if progress and n >= 512 and i % 256 == 255:
            sys.stderr.write(
                "  progress: %d/%d pairs\n" % (done, total)
            )
            sys.stderr.flush()
    if bad:
        return (False, "non-orthogonal row pairs: %s" % bad, done)
    return (True, "all %d row pairs orthogonal" % total, done)


def verify_file(path, progress=False):
    """Verify one file. Returns process exit code, printing the verdict."""
    try:
        with open(path, "r", encoding="ascii") as fh:
            text = fh.read()
    except OSError as exc:
        print("ERROR: cannot read %s: %s" % (path, exc))
        return 2
    try:
        rows = parse_matrix(text)
    except ValueError as exc:
        print("ERROR: parse: %s" % exc)
        return 2
    n = len(rows)
    ok, msg, pairs = check_hadamard(rows, progress=progress)
    digest = canonical_sha256(rows)
    if ok:
        print(
            "VERDICT: HADAMARD order=%d %s canonical_sha256=%s"
            % (n, msg, digest)
        )
        return 0
    print("VERDICT: FAIL order=%d %s canonical_sha256=%s" % (n, msg, digest))
    return 1


# ---------------------------------------------------------------- selftest


def _sylvester(order):
    """Sylvester Hadamard matrix; order must be a power of 2."""
    h = [[1]]
    while len(h) < order:
        m = len(h)
        nh = [row + row for row in h]
        nh += [row + [-v for v in row] for row in h]
        h = nh
    if len(h) != order:
        raise ValueError("order %d is not a power of 2" % order)
    return h


def selftest():
    """Prove the verifier on knowns and known-bads. Exit 0 iff all pass."""
    failures = []

    def expect(name, rows, want_ok):
        ok, msg, _ = check_hadamard(rows)
        verdict = "ok" if ok == want_ok else "UNEXPECTED"
        print(
            "  selftest %-34s -> %s (%s)"
            % (name, "PASS" if ok else "FAIL", verdict)
        )
        if ok != want_ok:
            failures.append(name)

    for k in range(0, 9):  # orders 1..256
        order = 2 ** k
        expect("sylvester-%d" % order, _sylvester(order), True)

    h = _sylvester(64)

    # Hadamard-preserving operations MUST still pass (verifier not over-strict)
    neg_row = [list(r) for r in h]
    neg_row[5] = [-v for v in neg_row[5]]
    expect("H64 negate row 5", neg_row, True)

    swap = [list(r) for r in h]
    swap[2], swap[7] = swap[7], swap[2]
    expect("H64 swap rows 2,7", swap, True)

    neg_col = [list(r) for r in h]
    for r in neg_col:
        r[3] = -r[3]
    expect("H64 negate col 3", neg_col, True)

    # Corruptions MUST fail
    flip1 = [list(r) for r in h]
    flip1[10][20] = -flip1[10][20]
    expect("H64 flip one entry", flip1, False)

    dup = [list(r) for r in h]
    dup[11] = list(dup[12])
    expect("H64 duplicate row", dup, False)

    allplus = [[1] * 64 for _ in range(64)]
    expect("64x64 all +1", allplus, False)

    expect("order 6 all-cases", [[1] * 6 for _ in range(6)], False)

    # Parse rejections
    for bad_text, why in [
        ("++\n+", "ragged"),
        ("+0\n-+", "bad symbol 0"),
        ("1 2\n1 1", "entry 2"),
        ("", "empty"),
    ]:
        try:
            parse_matrix(bad_text)
        except ValueError:
            print("  selftest parse-reject (%-12s)       -> rejected (ok)" % why)
        else:
            print("  selftest parse-reject (%-12s)       -> ACCEPTED (BAD)" % why)
            failures.append("parse:" + why)

    if failures:
        print("SELFTEST: FAIL (%d): %s" % (len(failures), failures))
        return 1
    print("SELFTEST: PASS")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("matrix", nargs="?", help="matrix file to verify")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument(
        "--progress", action="store_true", help="progress to stderr for big n"
    )
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()
    if not args.matrix:
        ap.print_usage()
        return 2
    return verify_file(args.matrix, progress=args.progress)


if __name__ == "__main__":
    sys.exit(main())
