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

The public API is `check_hadamard(rows)`.  It validates its own input --
nonempty, square, every entry the Python integer +1 or -1 -- and does not
rely on the file parser having done so.  Bools are refused explicitly
(`bool` is a subclass of `int`, so `True` would otherwise pack as +1), and
so is anything that is not an `int`: floats, strings, and non-`int`
integer types such as numpy scalars.  Callers holding such values must
convert them to plain `int` before calling.

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


def validate_rows(rows):
    """Return None if `rows` is a nonempty square matrix of the integers
    +1 and -1; otherwise a string naming the first violation.

    This is the input contract of `check_hadamard`.  It duplicates, for
    the direct-API caller, the guarantees `parse_matrix` already gives the
    CLI: without it, any value other than the integer 1 packs silently as
    -1, and `check_hadamard([[1, 1], [1, 0]])` returns True.

    Bools are rejected as a type error rather than coerced: `bool` is a
    subclass of `int`, so a value test alone would accept `True` as +1.
    Non-`int` types are rejected for the same reason -- `-1.0 == -1` is
    True in Python, and this file admits no floats.
    """
    if not isinstance(rows, (list, tuple)):
        return "matrix is not a list of rows"
    n = len(rows)
    if n == 0:
        return "empty matrix: no rows"
    for i, row in enumerate(rows):
        if not isinstance(row, (list, tuple)):
            return "row %d is not a list of entries" % i
        if len(row) != n:
            return (
                "row %d has length %d, expected %d (square matrix required)"
                % (i, len(row), n)
            )
        for j, v in enumerate(row):
            if isinstance(v, bool):
                return (
                    "row %d entry %d is %r: bool is not an accepted entry "
                    "type (entries must be the integers +1 or -1)" % (i, j, v)
                )
            if not isinstance(v, int):
                return (
                    "row %d entry %d is %r of type %s: entries must be the "
                    "integers +1 or -1" % (i, j, v, type(v).__name__)
                )
            if v != 1 and v != -1:
                return (
                    "row %d entry %d is %r: entries must be +1 or -1"
                    % (i, j, v)
                )
    return None


def check_hadamard(rows, progress=False):
    """Return (ok, message, pairs_checked). Exact arithmetic throughout.

    The input contract is enforced here, not assumed: see `validate_rows`.
    """
    bad_input = validate_rows(rows)
    if bad_input is not None:
        return (False, bad_input, 0)
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

    # Direct-API input contract.  These call check_hadamard(rows) with no
    # parser in front of it: the CLI's guarantees do not cover this path,
    # so the API must enforce them itself.  Without validate_rows, a value
    # test alone accepts every case below.
    expect("direct-API H2 accepted", [[1, 1], [1, -1]], True)
    expect("direct-API H4 accepted", _sylvester(4), True)
    expect("direct-API -H4 accepted",
           [[-v for v in r] for r in _sylvester(4)], True)
    expect("direct-API zero entry", [[1, 1], [1, 0]], False)
    expect("direct-API entry 2", [[1, 1], [1, 2]], False)
    expect("direct-API float entry -1.0", [[1, 1], [1, -1.0]], False)
    expect("direct-API all-float matrix", [[1.0, 1.0], [1.0, -1.0]], False)
    expect("direct-API bool True as +1", [[1, 1], [1, True]], False)
    expect("direct-API all-bool matrix",
           [[True, True], [True, False]], False)
    expect("direct-API str entry", [["+", "+"], ["+", "-"]], False)
    expect("direct-API non-square 2x3", [[1, 1, 1], [1, -1, 1]], False)
    expect("direct-API ragged rows", [[1, -1], [1]], False)
    expect("direct-API empty", [], False)
    expect("direct-API row not a list", [1, -1], False)

    # A rejection must be a rejection for the stated reason, not an
    # accidental non-orthogonality verdict.
    for name, rows_, needle in [
        ("zero entry", [[1, 1], [1, 0]], "entries must be +1 or -1"),
        ("float entry", [[1, 1], [1, -1.0]], "of type float"),
        ("bool entry", [[1, 1], [1, True]], "bool is not an accepted"),
        ("non-square", [[1, 1, 1], [1, -1, 1]], "square matrix required"),
        ("empty", [], "empty matrix"),
    ]:
        _, msg, _ = check_hadamard(rows_)
        good = needle in msg
        print("  selftest reason (%-11s)             -> %s"
              % (name, "ok" if good else "WRONG REASON: " + msg))
        if not good:
            failures.append("reason:" + name)

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
