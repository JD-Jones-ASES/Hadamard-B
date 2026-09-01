#!/usr/bin/env python3
"""bank_exact.py -- cash in the exact order-2060 4-profiles, all or nothing.

  WHAT THIS IS.  The producer of the exact |T4| 4-profile at order 2060 is
  the source laboratory's experiments/inequiv/exact_profile_big.py, run
  under the written pre-registration REGISTRATION-2060-exact.md.  This
  script converts that producer's output JSONs into the versioned bank
  cert 07 accepts, and REFUSES to write anything unless all four are
  present, mutually consistent, and bound to the matrices this repository
  rebuilds from data/sep2060-records.json.

  WHY FOUR FILES.  The registration fixes the order of runs: --impl blas
  first, on BOTH matrices; --impl bits on both matrices ONLY if the blas
  profiles differ, as the independent-arithmetic confirmation, "never as a
  rubber stamp on a null".  A PROVEN claim therefore requires the bits
  confirmation, so the complete bank is exactly

      data/sep2060-exact-blas-plain.json
      data/sep2060-exact-blas-gist.json
      data/sep2060-exact-bits-plain.json
      data/sep2060-exact-bits-gist.json

  and anything less must not upgrade the certificate.  Cert 07 enforces
  that on its side; this script refuses to create a partial bank on this
  one.

  WHAT THE BANK PROVES, AND WHAT IT DOES NOT.  Every check here is a check
  on the producer's SELF-DECLARED record: that its matrix_sha256 is the
  matrix this repository rebuilds, that its counts satisfy the order-2060
  identities, that two arithmetics agree.  None of that is proof that the
  counts were computed from that matrix.  The proof of computation is the
  registered run itself.  Say so wherever the bank is used.

  REFUSALS (nothing is written, exit 2):
    * any of the four producer files missing, unreadable, or not JSON;
    * a producer file whose n, impl, totals or second moment are wrong;
    * a producer file whose matrix_sha256 is not the canonical digest of
      the matrix rebuilt here for its tag;
    * the two blas profiles identical bin for bin -- registration S6: that
      outcome is "did not separate", the pair returns to UNKNOWN, and
      there is nothing to bank;
    * blas and bits disagreeing in any bin -- registration kill criterion
      5: hard stop, no claim in either direction;
    * an assembled bank that cert 07's own acceptance predicate rejects;
    * an output file that already exists (use --force deliberately).

Usage:
  python certs/07-2060-evidence/bank_exact.py \
      --blas-plain <producer json> --blas-gist <producer json> \
      --bits-plain <producer json> --bits-gist <producer json> \
      [--out-dir data] [--dry-run] [--force]

Standard library only.  Exact integer arithmetic only.  No numpy.
"""

import argparse
import datetime
import json
import os
import sys

sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import run as CERT                                          # noqa: E402

# Producer fields carried into the bank when present (cert 07's
# OPTIONAL_FIELDS; anything else the producer wrote is dropped, because
# the bank schema rejects unrecognised fields).
CARRY = ("seconds", "engine", "enumeration", "folded", "q", "threads",
         "geom", "peak_rss_mb")


class Refusal(Exception):
    pass


def rel(p):
    a = os.path.abspath(p)
    if a == ROOT or a.startswith(ROOT + os.sep):
        return os.path.relpath(a, ROOT).replace("\\", "/")
    return a


def read_producer(path, impl, tag, digest):
    """One exact_profile_big.py output, checked against the rebuilt matrix."""
    where = os.path.basename(path)
    if not os.path.isfile(path):
        raise Refusal("%s: no such file" % path)
    try:
        with open(path, "r", encoding="ascii") as fh:
            blob = json.load(fh)
    except Exception as exc:
        raise Refusal("%s: not readable as ASCII JSON (%s)" % (where, exc))
    if not isinstance(blob, dict):
        raise Refusal("%s: top level is not an object" % where)
    for f in ("n", "impl", "profile", "total", "second_moment",
              "matrix_sha256"):
        if f not in blob:
            raise Refusal("%s: the producer output has no %r field" % (where, f))
    if blob["n"] != CERT.N:
        raise Refusal("%s: n is %r, want %d" % (where, blob["n"], CERT.N))
    if blob["impl"] != impl:
        raise Refusal("%s: declares impl %r, but was passed as --%s-%s"
                      % (where, blob["impl"], impl, tag))
    if blob["matrix_sha256"] != digest:
        raise Refusal(
            "%s: the producer read a matrix with sha256 %s..., but the %s "
            "matrix this repository rebuilds is %s...  This is the binding "
            "that matters; it must not be overridden."
            % (where, str(blob["matrix_sha256"])[:16], tag, digest[:16]))
    try:
        prof = CERT.parse_profile(blob["profile"], where)
    except CERT.BankError as e:
        raise Refusal("%s: %s" % (e.reason, e.detail))
    tot = sum(prof.values())
    m2 = sum(k * k * v for k, v in prof.items())
    if not CERT.is_int(blob["total"]) or blob["total"] != tot:
        raise Refusal("%s: declares total %r, the bins sum to %d"
                      % (where, blob["total"], tot))
    if tot != CERT.c_n_4(CERT.N):
        raise Refusal("%s: total %d, want C(%d,4) = %d"
                      % (where, tot, CERT.N, CERT.c_n_4(CERT.N)))
    if not CERT.is_int(blob["second_moment"]) or blob["second_moment"] != m2:
        raise Refusal("%s: declares second moment %r, the bins give %d"
                      % (where, blob["second_moment"], m2))
    if m2 != CERT.second_moment_want(CERT.N):
        raise Refusal("%s: sum T4^2 = %d, want n^3(n-1)(n-2)/24 = %d"
                      % (where, m2, CERT.second_moment_want(CERT.N)))
    if "C_n_4" in blob and blob["C_n_4"] != CERT.c_n_4(CERT.N):
        raise Refusal("%s: declares C_n_4 = %r" % (where, blob["C_n_4"]))
    if ("second_moment_want" in blob
            and blob["second_moment_want"] != CERT.second_moment_want(CERT.N)):
        raise Refusal("%s: declares second_moment_want = %r"
                      % (where, blob["second_moment_want"]))
    return blob, prof


def make_bank_blob(src_path, blob, prof, impl, tag, digest, stamp):
    out = {"schema": CERT.SCHEMA,
           "n": CERT.N,
           "tag": tag,
           "impl": impl,
           "matrix_canonical_sha256": digest,
           "producer_matrix_sha256": blob["matrix_sha256"],
           "producer": "%s (%s, %s) -- Hadamard-2060 experiments/inequiv, "
                       "REGISTRATION-2060-exact.md, --impl %s"
                       % (blob.get("engine", "exact_profile_big.py"),
                          blob.get("enumeration", "canonical-split"),
                          blob.get("folded", "|T4|"), impl),
           "source_file": os.path.basename(src_path),
           "banked_utc": stamp,
           "profile": dict((str(k), v) for k, v in sorted(prof.items())),
           "total": sum(prof.values()),
           "second_moment": sum(k * k * v for k, v in prof.items())}
    for f in CARRY:
        if f in blob:
            out[f] = blob[f]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="convert the registered exact-profile run into cert 07's "
                    "bank; all four files or nothing")
    for impl in CERT.IMPLS:
        for tag in CERT.TAGS:
            ap.add_argument("--%s-%s" % (impl, tag), required=True,
                            metavar="JSON",
                            help="exact_profile_big.py --impl %s output for "
                                 "2060-%s" % (impl, tag))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "data"),
                    help="default: the repository's data/")
    ap.add_argument("--dry-run", action="store_true",
                    help="validate everything and write nothing")
    ap.add_argument("--force", action="store_true",
                    help="overwrite bank files that already exist")
    args = ap.parse_args(argv)
    src = dict(((impl, tag),
                getattr(args, "%s_%s" % (impl, tag)))
               for impl in CERT.IMPLS for tag in CERT.TAGS)

    print("=" * 74)
    print("bank_exact -- cashing in the exact order-2060 4-profiles")
    print("=" * 74)

    try:
        print("\n[1] rebuild both matrices for their canonical digests")
        with open(os.path.join(ROOT, "data", "sep2060-records.json"),
                  "r", encoding="ascii") as fh:
            recs = json.load(fh)
        plain, gist, _raw = CERT.build_2060(recs)
        digests = {"plain": CERT.rows_sha256(plain),
                   "gist": CERT.rows_sha256(gist)}
        del plain, gist
        if (digests["plain"] != CERT.SHA_PLAIN
                or digests["gist"] != CERT.SHA_GIST):
            raise Refusal("the rebuilt matrices do not reproduce the pins in "
                          "run.py; nothing may be banked against them")
        for tag in CERT.TAGS:
            print("      %-5s %s" % (tag, digests[tag]))

        print("\n[2] the four producer outputs")
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ")
        prof, bank = {}, {}
        for impl in CERT.IMPLS:
            for tag in CERT.TAGS:
                blob, p = read_producer(src[(impl, tag)], impl, tag,
                                        digests[tag])
                prof[(tag, impl)] = p
                bank[(impl, tag)] = make_bank_blob(src[(impl, tag)], blob, p,
                                                   impl, tag, digests[tag],
                                                   stamp)
                print("      %-5s %-4s  %3d bins  bound to %s...  <- %s"
                      % (tag, impl, len(p), digests[tag][:16],
                         os.path.basename(src[(impl, tag)])))

        print("\n[3] the registration's own gates")
        if prof[("plain", "blas")] == prof[("gist", "blas")]:
            raise Refusal(
                "the two blas profiles are identical bin for bin.  "
                "Registration S6: that is 'did not separate', NOT "
                "equivalence; the pair returns to UNKNOWN and the sampled "
                "separation is falsified.  There is nothing to bank -- write "
                "it up instead.")
        for tag in CERT.TAGS:
            if prof[(tag, "blas")] != prof[(tag, "bits")]:
                raise Refusal(
                    "%s: blas and bits disagree.  Registration kill criterion "
                    "5: hard stop, no claim in either direction until it is "
                    "explained." % tag)
        ks = sorted(set(prof[("plain", "blas")]) | set(prof[("gist", "blas")]))
        ndiff = sum(1 for k in ks
                    if prof[("plain", "blas")].get(k, 0)
                    != prof[("gist", "blas")].get(k, 0))
        print("      the blas profiles differ in %d of %d bins" % (ndiff,
                                                                   len(ks)))
        print("      blas == bits, bin for bin, on both matrices")

        print("\n[4] the assembled bank, through cert 07's own predicate")
        paths = dict((os.path.join(args.out_dir,
                                   CERT.exact_basename(impl, tag)),
                      bank[(impl, tag)])
                     for impl in CERT.IMPLS for tag in CERT.TAGS)
        try:
            CERT.accept_bank(sorted(paths), digests, lambda p: paths[p])
        except CERT.BankError as e:
            raise Refusal("cert 07 would reject this bank: %s -- %s"
                          % (e.reason, e.detail))
        print("      accepted by accept_bank(): 4 files, schema %s"
              % CERT.SCHEMA)

        if not args.dry_run:
            for p in sorted(paths):
                if os.path.exists(p) and not args.force:
                    raise Refusal("%s already exists; --force to overwrite"
                                  % rel(p))
    except Refusal as e:
        print("\n" + "=" * 74)
        print("REFUSED: %s" % e)
        print("Nothing was written.  Exact banking is all-or-nothing.")
        print("=" * 74)
        return 2

    print("\n[5] write")
    if args.dry_run:
        print("      --dry-run: nothing written.  All gates passed.")
        print("=" * 74)
        return 0
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    # All four are written to sibling .tmp files first and only then moved
    # into place, so a failure part-way through the writes -- a full disk,
    # an interrupt -- leaves the four bank names as they were rather than a
    # partial bank the docstring above promises cannot exist.
    tmps = dict((p, p + ".tmp") for p in sorted(paths))
    try:
        for p in sorted(paths):
            with open(tmps[p], "w", encoding="ascii", newline="\n") as fh:
                json.dump(paths[p], fh, indent=1)
                fh.write("\n")
    except BaseException:
        for t in tmps.values():
            if os.path.exists(t):
                os.remove(t)
        raise
    for p in sorted(paths):
        os.replace(tmps[p], p)
        print("      %s" % rel(p))
    print("\n[6] pin these digests in run.py EXACT_FILE_PINS, then re-run "
          "the cert")
    for p in sorted(paths):
        print("      %-40s %s" % (os.path.basename(p), CERT.file_sha256(p)))
    print("\n      wrote 4 files to %s" % rel(args.out_dir))
    print("      next: python -B certs/07-2060-evidence/run.py")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
