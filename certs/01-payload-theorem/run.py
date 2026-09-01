#!/usr/bin/env python3
"""cert 01 -- the twelve public records against the theorem.

Run from the repository root:

    python certs/01-payload-theorem/run.py

For each of the twelve parameter records banked in
`data/payload-records.json` this cert

  1. runs `tools/bordered_gs.py::check_record`, which checks the
     hypotheses (H0/H1/H2/H3/H4) of Theorem A in its house form
     (Theorem B) and the derived laws of Theorem C (D1, D3, D4, D5, the
     Sigma-bar law) in exact stdlib integer arithmetic;
  2. asserts, here, that the record lands on one of the three cells of
     the Theorem C classification, and that the twelve cells match the
     census stated in NOTE-B.md 2.1;
  3. assembles the matrix and hands the file to `verify/verify.py`, the
     repository's trust chain;
  4. compares the canonical SHA-256 the trust chain reports against the
     digest pinned in PINNED below (coded comparison, hard failure on
     mismatch);
  5. deletes the generated matrix.

Stdlib only.  Exact integers only.  No network.  Nothing is written
inside the repository.

Exit code 0 iff every record passed every check and every digest matched.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

DATA = os.path.join("data", "payload-records.json")
VERIFY = os.path.join("verify", "verify.py")
TOOLS = "tools"

# The pinned canonical SHA-256 of each assembled matrix.  Provenance of
# these twelve digests is in NOTES.md; they are the digests the lab's
# theorem-check report and its independent replay report both carry.
PINNED = {
    668:  "bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0",
    716:  "3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6",
    892:  "e77fc79ab287f5f5ba5bbdc10191bdc7593839052fe1015c1fb6a2e974ab54de",
    1132: "7d1c1e892149e90330d58bb0cf9ef2c888078df1b35fb55f8724d580ebf7b743",
    1244: "4cb747cf511eba1f203582b5121bdf6ab02671133e45579c1d023add8b2da143",
    1388: "a6b92584eb803b87026709d64fe892dec8f7182a120e13de9edd3065cf05bf0b",
    1436: "e4d745a4d44f39a5671f9cd86f5c1d0aef93504dcfb2e253451cadc9e4086728",
    1676: "8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99",
    1772: "1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2",
    1916: "be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9",
    1948: "fddc841ebf951f6e17e939551d058ea5df046251ea065b5f6e7ee2fd8d0f62ce",
    1964: "740b907cd442f1b7fd40dcc31f2b3aae9794842da6dc579f98dac1d0d9e1493d",
}

# The Theorem C classification census stated in NOTE-B.md 2.1:
# (s, i) -> how many of the twelve records must land in that cell.
CENSUS = {(0, 1): 5, (1, 1): 4, (3, 4): 1, (5, 6): 1, (7, 8): 1}


def die(msg):
    print("FAIL: %s" % msg)
    sys.exit(1)


def classification_cell_ok(s, i):
    """The three cells that survive Theorem C under w > 2s, n >= 3."""
    if (s, i) in ((0, 1), (1, 1)):
        return True
    return i == s + 1 and s % 2 == 1


def main():
    for p in (DATA, VERIFY, os.path.join(TOOLS, "bordered_gs.py")):
        if not os.path.isfile(p):
            die("missing %s -- run this from the repository root: "
                "python certs/01-payload-theorem/run.py" % p)
    sys.path.insert(0, TOOLS)
    import bordered_gs as T

    with open(DATA, "r", encoding="ascii") as fh:
        bank = json.load(fh)
    recs = bank["orders"]
    orders = [int(r["order"]) for r in recs]
    if sorted(orders) != sorted(PINNED):
        die("banked orders %s != pinned orders %s"
            % (sorted(orders), sorted(PINNED)))
    if len(recs) != 12:
        die("expected 12 records, found %d" % len(recs))

    tmp = tempfile.mkdtemp(prefix="cert01-")
    cells = {}
    bad = []
    t_all = time.time()
    try:
        for rec in recs:
            order = int(rec["order"])
            t0 = time.time()
            rep, rows = T.check_record(rec)
            if not rows:
                bad.append("%d: record refused: %s" % (order, rep["failures"]))
                continue
            if not rep["hypotheses_ok"]:
                bad.append("%d: hypotheses failed: %s"
                           % (order, rep["failures"]))
            s, i = rep["s"], rep["i"]
            cells[(s, i)] = cells.get((s, i), 0) + 1
            if not classification_cell_ok(s, i):
                bad.append("%d: cell (s,i)=(%d,%d) is outside the Theorem C "
                           "classification" % (order, s, i))
            # Derived laws of Theorem C, re-asserted here rather than only
            # read off the checker's report.  (For s = 0 the border is
            # empty and the border-side clauses are vacuous; the checker
            # does not populate them, so they are asserted arithmetically.)
            if not i <= s + 1:
                bad.append("%d: D1 (i <= s+1) failed" % order)
            if not rep["D5_rowsum_law"]:
                bad.append("%d: D5 row-sum law failed: %s"
                           % (order, rep["D5_value"]))
            if s >= 1 and i == s + 1 and s % 2 == 0:
                bad.append("%d: D2 (i = s+1 => s odd) failed" % order)
            if s > 0:
                if rep["w"] > 2 * s and not (rep["E_is_Hadamard"]
                                             and rep["PPt_eq_4i_I"]):
                    bad.append("%d: D3 forcing failed under w > 2s" % order)
                if not rep["SigmaBar_law"]:
                    bad.append("%d: Sigma-bar law failed" % order)
                if not rep["compression_lemma_crosscheck"]:
                    bad.append("%d: compression-lemma cross-check failed"
                               % order)

            path = os.path.join(tmp, "H%d.txt" % order)
            with open(path, "w", encoding="ascii", newline="\n") as fh:
                fh.write("\n".join(rows) + "\n")
            proc = subprocess.run([sys.executable, VERIFY, path],
                                  capture_output=True, text=True)
            line = ((proc.stdout or proc.stderr).strip().splitlines()
                    or ["(no output)"])[-1]
            os.remove(path)
            if proc.returncode != 0:
                bad.append("%d: verify.py rc=%d: %s"
                           % (order, proc.returncode, line))
                continue
            if "canonical_sha256=" not in line:
                bad.append("%d: no canonical_sha256 in verdict: %s"
                           % (order, line))
                continue
            got = line.rsplit("canonical_sha256=", 1)[1].strip()
            if got != PINNED[order]:
                bad.append("%d: DIGEST MISMATCH got=%s pinned=%s"
                           % (order, got, PINNED[order]))
                continue
            print("  order %-5d (s,i)=(%d,%d) w=%-4d hypotheses OK  "
                  "digest OK  (%.1fs)"
                  % (order, s, i, rep["w"], time.time() - t0))
            sys.stdout.flush()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if cells != CENSUS:
        bad.append("classification census %s != stated census %s"
                   % (sorted(cells.items()), sorted(CENSUS.items())))

    print("  cells: %s" % (sorted(cells.items()),))
    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 01: FAIL (%d problems)" % len(bad))
        return 1
    print("CERT 01: PASS -- 12/12 records satisfy every hypothesis of "
          "Theorem A/B and every derived law of Theorem C; 12/12 assembled "
          "matrices verified Hadamard and matched their pinned digests "
          "(%.1fs)" % (time.time() - t_all))
    return 0


if __name__ == "__main__":
    sys.exit(main())
