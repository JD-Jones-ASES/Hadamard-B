#!/usr/bin/env python3
"""cert 02 -- the four twisted i = 2 instances (Lemma T).

Run from the repository root:

    python certs/02-twisted-i2/run.py

`data/twisted-i2-records.json` banks a complete bordered-GS parameter
record, in the same format as `data/payload-records.json`, for the
psi-twisted i = 2 instance at each of the orders 668, 716, 1676, 1772.
The seeds are the psi-twist of the decoded s = 1 record at the same
order; psi is the character of Z_n whose kernel is the index-2 subgroup,
so psi(g) = (-1)^g and kappa(g) = g mod 2, i.e. `coset_divisors = [2]`.

This cert

  1. re-derives the twisted seeds from `data/payload-records.json` and
     checks them character-for-character against the banked record, so
     the bank cannot drift from its stated derivation;
  2. runs `tools/bordered_gs.py::check_record` -- every hypothesis of
     Theorem A/B and every derived law of Theorem C, in exact stdlib
     integer arithmetic -- and additionally asserts, here, the Lemma T
     signature: i = 2, w = n/2, the two-tier profile with the +4 tier
     actually populated (K is a proper subgroup), the doubling
     Q[2I+1] = -Q[2I] and the antisymmetry P[r][2J+1] = -P[r][2J] that
     Theorem D forces;
  3. assembles, hands the file to `verify/verify.py`, compares the
     canonical SHA-256 against the digest pinned in PINNED below (coded
     comparison, hard failure on mismatch), and deletes the matrix;
  4. asserts that each twisted digest DIFFERS from the digest of the
     decoded record at the same order -- i.e. that these are distinct
     artifacts, which is a statement about the files and NOT a claim of
     Hadamard inequivalence.

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

TWISTED = os.path.join("data", "twisted-i2-records.json")
PAYLOAD = os.path.join("data", "payload-records.json")
VERIFY = os.path.join("verify", "verify.py")
TOOLS = "tools"

# canonical SHA-256 of each assembled twisted i = 2 matrix
PINNED = {
    668:  "600849b038e7588d15b64e02794517dd6d95e5d62d505aa7d3d77078392008a3",
    716:  "6b20c6f63875b78adbb1221fda935cb3718918df8b4c779d5763e2e5052f18a7",
    1676: "6a4938371ddbe4ad8bd35f21d7e61dad683b15f8f2ec1c88e88ce579c4907405",
    1772: "82484769a28ac93201f208ca3256bfd491f8edc8bb5e0309764c4f609a113378",
}

# canonical SHA-256 of the DECODED record at the same order (cert 01), used
# only to assert that the twisted artifact is a different file.
DECODED = {
    668:  "bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0",
    716:  "3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6",
    1676: "8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99",
    1772: "1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2",
}

CH = {1: "+", -1: "-"}


def die(msg):
    print("FAIL: %s" % msg)
    sys.exit(1)


def main():
    for p in (TWISTED, PAYLOAD, VERIFY, os.path.join(TOOLS,
                                                     "bordered_gs.py")):
        if not os.path.isfile(p):
            die("missing %s -- run this from the repository root: "
                "python certs/02-twisted-i2/run.py" % p)
    sys.path.insert(0, TOOLS)
    import bordered_gs as T

    with open(TWISTED, "r", encoding="ascii") as fh:
        recs = json.load(fh)["orders"]
    with open(PAYLOAD, "r", encoding="ascii") as fh:
        payload = {int(r["order"]): r for r in json.load(fh)["orders"]}
    if sorted(int(r["order"]) for r in recs) != sorted(PINNED):
        die("banked orders != pinned orders")

    tmp = tempfile.mkdtemp(prefix="cert02-")
    bad = []
    t_all = time.time()
    try:
        for rec in recs:
            order = int(rec["order"])
            t0 = time.time()

            # --- 1. the twist is re-derived, not trusted ---
            src = payload[order]
            if list(src["group"]) != list(rec["group"]) or int(src["s"]) != 1:
                bad.append("%d: banked record is not the s=1 sibling of the "
                           "decoded record" % order)
                continue
            n = rec["group"][0]
            if len(rec["group"]) != 1 or n % 2:
                bad.append("%d: expected a cyclic group of even order" % order)
                continue
            psi = [1 if g % 2 == 0 else -1 for g in range(n)]
            base = [T.signs(x) for x in src["seeds"]]
            want = ["".join(CH[psi[g] * base[q][g]] for g in range(n))
                    for q in range(4)]
            if want != list(rec["seeds"]):
                bad.append("%d: banked seeds are NOT the psi-twist of the "
                           "decoded seeds" % order)
                continue
            if list(rec["coset_divisors"]) != [2]:
                bad.append("%d: coset_divisors != [2]" % order)
                continue
            if list(rec["r_shift"]) != list(src["r_shift"]):
                bad.append("%d: reflection shift differs from the decoded "
                           "record" % order)
                continue

            # --- 2. the theorem, exactly ---
            rep, rows = T.check_record(rec)
            if not rows or not rep["hypotheses_ok"]:
                bad.append("%d: hypotheses failed: %s"
                           % (order, rep.get("failures")))
                continue
            if rep["i"] != 2 or rep["w"] != n // 2:
                bad.append("%d: expected i=2, w=n/2; got i=%d w=%d"
                           % (order, rep["i"], rep["w"]))
            # the +4 tier is actually populated: K is proper, so some g has
            # kappa(g) != 0.  (check_record verified the profile; this
            # asserts the instance is genuinely i = 2 and not i = 1 in
            # disguise.)
            kappa, i = T.coset_map(T.AbelianGroup(rec["group"]),
                                   list(rec["coset_divisors"]))
            if i != 2 or sum(1 for c in kappa if c) != n // 2:
                bad.append("%d: kappa is not the index-2 quotient" % order)
            # Theorem D's forced border shape
            colT = [T.signs(r) for r in rec["col_table"]]
            Q = [[colT[r][k] for r in range(4)] for k in range(8)]
            P = [T.signs(r) for r in rec["row_table"]]
            if any(Q[2 * I + 1][c] != -Q[2 * I][c]
                   for I in range(4) for c in range(4)):
                bad.append("%d: (D-b) Q[2I+1] = -Q[2I] fails" % order)
            if any(P[r][2 * J + 1] != -P[r][2 * J]
                   for r in range(4) for J in range(4)):
                bad.append("%d: (D-c) P[r][2J+1] = -P[r][2J] fails" % order)
            if not (rep["E_is_Hadamard"] and rep["PPt_eq_4i_I"]):
                bad.append("%d: D3 forcing failed under w > 2s" % order)

            # --- 3. trust chain + pinned digest ---
            path = os.path.join(tmp, "H%d_i2.txt" % order)
            with open(path, "w", encoding="ascii", newline="\n") as fh:
                fh.write("\n".join(rows) + "\n")
            proc = subprocess.run([sys.executable, VERIFY, path],
                                  capture_output=True, text=True)
            line = ((proc.stdout or proc.stderr).strip().splitlines()
                    or ["(no output)"])[-1]
            os.remove(path)
            if proc.returncode != 0 or "canonical_sha256=" not in line:
                bad.append("%d: verify.py rc=%d: %s"
                           % (order, proc.returncode, line))
                continue
            got = line.rsplit("canonical_sha256=", 1)[1].strip()
            if got != PINNED[order]:
                bad.append("%d: DIGEST MISMATCH got=%s pinned=%s"
                           % (order, got, PINNED[order]))
                continue

            # --- 4. distinct artifact from the decoded record ---
            if got == DECODED[order]:
                bad.append("%d: the twisted artifact is byte-identical to "
                           "the decoded one" % order)

            print("  order %-5d twist OK  i=2 w=%-4d D-b/D-c border OK  "
                  "digest OK  distinct from decoded  (%.1fs)"
                  % (order, rep["w"], time.time() - t0))
            sys.stdout.flush()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if bad:
        for b in bad:
            print("FAIL: %s" % b)
        print("CERT 02: FAIL (%d problems)" % len(bad))
        return 1
    print("CERT 02: PASS -- 4/4 twisted i=2 instances re-derived from the "
          "decoded seeds, satisfying every hypothesis of Theorem A/B, the "
          "forced border shape of Theorem D, and their pinned digests "
          "(%.1fs)" % (time.time() - t_all))
    return 0


if __name__ == "__main__":
    sys.exit(main())
