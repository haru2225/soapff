"""Why the depth is 4 (paper: "the SOAP of the atoms within 4 bonds"): the test structures of test_typing.py typed with
oplsaa (depth 4) and with a library that extends oplsaa and keeps the atoms within 3 bonds only (depth 3), against Foyer.
With depth 3 the H of an ester O-CH3 (opls_469) and the H on a C bonded to an ether O (opls_185) are exchanged: the
carbonyl O that tells them apart is 4 bonds from the H.  Both have +0.03 e, so the charges stay right and soapff does not stop.

    python depth.py TESTDIR WORK results/depth.json      (soapff)

For each depth: the wrong types and charges in the test structures (TESTDIR/*.extxyz, written with strict=False), per
structure, the pairs of types exchanged (Foyer's -> soapff's), and the structures on which soapff.build would stop."""
import collections
import glob
import json
import os
import sys
import tempfile

import ase.io
import numpy as np

import soapff

T, W, OUT = sys.argv[1:4]
W = os.path.abspath(W)
os.makedirs(W, exist_ok=True)
STOPS = ("atoms are not in the library", "molecules have a charge that is not a whole number")   # the strict checks
tests = sorted(glob.glob(os.path.join(T, "*.extxyz")))
res = {}
for depth in (4, 3):
    if depth == 4:
        lib = "oplsaa"
    else:                                        # the [soap] of oplsaa with another depth
        soap = dict(soapff.load_library("oplsaa").soap, depth=depth)
        lib = os.path.join(W, f"oplsaa_depth{depth}.txt")
        open(lib, "w").write(f"[library]\nname = oplsaa_depth{depth}\nextends = oplsaa\ndescription = oplsaa with depth "
                             f"{depth}\n\n[soap]\n" + "".join(f"{k} = {v}\n" for k, v in soap.items()))
    wrong_t = wrong_q = n_atoms = 0
    per, pairs, stops, smin = {}, collections.Counter(), {}, 1.0
    for f in tests:
        at = ase.io.read(f)
        key = os.path.basename(f)[:-7]
        with tempfile.TemporaryDirectory(dir=W) as d:
            s = soapff.build(lib, f, d, strict=False, log=lambda *a: None)
            types, q = np.array(s.types), np.array(s.charges, float)
            why = [w for w in s.report.get("warnings", []) if any(x in w for x in STOPS)]
            if why:
                stops[key] = why[0][:300]
            smin = min(smin, min(float(r["similarity_soap"]) for r in s.table() if r["similarity_soap"]))
        ft, fq = np.array(at.arrays["fftype"]), np.array(at.arrays["ffcharge"], float)
        wt, wq = int((types != ft).sum()), int((np.abs(q - fq) > 1e-6).sum())
        pairs.update(f"{a}->{b}" for a, b in zip(ft[types != ft], types[types != ft]))
        wrong_t, wrong_q, n_atoms = wrong_t + wt, wrong_q + wq, n_atoms + len(at)
        if wt or wq:
            per[key] = [wt, wq]
    res[f"depth {depth}"] = dict(structures=len(tests), atoms=n_atoms, wrong_types=wrong_t, wrong_charges=wrong_q,
                                 structures_with_errors=per, types_exchanged=dict(pairs), structures_that_stop=stops,
                                 min_similarity=round(smin, 4))
    print(f"depth {depth}", json.dumps(res[f"depth {depth}"]), flush=True)
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=1)
