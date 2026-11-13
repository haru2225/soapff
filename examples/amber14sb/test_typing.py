"""soapff (library amber14sb) on the test structures, which are not references of the library, against the types and charges
OpenMM gives them (the columns fftype / ffcharge of test/*.extxyz, written by make_refs.py).

    python test_typing.py test WORKDIR results/typing.json"""
import glob
import json
import os
import sys
import time

import ase.io
import numpy as np

import soapff

T, W, OUT = sys.argv[1:4]
quiet = lambda *a: None
t0 = time.time()
lib = soapff.load_library("amber14sb")
lib.table_path(log=quiet)
res = dict(library=repr(lib), table_seconds=round(time.time() - t0, 1), structures={})
for f in sorted(glob.glob(os.path.join(T, "*.extxyz"))):
    name = os.path.basename(f)[:-7]
    ref = ase.io.read(f)
    t1 = time.time()
    s = soapff.build("amber14sb", f, os.path.join(W, name), md=False, log=quiet)
    tt, qq = np.array(s.types), s.charges
    r = dict(atoms=len(ref), wrong_type=int((tt != ref.arrays["fftype"]).sum()),
             wrong_charge=int((np.abs(qq - ref.arrays["ffcharge"]) > 1e-5).sum()),
             min_similarity_soap=round(s.report["min_similarity_soap"], 4),
             min_similarity_fragment=round(s.report["min_similarity_fragment"], 5), terms=s.report["terms"],
             seconds=round(time.time() - t1, 2))
    res["structures"][name] = r
    print(f"{name:22s} {r['atoms']:5d} atoms  wrong type {r['wrong_type']}  wrong charge {r['wrong_charge']}  "
          f"min similarity soap {r['min_similarity_soap']:.3f} fragment {r['min_similarity_fragment']:.4f}  {r['seconds']:.1f} s",
          flush=True)
json.dump(res, open(OUT, "w"), indent=1)
