"""soapff (library oplsaa, the SOAP of each atom only) on the test structures, which are not references of the library,
against the types and charges Foyer gives them (the columns fftype / ffcharge of test/*.extxyz, written by make_refs.py).

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
quiet = lambda *a: None                                                          # noqa: E731
t0 = time.time()
lib = soapff.load_library("oplsaa")
lib.table_path(log=quiet)
res = dict(library=f"{lib.name} ({len(lib.references)} references)", table_seconds=round(time.time() - t0, 1), structures={})
for f in sorted(glob.glob(os.path.join(T, "*.extxyz"))):
    name = os.path.basename(f)[:-7]
    ref = ase.io.read(f)
    t1 = time.time()
    try:
        s = soapff.build("oplsaa", f, os.path.join(W, name), log=quiet)
    except soapff.TypingError as e:                  # stopped (e.g. a bond stretched beyond the cutoff in a hot frame)
        res["structures"][name] = dict(atoms=len(ref), stopped=str(e)[:600])
        print(f"{name:16s} {len(ref):5d} atoms  stopped: {str(e)[:300]}", flush=True)
        continue
    tt, qq = np.array(s.types), s.charges
    wt = tt != ref.arrays["fftype"]
    r = dict(atoms=len(ref), wrong_type=int(wt.sum()), wrong_charge=int((np.abs(qq - ref.arrays["ffcharge"]) > 1e-6).sum()),
             types=len(set(tt)), min_similarity_soap=round(s.report["min_similarity_soap"], 4),
             net_charge=round(float(qq.sum()), 6), terms=s.report["terms"], seconds=round(time.time() - t1, 2))
    if wt.any():
        r["confused"] = sorted({f"{a}->{b}" for a, b in zip(ref.arrays["fftype"][wt], tt[wt])})
    res["structures"][name] = r
    print(f"{name:16s} {r['atoms']:5d} atoms  wrong type {r['wrong_type']}  wrong charge {r['wrong_charge']}  "
          f"types {r['types']:2d}  min similarity {r['min_similarity_soap']:.3f}  net {r['net_charge']:+.4f}  "
          f"{r['seconds']:.1f} s", flush=True)
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=1)
