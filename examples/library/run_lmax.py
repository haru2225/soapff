"""The mineral library (22 structures, as run.py (a)) rebuilt at l_max 3, 4 and 5, and its one-to-one check by charge and by type (the nearest
other row of the same element has the same charge / the same type).  The number of rows changes with l_max because the slab rows are pruned at
similarity 0.999.

    python run_lmax.py     (with a Python where soapff is installed; a few minutes)  ->  lmax.json, mineral_l<l_max>.npz (not kept)"""
import json
import os

import numpy as np

from soapff_cui import core
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = HERE
BULK = ["brucite_bulk", "kaolinite_bulk", "lizerdite_bulk", "muscovite_1M_bulk", "muscovite_2M1_bulk", "nanpingite_bulk",
        "phlogopite_1M_bulk", "pyrophyllite_2_bulk", "pyrophyllite_bulk", "talc_bulk", "montmorillonite_dry_Na"]
BASAL = ["brucite", "kaolinite", "lizerdite", "muscovite_1M", "muscovite_2M1", "nanpingite", "phlogopite_1M", "pyrophyllite_2",
         "pyrophyllite", "talc", "montmorillonite_dry_Na"]
res = {}
for l in (3, 4, 5):
    refs = [f"{LIB}/refs_bulk/{b}.data" for b in BULK] + [f"{LIB}/refs_basal/{m}_001.data" for m in BASAL]
    out = f"{HERE}/mineral_l{l}.npz"
    core.build_library(refs, out=out, water_in_soap=False, prune=0.999, protect=len(BULK), soap=dict(l_max=l), log=lambda *a: None)
    d = np.load(out, allow_pickle=False)
    X, q, el, ty = d["X"].astype(np.float64), d["q"], d["element"], d["type"]
    wq = wt = 0
    for e in sorted(set(el)):
        m = np.nonzero(el == e)[0]
        S = X[m] @ X[m].T
        np.fill_diagonal(S, -np.inf)
        j = m[S.argmax(1)]
        wq += int((np.abs(q[j] - q[m]) > 1e-6).sum()); wt += int((ty[j] != ty[m]).sum())
    res[l] = dict(rows=int(len(q)), dim=int(X.shape[1]), wrong_charge=wq, wrong_type=wt)
    print(l, res[l], flush=True)
json.dump(res, open(f"{HERE}/lmax.json", "w"), indent=1)
