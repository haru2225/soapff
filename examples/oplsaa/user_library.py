"""A user library: oplsaa with one more reference, the junction of a PE-b-PEO block copolymer, which is in no reference of
oplsaa (make_refs.py writes extra/PE4_PEO3.extxyz, a short PE-b-PEO typed by Foyer).

    python user_library.py extra/block_PE_PEO.extxyz extra/PE4_PEO3.extxyz WORK results/user_library.json

The block copolymer (PE 20 C - PEO 10 units, extra/) is typed with the user library and compared with the types and charges
of Foyer.  The library file is written in WORK (soapff.new_library(..., folder=WORK)), not in ~/.soapff."""
import json
import os
import sys

import ase.io
import numpy as np

import soapff

TEST, EXTRA, W, OUT = sys.argv[1:5]
os.makedirs(W, exist_ok=True)
lib = os.path.join(W, "pe_peo.txt")
if not os.path.isfile(lib):
    soapff.new_library("pe_peo", extends="oplsaa", folder=W)
    soapff.add_reference(lib, EXTRA)
ref = ase.io.read(TEST)
key = "oplsaa + PE4_PEO3"
s = soapff.build(lib, TEST, os.path.join(W, "pe_peo"), log=lambda *a: None)
tt = np.array(s.types)
rows = s.table()
wrong = tt != ref.arrays["fftype"]
res = {key: dict(atoms=len(ref), wrong_type=int(wrong.sum()),
                 wrong_charge=int((np.abs(s.charges - ref.arrays["ffcharge"]) > 1e-6).sum()),
                 net_charge=round(float(s.charges.sum()), 6), min_similarity_soap=round(s.report["min_similarity_soap"], 4),
                 wrong_atoms=[dict(atom=int(i) + 1, element=ref.get_chemical_symbols()[i], foyer=str(ref.arrays["fftype"][i]),
                                   soapff=str(tt[i]), similarity=float(rows[i]["similarity_soap"]),
                                   nearest=rows[i]["nearest_reference"]) for i in np.nonzero(wrong)[0]])}
r = res[key]
print(f"{key:18s} wrong type {r['wrong_type']}  wrong charge {r['wrong_charge']}  net {r['net_charge']:+.4f}  "
      f"min similarity {r['min_similarity_soap']:.3f}  {r['wrong_atoms']}")
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=1)
