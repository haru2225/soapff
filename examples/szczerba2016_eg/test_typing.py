"""The typing of the SWy-1 + EG + water models (build.py) with the bundled clayff_oplsaa, from the elements, positions and
cell alone, compared with the right answer of build.py.  EG is not among the references; its atoms get the types of the
same atoms at the HO-CH2- end of the PVA references.

    python test_typing.py      (soapff)  ->  results/typing_start.json"""
import json
import os
import tempfile

import numpy as np
from ase.io import read

import soapff

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = "clayff_oplsaa"
OPLS = set(soapff.library.load("oplsaa").types)


def main():
    res = {}
    for eg in (1.7, 1.5):
        at = read(os.path.join(HERE, "structures", f"SWy-1_EG{eg}_start.extxyz"))
        with tempfile.TemporaryDirectory() as d:
            s = soapff.build(LIB, at, os.path.join(d, "ff"), log=lambda *a: None)
            t, q = np.array(s.types), s.charges
            sim = np.array([float(r["similarity_soap"]) for r in s.table()])
        truth = at.arrays["fftype"]
        wrong = (t != truth) | (np.abs(q - at.arrays["ffcharge"]) > 1e-4)
        org = np.isin(truth, list(OPLS))
        res[f"EG{eg}_start"] = dict(atoms=len(at), wrong=int(wrong.sum()),
                                    min_similarity=dict(clayff=float(sim[~org].min()), oplsaa=float(sim[org].min())),
                                    terms=s.report["terms"], frameworks=s.report.get("frameworks"))
        print(f"EG {eg}: {res[f'EG{eg}_start']}")
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, "results", "typing_start.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
