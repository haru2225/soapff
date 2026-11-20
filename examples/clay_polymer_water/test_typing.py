"""The typing of the clay + polymer + water tests with clayff_oplsaa: the start structures (structures/*_start.extxyz), from
the elements, positions and cell alone, compared with the right answer of build.py.

    python test_typing.py      (soapff)  ->  results/typing.json"""
import json
import os
import tempfile

import numpy as np
from ase.io import read

import soapff

HERE = os.path.dirname(os.path.abspath(__file__))
CLAYFF = set(soapff.library.load("clayff_oplsaa").types) - set(soapff.library.load("oplsaa").types)   # mineral, water, ions


def score(at, s, truth_t, truth_q):
    t = np.array(s.types)
    q = s.charges
    sim = np.array([float(r["similarity_soap"]) for r in s.table()])
    wrong = (t != truth_t) | (np.abs(q - truth_q) > 1e-4)
    clay = np.isin(truth_t, list(CLAYFF))
    rep = s.report
    return dict(atoms=len(at), wrong=int(wrong.sum()),
                wrong_atoms=[dict(id=int(i + 1), right=str(truth_t[i]), got=str(t[i])) for i in np.nonzero(wrong)[0][:10]],
                atoms_clayff=int(clay.sum()), atoms_oplsaa=int((~clay).sum()),
                min_similarity=dict(clayff=float(sim[clay].min()), oplsaa=float(sim[~clay].min())),
                terms=rep["terms"], bonds_not_in_force_field=rep.get("bonds_not_in_force_field"),
                frameworks=rep.get("frameworks"), warnings=rep["warnings"])


def main():
    res = {}
    for name in ("PEO", "PVA"):
        start = read(os.path.join(HERE, "structures", f"{name}_start.extxyz"))
        T, Q = start.arrays["fftype"], start.arrays["ffcharge"]
        with tempfile.TemporaryDirectory() as d:
            s = soapff.build("clayff_oplsaa", start, os.path.join(d, "start"), log=lambda *a: None)
            res[name] = score(start, s, T, Q)
        v = res[name]
        print(name, v["atoms"], "wrong", v["wrong"], "min similarity",
              {a: round(b, 3) for a, b in v["min_similarity"].items()}, v["terms"])
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, "results", "typing.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
