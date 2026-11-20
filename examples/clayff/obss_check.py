"""A type that clayff has no reference of (paper: the obss example).  In the dry Na montmorillonite reference, one Si next to
an obos O (an O next to an octahedral Mg) is changed to Al.  That O then touches both a tetrahedral and an octahedral
substitution, which ClayFF types obss; no clayff reference has an obss.  soapff finds the obts rows nearest, but the bonded
elements of the O (Al2 Mg1) are those of no obts reference atom, so the build stops with a TypingError and names the atoms.

    python obss_check.py WORK          (soapff)  ->  ../results/clayff_obss_check.json"""
import json
import os
import sys

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list

import soapff
from soapff_cui import core

EX = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(os.path.dirname(soapff.__file__), "libraries", "clayff", "bulk", "montmorillonite_dry_Na.data")
OUT = os.path.join(EX, "results", "clayff_obss_check.json")


def main(work):
    at, _, t = core.read_lammps_data(REF)
    t = np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t])
    el = np.array(at.get_chemical_symbols())
    i, j = neighbor_list("ij", at, 2.0)
    o = int(np.nonzero(t == "obos")[0][0])                 # the first obos O
    si = int([b for a, b in zip(i, j) if a == o and el[b] == "Si"][0])
    el[si] = "Al"                                          # its Si -> Al: the O touches both substitutions (obss)
    new = Atoms(el, positions=at.positions, cell=at.cell, pbc=True)
    rule, _ = core.clayff_answer(new)                      # the bond-graph rule of ClayFF: obss
    res = dict(reference=os.path.basename(REF), changed=f"Si{si + 1} -> Al (next to O{o + 1}, obos)",
               rule_type_of_O=str(rule[o]))
    try:
        soapff.build("clayff", new, work, log=lambda *a: None)
        res["stopped"] = False
    except soapff.TypingError as e:
        res["stopped"] = True
        res["message"] = str(e)
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
