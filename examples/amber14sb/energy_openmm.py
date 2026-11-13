"""OpenMM energies (kcal/mol) of the test structures: ff14SB + TIP3P, no cutoff, no constraints, every force separately.
The topology is built again from the PDB file as make_refs.py does; the coordinates are those of the extxyz file.
    python energy_openmm.py TESTDIR OUT.json"""
import json, os, sys
import numpy as np
from omm_common import FF, atom_lines, fragment, prepare, system_vacuum, energies, typed
T = sys.argv[1]
cases = {"1UBQ": ("1UBQ", None), "1UBQ_md300K_5ps": ("1UBQ", None), "1UBQ_md300K_15ps": ("1UBQ", None), "5PTI": ("5PTI", None),
         "1CRN": ("1CRN", None), "1L2Y": ("1L2Y", None), "1UBQ_ACE_10-15_NME": ("1UBQ", (10, 15, True)),
         "1UBQ_ACE_50-55_NME": ("1UBQ", (50, 55, True)), "1UBQ_1-5": ("1UBQ", (1, 5, False)), "1UBQ_70-76": ("1UBQ", (70, 76, False))}
out = {}
for name, (pdbid, frag) in cases.items():
    L = atom_lines(pdbid)
    txt = "".join(L) + "END\n" if frag is None else fragment(L, frag[0], frag[1], ace=frag[2], nme=frag[2])
    top, pos = prepare(txt)
    lines = open(os.path.join(T, name + ".extxyz")).read().splitlines()[2:]
    sym = [l.split()[0] for l in lines]
    xyz = np.array([[float(x) for x in l.split()[1:4]] for l in lines])
    ftype = [l.split()[4] for l in lines]
    d = typed(top, pos)
    assert d["symbols"] == sym and list(d["fftype"]) == ftype, name
    out[name] = energies(system_vacuum(top), xyz)
    print(name, {k: round(v, 4) for k, v in out[name].items()}, flush=True)
json.dump(out, open(sys.argv[2], "w"), indent=1)
