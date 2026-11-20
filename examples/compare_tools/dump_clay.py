"""The 16 clayff test structures in the reference arrangement (the 10 (001) slabs with water and the 6 water shapes) ->
IN/NAME.npz (elements, positions, cell, soapff's types) and IN/NAME.txt, the input of run_atom.m: line 1 Box_dim
(GROMACS order, Angstrom), then one atom per line: molid resname type x y z.  As in a .gro file of a solvated clay, the
mineral is one molecule (MIN, element names), and each water molecule (SOL: Ow Hw Hw) and each ion (ION) its own.
usage (soapff's Python): python dump_clay.py IN"""
import glob, os, sys
import numpy as np
from soapff_cui import core
EX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = sys.argv[1]
IONS = {"Na", "K", "Cs", "Ca", "Ba", "Cl"}
os.makedirs(OUT, exist_ok=True)


def save(group, name, at, truth):
    sym, pos, cell = np.array(at.get_chemical_symbols()), at.get_positions(), np.array(at.cell[:])
    np.savez_compressed(os.path.join(OUT, f"{group}__{name}.npz"), symbols=sym, pos=pos, cell=cell, truth=truth)
    pos = (np.linalg.solve(cell.T, pos.T).T % 1.0) @ cell
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    box = [cell[0, 0], cell[1, 1], cell[2, 2]] + ([0, 0, xy, 0, xz, yz] if max(abs(xy), abs(xz), abs(yz)) > 1e-6 else [])
    rows, mol = [], 1
    for s, p, t in zip(sym, pos, truth):
        if t == "ow":
            mol += 1; res, typ, m = "SOL", "Ow", mol
        elif t == "hw":
            res, typ, m = "SOL", "Hw", mol
        elif t in IONS:
            mol += 1; res, typ, m = "ION", str(s), mol
        else:
            res, typ, m = "MIN", str(s), 1
        rows.append(f"{m} {res} {typ} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
    with open(os.path.join(OUT, f"{group}__{name}.txt"), "w") as o:
        o.write(" ".join(f"{b:.6f}" for b in box) + "\n" + "\n".join(rows) + "\n")
    print(group, name, len(sym), flush=True)


for g, sub in (("interfaces", "interfaces_001_water"), ("shapes", "interface_shapes")):
    for d in sorted(glob.glob(os.path.join(EX, sub, "*", "system.data"))):
        at, q, t = core.read_lammps_data(d)
        save(g, d.split("/")[-2], at, np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t]))
