"""The clay + polymer + water test structures: a smectite layer built here (not a reference of clayff), with water, Na+ and
two polymer chains in its interlayer (packmol), and the right types and charges of every atom.

- the layer: the 2:1 layer of the pyrophyllite crystal of clayff (its z = 9.2 A layer, 2 x 2 cells), repeated 3 x 2
  (6 x 4 cells, 31.0 x 35.9 A), with substitutions placed at random apart from each other:
    PEO: 16 of the 96 octahedral Al by Mg, no two Mg nearer than 5 A (a montmorillonite, 0.33 per O10(OH)2);
    PVA: 18 of the 192 tetrahedral Si by Al, no two Al nearer than 4 A (Loewenstein; a beidellite, 0.375 per O10(OH)2);
  the clayff library has one montmorillonite (Mg in 1 of 8 octahedral sites) and tetrahedral Al only in micas (1 in 4).
- the interlayer: 22 A between the O planes, with as many Na+ as substitutions, two chains of chains/PEO or chains/PVA and
  water at about 1 g/cm3 (packmol, periodic in x and y).
- the right answer: the mineral, the ClayFF type of the bond-graph rule of ClayFF (soapff_cui.core.rule_types) with the
  charges of ClayFF Table 1; the water ow/hw (-0.82, +0.41), Na+ (+1); the chains, Foyer's (chains.py).

    python build.py        (soapff)  ->  structures/{PEO,PVA}_start.extxyz (fftype, ffcharge: the right answer)"""
import os
import subprocess
import tempfile

import numpy as np
from ase import Atoms
from ase.io import read, write

from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff")
GAP = 22.0                   # A between the O planes of the interlayer
SYSTEMS = {"PEO": dict(kind="oct", n=16, dmin=5.0, seed=1), "PVA": dict(kind="tet", n=18, dmin=4.0, seed=2)}
WATER = 0.0334               # molecules per A3 (1 g/cm3)


def layer():
    """the 2:1 layer of pyrophyllite (2 x 2 cells) repeated 3 x 2, with c along z"""
    at, _, _ = core.read_lammps_data(os.path.join(LIB, "bulk", "pyrophyllite_bulk.data"))
    z = at.positions[:, 2]
    one = at[(z > 5.0) & (z < 13.5)]
    cell = np.array(at.cell)
    one.set_cell([cell[0], cell[1], [0.0, 0.0, 30.0]])
    big = one.repeat((3, 2, 1))
    big.positions[:, 2] -= big.positions[:, 2].min() - 0.5
    return big


def substitute(at, kind, n, dmin, seed):
    """n substitutions (oct: Al -> Mg, tet: Si -> Al), at random, none nearer than dmin to another"""
    rng = np.random.default_rng(seed)
    el = np.array(at.get_chemical_symbols())
    site = np.nonzero(el == ("Al" if kind == "oct" else "Si"))[0]
    for _ in range(2000):
        pick = []
        for i in rng.permutation(site):
            if all(at.get_distance(i, j, mic=True) > dmin for j in pick):
                pick.append(int(i))
            if len(pick) == n:
                break
        if len(pick) == n:
            break
    assert len(pick) == n, "no arrangement found"
    el[pick] = "Mg" if kind == "oct" else "Al"
    at.set_chemical_symbols(el)
    return at


def water():
    """an SPC water molecule (O-H 1.0 A, H-O-H 109.47 deg)"""
    t = np.radians(109.47 / 2)
    return Atoms("OH2", positions=[[0, 0, 0], [np.sin(t), np.cos(t), 0], [-np.sin(t), np.cos(t), 0]])


def pack(cell, zlo, zhi, items, seed):
    """packmol: items [(Atoms, number)] in x, y periodic (0..Lx, 0..Ly) and zlo..zhi -> the packed atoms, in order"""
    Lx, Ly = cell[0, 0], cell[1, 1]
    with tempfile.TemporaryDirectory() as d:
        L = [f"tolerance 2.0\nfiletype xyz\noutput {d}/out.xyz\nseed {seed}\n"
             f"pbc 0.0 0.0 {zlo - 3:.3f} {Lx:.3f} {Ly:.3f} {zhi + 3:.3f}\n"]
        for k, (a, num) in enumerate(items):
            write(f"{d}/m{k}.xyz", a)
            L.append(f"structure {d}/m{k}.xyz\n  number {num}\n  inside box 0.0 0.0 {zlo:.3f} {Lx:.3f} {Ly:.3f} {zhi:.3f}\n"
                     f"end structure\n")
        open(f"{d}/in.inp", "w").write("".join(L))
        r = subprocess.run(f"packmol < {d}/in.inp", shell=True, capture_output=True, text=True)
        if "Success" not in r.stdout:
            raise RuntimeError(r.stdout[-2000:])
        return read(f"{d}/out.xyz")


def main():
    os.makedirs(os.path.join(HERE, "structures"), exist_ok=True)
    for name, opt in SYSTEMS.items():
        lay = substitute(layer(), **opt)
        el = np.array(lay.get_chemical_symbols())
        zO = lay.positions[el == "O", 2]
        Lz = (zO.max() - zO.min()) + GAP
        cell = np.array(lay.cell)
        cell[2] = [0.0, 0.0, Lz]
        lay.set_cell(cell)
        lay.pbc = True
        # the right answer of the mineral: the rule of ClayFF on the layer alone
        qm, tm = core.rule_types(lay, np.array(["mineral"] * len(lay)))
        tm = [core.TYPE_ALIAS.get(str(t), str(t)) for t in tm]
        chain = read(os.path.join(HERE, "chains", f"{name}.extxyz"))
        nion = opt["n"]
        area = abs(np.linalg.det(cell[:2, :2]))
        vpoly = 2 * chain.get_masses().sum() / (1.1 * 0.60221)          # A3 at 1.1 g/cm3
        nw = int(round(WATER * (area * (GAP - 2.6) - vpoly))) - nion
        zlo, zhi = zO.max() + 1.8, zO.min() + Lz - 1.8
        packed = pack(cell, zlo, zhi, [(chain, 2), (Atoms("Na"), nion), (water(), nw)], seed=opt["seed"])
        at = lay + packed
        at.set_cell(cell)
        at.pbc = True
        types = list(tm) + list(chain.arrays["fftype"]) * 2 + ["Na"] * nion + ["ow", "hw", "hw"] * nw
        q = np.concatenate([qm, np.tile(chain.arrays["ffcharge"], 2), np.ones(nion), np.tile([-0.82, 0.41, 0.41], nw)])
        at.set_array("fftype", np.array(types))
        at.set_array("ffcharge", q)
        out = os.path.join(HERE, "structures", f"{name}_start.extxyz")
        write(out, at)
        sub = {t: types[:len(lay)].count(t) for t in sorted(set(types[:len(lay)]))}
        print(f"{name}: {len(at)} atoms (layer {len(lay)}: {sub}; 2 chains of {len(chain)}; {nion} Na; {nw} water); "
              f"cell {np.round(at.cell.lengths(), 2)}; net charge {q.sum():+.4f}; -> {out}")


if __name__ == "__main__":
    main()
