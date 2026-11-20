"""The SWy-1 model of Szczerba & Kalinichev (2016) with ethylene glycol (EG) and water, built here, with the right types and
charges of every atom.

As theirs: the pyrophyllite crystal (the one of clayff, 2 layers) repeated to 8 x 4 x 2 unit cells (41.3 x 35.9 A, 2560
atoms), Al replaced by Mg in the octahedral sheets with the Mg as far apart as can be (the best of many random
arrangements), as many as the octahedral charge of SWy-1 (0.28 per half unit cell, phuc: 18 Mg in each layer of 64 half
cells), Ca2+ in the interlayers (9 per interlayer), 0.8 H2O phuc (51 per interlayer) and 1.7 or 1.5 EG phuc (109 or 96 per
interlayer).  The interlayers are opened to 21 A for packmol (2 A kept between the packed atoms and the O planes).

The right answer: the mineral, the ClayFF type of the bond-graph rule of ClayFF (soapff_cui.core.rule_types) with the
charges of ClayFF Table 1; water ow/hw (-0.82, +0.41); Ca (+2); EG, Foyer's (refs/EG.extxyz).

    python build.py            (soapff)  ->  structures/SWy-1_EG{1.7,1.5}_start.extxyz (fftype, ffcharge)"""
import os
import subprocess
import tempfile

import numpy as np
from ase import Atoms
from ase.io import read, write

from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
CRYSTAL = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "bulk", "pyrophyllite_bulk.data")
HALF_CELLS = 64                 # half unit cells (O10(OH)2) per layer: 8 x 4 unit cells
OCT_CHARGE = 0.28               # SWy-1, Table 1 of Szczerba & Kalinichev (tetrahedral charge 0)
WATER_PHUC = 0.8
D0 = 21.0                       # A: the interlayer spacing for packmol


def crystal():
    """the 2-layer pyrophyllite crystal repeated 4 x 2 x 1 (8 x 4 x 2 unit cells), z shifted so that no layer is cut by
    the cell"""
    at, _, _ = core.read_lammps_data(CRYSTAL)
    at = at.repeat((4, 2, 1))
    sp = at.get_scaled_positions()
    sp[:, 2] = (sp[:, 2] + 0.25) % 1.0
    at.set_scaled_positions(sp)
    return at


def spread(at, sites, n, seed, tries=400):
    """n of the sites, as far apart as can be: the best (largest smallest distance) of many random greedy arrangements"""
    rng = np.random.default_rng(seed)
    D = at.get_all_distances(mic=True)[np.ix_(sites, sites)]
    best, bd = None, -1.0
    for _ in range(tries):
        order = rng.permutation(len(sites))
        pick = [order[0]]
        for k in order[1:]:
            if len(pick) == n:
                break
            if D[k, pick].min() > bd:
                pick.append(k)
        if len(pick) < n:
            continue
        d = D[np.ix_(pick, pick)][np.triu_indices(n, 1)].min()
        if d > bd:
            best, bd = pick, d
    return [int(sites[k]) for k in best], bd


def water():
    t = np.radians(109.47 / 2)
    return Atoms("OH2", positions=[[0, 0, 0], [np.sin(t), np.cos(t), 0], [-np.sin(t), np.cos(t), 0]])


def pack(cell, zlo, zhi, items, seed):
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


def model(eg_phuc, seed=1):
    at = crystal()
    el = np.array(at.get_chemical_symbols())
    z = at.positions[:, 2]
    Lz0 = at.cell[2, 2]
    lower = z < Lz0 / 2                                      # the two layers
    n_mg = int(round(OCT_CHARGE * HALF_CELLS))
    dmin = []
    for k, lay in enumerate((lower, ~lower)):
        sites = np.nonzero(lay & (el == "Al"))[0]
        pick, d = spread(at, sites, n_mg, seed + k)
        el[pick] = "Mg"
        dmin.append(round(d, 2))
    at.set_chemical_symbols(el)
    # open the interlayers to D0: move the upper layer up, lengthen c along z
    dz = D0 - Lz0 / 2
    pos = at.positions
    pos[~lower, 2] += dz
    at.positions = pos
    cell = np.array(at.cell)
    cell[2, 2] = 2 * D0
    at.set_cell(cell)
    at.pbc = True
    qm, tm = core.rule_types(at, np.array(["mineral"] * len(at)))
    tm = [core.TYPE_ALIAS.get(str(t), str(t)) for t in tm]
    eg = read(os.path.join(HERE, "refs", "EG.extxyz"))
    n_eg = int(round(eg_phuc * HALF_CELLS))
    n_w = int(round(WATER_PHUC * HALF_CELLS))
    n_ca = n_mg // 2
    oz = at.positions[el == "O", 2]
    lowO, upO = oz[lower[el == "O"]], oz[~lower[el == "O"]]
    gaps = [(lowO.max() + 2.0, upO.min() - 2.0), (upO.max() + 2.0, lowO.min() + 2 * D0 - 2.0)]
    out = at
    types, q = list(tm), list(qm)
    for k, (zlo, zhi) in enumerate(gaps):
        p = pack(cell, zlo, zhi, [(eg, n_eg), (Atoms("Ca"), n_ca), (water(), n_w)], seed=seed + 10 + k)
        out = out + p
        types += list(eg.arrays["fftype"]) * n_eg + ["Ca"] * n_ca + ["ow", "hw", "hw"] * n_w
        q += list(np.tile(eg.arrays["ffcharge"], n_eg)) + [2.0] * n_ca + [-0.82, 0.41, 0.41] * n_w
    out.set_cell(cell)
    out.pbc = True
    out.set_array("fftype", np.array(types))
    out.set_array("ffcharge", np.array(q))
    return out, dict(mg_per_layer=n_mg, mg_min_distance=dmin, eg_per_interlayer=n_eg, water_per_interlayer=n_w,
                     ca_per_interlayer=n_ca, mineral={t: tm.count(t) for t in sorted(set(tm))})


def main():
    os.makedirs(os.path.join(HERE, "structures"), exist_ok=True)
    for eg in (1.7, 1.5):
        at, info = model(eg)
        path = os.path.join(HERE, "structures", f"SWy-1_EG{eg}_start.extxyz")
        write(path, at)
        print(f"SWy-1, {eg} EG phuc: {len(at)} atoms, net charge {at.arrays['ffcharge'].sum():+.4f}; {info}")


if __name__ == "__main__":
    main()
