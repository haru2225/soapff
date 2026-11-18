"""Reference systems of water and aqueous ions in vacuum for the clayff libraries (src/soapff/libraries/clayff/water/vac_*.data,
the refs/vac_*/system.data written here): the environment of a water molecule or an ion at a free water surface (a water
molecule is typed on its own bonds, taken apart from the mineral in step 1, so a water at the surface of a film or a droplet on
a mineral sees this environment).  packmol, flexible SPC water, ClayFF Table 1 ions, not relaxed, no MD.  The stored files
were typed by an earlier soapff, with the water and the ions by rule; the bundled clayff library gives them the same types.

  vac_water_film       water film 20 A thick in a 24.83 x 24.83 A cell (periodic in x and y), vacuum above and below (cell 50 A)
  vac_water_film_thin  the same, 10 A thick
  vac_water_droplet    water sphere of radius 13 A in a 50 A cube
  vac_<salt>_film      the 20 A film with NaCl, KCl, CsCl (7 + 7 ions), CaCl2, BaCl2 (5 cations + 10 Cl-); the concentration of the liquids of build_refs.py

    python build_refs_vacuum.py     (with a Python where soapff is installed)  ->  refs/vac_*/ (structure.extxyz, system.data)"""
import os
import subprocess
import tempfile

import ase.io
import numpy as np
from ase import Atoms

from soapff_cui import core, water

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "prebuilt", "library_clayff_bulk.npz")
L = 24.83
H = 50.0                                                     # cell height of the films and the droplet cube (A)
RHO = 0.03343                                                # molecules per A^3 at 1 g/cm3
quiet = lambda *a: None
SYSTEMS = {"vac_water_film": ("film", 20.0, 0, []), "vac_water_film_thin": ("film", 10.0, 0, []),
           "vac_water_droplet": ("sphere", 13.0, 0, []),
           "vac_NaCl_film": ("film", 20.0, 0, [("Na", 7), ("Cl", 7)]), "vac_KCl_film": ("film", 20.0, 0, [("K", 7), ("Cl", 7)]),
           "vac_CsCl_film": ("film", 20.0, 0, [("Cs", 7), ("Cl", 7)]), "vac_CaCl2_film": ("film", 20.0, 0, [("Ca", 5), ("Cl", 10)]),
           "vac_BaCl2_film": ("film", 20.0, 0, [("Ba", 5), ("Cl", 10)])}


def build(kind, size, ions, seed=2, tolerance=2.0):
    h = np.radians(water.HOH / 2)
    mol = {"water": (["O", "H", "H"], [[0, 0, 0], [water.OH * np.sin(h), 0, water.OH * np.cos(h)],
                                       [-water.OH * np.sin(h), 0, water.OH * np.cos(h)]])}
    if kind == "film":
        cell, region, nw = [L, L, H], f"inside box 0. 0. {(H - size) / 2} {L} {L} {(H + size) / 2}", int(round(RHO * L * L * size))
    else:
        cell, region, nw = [H, H, H], f"inside sphere {H / 2} {H / 2} {H / 2} {size}", int(round(RHO * 4 / 3 * np.pi * size ** 3))
    with tempfile.TemporaryDirectory() as d:
        inp = [f"tolerance {tolerance}", "filetype xyz", "output packed.xyz", f"seed {seed}", f"pbc 0. 0. 0. {cell[0]} {cell[1]} {cell[2]}"]
        for m, n in [("water", nw)] + list(ions):
            sym, pos = mol.get(m, ([m], [[0, 0, 0]]))
            water._xyz(os.path.join(d, m + ".xyz"), sym, np.array(pos, float))
            inp += ["", f"structure {m}.xyz", f"  number {n}", f"  {region}", "end structure"]
        open(os.path.join(d, "pack.inp"), "w").write("\n".join(inp) + "\n")
        subprocess.run("packmol < pack.inp > packmol.log", shell=True, cwd=d, check=True)
        if "Success!" not in open(os.path.join(d, "packmol.log")).read():
            raise RuntimeError(f"packmol did not converge: {kind} {size} {ions}")
        sym, pos = water._read_xyz(os.path.join(d, "packed.xyz"))
    return Atoms(sym, positions=pos, cell=cell, pbc=True)


if __name__ == "__main__":
    for name, (kind, size, _, ions) in SYSTEMS.items():
        d = os.path.join(HERE, "refs", name)
        os.makedirs(d, exist_ok=True)
        s = build(kind, size, ions)
        ase.io.write(os.path.join(d, "structure.extxyz"), s, format="extxyz")
        r = core.make(os.path.join(d, "structure.extxyz"), d, LIB, strict=True, log=quiet)
        print(name, r["atoms"], "atoms", flush=True)
