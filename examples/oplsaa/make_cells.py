"""Amorphous cells of polymer chains (test structures, and the start of the melt MD), typed by Foyer (OPLS-AA).

    python make_cells.py TESTDIR [P ...]       (rdkit, openbabel, foyer, openmm, packmol) -> TESTDIR/cell_<p>.extxyz (+ .json)

Chains of PE (9 x C100), atactic PP (6 x 50 units), PIB (6 x 40), PEO (8 x 20) and atactic PMMA (6 x 15): every chain built
as in make_refs.py (Open Babel, minimised with Foyer's force field), made compact by 15 ps of MD at 600 K in a sphere
(polymers.cell), packed by packmol with periodic boundaries at 0.25 g/cm3, and the cell minimised with Foyer's force field
(PME; the bonds and angles of the 600 K MD relax, the compact shapes stay).  md_prepare.py compresses the cell to the melt.
A cell that exists is kept (a cell takes minutes).  The chains are whole (they may stick out of the box)."""
import os
import sys

import numpy as np

import polymers as P

T = os.path.abspath(sys.argv[1])
CELLS = {"PE": (50, 9), "PP": (50, 6), "PIB": (40, 6), "PEO": (20, 8), "PMMA": (15, 6)}     # repeat units, chains
for p in sys.argv[2:] or CELLS:
    n, k = CELLS[p]
    name = f"cell_{p}"
    if os.path.isfile(os.path.join(T, name + ".extxyz")):
        continue
    mols, xyz, box, smi3 = P.cell([P.chain_smiles(p, n, seed=20 + i) for i in range(k)], 0.25, seed=20)
    s = P.save(os.path.join(T, name + ".extxyz"), mols, xyz, box=box, comment=f"{k} chains of {n} repeat units", smiles=smi3)
    q = np.array([a.charge for a in s.atoms])
    print(f"test {name:18s} {len(q):5d} atoms  net {q.sum():+.4f}  box {box[0]:.2f} A", flush=True)
