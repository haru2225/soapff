"""One ethylene glycol (EG, HOCH2CH2OH) molecule minimised with Foyer's OPLS-AA (OpenMM, vacuum) and typed by Foyer, as the
references of oplsaa (../oplsaa/polymers.py): the molecule that build.py and eg_liquid.py pack, and its right answer.  It is
not a reference of a library: oplsaa and clayff_oplsaa type EG as Foyer without it.  Needs foyer, openmm and parmed (the
environment of ../oplsaa).

    python eg_refs.py          -> refs/EG.extxyz (+ .json)"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "oplsaa"))
import polymers  # noqa: E402

# HO-CH2-CH2-OH, trans O-C-C-O, from bond lengths and tetrahedral angles (minimised below)
EL = ["O", "C", "C", "O", "H", "H", "H", "H", "H", "H"]
BONDS = [(0, 1), (1, 2), (2, 3), (0, 4), (1, 5), (1, 6), (2, 7), (2, 8), (3, 9)]


def start():
    x = np.zeros((10, 3))
    x[1] = [0.0, 0.0, 0.0]
    x[2] = [1.53, 0.0, 0.0]
    x[0] = [-0.48, 1.35, 0.0]
    x[3] = [2.01, -1.35, 0.0]
    x[4] = [-1.44, 1.35, 0.0]
    x[9] = [2.97, -1.35, 0.0]
    x[5] = [-0.36, -0.51, 0.89]
    x[6] = [-0.36, -0.51, -0.89]
    x[7] = [1.89, 0.51, 0.89]
    x[8] = [1.89, 0.51, -0.89]
    return x


def dihedral(x, a, b, c, d):
    b0, b1, b2 = x[a] - x[b], x[c] - x[b], x[d] - x[c]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def main():
    os.makedirs(os.path.join(HERE, "refs"), exist_ok=True)
    mol = (EL, BONDS)
    x = polymers.relax(mol, start())
    s = polymers.save(os.path.join(HERE, "refs", "EG.extxyz"), [mol], x - x.mean(0),
                      comment="ethylene glycol, minimised (Foyer OPLS-AA); types and charges by Foyer", smiles=["OCCO"])
    print("EG", [a.type for a in s.atoms], [round(a.charge, 3) for a in s.atoms], f"O-C-C-O {dihedral(x, 0, 1, 2, 3):.0f}")


if __name__ == "__main__":
    main()
