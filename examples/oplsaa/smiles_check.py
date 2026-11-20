"""The check behind examples/polymer_smiles.py: the PMMA 10-mer that RDKit makes from its SMILES gets from soapff (library
oplsaa, elements and positions only) the same OPLS-AA types and charges as from Foyer.

    python smiles_check.py foyer   RDKit and Foyer (the environment of polymers.py): writes extra/PMMA_10mer_rdkit.extxyz
                                   and its .json, the chain of polymer_smiles.py typed by Foyer
    python smiles_check.py WORK    soapff: types that chain in WORK and writes results/smiles_check.json"""
import json
import os
import sys

import ase.io
import numpy as np

import polymers as P

HERE = os.path.dirname(os.path.abspath(__file__))
SMILES = "CC(C)(C(=O)OC)" * 10                    # the SMILES of polymer_smiles.py
PATH = os.path.join(HERE, "extra", "PMMA_10mer_rdkit.extxyz")
OUT = os.path.join(HERE, "results", "smiles_check.json")


def foyer():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))      # as in polymer_smiles.py
    AllChem.EmbedMolecule(mol, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(mol)
    m = ([a.GetSymbol() for a in mol.GetAtoms()], [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])
    P.save(PATH, [m], mol.GetConformer().GetPositions(), comment="PMMA 10-mer from its SMILES (RDKit), typed by Foyer",
           smiles=[SMILES])
    print(len(m[0]), "atoms typed by Foyer:", PATH)


def check(work):
    import soapff
    ref = ase.io.read(PATH)
    s = soapff.build("oplsaa", PATH, work, log=lambda *a: None)
    tt, qq = np.array(s.types), s.charges
    r = dict(structure=os.path.relpath(PATH, HERE), smiles=SMILES, atoms=len(ref),
             wrong_type=int((tt != ref.arrays["fftype"]).sum()),
             wrong_charge=int((np.abs(qq - ref.arrays["ffcharge"]) > 1e-6).sum()),
             types=len(set(tt)), min_similarity_soap=round(s.report["min_similarity_soap"], 4),
             net_charge=round(float(qq.sum()), 6), terms=s.report["terms"])
    json.dump(r, open(OUT, "w"), indent=1)
    print(r)


if __name__ == "__main__":
    foyer() if sys.argv[1] == "foyer" else check(sys.argv[1])
