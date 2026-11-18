"""A polymer from its SMILES in a few lines: a PMMA 10-mer (H at both ends) made in 3D with hydrogens by RDKit, its OPLS-AA
force field from soapff (library oplsaa), and MD at 413 K in LAMMPS (its Python module).  The SMILES gives only the
elements and the positions to soapff: no types, charges or bond orders are taken from RDKit.

    python examples/polymer_smiles.py [ps]        (soapff, pip install rdkit "lammps[mpi]"; default 100 ps)"""
import sys

from ase import Atoms
from lammps import lammps
from rdkit import Chem
from rdkit.Chem import AllChem

import soapff

ps = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0

# a PMMA 10-mer from its SMILES, 3D with hydrogens (RDKit)
mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)(C(=O)OC)" * 10))
AllChem.EmbedMolecule(mol, randomSeed=1)
AllChem.MMFFOptimizeMolecule(mol)
chain = Atoms([a.GetSymbol() for a in mol.GetAtoms()], positions=mol.GetConformer().GetPositions())

# the force field: no cell, so the chain is put in a periodic box with 10 A of padding
pmma = soapff.build("oplsaa", chain, "polymer_smiles_out")
print(pmma.report["molecule_kinds"], pmma.report["terms"])
for row in pmma.table()[:6]:                               # charges.csv: the type, charge and similarity of every atom
    print(row["element"], row["type"], row["charge"], row["similarity_soap"])

# the MD is written in LAMMPS: minimisation, then NVT at 413 K
lmp = lammps(cmdargs=["-log", "polymer_smiles_out/log.lammps", "-screen", "none"])
lmp.commands_list(pmma.commands())
lmp.commands_string(f"""
    thermo_style custom step temp pe ebond eangle edihed evdwl ecoul elong
    thermo 1000
    minimize 1.0e-4 1.0e-6 1000 10000
    reset_timestep 0
    velocity all create 413 1 mom yes rot yes dist gaussian
    timestep 1.0
    fix 1 all nvt temp 413 413 100
    compute rg all gyration
    fix avg all ave/time 10 {int(ps * 50)} {int(ps * 1000)} c_thermo_temp c_rg
    run {int(ps * 1000)}
""")
T, rg = (lmp.extract_fix("avg", 0, 1, i) for i in (0, 1))      # averages over the second half of the run
print(f"{ps:g} ps of NVT: mean temperature {T:.0f} K, radius of gyration {rg:.2f} A (second half)")
lmp.close()
