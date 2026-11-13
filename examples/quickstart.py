"""soapff in a few lines: a structure file of any format + a library -> the LAMMPS input (and the run).

    python examples/quickstart.py [--run]        (from the root of the repository)"""
import os
import sys

import ase.io

import soapff

for lib in soapff.libraries():                        # the libraries that can be given by name
    print(f"{lib['name']:12s} {lib['engine']:9s} {lib['description'][:90]}")

# ClayFF: kaolinite with water on its basal plane (any format ASE reads: here extxyz; cif, POSCAR, xyz with a cell, ... alike)
clay = soapff.build("clayff", "examples/clayff/kaolinite_water.extxyz", "quickstart_out/kaolinite_water",
                    relax_ps=1.0, nvt_ps=5.0)
print(clay)

# a protein: Trp-cage (PDB 1L2Y with hydrogens) from an xyz file without types; no cell -> put in a box with 10 A of padding
os.makedirs("quickstart_out", exist_ok=True)
ase.io.write("quickstart_out/trpcage.xyz", ase.io.read("examples/amber14sb/test/1L2Y.extxyz"))
prot = soapff.build("amber14sb", "quickstart_out/trpcage.xyz", "quickstart_out/trpcage", md=False)
print(prot, prot.report["terms"], prot.report["molecule_kinds"])
print(prot.table()[:3])                                # charges.csv: id, element, type, charge, molecule, similarities, ...

if "--run" in sys.argv:
    prot.run(np=1)                                     # lmp (and mpirun when np > 1); $SOAPFF_LMP picks the executable
    print("final energy", prot.thermo()["PotEng"][-1], "kcal/mol")
