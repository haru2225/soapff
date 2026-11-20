"""soapff in a few lines.  soapff types every atom of a structure file with a library and writes the force field for
LAMMPS; what is done with it (minimisation, MD) is written in LAMMPS, here with its Python module.

    python examples/quickstart.py [--run]        (from the root of the repository; --run needs pip install "lammps[mpi]")"""
import os
import sys

import ase.io

import soapff

for lib in soapff.libraries():                        # the libraries that can be given by name
    print(f"{lib['name']:14s} {lib['reference'][:90]}")

# ClayFF: kaolinite with water on its basal plane (any format ASE reads: here extxyz; cif, POSCAR, xyz with a cell, ... alike)
clay = soapff.build("clayff", "examples/clayff/kaolinite_water.extxyz", "quickstart_out/kaolinite_water")
print(clay, clay.type_names)

# a PMMA chain from a plain xyz file (elements and positions only, no cell -> put in a box with 10 A of padding)
os.makedirs("quickstart_out", exist_ok=True)
ase.io.write("quickstart_out/pmma.xyz", ase.io.read("examples/oplsaa/test/chain_PMMA.extxyz"), format="xyz")
pmma = soapff.build("oplsaa", "quickstart_out/pmma.xyz", "quickstart_out/pmma")
print(pmma, pmma.report["terms"], pmma.report["molecule_kinds"])
print(pmma.table()[:3])                                # charges.csv: id, element, type, charge, molecule, similarities, ...

# a polyethylene melt: 9 chains of C100 in a periodic box, equilibrated at 413 K (examples/oplsaa/md_prepare.py)
pe = soapff.build("oplsaa", "examples/oplsaa/md/start_PE.extxyz", "quickstart_out/pe")
print(pe, pe.types[:3], pe.charges[:3])
print(pe.commands()[:6], "...")                        # the force field as LAMMPS commands (system.settings)

if "--run" in sys.argv:
    from lammps import lammps
    # the MD is written in LAMMPS: here 10 ps of NPT at 413 K and 1 atm after a minimisation
    lmp = lammps(cmdargs=["-log", "quickstart_out/pe/log.lammps", "-screen", "none"])
    lmp.commands_list(pe.commands())                   # the force field: styles, read_data, coefficients, kspace
    lmp.commands_string("""
        thermo_style custom step temp press pe density
        thermo 1000
        minimize 1.0e-4 1.0e-6 1000 10000
        reset_timestep 0
        velocity all create 413 1 dist gaussian
        timestep 1.0
        fix 1 all npt temp 413 413 100 iso 1 1 1000
        run 10000
    """)
    print(f"PE melt after 10 ps of NPT at 413 K: {lmp.get_thermo('density'):.3f} g/cm3, {lmp.get_thermo('temp'):.0f} K")
    lmp.close()
