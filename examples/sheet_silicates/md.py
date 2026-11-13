"""MD of the 8 sheet-silicate supercells (paper: Table "sheet silicates", "MD with LAMMPS"): soapff's ClayFF force field (bundled
library) and its MD input (minimisation, 10 ps NVT relaxation, 10 ps NVT, 40 ps NPT with the full cell and tilt at 1 atm;
300 K, 0.5 fs, frames every 1 ps), run with LAMMPS in md_<mineral>/ (supercell.json, the repeats of the unit cell and the
input cell, is kept).  run.py then types the last frame and averages the lattice.

    python md.py [mpirun lmp np] [mineral ...]   (with a Python where soapff is installed; default: mpirun lmp 2, all 8)  ->  md_<mineral>/"""
import os
import subprocess
import sys

from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "prebuilt", "library_clayff_bulk.npz")
MINERALS = ["muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "kaolinite", "lizerdite", "talc", "pyrophyllite_2"]
MPIRUN, LMP, NP = (sys.argv[1:4] + ["mpirun", "lmp", "2"][len(sys.argv[1:4]):])

for m in sys.argv[4:] or MINERALS:
    d = os.path.join(HERE, "md_" + m)
    core.make(os.path.join(HERE, m + ".extxyz"), d, LIB, nvt_ps=10, npt_ps=40, baro="tri", relax_ps=10, dump_ps=1.0,
              stop_on_fallback=True, log=lambda *a: None)
    subprocess.run(f"{MPIRUN} -np {NP} {LMP} -in in.lammps -log log.lammps -screen none", shell=True, cwd=d, check=True)
    print(m, "done", flush=True)
