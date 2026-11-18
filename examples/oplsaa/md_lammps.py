"""NPT MD of a polymer melt with soapff's force field (library oplsaa) in LAMMPS, run from Python with the LAMMPS Python
module, for the comparison with Foyer's OPLS-AA in OpenMM (md_openmm.py): the same start (md_prepare.py), 413 K and 1 bar.

    python md_lammps.py WORK PS P1 [P2 ...]      (soapff, LAMMPS Python module) -> WORK/lammps_<p>/log.lammps, WORK/lammps_<p>.csv

soapff writes the force field (system.data, system.settings: PPPM 1e-5, LJ cut at 12 A with the tail correction,
geometric mixing, 1-4 pairs scaled by 0.5).  The MD is written here: velocities at 413 K, Nose-Hoover NPT (isotropic,
1 bar = 0.986923 atm; damping 100 fs and 1000 fs), 1 fs, no minimisation (the start is minimised and equilibrated by
md_prepare.py); the density every 100 steps from the thermo output.

In parallel: mpirun -np 4 python md_lammps.py WORK PS P.  Every process runs this script; soapff runs on the first one and
writes the force field as LAMMPS commands (s.commands()) to forcefield.lammps, which LAMMPS reads on the first process and
passes to the others.  The stored results (md/lammps_*.csv) were run with the same commands in an input file, with the
lmp_mpi executable on 3 processes."""
import os
import sys

import numpy as np
from lammps import lammps
from lammps.formats import LogFile

import soapff

W, PS = os.path.abspath(sys.argv[1]), float(sys.argv[2])
for p in sys.argv[3:]:
    out = os.path.join(W, f"lammps_{p}")
    os.makedirs(out, exist_ok=True)
    log, ff = os.path.join(out, "log.lammps"), os.path.join(out, "forcefield.lammps")
    lmp = lammps(cmdargs=["-log", log, "-screen", "none"])
    first = lmp.extract_setting("world_rank") == 0
    if first:
        s = soapff.build("oplsaa", os.path.join(W, f"start_{p}.extxyz"), out, log=lambda *a: None)
        open(ff, "w").write("\n".join(s.commands()) + "\n")
    lmp.file(ff)                                                      # the force field
    s = soapff.System(out)                                            # written by the first process
    lmp.commands_string(f"""
        thermo_style custom step temp press pe ebond eangle edihed eimp evdwl ecoul elong vol density
        thermo 100
        velocity all create 413 4928459 mom yes rot no dist gaussian
        timestep 1.0
        fix npt all npt temp 413 413 100 iso 0.986923 0.986923 1000
        dump d all custom 10000 {os.path.join(s.out_dir, "traj.lammpstrj")} id mol type element x y z ix iy iz q
        dump_modify d sort id element {" ".join(s.elements)}
        run {int(round(PS * 1000))}
        write_data {os.path.join(s.out_dir, "final.data")}
    """)
    lmp.close()
    if not first:
        continue
    run = LogFile(log).runs[-1]
    rows = np.column_stack([np.array(run["Step"]) / 1000.0, run["Density"], run["Temp"]])
    np.savetxt(os.path.join(W, f"lammps_{p}.csv"), rows, delimiter=",", fmt=["%.1f", "%.5f", "%.1f"],
               header="time_ps,density,temperature", comments="")
    print(f"{p}: {len(rows)} thermo rows, last {rows[-1, 0]:g} ps, density {np.mean(rows[len(rows) // 5:, 1]):.4f} g/cm3",
          flush=True)
