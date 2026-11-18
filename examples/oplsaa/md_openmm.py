"""NPT MD of a polymer melt with Foyer's OPLS-AA in OpenMM, for the comparison with the MD of the files soapff writes
(md_lammps.py): the same start (md_prepare.py), 413 K and 1 bar.

    python md_openmm.py WORK PS P1 [P2 ...]       (foyer, openmm) -> WORK/openmm_<p>.csv (time ps, density g/cm3, T K)

Langevin (1/ps) at 413 K, Monte Carlo barostat (every 25 steps), 1 fs, flexible bonds, PME 1e-5, LJ cut at 12 A with the
long-range correction; the density every ps.  Platform: polymers.platform() (OpenCL in single precision when there is a
GPU; $OPENMM_PLATFORM chooses another)."""
import os
import sys

import numpy as np
import openmm as mm
from openmm import unit

import polymers as P

W, PS = sys.argv[1], float(sys.argv[2])
for p in sys.argv[3:]:
    mols, xyz, box, _, _ = P.load(os.path.join(W, f"start_{p}.extxyz"))
    s, sysm = P.system(mols, xyz, box, pme_tol=1e-5, dispersion=True)
    sysm.addForce(mm.MonteCarloBarostat(1 * unit.bar, 413 * unit.kelvin, 25))
    integ = mm.LangevinMiddleIntegrator(413 * unit.kelvin, 1.0 / unit.picosecond, 1.0 * unit.femtoseconds)
    integ.setRandomNumberSeed(3)
    ctx = mm.Context(sysm, integ, *P.platform())
    ctx.setPeriodicBoxVectors(*(mm.Vec3(*(np.eye(3)[i] * box[i] * 0.1)) for i in range(3)))
    ctx.setPositions(xyz * 0.1)
    ctx.setVelocitiesToTemperature(413 * unit.kelvin, 3)
    m = P.mass(mols) / 0.60221408
    ndf = 3 * len(xyz)
    out = open(os.path.join(W, f"openmm_{p}.csv"), "w")
    out.write("time_ps,density,temperature\n")
    for k in range(int(PS)):
        integ.step(1000)
        st = ctx.getState(getEnergy=True)
        v = st.getPeriodicBoxVolume().value_in_unit(unit.angstrom ** 3)
        temp = 2 * st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole) / (ndf * 0.0083144626)
        out.write(f"{k + 1},{m / v:.5f},{temp:.1f}\n")
        out.flush()
    print(f"{p}: {PS:g} ps done", flush=True)
