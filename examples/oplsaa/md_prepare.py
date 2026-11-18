"""The start of the melt MD (md_lammps.py, md_openmm.py): an amorphous cell of make_cells.py (packmol, 0.25 g/cm3)
compressed and equilibrated with Foyer's OPLS-AA in OpenMM.

    python md_prepare.py TESTDIR WORK PE PP PIB          (foyer, openmm) -> WORK/start_<p>.extxyz (+ .json)

Minimisation, 50 ps NVT at 600 K, 100 ps NPT at 413 K and 1000 bar (compression), 200 ps NPT at 413 K and 1 bar; Langevin
(1/ps), Monte Carlo barostat, 1 fs, flexible bonds, PME 1e-5, LJ cut at 12 A with the long-range correction (as the files
of soapff: pair_modify tail yes).  Platform: polymers.platform() (OpenCL in single precision when there is a GPU)."""
import os
import sys

import numpy as np
import openmm as mm
from openmm import unit

import polymers as P

T, W = sys.argv[1:3]
os.makedirs(W, exist_ok=True)
for p in sys.argv[3:]:
    mols, xyz, box, _, _ = P.load(os.path.join(T, f"cell_{p}.extxyz"))
    s, sysm = P.system(mols, xyz, box, pme_tol=1e-5, dispersion=True)
    baro = mm.MonteCarloBarostat(1000 * unit.bar, 413 * unit.kelvin, 0)       # off during the NVT at 600 K
    sysm.addForce(baro)
    integ = mm.LangevinMiddleIntegrator(600 * unit.kelvin, 1.0 / unit.picosecond, 1.0 * unit.femtoseconds)
    integ.setRandomNumberSeed(11)
    ctx = mm.Context(sysm, integ, *P.platform())
    ctx.setPeriodicBoxVectors(*(mm.Vec3(*(np.eye(3)[i] * box[i] * 0.1)) for i in range(3)))
    ctx.setPositions(xyz * 0.1)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 2000)
    ctx.setVelocitiesToTemperature(600 * unit.kelvin, 11)
    m = P.mass(mols) / 0.60221408                         # g/mol -> 1e-24 g (density in g/cm3 = m / V[A^3])

    def dens():
        v = ctx.getState().getPeriodicBoxVolume().value_in_unit(unit.angstrom ** 3)
        return m / v
    integ.step(50000)
    print(f"{p}: after 50 ps at 600 K (NVT) {dens():.3f} g/cm3", flush=True)
    baro.setFrequency(25)
    ctx.reinitialize(preserveState=True)
    integ.setTemperature(413 * unit.kelvin)
    for k in range(10):
        integ.step(10000)
    print(f"{p}: after 100 ps at 413 K, 1000 bar {dens():.3f} g/cm3", flush=True)
    ctx.setParameter(mm.MonteCarloBarostat.Pressure(), 1.0)                     # bar
    for k in range(20):
        integ.step(10000)
    st = ctx.getState(getPositions=True)
    x = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    a, b, c = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
    box = np.array([a[0], b[1], c[2]])
    print(f"{p}: after 200 ps at 413 K, 1 bar {dens():.3f} g/cm3, box {box[0]:.2f} A", flush=True)
    P.save(os.path.join(W, f"start_{p}.extxyz"), mols, x, box=box, comment=f"cell_{p} compressed and equilibrated at 413 K")
