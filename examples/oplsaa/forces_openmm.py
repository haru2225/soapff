"""OpenMM forces (kcal/mol/A) on every atom with Foyer's OPLS-AA, for the comparison with the forces of the files soapff
writes (forces_lammps.py).

    python forces_openmm.py TESTDIR OUT.npz          (foyer, parmed, openmm)

Every test structure (TESTDIR/*.extxyz and its .json): the topology is rebuilt from the molecules of the .json, typed by
Foyer (the types are checked against the extxyz) and made an OpenMM System by ParmEd (geometric mixing as a
CustomNonbondedForce, 1-4 x 0.5).  Chains in vacuum: no cutoff; amorphous cells: PME (error tolerance 1e-8: the
difference to LAMMPS is then 100 times smaller than with 1e-6, so what is left is the k-space accuracy) and LJ cut at 12 A,
no switching, no dispersion correction.  Reference platform (double precision); flexible bonds."""
import glob
import os
import sys

import numpy as np
import openmm as mm
from openmm import unit

import polymers as P

T, OUT = sys.argv[1:3]
KCAL_A = unit.kilocalorie_per_mole / unit.angstrom
out = {}
for f in sorted(glob.glob(os.path.join(T, "*.extxyz"))):
    name = os.path.basename(f)[:-7]
    mols, xyz, box, types, _ = P.load(f)
    s, sysm = P.system(mols, xyz, box, pme_tol=1e-8)
    assert [a.type for a in s.atoms] == types, name
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    if box is not None:
        ctx.setPeriodicBoxVectors(*(mm.Vec3(*(np.eye(3)[i] * box[i] * 0.1)) for i in range(3)))
    ctx.setPositions(xyz * 0.1)
    st = ctx.getState(getForces=True, getEnergy=True)
    out[name] = st.getForces(asNumpy=True).value_in_unit(KCAL_A)
    print(f"{name:16s} {len(xyz):6d} atoms  E {st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole):.4f} "
          f"kcal/mol", flush=True)
np.savez_compressed(OUT, **out)
