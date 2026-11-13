"""ubiquitin in TIP3P water with 0.15 M NaCl (OpenMM Modeller) -> ubq_water.pdb (with CRYST1), ubq_water.extxyz (types, charges)
and the OpenMM energy of these files with PME (ewald error 1e-6, cutoff 10 A, no dispersion correction, flexible water):
results/energies_solvated.json      (run in this folder: python solvated.py; OpenMM, pdbfixer)"""
import json, os, sys
import numpy as np, openmm as mm
from openmm import app, unit
from omm_common import FF, atom_lines, prepare, typed, write_extxyz, energies
top, pos = prepare("".join(atom_lines("1UBQ")) + "END\n")
mod = app.Modeller(top, pos)
mod.addSolvent(FF, padding=1.0 * unit.nanometer, ionicStrength=0.15 * unit.molar, boxShape="cube")
# addSolvent can leave water molecules overlapping across the periodic boundary (an H 1 A from the O of another water): one
# water of each such pair is removed
xyz = np.array(mod.positions.value_in_unit(unit.angstrom))
L = np.array([mod.topology.getPeriodicBoxVectors()[i][i].value_in_unit(unit.angstrom) for i in range(3)])
wat = [r for r in mod.topology.residues() if r.name == "HOH"]
O = np.array([[a.index for a in r.atoms() if a.element.symbol == "O"][0] for r in wat])
H = np.array([a.index for r in wat for a in r.atoms() if a.element.symbol == "H"])
Hres = np.array([k for k, r in enumerate(wat) for a in r.atoms() if a.element.symbol == "H"])
drop = set()
for k, o in enumerate(O):
    d = xyz[H] - xyz[o]
    d -= L * np.round(d / L)
    close = np.nonzero((np.linalg.norm(d, axis=1) < 1.5) & (Hres != k))[0]
    for c in close:
        if Hres[c] not in drop and k not in drop:
            drop.add(max(k, int(Hres[c])))
print("overlapping water molecules removed:", len(drop))
mod.delete([wat[k] for k in sorted(drop)])
app.PDBFile.writeFile(mod.topology, mod.positions, open("ubq_water.pdb", "w"))
d = typed(mod.topology, mod.positions)
write_extxyz("ubq_water.extxyz", d, "ubiquitin in TIP3P, 0.15 M NaCl")
# the energy of what the files hold: the topology and box of the PDB file (CRYST1 rounds the box to 0.001 A) and the coordinates
# of the extxyz file. The box of the files matters: the water box of addSolvent is not relaxed at the boundary, and the box
# of the Modeller (56.7416 A) against that of the PDB file (56.742 A) changes the nonbonded energy by about 20 kcal/mol
pdb = app.PDBFile("ubq_water.pdb")
lines = open("ubq_water.extxyz").read().splitlines()[2:]
xyz = np.array([[float(x) for x in l.split()[1:4]] for l in lines])
sysm = FF.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=None,
                       rigidWater=False, ewaldErrorTolerance=1e-6, removeCMMotion=False)
for f in sysm.getForces():
    if isinstance(f, mm.NonbondedForce):
        f.setUseDispersionCorrection(False)
e = energies(sysm, xyz)
box = pdb.topology.getPeriodicBoxVectors()
print(len(d["symbols"]), "atoms, box", [round(box[i][i].value_in_unit(unit.angstrom), 3) for i in range(3)], e)
json.dump(e, open(os.path.join("results", "energies_solvated.json"), "w"), indent=1)
