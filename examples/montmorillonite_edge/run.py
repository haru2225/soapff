"""The optional library "Montmorillonite edge (Underwood & Bourg 2020)" on the montmorillonite particle with edges on every
side (paper: test (7)).  The particle files are the references of the clayff_edge library,
src/soapff/libraries/clayff/edge/ (clay_montmorillonite.data, .settings and clay_montmorillonite_2.data, the particle twice
along c).
  1. the particle (2,292 atoms) and the particle twice (4,584 atoms) as elements, coordinates and box only -> soapff's force
     field with the edge library (src/soapff/libraries/clayff/prebuilt/library_clayff_edge.npz); types and charges against
     the files, the pair, bond and angle coefficients against the settings, the numbers of bonds and angles;
  2. a LAMMPS single point (Ewald 1e-5, cutoff 10 A, special_bonds 0 0 0) of the written system.data / system.settings against
     the original files (the original data file needs a blank line after "Angle Coeffs" to be read by LAMMPS: a copy is fixed);
  3. minimisations (Ewald 1e-4): the particle, the particle twice, and the particle with a water film (1 g/cm3, 2-12 A above it)
     and 4 Na+ + 4 Cl- (packmol): the water and Cl- take ClayFF (flexible SPC), the particle and Na+ the edge model.  The
     particle, the water and the ions are typed together with the edge library (its references include water and ions; as
     soapff.build("clayff_edge", ...) does), the water molecules taken apart from the particle in step 1.
LAMMPS runs in this Python, with its Python module (the minimisation: ../clay_md.py).

    python run.py            (with a Python where soapff and LAMMPS are installed: pip install "lammps[mpi]"; packmol)
                              ->  out_*/ (structure.extxyz; soapff's system.data, system.settings, charges.csv, report.json;
                                  the minimisation: its commands in.min, its log log.min), single_point/, results.json"""
import csv
import json
import os
import shutil
import sys

import ase.io
from ase import Atoms
from lammps import lammps

import soapff
from soapff_cui import core, steps, water as W

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import clay_md  # noqa: E402  (examples/clay_md.py: the LAMMPS commands of the clay examples)

REFS = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "edge")
DATA = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "prebuilt")
EDGE = os.path.join(DATA, "library_clayff_edge.npz")
EWALD = "kspace_style ewald 1.0e-4"          # this test: Ewald 1e-4 instead of the PPPM of system.settings
quiet = lambda *a: None


def bare(name):
    """the elements, coordinates and box of a particle file -> <name>.extxyz; returns (path, charges, types)"""
    at, q, t = core.read_lammps_data(os.path.join(REFS, name + ".data"))
    f = os.path.join(HERE, name + ".extxyz")
    ase.io.write(f, Atoms(at.get_chemical_symbols(), positions=at.positions, cell=at.cell, pbc=True), format="extxyz")
    return f, q, t


def minimize(d):
    """the two-stage minimisation (clay_md.minimize) of the force field in `d` with Ewald 1e-4 after system.settings, with the
    LAMMPS Python module (log: log.min); the same commands are written to in.min (lmp -in in.min in the folder)"""
    s = soapff.System(d)
    cmds = [EWALD] + clay_md.minimize()
    clay_md.input_file(s, cmds, "in.min")
    lmp = lammps(cmdargs=["-log", os.path.join(d, "log.min"), "-screen", "none"])
    lmp.commands_list(s.commands() + cmds)
    lmp.close()
    L = open(os.path.join(d, "log.min")).read().splitlines()
    k = [i for i, l in enumerate(L) if "Energy initial, next-to-last, final" in l]
    return dict(stop=[l.split("=")[1].strip() for l in L if "Stopping criterion" in l],
                energy_initial_final=[float(L[k[0] + 1].split()[0]), float(L[k[-1] + 1].split()[2])] if k else None)


def single_point(commands, log):
    """the energy terms (kcal/mol) of a LAMMPS run 0 after `commands`"""
    lmp = lammps(cmdargs=["-log", log, "-screen", "none"])
    lmp.commands_list(commands + ["thermo_style custom pe ebond eangle evdwl ecoul elong", "thermo_modify format float %.6f",
                                  "run 0"])
    e = {k: round(lmp.get_thermo(k), 6) for k in ("pe", "ebond", "eangle", "evdwl", "ecoul", "elong")}
    lmp.close()
    return e


res = {}
ff = steps.load_table(EDGE)["forcefield"]
# 1. the particle and the particle twice with the edge library
for name in ("clay_montmorillonite", "clay_montmorillonite_2"):
    f, q, t = bare(name)
    d = os.path.join(HERE, "out_" + name)
    shutil.rmtree(d, ignore_errors=True)
    rep = core.make(f, d, EDGE, strict=True, log=quiet)
    shutil.copy(f, os.path.join(d, "structure.extxyz"))
    rows = list(csv.DictReader(open(os.path.join(d, "charges.csv"))))
    D = open(os.path.join(d, "system.data")).read().splitlines()
    S = open(os.path.join(d, "system.settings")).read().splitlines()
    tid = {l.split("#")[1].strip(): int(l.split()[0]) for l in core._section(D, "Masses", raw=True)}
    inv = {v: k for k, v in tid.items()}
    lj_differ = [l for l in S if l.startswith("pair_coeff") and [float(x) for x in l.split()[3:5]] != ff["lj"][inv[int(l.split()[1])]]]
    res[name] = dict(atoms=len(rows), wrong_type=sum(r["clayff_type"] != x for r, x in zip(rows, t)),
                     wrong_charge=int(sum(abs(float(r["charge"]) - x) > 1e-6 for r, x in zip(rows, q))),
                     min_similarity=round(min(float(r["similarity"]) for r in rows), 4), type_check_differ=rep["type_check"]["differ"],
                     pair_coeff_differ=lj_differ, bonds=next(l for l in D if l.endswith(" bonds")),
                     angles=next(l for l in D if l.endswith(" angles")),
                     bonded_coeffs=[l for l in S if l.startswith(("bond_coeff", "angle_coeff"))])
    res[name]["minimization"] = minimize(d)
    print(name, res[name], flush=True)
# 2. single point: the original files against the written ones (system.settings with Ewald 1e-5 instead of PPPM)
sp = os.path.join(HERE, "single_point")
shutil.rmtree(sp, ignore_errors=True)
os.makedirs(sp)
s = open(os.path.join(REFS, "clay_montmorillonite.data")).read().replace("Angle Coeffs # harmonic\n1 ", "Angle Coeffs # harmonic\n\n1 ", 1)
open(os.path.join(sp, "original.data"), "w").write(s)
shutil.copy(os.path.join(REFS, "clay_montmorillonite.settings"), os.path.join(sp, "original.settings"))
original = ["units real", "atom_style full", "bond_style harmonic", "angle_style harmonic", "pair_style lj/cut/coul/long 10.0",
            "special_bonds lj/coul 0.0 0.0 0.0", f'read_data "{os.path.join(sp, "original.data")}"',
            f'include "{os.path.join(sp, "original.settings")}"']
written = [c.replace("kspace_style pppm 1.0e-5", "kspace_style ewald 1.0e-5")
           for c in soapff.System(os.path.join(HERE, "out_clay_montmorillonite")).commands()]
res["single_point"] = {x: single_point(c, os.path.join(sp, "log." + x)) for x, c in (("original", original), ("written", written))}
print("single point", res["single_point"], flush=True)
# 3. the particle with a water film and NaCl
at = ase.io.read(os.path.join(HERE, "clay_montmorillonite.extxyz"))
fr = W.frame_of(at)
cons = [("inside", "box", [0.0, 0.0, fr["z_atoms"][1] + 2.0, fr["ax"], fr["by"], fr["z_atoms"][1] + 12.0])]
items = [dict(molecule="water", number=W.estimate(at, cons, density=1.0)["n_water"], constraints=cons),
         dict(molecule="Na", number=4, constraints=cons), dict(molecule="Cl", number=4, constraints=cons)]
s, rp = W.pack(at, items, seed=1, log=quiet)
f = os.path.join(HERE, "particle_water.extxyz")
ase.io.write(f, s, format="extxyz")
d = os.path.join(HERE, "out_particle_water")
shutil.rmtree(d, ignore_errors=True)
rep = core.make(f, d, EDGE, strict=True, log=quiet)
shutil.copy(f, os.path.join(d, "structure.extxyz"))
res["particle_water"] = dict(atoms=rep["atoms"], placed={x["molecule"]: x["n"] for x in rp["items"]}, roles=rep["roles"],
                             type_check_differ=rep["type_check"]["differ"], water_apart=rep["water_apart"],
                             min_similarity=round(rep["min_similarity"], 4),
                             charge_sum=round(rep["charge_sum"], 4), minimization=minimize(d))
print("particle + water", res["particle_water"], flush=True)
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
