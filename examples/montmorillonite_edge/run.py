"""The optional library "Montmorillonite edge (Underwood & Bourg 2020)" on the montmorillonite particle with edges on every
side (paper: test (7)).  The particle files are the references of the clayff_edge library,
src/soapff/libraries/clayff/edge/ (clay_montmorillonite.data, .settings and clay_montmorillonite_2.data, the particle twice
along c).
  1. the particle (2,292 atoms) and the particle twice (4,584 atoms) as elements, coordinates and box only -> soapff's force
     field (soapff.build("clayff_edge", ...)); types and charges against the files, the pair coefficients against the
     library, the numbers of bonds and angles;
  2. a LAMMPS single point (Ewald 1e-5, cutoff 10 A, special_bonds 0 0 0) of the written system.data / system.settings against
     the original files (the original data file needs a blank line after "Angle Coeffs" to be read by LAMMPS: a copy is fixed);
  3. minimisations (Ewald 1e-4): the particle, the particle twice, and the particle with a water film (1 g/cm3, 2-12 A above it)
     and 4 Na+ + 4 Cl- (packmol): the water and Cl- take ClayFF (flexible SPC), the particle and Na+ the edge model.  The
     particle, the water and the ions are typed together (the references of the library include water and ions), the water
     molecules taken apart from the particle in step 1; the right answer: the particle file, ow/hw, Na and Cl.  A
     particle_water.extxyz kept from packmol is used again.
LAMMPS runs with its Python module, or the lmp executable of $LMP (../clay_md.py).

    python run.py            (soapff and the LAMMPS Python module, pip install "lammps[mpi]", or LMP=/path/to/lmp; packmol)
                              ->  out_*/ (structure.extxyz; soapff's system.data, system.settings, charges.csv, report.json;
                                  the minimisation: its commands in.min, its log log.min), single_point/, results.json"""
import json
import os
import shutil
import sys

import subprocess

import ase.io
import numpy as np
from ase import Atoms

import soapff
from soapff_cui import core, water as W

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import clay_md  # noqa: E402  (examples/clay_md.py: the LAMMPS commands of the clay examples)

REFS = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "edge")
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
    L = open(clay_md.run(s, [EWALD] + clay_md.minimize(), "in.min", "log.min")).read().splitlines()
    k = [i for i, l in enumerate(L) if "Energy initial, next-to-last, final" in l]
    return dict(stop=[l.split("=")[1].strip() for l in L if "Stopping criterion" in l],
                energy_initial_final=[float(L[k[0] + 1].split()[0]), float(L[k[-1] + 1].split()[2])] if k else None)


def single_point(commands, log):
    """the energy terms (kcal/mol) of a LAMMPS run 0 after `commands` (the LAMMPS Python module, or the lmp of $LMP)"""
    keys = ("pe", "ebond", "eangle", "evdwl", "ecoul", "elong")
    cmds = commands + ["thermo_style custom " + " ".join(keys), "thermo_modify format float %.6f", "run 0"]
    try:
        from lammps import lammps
    except ImportError:
        lammps = None
    if lammps is not None:
        lmp = lammps(cmdargs=["-log", log, "-screen", "none"])
        lmp.commands_list(cmds)
        e = {k: round(lmp.get_thermo(k), 6) for k in keys}
        lmp.close()
        return e
    open(log + ".in", "w").write("\n".join(cmds) + "\n")
    subprocess.run([os.environ["LMP"], "-in", log + ".in", "-log", log, "-screen", "none"], check=True)
    L = open(log).read().splitlines()
    k = next(i for i, l in enumerate(L) if l.split()[:1] == ["PotEng"])
    return {key: round(float(v), 6) for key, v in zip(keys, L[k + 1].split())}


res = {}
lib = soapff.load_library("clayff_edge")
lj = {t: [v["epsilon"], v["sigma"]] for t, v in lib.types.items()}
# 1. the particle and the particle twice with the edge library
for name in ("clay_montmorillonite", "clay_montmorillonite_2"):
    f, q, t = bare(name)
    d = os.path.join(HERE, "out_" + name)
    shutil.rmtree(d, ignore_errors=True)
    sy = soapff.build("clayff_edge", f, d, log=quiet)
    shutil.copy(f, os.path.join(d, "structure.extxyz"))
    rows = sy.table()
    D = open(os.path.join(d, "system.data")).read().splitlines()
    S = open(os.path.join(d, "system.settings")).read().splitlines()
    inv = {v: k for k, v in sy.type_ids.items()}
    lj_differ = [l for l in S if l.startswith("pair_coeff") and
                 not np.allclose([float(x) for x in l.split()[3:5]], lj[inv[int(l.split()[1])]], rtol=1e-9)]
    res[name] = dict(atoms=len(rows), wrong_type=sum(r["type"] != x for r, x in zip(rows, t)),
                     wrong_charge=int(sum(abs(float(r["charge"]) - x) > 1e-6 for r, x in zip(rows, q))),
                     min_similarity=round(sy.report["min_similarity_soap"], 4),
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
old = json.load(open(os.path.join(HERE, "results.json"))) if os.path.isfile(os.path.join(HERE, "results.json")) else {}
at = ase.io.read(os.path.join(HERE, "clay_montmorillonite.extxyz"))
f = os.path.join(HERE, "particle_water.extxyz")
placed = old.get("particle_water", {}).get("placed")
if not os.path.isfile(f):
    fr = W.frame_of(at)
    cons = [("inside", "box", [0.0, 0.0, fr["z_atoms"][1] + 2.0, fr["ax"], fr["by"], fr["z_atoms"][1] + 12.0])]
    items = [dict(molecule="water", number=W.estimate(at, cons, density=1.0)["n_water"], constraints=cons),
             dict(molecule="Na", number=4, constraints=cons), dict(molecule="Cl", number=4, constraints=cons)]
    st, rp = W.pack(at, items, seed=1, log=quiet)
    ase.io.write(f, st, format="extxyz")
    placed = {x["molecule"]: x["n"] for x in rp["items"]}
st = ase.io.read(f)
d = os.path.join(HERE, "out_particle_water")
shutil.rmtree(d, ignore_errors=True)
sy = soapff.build("clayff_edge", f, d, log=quiet)
shutil.copy(f, os.path.join(d, "structure.extxyz"))
# the right answer: the particle file for the particle; the molecules added: ow, hw, Na, Cl
_, q0, t0 = bare("clay_montmorillonite")
el = st.get_chemical_symbols()[len(t0):]
added_t = {"O": "ow", "H": "hw", "Na": "Na", "Cl": "Cl"}
added_q = {"O": -0.82, "H": 0.41, "Na": 1.0, "Cl": -1.0}
ans_t = list(t0) + [added_t[e] for e in el]
ans_q = list(q0) + [added_q[e] for e in el]
res["particle_water"] = dict(atoms=len(st), placed=placed, **core.score(sy.types, sy.charges, ans_t, ans_q),
                             water_apart=sy.report["water_apart"], min_similarity=round(sy.report["min_similarity_soap"], 4),
                             charge_sum=round(sy.report["charge_sum"], 4), minimization=minimize(d))
print("particle + water", res["particle_water"], flush=True)
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
