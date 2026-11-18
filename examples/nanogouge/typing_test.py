"""Many edge-bearing particles in one cell (paper: Table "clay test cases" and Fig. "nanogouge"): the initial state of a
nano-gouge (Underwood & Bourg 2020), typed by soapff from its elements, coordinates and box alone.

build_gouge.py places N montmorillonite particles with edges on every side (the particle of ../montmorillonite_edge,
2,292 atoms) at random positions and orientations in a periodic box (no two atoms of two particles nearer than 3.5 A; of
300 random configurations the most random one, the lowest order parameter S) and fills the box with flexible SPC water at
a water/clay mass ratio of 3 (Underwood & Bourg, J. Phys. Chem. C 124, 3702 (2020)).  Its data file carries the types of the
particle (those of the particle's files) and ow, hw.  Typed here:
  wet      the initial state as built (particles, their Na+, water).  The water is put on a jittered grid with random
           orientations, and a few molecules overlap (an H 1.1-1.25 A from the O of another molecule): soapff bonds them
           (H-O bond length 1.25 A) and stops on these atoms, whose bonds are not those of their type;
  wet_min  the same after a minimisation with the builder's force field (../clay_md.py: minimize; PPPM 1e-4, cutoff
           10 A), the first step of an MD of the nano-gouge.  Water O come nearer than the Al-O and Mg-O
           bond lengths to edge Al and Mg (from 1.92 and 2.19 A): soapff takes the water molecules apart from the particles
           (water_apart in the results: the bonds taken out), so that they are not bonded to these metals;
  dry      the same particles without the water: the edges in vacuum, next to the edges and faces of other particles.
Each is typed with the edge library (soapff.build("clayff_edge", ...), strict) and compared with the types and charges of the
data file.  For the smallest N also a LAMMPS single point (Ewald 1e-5, cutoff 10 A, special_bonds 0 0 0) of soapff's force
field against the builder's data and settings files (wet_min: the minimised data file; dry: the water deleted in LAMMPS).

    python typing_test.py WORK [N ...]      (default N: 3 12; soapff and the LAMMPS Python module)
        -> typing_results.json; WORK/gouge_N<N>.*: the builder's files, WORK/out_N<N>_{wet,wet_min,dry}/: soapff's files"""
import csv
import json
import os
import subprocess
import sys
import time

import ase.io
import numpy as np
from ase import Atoms
from lammps import lammps

import soapff
from soapff_cui import core

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import clay_md  # noqa: E402  (examples/clay_md.py: the minimisation of the clay examples)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "typing_results.json")
WORK = os.path.abspath(sys.argv[1])
NS = [int(x) for x in sys.argv[2:]] or [3, 12]
os.makedirs(WORK, exist_ok=True)
alias = lambda t: np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t])  # noqa: E731
quiet = lambda *a: None  # noqa: E731
WATER = ("ow", "hw")


def build(n):
    """build_gouge.py -n N --seed 1 (300 configurations, m_w/m_c 3) -> WORK/gouge_N<n> (.data, .settings, .json)"""
    out = os.path.join(WORK, f"gouge_N{n}")
    if not os.path.isfile(out + ".data"):
        subprocess.run([sys.executable, os.path.join(HERE, "build_gouge.py"), "-n", str(n), "--seed", "1", "-o", out],
                       check=True)
    return out


def typed(at, truth, q0, name):
    """at (elements, coordinates, box) -> soapff.build('clayff_edge') in WORK/out_<name>; compared with truth and q0"""
    f = os.path.join(WORK, name + ".extxyz")
    a = Atoms(at.get_chemical_symbols(), positions=at.positions, cell=at.cell, pbc=True)
    a.wrap()
    ase.io.write(f, a, format="extxyz")
    d = os.path.join(WORK, "out_" + name)
    t0 = time.time()
    try:
        soapff.build("clayff_edge", f, d, log=quiet)
    except soapff.TypingError as e:
        return dict(atoms=len(at), stopped=str(e)[:600], seconds=round(time.time() - t0, 1))
    sec = time.time() - t0
    rows = list(csv.DictReader(open(os.path.join(d, "charges.csv"))))
    got = alias([r["clayff_type"] for r in rows])
    q = np.array([float(r["charge"]) for r in rows])
    sim = np.array([float(r["similarity"]) for r in rows])
    wt, wq = got != truth, np.abs(q - q0) > 1e-6
    water = np.isin(truth, WATER)
    res = dict(atoms=len(at), water_molecules=int((truth == "ow").sum()), wrong_type=int(wt.sum()), wrong_charge=int(wq.sum()),
               min_similarity=round(float(sim.min()), 4),
               min_similarity_particles=round(float(sim[~water].min()), 4),
               min_similarity_water=round(float(sim[water].min()), 4) if water.any() else None,
               charge_sum=round(float(q.sum()), 4), seconds=round(sec, 1),
               water_apart=json.load(open(os.path.join(d, "report.json"))).get("water_apart"))
    if wt.any():
        res["confused"] = sorted({f"{x}->{y}" for x, y in zip(truth[wt], got[wt])})
    return res


def builder_commands(data, settings):
    """the builder's force field (as the nano-gouge's LAMMPS input, cutoff 10 A) without the k-space"""
    return ["units real", "atom_style full", "boundary p p p", "bond_style harmonic", "angle_style harmonic",
            "pair_style lj/cut/coul/long 10.0", "pair_modify mix arithmetic tail no", "special_bonds lj/coul 0.0 0.0 0.0",
            f'read_data "{data}"', f'include "{settings}"']


def minimised(g, at):
    """the builder's structure minimised with its force field (PPPM 1e-4) -> (Atoms: the minimised positions in the order of
    the atom ids, the same box; WORK/gouge_N<n>_min.data)"""
    out = g + "_min.data"
    if not os.path.isfile(out):
        lmp = lammps(cmdargs=["-log", g + "_min.log", "-screen", "none"])
        lmp.commands_list(builder_commands(g + ".data", g + ".settings") + ["kspace_style pppm 1.0e-4", "neighbor 2.0 bin",
                          "neigh_modify delay 0 every 1 check yes"] + clay_md.minimize() + [f'write_data "{out}" nocoeff'])
        lmp.close()
    m, _, _ = core.read_lammps_data(out)
    a = at.copy()
    a.positions = m.positions
    return a


def single_point(commands, log):
    """the energy terms (kcal/mol) of a LAMMPS run 0 after `commands`"""
    lmp = lammps(cmdargs=["-log", log, "-screen", "none"])
    lmp.commands_list(commands + ["thermo_style custom pe ebond eangle evdwl ecoul elong", "thermo_modify format float %.6f",
                                  "run 0"])
    e = {k: round(lmp.get_thermo(k), 6) for k in ("pe", "ebond", "eangle", "evdwl", "ecoul", "elong")}
    lmp.close()
    return e


res = json.load(open(OUT)) if os.path.isfile(OUT) else {}
for n in NS:
    g = build(n)
    at, q0, t0 = core.read_lammps_data(g + ".data")
    truth, q0 = alias(t0), np.asarray(q0, float)
    meta = json.load(open(g + ".json"))
    r = dict(particles=n, S=round(meta["S_initial"], 3), box=[round(x, 2) for x in meta["box"]], mw_mc=meta["mwmc"])
    r["wet"] = typed(at, truth, q0, f"N{n}_wet")
    print(n, "wet", r["wet"], flush=True)
    r["wet_min"] = typed(minimised(g, at), truth, q0, f"N{n}_wet_min")
    print(n, "wet_min", r["wet_min"], flush=True)
    dry = ~np.isin(truth, WATER)
    r["dry"] = typed(at[dry], truth[dry], q0[dry], f"N{n}_dry")
    print(n, "dry", r["dry"], flush=True)
    if n == min(NS):
        sp = {}
        for kind, built in (("wet_min", builder_commands(g + "_min.data", g + ".settings")),
                            ("dry", builder_commands(g + ".data", g + ".settings")
                             + [f"group water type {meta['type_hw']} {meta['type_ow']}", "delete_atoms group water bond yes"])):
            if "stopped" not in r[kind]:
                written = [c.replace("kspace_style pppm 1.0e-5", "kspace_style ewald 1.0e-5")
                           for c in soapff.System(os.path.join(WORK, f"out_N{n}_{kind}")).commands()]
                sp[kind] = {x: single_point(c, os.path.join(WORK, f"log.sp_N{n}_{kind}_{x}"))
                            for x, c in (("builder", built + ["kspace_style ewald 1.0e-5"]), ("soapff", written))}
                print(n, kind, "single point", sp[kind], flush=True)
        r["single_point"] = sp
    res[f"N{n}"] = r
    json.dump(res, open(OUT, "w"), indent=1)
