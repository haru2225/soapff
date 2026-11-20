"""Systems with water and aqueous ions: structures made with soapff (packmol, soapff_cui.water) and its force field (flexible
SPC water, ClayFF Table 1 ions), then a 300 K MD with LAMMPS.  The packmol liquids (refs/<liquid>/system.data) are references
of the clayff library (src/soapff/libraries/clayff/water/); the last MD frames (refs/<system>/final.data) are the test
frames of ../clayff/water_md_frames.py (none is a reference).

  water/           liquid water, 512 molecules in a 24.83 A cube
  NaCl/ KCl/ CsCl/ 8 cations + 8 Cl- in 496 water molecules (about 0.9 mol/L)
  CaCl2/ BaCl2/    6 cations + 12 Cl- in 494 water molecules (about 0.7 mol/L)
  if_<mineral>/    the library's basal (0 0 1) slab (src/soapff/libraries/clayff/basal/) with the vacuum above it filled with
                   water (1 g/cm3, packmol seed 2; the test structures of ../interfaces_001_water/ use seed 1)
  hydrate_2W/      the 2-layer Na-montmorillonite hydrate of ../montmorillonite_2W/ (hydrate_2W_start.data, the
                   system.data of its diffusion/start/ run)

MD (the commands of ../clay_md.py, the same as the in.lammps soapff wrote before): minimisation, 5 ps NVT relaxation, then
the liquids and the hydrate 5 ps NVT + 20 ps NPT at 1 atm (liquids isotropic, hydrate full cell), the interfaces 25 ps NVT;
0.5 fs, frames every 1 ps.  The same commands are written to refs/<system>/in.lammps (include system.settings), for an lmp
executable.  The stored systems were typed by an earlier soapff, with the water and the ions by rule; the bundled clayff
library (whose references include these liquids) gives them the same types and charges.

    python build_refs.py [--jobs 3] [system ...]       the force field and the MD here, with the LAMMPS Python module
                                                       (pip install "lammps[mpi]"), `jobs` systems at a time
    python build_refs.py --inputs-only [system ...]    only the inputs (system.data, system.settings, in.lammps), for a run
                                                       with an lmp executable"""
import argparse
import concurrent.futures
import os
import subprocess
import sys
import tempfile

import numpy as np

import ase.io
from ase import Atoms

import soapff
from soapff_cui import core, water

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))                 # clay_md.py
import clay_md  # noqa: E402

REFS = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff")
L = 24.83                                                    # cube edge (A): 512 water molecules at 1 g/cm3
SOLUTIONS = {"water": (512, []), "NaCl": (496, [("Na", 8), ("Cl", 8)]), "KCl": (496, [("K", 8), ("Cl", 8)]),
             "CsCl": (496, [("Cs", 8), ("Cl", 8)]), "CaCl2": (494, [("Ca", 6), ("Cl", 12)]),
             "BaCl2": (494, [("Ba", 6), ("Cl", 12)])}
MINERALS = ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
            "pyrophyllite_2", "kaolinite", "lizerdite", "talc", "brucite"]
quiet = lambda *a: None


def liquid(nw, ions, seed=2, tolerance=2.0):
    """water (+ ions) filling the periodic cube: packmol directly, periodic in x, y and z (water.pack needs a structure to place
    them around); the water geometry of soapff_cui.water (O-H 1.0 A, H-O-H 109.47 deg)"""
    h = np.radians(water.HOH / 2)
    mol = {"water": (["O", "H", "H"], [[0, 0, 0], [water.OH * np.sin(h), 0, water.OH * np.cos(h)],
                                       [-water.OH * np.sin(h), 0, water.OH * np.cos(h)]])}
    with tempfile.TemporaryDirectory() as d:
        inp = [f"tolerance {tolerance}", "filetype xyz", "output packed.xyz", f"seed {seed}", f"pbc 0. 0. 0. {L} {L} {L}"]
        for m, n in [("water", nw)] + list(ions):
            sym, pos = mol.get(m, ([m], [[0, 0, 0]]))
            water._xyz(os.path.join(d, m + ".xyz"), sym, np.array(pos, float))
            inp += ["", f"structure {m}.xyz", f"  number {n}", f"  inside box 0. 0. 0. {L} {L} {L}", "end structure"]
        open(os.path.join(d, "pack.inp"), "w").write("\n".join(inp) + "\n")
        subprocess.run("packmol < pack.inp > packmol.log", shell=True, cwd=d, check=True)
        if "Success!" not in open(os.path.join(d, "packmol.log")).read():
            raise RuntimeError(f"packmol did not converge for {ions}")
        sym, pos = water._read_xyz(os.path.join(d, "packed.xyz"))
    return Atoms(sym, positions=pos, cell=[L, L, L], pbc=True)


def structure(name):
    """the initial structure of a reference system and its MD settings (clay_md.md)"""
    if name in SOLUTIONS:
        return liquid(*SOLUTIONS[name]), dict(relax_ps=5, nvt_ps=5, npt_ps=20, baro="iso")
    if name == "hydrate_2W":
        return core.read_structure(os.path.join(HERE, "hydrate_2W_start.data")), \
            dict(relax_ps=5, nvt_ps=5, npt_ps=20, baro="tri")
    at, q, types = core.read_lammps_data(os.path.join(REFS, "basal", name[3:] + "_001.data"))
    gap = water.gaps(at, min_gap=8.0, density=1.0)[0]                 # the vacuum above the slab
    s, rep = water.add_water(at, [(gap, gap["n_water"])], seed=2, log=quiet)
    return s, dict(relax_ps=5, nvt_ps=25, npt_ps=0)


def run(name, inputs_only=False):
    """refs/<name>/: the structure, its force field (the bundled clayff library), in.lammps, and the MD with the LAMMPS Python
    module (log.lammps, traj.lammpstrj, msd_*.dat, final.data)"""
    d = os.path.join(HERE, "refs", name)
    os.makedirs(d, exist_ok=True)
    at, md = structure(name)
    ase.io.write(os.path.join(d, "structure.extxyz"), at, format="extxyz")
    s = soapff.build("clayff", os.path.join(d, "structure.extxyz"), d, log=quiet)
    r = s.report
    commands = clay_md.minimize() + clay_md.md(s, temp=300.0, dump_ps=1.0, **md)
    clay_md.input_file(s, commands)                           # the same commands for lmp -in in.lammps
    if inputs_only:
        return name, r["atoms"], "inputs written"
    lmp = clay_md.start(s)                                    # the force field; log in refs/<name>/log.lammps
    try:
        lmp.commands_list(commands)
        status = "done"
    except Exception as e:                                    # a LAMMPS error stops this system only
        status = f"LAMMPS error: {e}"
    finally:
        lmp.close()
    return name, r["atoms"], status


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="the reference systems with water and their 300 K MD")
    ap.add_argument("systems", nargs="*", help="default: all")
    ap.add_argument("--jobs", type=int, default=3, help="systems run at the same time (one LAMMPS process each)")
    ap.add_argument("--inputs-only", action="store_true", help="only the inputs, for a run with an lmp executable")
    a = ap.parse_args()
    names = a.systems or list(SOLUTIONS) + ["hydrate_2W"] + ["if_" + m for m in MINERALS]
    with concurrent.futures.ProcessPoolExecutor(1 if a.inputs_only else a.jobs) as ex:
        for name, n, status in ex.map(run, names, [a.inputs_only] * len(names)):
            print(name, n, "atoms,", status, flush=True)
