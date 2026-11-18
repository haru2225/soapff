"""Short MD with the force fields soapff writes, run in LAMMPS from Python (the LAMMPS Python module): TIP3P and SPC/E
boxes of 216 waters (rigid: SHAKE, part of the force field), NVT then NPT; and kaolinite with water (ClayFF), NVT.

    python examples/md_check.py WORKDIR      (from the root of the repository) -> WORKDIR/md.json"""
import json
import os
import sys
import time

import numpy as np
from ase import Atoms
from ase.build import molecule
from lammps import lammps
from lammps.formats import LogFile

import soapff

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clay_md  # noqa: E402

W = os.path.abspath(sys.argv[1])
os.makedirs(W, exist_ok=True)
out = {}
rng = np.random.default_rng(7)
wat = molecule("H2O")
box = Atoms(cell=[18.6] * 3, pbc=True)
for ix in range(6):
    for iy in range(6):
        for iz in range(6):
            m = wat.copy()
            m.rotate(rng.uniform(0, 360), rng.normal(size=3))
            m.translate((np.array([ix, iy, iz]) + 0.5) * 3.1 - m.positions.mean(0))
            box += m
box.write(os.path.join(W, "water216.extxyz"))

for lib in ("tip3p", "spce"):
    t0 = time.time()
    s = soapff.build(lib, os.path.join(W, "water216.extxyz"), os.path.join(W, lib), log=lambda *a: None)
    log = os.path.join(s.out_dir, "log.lammps")
    lmp = lammps(cmdargs=["-log", log, "-screen", "none"])
    lmp.commands_list(s.commands())          # the force field, with SHAKE on the water (fix rigid)
    lmp.commands_string("""
        thermo_style custom step temp press pe evdwl ecoul elong vol density
        thermo 100
        minimize 1.0e-4 1.0e-6 1000 10000
        reset_timestep 0
        velocity all create 300 4928459 mom yes rot no dist gaussian
        timestep 1.0
        fix 1 all nvt temp 300 300 100
        run 2000
        reset_timestep 0
        run 5000
        unfix 1
        fix 1 all npt temp 300 300 100 iso 1 1 1000
        run 20000
    """)
    lmp.close()
    npt = LogFile(log).runs[-1]                # the NPT: 20 ps, thermo every 100 fs
    half = len(npt["Step"]) // 2
    out[lib] = dict(seconds=round(time.time() - t0, 1), types=sorted(set(s.types)), steps=int(npt["Step"][-1]),
                    temp_mean_npt=round(float(np.mean(npt["Temp"][half:])), 1),
                    density_g_cm3=round(float(np.mean(npt["Density"][half:])), 3))
    print(lib, out[lib], flush=True)

t0 = time.time()
s = soapff.build("clayff", "examples/clayff/kaolinite_water.extxyz", os.path.join(W, "kaolinite_water"), log=lambda *a: None)
lmp = clay_md.start(s)                           # LAMMPS with the force field
lmp.commands_list(clay_md.minimize() + clay_md.md(s, temp=300.0, relax_ps=0.5, nvt_ps=2.0))
lmp.close()
nvt = LogFile(os.path.join(s.out_dir, "log.lammps")).runs[-1]
out["clayff_kaolinite_water"] = dict(seconds=round(time.time() - t0, 1), atoms=s.report.get("atoms"), steps=int(nvt["Step"][-1]),
                                     temp_mean=round(float(np.mean(nvt["Temp"][len(nvt["Temp"]) // 2:])), 1),
                                     final_data=os.path.isfile(os.path.join(s.out_dir, "final.data")))
print("clayff", out["clayff_kaolinite_water"], flush=True)
json.dump(out, open(os.path.join(W, "md.json"), "w"), indent=1)
