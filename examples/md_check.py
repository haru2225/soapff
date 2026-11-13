"""short LAMMPS runs of what soapff writes: TIP3P and SPC/E water boxes (NVT then NPT, SHAKE) and kaolinite / water (ClayFF).
Ewald instead of PPPM: the local Homebrew LAMMPS crashes when PPPM is set up a second time (FFTW3 with threads).
    python examples/md_check.py WORKDIR      (from the root of the repository)"""
import json, os, sys, time
import numpy as np
from ase import Atoms
from ase.build import molecule
import soapff
W = sys.argv[1]
out = {}
rng = np.random.default_rng(7)
wat = molecule("H2O")
box = Atoms(cell=[18.6] * 3, pbc=True)
for ix in range(6):
    for iy in range(6):
        for iz in range(6):
            m = wat.copy(); m.rotate(rng.uniform(0, 360), rng.normal(size=3))
            m.translate((np.array([ix, iy, iz]) + 0.5) * 3.1 - m.positions.mean(0)); box += m
box.write(os.path.join(W, "water216.extxyz"))
for lib in ("tip3p", "spce"):
    t0 = time.time()
    s = soapff.build(lib, os.path.join(W, "water216.extxyz"), os.path.join(W, lib), relax_ps=2.0, nvt_ps=5.0, npt_ps=20.0,
                     baro="iso", dump_ps=1.0, kspace="ewald 1.0e-5", log=lambda *a: None)
    s.run(lmp="lmp_serial", echo=False)
    th = s.thermo()
    vol = th.get("Volume")
    n = len(s.types) // 3
    dens = (n * 18.015 / 6.02214076e23) / (vol[len(vol) // 2:].mean() * 1e-24) if vol is not None else None
    out[lib] = dict(seconds=round(time.time() - t0, 1), types=sorted(set(s.types)), steps=int(th["Step"][-1]),
                    temp_mean_npt=round(float(th["Temp"][len(th["Temp"]) // 2:].mean()), 1),
                    density_g_cm3=round(float(dens), 3) if dens else None, columns=list(th))
    print(lib, out[lib], flush=True)
t0 = time.time()
s = soapff.build("clayff", "examples/clayff/kaolinite_water.extxyz", os.path.join(W, "kaolinite_water"),
                 relax_ps=0.5, nvt_ps=2.0, kspace="ewald 1.0e-5", log=lambda *a: None)
s.run(lmp="lmp_serial", echo=False)
th = s.thermo()
out["clayff_kaolinite_water"] = dict(seconds=round(time.time() - t0, 1), atoms=s.report.get("atoms"), steps=int(th["Step"][-1]),
                                     temp_mean=round(float(th["Temp"][len(th["Temp"]) // 2:].mean()), 1),
                                     final_data=os.path.isfile(os.path.join(s.out_dir, "final.data")))
print("clayff", out["clayff_kaolinite_water"], flush=True)
json.dump(out, open(os.path.join(W, "md.json"), "w"), indent=1)
