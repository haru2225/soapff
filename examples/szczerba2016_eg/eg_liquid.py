"""Liquid ethylene glycol with the OPLS-AA of soapff (oplsaa_eg), as Table 2 of Szczerba & Kalinichev (2016): a box of 230
EG molecules (packmol, 1.0 g/cm3), typed from its elements, positions and cell, then in LAMMPS: minimisation, 0.1 ns NVT at
293 K and 0.4 ns NPT at 293 K and 1 bar; the density is the mean of the last 0.3 ns (it levels off within
0.05 ns of NPT).

Two settings of the interactions: liquid/, those soapff writes for oplsaa (LJ cutoff 12 A with the tail correction);
liquid_8.5/, those of Szczerba & Kalinichev (cutoff 8.5 A, Ewald for the Coulomb; no tail correction, the default of
LAMMPS; clayff_oplsaa_eg.txt has the same cutoff, with Lorentz-Berthelot mixing).

    python eg_liquid.py build  (soapff)  ->  liquid/structure.extxyz, liquid/ff, liquid/in.md
    python eg_liquid.py build85          ->  liquid_8.5/ff, liquid_8.5/in.md (liquid/ff with cutoff 8.5 A, no tail)
    cd liquid && mpirun -np 2 lmp -in in.md          ->  density.dat (the same in liquid_8.5)
    python eg_liquid.py result                      ->  results/eg_liquid.json"""
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from ase.io import read, write

HERE = os.path.dirname(os.path.abspath(__file__))
N = 230
D = os.path.join(HERE, "liquid")
D85 = os.path.join(HERE, "liquid_8.5")

IN = """# liquid ethylene glycol, OPLS-AA (soapff oplsaa_eg), 293 K, 1 bar
include ff/system.settings
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 1000
thermo_style custom step temp press pe density
min_style cg
minimize 1.0e-4 1.0e-6 2000 20000
reset_timestep 0
velocity all create 293 5 dist gaussian mom yes rot yes
timestep 1.0
fix 1 all nvt temp 293 293 100
run 100000
unfix 1
reset_timestep 0
fix 2 all npt temp 293 293 100 iso 1 1 1000
variable rho equal density
fix 3 all ave/time 100 10 1000 v_rho file density.dat
run 400000
write_data final.data
"""


def build():
    import soapff
    eg = read(os.path.join(HERE, "refs", "EG.extxyz"))
    L = (N * 62.068 / 0.60221 / 1.0) ** (1 / 3)
    os.makedirs(D, exist_ok=True)
    with tempfile.TemporaryDirectory() as d:
        write(f"{d}/eg.xyz", eg)
        open(f"{d}/in.inp", "w").write(f"tolerance 2.0\nfiletype xyz\noutput {d}/out.xyz\nseed 3\npbc 0 0 0 {L:.3f} {L:.3f} {L:.3f}\n"
                                       f"structure {d}/eg.xyz\n  number {N}\n  inside box 0 0 0 {L:.3f} {L:.3f} {L:.3f}\nend structure\n")
        r = subprocess.run(f"packmol < {d}/in.inp", shell=True, capture_output=True, text=True)
        assert "Success" in r.stdout, r.stdout[-1500:]
        at = read(f"{d}/out.xyz")
    at.set_cell([L, L, L])
    at.pbc = True
    at.set_array("fftype", np.array(list(eg.arrays["fftype"]) * N))
    at.set_array("ffcharge", np.tile(eg.arrays["ffcharge"], N))
    write(os.path.join(D, "structure.extxyz"), at)
    s = soapff.build(os.path.join(HERE, "oplsaa_eg.txt"), at, os.path.join(D, "ff"), log=lambda *a: None)
    right = bool((np.array(s.types) == at.arrays["fftype"]).all())
    txt = open(os.path.join(D, "ff", "system.settings")).read().replace("read_data system.data", "read_data ff/system.data")
    open(os.path.join(D, "ff", "system.settings"), "w").write(txt)
    open(os.path.join(D, "in.md"), "w").write(IN)
    print(f"{len(at)} atoms, box {L:.2f} A; every atom typed as Foyer: {right}; -> {D}/in.md")


def build85():
    """liquid_8.5: the force field of liquid/ff with the cutoff and tail of Szczerba & Kalinichev"""
    import shutil
    os.makedirs(D85, exist_ok=True)
    if os.path.isdir(os.path.join(D85, "ff")):
        shutil.rmtree(os.path.join(D85, "ff"))
    shutil.copytree(os.path.join(D, "ff"), os.path.join(D85, "ff"))
    path = os.path.join(D85, "ff", "system.settings")
    txt = open(path).read()
    for a, b in (("pair_style lj/cut/coul/long 12.0", "pair_style lj/cut/coul/long 8.5"),
                 ("pair_modify mix geometric tail yes", "pair_modify mix geometric tail no")):
        assert txt.count(a) == 1, a
        txt = txt.replace(a, b)
    open(path, "w").write(txt)
    open(os.path.join(D85, "in.md"), "w").write(IN)
    print(f"-> {D85}/in.md")


def result():
    out = dict(molecules=N, temperature_K=293, pressure_bar=1, average_over_ns=0.3,
               szczerba2016_table2=dict(oplsaa=1.060, experiment=1.115))
    for key, d in (("cutoff12_tail", D), ("cutoff8.5_notail", D85)):
        if not os.path.isfile(os.path.join(d, "density.dat")):
            continue
        t, rho = np.loadtxt(os.path.join(d, "density.dat"), unpack=True)
        last = rho[t > t.max() - 300000]
        m = [b.mean() for b in np.array_split(last, 5)]
        out[key] = dict(density_g_cm3=float(last.mean()), stderr=float(np.std(m, ddof=1) / np.sqrt(5)), npt_ns=t.max() / 1e6)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "results", "eg_liquid.json"), "w"), indent=1)
    print(out)


if __name__ == "__main__":
    {"build": build, "build85": build85, "result": result}[sys.argv[1]]()
