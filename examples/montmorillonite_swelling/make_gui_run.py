"""The inputs of gui_run/ without a window: the force field of dry.extxyz with the bundled ClayFF library and its solution
library (soapff.clayff.prepare, md=False), then the swelling-curve block with the water counts 4 8 12 16 20 25 32 40 48 59
(NPT 30 + 20 ps each): the water types are reserved in system.data / system.settings (soapff_cui.core.reserve_water_types),
water.mol is written and the block of soapff_cui.core.swelling_lines is appended to in.lammps.

    python make_gui_run.py      (with a Python where soapff is installed)  ->  gui_run/{in.lammps, system.data, system.settings, water.mol, ...}"""
import json
import os
import shutil

from soapff_cui import core
from soapff import clayff as runner

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "gui_run")
WATERS = [4, 8, 12, 16, 20, 25, 32, 40, 48, 59]
EQ_PS, AVG_PS = 30.0, 20.0


def type_ids(out):
    """{ClayFF type name: type number} from the Masses of system.data"""
    L = open(os.path.join(out, "system.data")).read().splitlines()
    i = next(k for k, l in enumerate(L) if l.strip() == "Masses") + 2
    ids = {}
    while i < len(L) and L[i].strip():
        if "#" in L[i]:
            ids[L[i].split("#")[1].strip()] = int(L[i].split()[0])
        i += 1
    return ids


if os.path.isdir(OUT):
    shutil.rmtree(OUT)
# the default settings of the run (temperature K, MD lengths ps, pressure atm, barostat, boundary, SOAP threshold, cutoff A,
# dump interval ps, ClayFF charge scheme, water model)
s = dict(temp=300.0, relax_ps=10.0, nvt_ps=50.0, npt_ps=0.0, press=1.0, baro="aniso", boundary="p p p", threshold=0.9,
         cutoff=10.0, dump_ps=1.0, charges="clayff2004", water="spc")
md = {k: s[k] for k in ("temp", "relax_ps", "nvt_ps", "npt_ps", "press", "baro", "boundary", "threshold", "cutoff", "dump_ps",
                        "charges", "water")}
rep = runner.prepare(os.path.join(HERE, "dry.extxyz"), OUT, md=False, **md)
print("fallback", rep.get("fallback"), "type check differ", rep["type_check"]["differ"], "charge sum", round(rep["charge_sum"], 4))
if "ow" not in type_ids(OUT):
    core.reserve_water_types(OUT)
tid = type_ids(OUT)
D = open(os.path.join(OUT, "system.data")).read().splitlines()
nat = int(D[2].split()[0])
k = D.index(next(l for l in D if l.startswith("Atoms")))
atoms = [l.split() for l in D[k + 2:k + 2 + nat]]
n_water0 = sum(1 for a in atoms if int(a[2]) == tid["ow"])
mol_max = max(int(a[1]) for a in atoms)
tilted = any(l.strip().endswith("xy xz yz") and any(abs(float(x)) > 1e-8 for x in l.split()[:3]) for l in D[:20])
S = open(os.path.join(OUT, "system.settings")).read().splitlines()
bond = int(next(l for l in S if l.startswith("bond_coeff") and "water" in l).split()[1])
angle = int(next(l for l in S if l.startswith("angle_coeff") and "HOHw" in l).split()[1])
water = json.load(open(os.path.join(OUT, "report.json"))).get("water_model", "spc")
open(os.path.join(OUT, "water.mol"), "w").write(core.water_molecule(tid, bond, angle, water=water))
dt = float(getattr(core, "TIMESTEP_FS", 0.5))
block = core.swelling_lines(tid, nat - 3 * n_water0, n_water0, mol_max, s["temp"], s["press"], dt, tilted, WATERS, [],
                            eq_ps=EQ_PS, avg_ps=AVG_PS)
inp = os.path.join(OUT, "in.lammps")
txt = open(inp).read().rstrip("\n").splitlines()
open(inp, "w").write("\n".join(txt + block) + "\n")
print("wrote", OUT, "atoms", nat, "types", tid)
