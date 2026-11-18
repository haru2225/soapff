"""The run of run/: the force field of dry.extxyz with the bundled clayff library (soapff.clayff.prepare), the water types
reserved in system.data / system.settings (soapff_cui.core.reserve_water_types), water.mol, then the minimisation of ../clay_md.py and the swelling curve of swelling.py with the water counts
4 8 12 16 20 25 32 40 48 59 (NPT 30 + 20 ps each).  The commands are written to in.lammps (include system.settings; the
input of the stored run, made on another computer with lmp -in in.lammps) and, with --run, given to the LAMMPS Python
module here (pip install "lammps[mpi]").

    python make_run.py [--out run] [--run]
        -> run/{system.data, system.settings, charges.csv, report.json, water.mol, in.lammps}  (the folder is replaced)
           --run: also log.lammps, swelling.dat, swelling_states.lammpstrj"""
import argparse
import json
import os
import shutil
import sys

import soapff
from soapff import clayff as runner
from soapff_cui import core
from soapff_cui._clayff.ff import ELEMENT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))                 # clay_md.py
import clay_md  # noqa: E402
import swelling  # noqa: E402

WATERS = [4, 8, 12, 16, 20, 25, 32, 40, 48, 59]
EQ_PS, AVG_PS = 30.0, 20.0
TEMP, PRESS = 300.0, 1.0                                     # K, atm


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


ap = argparse.ArgumentParser(description="the swelling curve of the Na montmorillonite of Cygan et al. (2004)")
ap.add_argument("--out", default=os.path.join(HERE, "run"), help="output folder (replaced; default run)")
ap.add_argument("--run", action="store_true", help="run it here with the LAMMPS Python module")
a = ap.parse_args()
OUT = os.path.abspath(a.out)
if os.path.isdir(OUT):
    shutil.rmtree(OUT)
# the force field: boundary, SOAP threshold, cutoff (A), ClayFF charge scheme, water model
rep = runner.prepare(os.path.join(HERE, "dry.extxyz"), OUT, boundary="p p p", threshold=0.9, cutoff=10.0,
                     charges="clayff2004", water="spc")
print("min similarity", round(rep["min_similarity"], 4), "type check differ", rep["type_check"]["differ"], "charge sum",
      round(rep["charge_sum"], 4))
if "ow" not in type_ids(OUT):
    core.reserve_water_types(OUT)
tid = type_ids(OUT)
# report.json: the type numbers after the water types were reserved (read by soapff.System)
rep = json.load(open(os.path.join(OUT, "report.json")))
rep["type_ids"] = tid
rep["type_elements"] = [ELEMENT[t] for t, _ in sorted(tid.items(), key=lambda kv: kv[1])]
json.dump(rep, open(os.path.join(OUT, "report.json"), "w"), indent=1)
D = open(os.path.join(OUT, "system.data")).read().splitlines()
nat = int(D[2].split()[0])
k = D.index(next(l for l in D if l.startswith("Atoms")))
atoms = [l.split() for l in D[k + 2:k + 2 + nat]]
n_water0 = sum(1 for x in atoms if int(x[2]) == tid["ow"])
mol_max = max(int(x[1]) for x in atoms)
tilted = any(l.strip().endswith("xy xz yz") and any(abs(float(x)) > 1e-8 for x in l.split()[:3]) for l in D[:20])
S = open(os.path.join(OUT, "system.settings")).read().splitlines()
bond = int(next(l for l in S if l.startswith("bond_coeff") and "water" in l).split()[1])
angle = int(next(l for l in S if l.startswith("angle_coeff") and "HOHw" in l).split()[1])
open(os.path.join(OUT, "water.mol"), "w").write(core.water_molecule(tid, bond, angle, water=rep.get("water_model", "spc")))

s = soapff.System(OUT)
dt = float(s.report.get("timestep_fs", clay_md.TIMESTEP_FS))
args = (tid, nat - 3 * n_water0, n_water0, mol_max, TEMP, PRESS, dt, tilted, WATERS, [])
clay_md.input_file(s, clay_md.minimize() + swelling.swelling_lines(*args, eq_ps=EQ_PS, avg_ps=AVG_PS, out=OUT))
print("wrote", OUT, "atoms", nat, "types", tid)
if a.run:
    lmp = clay_md.start(s)                                    # the force field; log in OUT/log.lammps
    lmp.commands_list(clay_md.minimize())
    swelling.run(lmp, *args, eq_ps=EQ_PS, avg_ps=AVG_PS, out=OUT)
    lmp.close()
    print("swelling curve in", os.path.join(OUT, "swelling.dat"), "(python analyze_swelling.py)")
