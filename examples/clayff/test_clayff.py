"""Test case of the clayff library (the clay minerals of the paper): structures that are not the references get, atom by atom,
the ClayFF type and charge of the reference they were made from.

    python examples/clayff/test_clayff.py [WORKDIR]        (from the root of the repository; WORKDIR default clayff_test)

1. the 11 bulk crystals and the 11 basal (001) slabs of the library, as 2 x 2 x 1 supercells (not the cells of the references)
2. the same supercells with every atom displaced at random (0.05 A per coordinate; an H moves with its O, plus 0.02 A),
   like a frame of an MD at room temperature
3. kaolinite written to a CIF file (no types, no charges in the file)
4. kaolinite with water on its basal plane (examples/clayff/kaolinite_water.extxyz): mineral first, then water by the
   solution library; no atom may fall back to a rule
The truth is the type (Masses comment) and the charge of each atom in the LAMMPS data file of the reference.
Results: examples/results/clayff.json"""
import glob
import json
import os
import sys
import time

import ase.io
import numpy as np

import soapff

W = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "clayff_test")
LIB = os.path.join(os.path.dirname(soapff.__file__), "libraries", "clayff")
os.makedirs(W, exist_ok=True)
quiet = dict(md=False, log=lambda *a: None)


def reference(path):
    """ase.Atoms of a reference data file with the arrays 'fftype' (type name) and 'ffcharge'"""
    names, sec = {}, None
    for line in open(path):
        s = line.strip()
        if s in ("Masses", "Atoms", "Atoms  # full") or s.startswith(("Atoms", "Bonds", "Angles", "Velocities")):
            sec = s.split()[0]
            continue
        if sec == "Masses" and "#" in s:
            names[int(s.split()[0])] = s.split("#")[1].split()[0]
    at = ase.io.read(path, format="lammps-data", atom_style="full", units="real")
    at.new_array("fftype", np.array([names[int(t)] for t in at.arrays["type"]]))
    at.new_array("ffcharge", np.asarray(at.arrays["initial_charges"], float))
    for k in ("bonds", "angles", "mmcharges", "initial_charges", "id", "type", "mol-id", "masses"):
        if k in at.arrays:
            del at.arrays[k]
    return at


def check(name, at, structure=None):
    """type `structure` (default: at) with clayff; compare with the fftype / ffcharge of `at`"""
    t0 = time.time()
    s = soapff.build("clayff", at if structure is None else structure, os.path.join(W, name), **quiet)
    rows = s.table()
    sim = np.array([float(r["similarity"]) for r in rows])
    wrong_t = int((np.array(s.types) != at.arrays["fftype"]).sum())
    wrong_q = int((np.abs(np.asarray(s.charges) - at.arrays["ffcharge"]) > 0.01).sum())
    out = dict(atoms=len(at), wrong_type=wrong_t, wrong_charge=wrong_q, min_similarity=round(float(sim.min()), 4),
               seconds=round(time.time() - t0, 1))
    print(f"{name:42s} {len(at):6d} atoms  wrong type {wrong_t}  wrong charge {wrong_q}  min similarity {sim.min():.4f}")
    return out


def displacement(at, heavy=0.05, h=0.02):
    """random displacements like a frame of an MD: 0.05 A (Gaussian, per coordinate) for the heavy atoms; an H moves with the O
    it is bonded to, plus 0.02 A (the stiff O-H bond is not stretched as a free random move of both atoms would)"""
    d = rng.normal(0.0, heavy, at.positions.shape)
    iH = np.nonzero(at.numbers == 1)[0]
    if len(iH):
        D = at.get_all_distances(mic=True)[iH][:, at.numbers == 8]
        iO = np.nonzero(at.numbers == 8)[0][D.argmin(axis=1)]
        d[iH] = d[iO] + rng.normal(0.0, h, (len(iH), 3))
    return d


res = {}
rng = np.random.default_rng(1)
for kind in ("bulk", "basal"):
    for path in sorted(glob.glob(os.path.join(LIB, kind, "*.data"))):
        name = os.path.splitext(os.path.basename(path))[0]
        at = reference(path)
        sc = at.repeat((2, 2, 1))
        res[f"{name} 2x2x1"] = check(f"{name}_2x2x1", sc)
        rat = sc.copy()
        rat.positions += displacement(rat)
        res[f"{name} 2x2x1 displaced"] = check(f"{name}_2x2x1_displaced", rat)

# kaolinite from a CIF file
kao = reference(os.path.join(LIB, "bulk", "kaolinite_bulk.data"))
cif = os.path.join(W, "kaolinite.cif")
ase.io.write(cif, kao)
back = ase.io.read(cif)
assert (back.numbers == kao.numbers).all() and np.allclose(back.positions, kao.positions, atol=1e-3)
res["kaolinite from CIF"] = check("kaolinite_cif", kao, structure=cif)

# kaolinite with water: no truth per atom in the file; the build must not fall back for any atom
s = soapff.build("clayff", "examples/clayff/kaolinite_water.extxyz", os.path.join(W, "kaolinite_water"), **quiet)
rows = s.table()
src = {}
for r in rows:
    src[r["source"]] = src.get(r["source"], 0) + 1
roles = {}
for r in rows:
    roles[r["role"]] = roles.get(r["role"], 0) + 1
res["kaolinite with water"] = dict(atoms=len(rows), sources=src, roles=roles, net_charge=round(float(np.sum(s.charges)), 6),
                                   types=sorted(set(s.types)))
print("kaolinite with water", res["kaolinite with water"])

bad = {k: v for k, v in res.items() if v.get("wrong_type") or v.get("wrong_charge")}
print(f"{len(res) - 1} structures with a truth: {len(bad)} with errors")
os.makedirs("examples/results", exist_ok=True)
json.dump(res, open("examples/results/clayff.json", "w"), indent=1)
