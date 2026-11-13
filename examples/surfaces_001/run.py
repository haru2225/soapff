"""Basal (0 0 1) surfaces of the 11 layered minerals (paper: "basal planes and interfaces", tests (2) and (3)).
  (2) every face |h|, |k|, |l| <= 1 (13 faces) of each bulk reference (../library/refs_bulk/) is cut with pymatgen
      (soapff_cui.surfaces.make_slab): only the faces that cut no cation-O bond are made;
  (3) each (0 0 1) slab is cut again from the bulk doubled along a, with 30 A of vacuum (the default), written as
      elements, coordinates and box only (slabs/<mineral>_001.extxyz), and typed with the bundled library; the truth is the
      type each atom had in the bulk and its ClayFF Table 1 charge.

    python run.py            (with a Python where soapff is installed)  ->  slabs/, out/<mineral>/ (structure.extxyz, the slab, and the
                                                             files soapff writes for it), results.json"""
import csv
import json
import os
import shutil

import ase.io
import numpy as np
from ase import Atoms

from soapff_cui import core, surfaces
from soapff_cui.clayff import type_charges

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(HERE, "..", "library", "refs_bulk")
LIB = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "prebuilt", "library_clayff_bulk.npz")
BULK = ["montmorillonite_dry_Na", "muscovite_2M1_bulk", "muscovite_1M_bulk", "phlogopite_1M_bulk", "nanpingite_bulk",
        "pyrophyllite_bulk", "pyrophyllite_2_bulk", "kaolinite_bulk", "lizerdite_bulk", "talc_bulk", "brucite_bulk"]
Q1 = type_charges("clayff2004")

res = dict(faces={}, slabs={})
for b in BULK:
    name = b.replace("_bulk", "")
    at, q, types = core.read_lammps_data(os.path.join(REFS, b + ".data"))
    keep = ~np.isin(types, ["ow", "hw"])
    at, types = at[keep], [t for t, k in zip(types, keep) if k]
    made, refused = [], []
    for m in surfaces.faces(1):                                  # (2) the 13 faces
        try:
            surfaces.make_slab(at, types, m)
            made.append(surfaces.hkl_name(m))
        except ValueError:
            refused.append(surfaces.hkl_name(m))
    res["faces"][name] = dict(made=made, refused=refused)
    # (3) the (0 0 1) slab from the bulk doubled along a, 30 A of vacuum
    slab, typ, info = surfaces.make_slab(at.repeat((2, 1, 1)), types * 2, (0, 0, 1), vacuum=30.0)
    f = os.path.join(HERE, "slabs", name + "_001.extxyz")
    ase.io.write(f, Atoms(slab.get_chemical_symbols(), positions=slab.positions, cell=slab.cell, pbc=True), format="extxyz")
    out = os.path.join(HERE, "out", name)
    rep = core.make(f, out, LIB, md=False, stop_on_fallback=True, log=lambda *a: None)
    shutil.copy(f, os.path.join(out, "structure.extxyz"))
    rows = list(csv.DictReader(open(os.path.join(out, "charges.csv"))))
    wrong = sum(r["clayff_type"] != t or abs(float(r["charge"]) - Q1[t]) > 1e-6 for r, t in zip(rows, typ))
    res["slabs"][name] = dict(atoms=len(rows), wrong=wrong, fallback=rep["fallback"], type_check_differ=rep["type_check"]["differ"],
                              min_similarity=round(min(float(r["similarity"]) for r in rows), 4),
                              charge_sum=round(rep["charge_sum"], 4), formal_charge=info["formal_charge"], ions=info["ions"],
                              ions_face=info["ions_face"], thickness=info["thickness"], vacuum=info["vacuum"])
    print(name, "faces made", made, "| slab", res["slabs"][name], flush=True)
n_made = sum(len(v["made"]) for v in res["faces"].values())
n_ref = sum(len(v["refused"]) for v in res["faces"].values())
res["summary"] = dict(faces_tried=n_made + n_ref, made=n_made, refused=n_ref)
print(res["summary"])
json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
