"""The solution library of the edge library "Montmorillonite edge (Underwood & Bourg 2020)" (library_clayff_edge_water.npz): the two particle files of the
edge library (../library/refs_edge/), and the structures that packmol builds with water and cations, not relaxed, no MD frame: liquid water and the five
chloride solutions (NaCl, KCl, CsCl, CaCl2, BaCl2; refs/<name>/<name>.data of build_refs.py, ClayFF water and ions) and the same in vacuum (build_refs_vacuum.py).  The particle is typed first, on the
structure without the water, so the library has no particle-water rows.  The water and the ions are ClayFF in every system (the model of the
particle has neither Cl- nor water).  Rows of the systems with water are pruned at 0.99, the second particle file at 0.999 (as the edge library).

    python build_library_clayff_edge_water.py     (with a Python where soapff is installed)  ->  library_clayff_edge_water.npz, build_library_clayff_edge_water.json"""
import json
import os

from soapff_cui import core
from build_library import HERE, reference

LIBDIR = os.path.join(HERE, "..", "library", "refs_edge")
SYSTEMS = ["water", "NaCl", "KCl", "CsCl", "CaCl2", "BaCl2"]
VACUUM = ["vac_water_film", "vac_water_film_thin", "vac_water_droplet", "vac_NaCl_film", "vac_KCl_film", "vac_CsCl_film", "vac_CaCl2_film", "vac_BaCl2_film"]
SYSTEMS = SYSTEMS + VACUUM
refs = [os.path.join(LIBDIR, "clay_montmorillonite.data"), os.path.join(LIBDIR, "clay_montmorillonite_2.data")] + \
       [reference(s) for s in SYSTEMS]
prune = [None, 0.999] + [0.99] * len(SYSTEMS)
out = os.path.join(HERE, "library_clayff_edge_water.npz")
core.build_library(refs, [], out=out, water_in_soap=True, prune=prune, protect=1)
lib = core.load_library(out)
r = core.check_library(out)
info = dict(rows=int(len(lib["q"])), species=list(lib["species"]), rows_per_type={t: int((lib["type"] == t).sum()) for t in sorted(set(lib["type"]))},
            one_to_one=dict(rows=r["total"], same_charge=r["total"] - r["wrong"], confusions=r["confusions"]))
json.dump(info, open(os.path.join(HERE, "build_library_clayff_edge_water.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
