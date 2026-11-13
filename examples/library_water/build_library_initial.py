"""The library with water and aqueous ions from the initial structures only (library_clayff_water.npz): the 11 bulk layered minerals and their
basal (0 0 1) slabs of the water-free library, plus the structures that packmol builds (build_refs.py; not relaxed, no MD frame): liquid water and the
five chloride solutions (NaCl, KCl, CsCl, CaCl2, BaCl2), and the same water and solutions in vacuum (films and a droplet: build_refs_vacuum.py).  The mineral is typed first, on the structure without the water (core.make), so the library
has no mineral-water interface rows; the water and the ions are typed in the second step on the structure of the water and the ions alone.  The 2-layer hydrate is not included (its structure
comes from an NPT run).  Rows of the systems with water are pruned at 0.99, the slab rows at 0.999 (as in build_library.py).

    python build_library_initial.py     (with a Python where soapff is installed)  ->  library_clayff_water.npz, build_library_initial.json"""
import json
import os

from soapff_cui import core
from build_library import BASAL, BULK, LIBDIR, HERE, reference

SYSTEMS = ["water", "NaCl", "KCl", "CsCl", "CaCl2", "BaCl2"]
VACUUM = ["vac_water_film", "vac_water_film_thin", "vac_water_droplet", "vac_NaCl_film", "vac_KCl_film", "vac_CsCl_film", "vac_CaCl2_film", "vac_BaCl2_film"]
SYSTEMS = SYSTEMS + VACUUM
refs = [os.path.join(LIBDIR, "refs_bulk", b + ".data") for b in BULK] + \
       [os.path.join(LIBDIR, "refs_basal", m + "_001.data") for m in BASAL] + [reference(s) for s in SYSTEMS]
prune = [None] * len(BULK) + [0.999] * len(BASAL) + [0.99] * len(SYSTEMS)
out = os.path.join(HERE, "library_clayff_water.npz")
core.build_library(refs, [], out=out, water_in_soap=True, prune=prune, protect=len(BULK))
lib = core.load_library(out)
r = core.check_library(out)
info = dict(rows=int(len(lib["q"])), species=list(lib["species"]), rows_per_type={t: int((lib["type"] == t).sum()) for t in sorted(set(lib["type"]))},
            one_to_one=dict(rows=r["total"], same_charge=r["total"] - r["wrong"], confusions=r["confusions"]))
json.dump(info, open(os.path.join(HERE, "build_library_initial.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
