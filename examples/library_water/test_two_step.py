"""soapff's two steps on the 16 initial structures (10 packmol interfaces, 6 water shapes; not relaxed, no MD frame) with library_clayff_water.npz
(the mineral first, on the structure without the water; then the water and the ions with the solution library: runner.prepare).  Also the SOAP of the water and the ions alone (the mineral taken away) and of the
whole structure against the solution library alone (what a structure with water would give without the first step).
Truth: the bond-graph ClayFF types of the mineral, ow/hw for the water, the element for the ions (water_common.truth_rule).
Counted per structure: wrong types and fallback of runner.prepare, atoms typed by the rule (water by shape, ions formal) instead of by SOAP;
wrong types, atoms below 0.9 and the minimum similarity per kind of the whole-structure SOAP.
    python test_two_step.py     (with a Python where soapff is installed)  ->  test_two_step.json"""
import collections
import csv
import json
import os
import tempfile

import numpy as np

from soapff_cui import core
from soapff import clayff as runner
from water_common import nearest, truth_rule, alias

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..")
THR = 0.9
LIB = core.load_library(runner.LIBRARY_WATER)
MIN = ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite", "kaolinite", "lizerdite", "talc", "brucite"]
cases = {f"interface_{m}": f"interfaces_001_water/{m}/structure.extxyz" for m in MIN}
cases.update({f"shape_{c}": f"interface_shapes/{c}/structure.extxyz" for c in ["sphere", "hemisphere", "cylinder", "cavity", "xygauss", "ions"]})
rows = {}
for name, rel in cases.items():
    at = core.read_structure(os.path.join(EX, rel))
    tr = truth_rule(at)
    d = tempfile.mkdtemp()
    rep = runner.prepare(at, d, runner.LIBRARY, log=lambda *a: None, md=False)
    ch = list(csv.DictReader(open(os.path.join(d, "charges.csv"))))
    t1 = alias(np.array([x["clayff_type"] for x in ch]))
    src = collections.Counter(x.get("source") for x in ch)
    t3, s3 = nearest(at, LIB)
    n0 = int(np.nonzero(np.asarray(at.arrays["added"], bool) if "added" in at.arrays else core.roles(at)[0] == "water")[0][0])
    tw_, sw_ = nearest(at[n0:], LIB)                          # the water and the ions alone (the mineral taken away)
    water = np.isin(tr, ["ow", "hw"]); ion = np.isin(tr, ["Na", "K", "Cs", "Ca", "Ba", "Cl"]); mine = ~water & ~ion
    rows[name] = dict(atoms=len(at), water=int((tr == "ow").sum()), ions=int(ion.sum()),
                      two_step=dict(wrong=int((t1 != tr).sum()), fallback=int(rep["fallback"]), by_rule=int(src.get("rule", 0)),
                                    type_check_differ=int(rep["type_check"]["differ"])),
                      water_only=dict(wrong=int((tw_ != tr[n0:]).sum()), low=int((sw_ < THR).sum()), min=float(sw_.min())),
                      whole=dict(wrong=int((t3 != tr).sum()), low=int((s3 < THR).sum()), low_mineral=int((s3[mine] < THR).sum()),
                                 min_water=float(s3[water].min()), min_ion=float(s3[ion].min()) if ion.any() else None, min_mineral=float(s3[mine].min())))
    print(name, rows[name], flush=True)
tot = dict(structures=len(rows), atoms=sum(v["atoms"] for v in rows.values()), water=sum(v["water"] for v in rows.values()),
           ions=sum(v["ions"] for v in rows.values()),
           two_step={k: sum(v["two_step"][k] for v in rows.values()) for k in ("wrong", "fallback", "by_rule", "type_check_differ")},
           water_only={k: sum(v["water_only"][k] for v in rows.values()) for k in ("wrong", "low")},
           whole={k: sum(v["whole"][k] for v in rows.values()) for k in ("wrong", "low", "low_mineral")})
print(tot)
json.dump(dict(cases=rows, total=tot), open(os.path.join(HERE, "test_two_step.json"), "w"), indent=1)
