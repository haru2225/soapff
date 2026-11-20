"""atom's types and charges against soapff's (the right ones of the validation table).  atom's type names are mapped to
Cygan's with the correspondence atom itself gives in forcefield_functions/clayff_param.m (O -> ob, Omg -> obos,
Oalt -> obts, Odsub -> obss, Ohmg -> ohs, Oh -> oh, Si -> st, Al -> ao, Alt -> at, Mgo -> mgo, Mgh -> mgh, H -> ho,
Ow -> ow, Hw -> hw), plus its newer names for the O of the basal plane and of the apex (Ob, Op -> ob; Obs, Ops -> obts).
An atom's type differs when its mapped name is not soapff's; its charge differs when it is off by more than 1e-3 e.
usage: python score_atom.py IN OUT      (IN from dump_clay.py, OUT from run_atom.m; writes OUT/score.json)"""
import collections, glob, json, os, re, sys
import numpy as np

IN, OUT = sys.argv[1:3]
MAP = {"O": "ob", "Ob": "ob", "Op": "ob", "Omg": "obos", "Oalt": "obts", "Obs": "obts", "Ops": "obts", "Odsub": "obss",
       "Oh": "oh", "Ohmg": "ohs", "Si": "st", "Al": "ao", "Alt": "at", "Mgo": "mgo", "Mgh": "mgh", "H": "ho",
       "Ow": "ow", "Hw": "hw", "Na": "Na", "K": "K", "Cs": "Cs", "Ca": "Ca", "Ba": "Ba", "Cl": "Cl"}
Q = {"ao": 1.575, "at": 1.575, "st": 2.1, "mgo": 1.36, "mgh": 1.05, "ho": 0.425, "ob": -1.05, "obos": -1.1808,
     "obts": -1.1688, "obss": -1.2996, "oh": -0.95, "ohs": -1.0808, "ow": -0.82, "hw": 0.41,
     "Na": 1, "K": 1, "Cs": 1, "Ca": 2, "Ba": 2, "Cl": -1}
total = dict(structures=0, atoms=0, same_type=0, same_charge=0, same_both=0, seconds=0.0)
detail = {}
for f in sorted(glob.glob(os.path.join(OUT, "*.txt"))):
    name = os.path.basename(f)[:-4]
    truth = np.load(os.path.join(IN, name + ".npz"))["truth"]
    lines = open(f).read().splitlines()
    sec = float(re.search(r"seconds ([\d.]+)", lines[0]).group(1))
    rows = [l.split() for l in lines[1:]]
    assert [int(r[0]) for r in rows] == list(range(1, len(truth) + 1)), name  # atom added or removed no atom
    at = np.array([r[1] for r in rows]); q = np.array([float(r[2]) for r in rows])
    dq = np.abs(q - np.array([Q[t] for t in truth]))
    same_t = np.array([MAP.get(a) == t for a, t in zip(at, truth)])
    same_q = dq <= 1e-3
    diff = collections.Counter(f"{t} -> {a}" for a, t, s in zip(at, truth, same_t) if not s)
    detail[name] = dict(atoms=len(truth), same_type=int(same_t.sum()), same_charge=int(same_q.sum()),
                        same_both=int((same_t & same_q).sum()), max_charge_difference=float(dq.max()),
                        seconds=sec, other_names=dict(diff))
    for k in ("same_type", "same_charge", "same_both"):
        total[k] += detail[name][k]
    total["structures"] += 1; total["atoms"] += len(truth); total["seconds"] += sec
for k, v in detail.items():
    print(f"{k:36s} {v['atoms']:6d} type {v['same_type']:6d} charge {v['same_charge']:6d} {v['seconds']:6.0f} s  {v['other_names']}")
print("total", total)
json.dump(dict(total=total, detail=detail), open(os.path.join(OUT, "score.json"), "w"), indent=1)
