"""clayff on 300 K MD frames with water (paper: Table "clay test cases", the bond length of H-O): the last frame of the MD of
the systems of ../library_water/build_refs.py (liquid water, the NaCl, KCl, CsCl, CaCl2 and BaCl2 solutions, the (0 0 1)
slabs of the 11 minerals with water, the 2-layer Na-montmorillonite hydrate; 300 K, 25-30 ps, flexible SPC water) and the
2-layer hydrate after its 100 ps NPT (../montmorillonite_2W/diffusion/start/final.data).  No frame is a reference of clayff.

Each frame (elements, coordinates, box) is typed in the two steps of soapff (soapff.molecular.type_atoms with the table of
the bundled clayff library) and compared with the types of its MD input (system.data: same atom order).  Also reported: whether soapff would
stop (an atom below 0.9, an atom whose bonds are not those of its reference atom), and the distances that decide the bond
length of H-O: the longest O-H bond and the shortest H...O of an H to an O it is not bonded to (hydrogen bonds), with the
bond length of the library (H-O 1.25 A) and the one of the default rule (halfway in the references: 1.54 A).

    python water_md_frames.py          (soapff)  ->  ../results/clayff_water_md_frames.json"""
import glob
import json
import os

import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list

import soapff
from soapff import molecular as M
from soapff_cui import core

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..")
OUT = os.path.join(EX, "results", "clayff_water_md_frames.json")
alias = lambda t: np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t])  # noqa: E731


def frames():
    for d in sorted(glob.glob(os.path.join(EX, "library_water", "refs", "*"))):
        if os.path.isfile(os.path.join(d, "final.data")) and os.path.isfile(os.path.join(d, "system.data")):
            yield os.path.basename(d), d
    yield "hydrate_2W_after_npt", os.path.join(EX, "montmorillonite_2W", "diffusion", "start")


def oh_distances(at):
    """the O-H bonds (each H to its nearest heavy atom, when an O) and, for every H, its nearest O that is not that atom"""
    el = np.array(at.get_chemical_symbols())
    I, J, d = neighbor_list("ijd", at, 3.0)
    heavy = el[J] != "H"
    bond, other = [], []
    for i in np.nonzero(el == "H")[0]:
        k = (I == i) & heavy
        jj, dd = J[k], d[k]
        p = np.argmin(dd)
        if el[jj[p]] == "O":
            bond.append(dd[p])
        o = (el[jj] == "O") & (np.arange(len(jj)) != p)
        if o.any():
            other.append(dd[o].min())
    return np.array(bond), np.array(other)


def main():
    lib = soapff.load_library("clayff")
    T = M.load_table(lib.table_path(log=lambda *a: None))
    res = dict(bond_lengths={k: v["r"] for k, v in T["bond_lengths"].items()}, frames={})
    for name, d in frames():
        at0, q0, t0 = core.read_lammps_data(os.path.join(d, "system.data"))
        fin = read(os.path.join(d, "final.data"), format="lammps-data", atom_style="full")
        fin = fin[np.argsort(fin.arrays["id"])] if "id" in fin.arrays else fin
        fin.pbc = True
        fin.set_chemical_symbols(at0.get_chemical_symbols())
        truth = alias(t0)
        r = M.type_atoms(lib, T, fin, log=lambda *a: None)
        got = alias(r["type"])
        wrong = got != truth
        bond, other = oh_distances(fin)
        soap = r["source"] == "soap"
        row = dict(atoms=len(fin), water=int((truth == "ow").sum()), wrong=int(wrong.sum()),
                   below_0_9=int((r["sim1"][soap] < 0.9).sum()), bonds_not_of_type=len(r["odd"]),
                   min_similarity=round(float(r["sim1"][soap].min()), 4),
                   longest_OH_bond=round(float(bond.max()), 3), shortest_HO_not_bonded=round(float(other.min()), 3))
        if wrong.any():
            row["confused"] = sorted({f"{a}->{b}" for a, b in zip(truth[wrong], got[wrong])})
        res["frames"][name] = row
        print(name, row, flush=True)
    F = res["frames"].values()
    res["total"] = dict(frames=len(F), atoms=sum(r["atoms"] for r in F), wrong=sum(r["wrong"] for r in F),
                        frames_that_stop=sum(bool(r["below_0_9"] or r["bonds_not_of_type"]) for r in F),
                        min_similarity=min(r["min_similarity"] for r in F),
                        longest_OH_bond=max(r["longest_OH_bond"] for r in F),
                        shortest_HO_not_bonded=min(r["shortest_HO_not_bonded"] for r in F))
    print(res["total"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
