"""The water molecules taken apart in step 1 (paper: "水との界面"): soapff bonds two atoms nearer than the bond length of their
elements and then takes the water molecules (an O bonded to two H only) apart from the rest (soapff.molecular.water_apart),
so that a water O nearer than the Al-O or Mg-O bond length to an edge Al or Mg is not bonded to it.  This script runs step 1
on the clay test structures with water and counts the bonds between water and the rest that are taken out: the packmol
interfaces and shapes, the NaCl droplet, the last frames of the MD with water (water_md_frames.py) and the edge particle with a
water film.  None is taken out in any of them (the water O stay farther from the cations than the bond lengths), so their types
are those of soapff before the water was taken apart.  The nano-gouge (../nanogouge/typing_test.py) is the case where they are.

    python water_apart_check.py          (soapff)  ->  ../results/clayff_water_apart_check.json"""
import glob
import json
import os

from soapff import library
from soapff import molecular as M
from soapff_cui import core

EX = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(EX, "results", "clayff_water_apart_check.json")


def jobs():
    for p in sorted(glob.glob(f"{EX}/interfaces_001_water/*/system.data") + glob.glob(f"{EX}/interface_shapes/*/system.data")
                    + [f"{EX}/meng2022_droplet/paper/system.data"]):
        yield p, "clayff"
    for p in sorted(glob.glob(f"{EX}/library_water/refs/*/final.data")) + [f"{EX}/montmorillonite_2W/diffusion/start/final.data"]:
        yield p, "clayff"                                  # the MD frames of water_md_frames.py
    yield f"{EX}/montmorillonite_edge/out_particle_water/system.data", "clayff_edge"


def main():
    T = {n: M.load_table(library.load(n).table_path(log=lambda *a: None)) for n in ("clayff", "clayff_edge")}
    out = {}
    for p, lib in jobs():
        if not os.path.isfile(p):
            print("missing (run its example first):", p)
            continue
        at, _, _ = core.read_lammps_data(p)
        at.pbc = True
        nb = M.bonds_by_distance(at, T[lib]["bond_lengths"])
        _, ow, cut = M.water_apart(at, nb)
        name = os.path.relpath(p, EX)
        out[name] = dict(library=lib, atoms=len(at), water_molecules=len(ow), bonds_taken_out=cut)
        print(f"{name}: {len(at)} atoms, {len(ow)} water molecules, bonds taken out {cut or 'none'}", flush=True)
    res = dict(structures=out, total=dict(structures=len(out), with_bonds_taken_out=sum(bool(v["bonds_taken_out"])
                                                                                       for v in out.values())))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(res["total"], "->", OUT)


if __name__ == "__main__":
    main()
