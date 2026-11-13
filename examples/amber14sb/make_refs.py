"""The references of the amber14sb library and a test set kept out of it, typed by OpenMM (ff14SB + TIP3P + ions).

references: hen egg-white lysozyme (1AKI) and T4 lysozyme (2LZM), whole, with every histidine as HID, HIE, HIP (and T4 lysozyme with
ASH, GLH, LYN); for each of the 20 amino acids (and HID, HIE, HIP) a fragment that starts with it (N terminus) and one that ends
with it (C terminus); fragments capped with ACE and NME; TIP3P water with Na+ Cl-; one ion of each kind.
test (not in the references): ubiquitin (1UBQ), BPTI (5PTI), crambin (1CRN), Trp-cage (1L2Y), fragments of ubiquitin, and three
frames of a 300 K MD of ubiquitin in vacuum (flexible bonds).

    python make_refs.py OUTDIR      (OpenMM, pdbfixer)"""
import json
import os
import sys

import numpy as np
import openmm as mm
from openmm import app, unit

from omm_common import FF, atom_lines, fragment, prepare, residue_numbers, typed, write_extxyz

OUT = sys.argv[1]
os.makedirs(os.path.join(OUT, "refs"), exist_ok=True)
os.makedirs(os.path.join(OUT, "test"), exist_ok=True)
manifest = {}


def save(kind, name, top, pos, comment=""):
    d = typed(top, pos)
    path = os.path.join(OUT, kind, name + ".extxyz")
    write_extxyz(path, d, comment)
    manifest[f"{kind}/{name}"] = dict(atoms=len(d["symbols"]), charge=round(sum(d["ffcharge"]), 6),
                                      templates=sorted(set(d["template"])))
    print(kind, name, len(d["symbols"]), round(sum(d["ffcharge"]), 4), flush=True)
    return d


def whole(pdbid, variant=None):
    return prepare("".join(atom_lines(pdbid)) + "END\n", variant)


# ---------------------------------------------------------------- references
save("refs", "1AKI", *whole("1AKI"), "hen egg-white lysozyme 1AKI (pH 7, OpenMM protonation)")
save("refs", "2LZM", *whole("2LZM"), "T4 lysozyme 2LZM (pH 7, OpenMM protonation)")
# protonation variants, each inside a capped fragment (ACE ... NME) of five residues around it
for pdbid, center, res, variants in (("1AKI", 15, "HIS", ("HID", "HIE", "HIP")), ("2LZM", 31, "HIS", ("HID", "HIE", "HIP")),
                                     ("2LZM", 10, "ASP", ("ASH",)), ("2LZM", 11, "GLU", ("GLH",)), ("2LZM", 16, "LYS", ("LYN",)),
                                     ("2LZM", 54, "CYS", ("CYS",)), ("2LZM", 97, "CYS", ("CYS",))):
    L = atom_lines(pdbid)
    for v in variants:
        save("refs", f"{pdbid}_{v}{center}_capped", *prepare(fragment(L, center - 2, center + 2, ace=True, nme=True), {res: v}),
             f"{pdbid} {center - 2}-{center + 2} capped, {res} {center} as {v}")

# termini: every amino acid at the start and at the end of a fragment of four residues
LEN = 4
done_n, done_c = set(), set()
for pdbid in ("2LZM", "1AKI"):
    L = atom_lines(pdbid)
    rn = residue_numbers(L)
    nums = [r for r, _ in rn]
    for k, (r, name) in enumerate(rn):
        variants = ("HID", "HIE", "HIP") if name == "HIS" else (None,)
        for v in variants:
            key = v or name
            var = {"HIS": v} if v else None
            if key not in done_n and k + LEN <= len(rn) and nums[k + LEN - 1] - r == LEN - 1:
                save("refs", f"{pdbid}_N{key}_{r}", *prepare(fragment(L, r, nums[k + LEN - 1]), var), f"{pdbid} {r}-{nums[k + LEN - 1]}")
                done_n.add(key)
            if key not in done_c and k - LEN + 1 >= 0 and r - nums[k - LEN + 1] == LEN - 1:
                save("refs", f"{pdbid}_C{key}_{r}", *prepare(fragment(L, nums[k - LEN + 1], r), var), f"{pdbid} {nums[k - LEN + 1]}-{r}")
                done_c.add(key)
print("N termini", len(done_n), sorted(done_n))
print("C termini", len(done_c), sorted(done_c))

# ACE / NME caps
L = atom_lines("1AKI")
for a, b in ((2, 6), (20, 24), (40, 44), (60, 64), (80, 84), (100, 104), (118, 122)):
    save("refs", f"1AKI_ACE_{a}-{b}_NME", *prepare(fragment(L, a, b, ace=True, nme=True)), f"1AKI {a}-{b} capped")

# TIP3P water with NaCl (OpenMM addSolvent)
mod = app.Modeller(app.Topology(), [])
mod.addSolvent(FF, boxSize=mm.Vec3(2.0, 2.0, 2.0) * unit.nanometer, positiveIon="Na+", negativeIon="Cl-",
               ionicStrength=0.5 * unit.molar)
save("refs", "tip3p_NaCl", mod.topology, mod.positions, "TIP3P water, NaCl 0.5 M")

# one ion of each kind (Joung-Cheatham monovalent, Li/Merz divalent 12-6 of OpenMM amber14/tip3p.xml)
for el, name in (("Li", "LI"), ("Na", "NA"), ("K", "K"), ("Rb", "RB"), ("Cs", "CS"), ("F", "F"), ("Cl", "CL"), ("Br", "BR"),
                 ("I", "I"), ("Mg", "MG"), ("Ca", "CA"), ("Zn", "ZN"), ("Sr", "SR"), ("Ba", "BA")):
    top = app.Topology()
    ch = top.addChain()
    res = top.addResidue(name, ch)
    top.addAtom(name, app.Element.getBySymbol(el), res)
    save("refs", f"ion_{el}", top, [mm.Vec3(0, 0, 0)] * unit.nanometer, f"{el} ion")

# ---------------------------------------------------------------- test set (not in the references)
for pdbid in ("1UBQ", "5PTI", "1CRN", "1L2Y"):
    save("test", pdbid, *whole(pdbid), pdbid)
L = atom_lines("1UBQ")
for a, b in ((1, 5), (30, 36), (70, 76)):
    save("test", f"1UBQ_{a}-{b}", *prepare(fragment(L, a, b)), f"1UBQ {a}-{b}")
for a, b in ((10, 15), (50, 55)):
    save("test", f"1UBQ_ACE_{a}-{b}_NME", *prepare(fragment(L, a, b, ace=True, nme=True)), f"1UBQ {a}-{b} capped")

# 300 K MD of ubiquitin in vacuum, flexible bonds (Langevin, 0.5 fs), three frames
top, pos = whole("1UBQ")
sysm = FF.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.5 * unit.femtosecond)
sim = app.Simulation(top, sysm, integ, mm.Platform.getPlatformByName("CPU"))
sim.context.setPositions(pos)
sim.minimizeEnergy(maxIterations=500)
sim.context.setVelocitiesToTemperature(300 * unit.kelvin, 1)
for k in range(3):
    sim.step(10000)                                           # 5 ps
    p = sim.context.getState(getPositions=True).getPositions()
    save("test", f"1UBQ_md300K_{(k + 1) * 5}ps", top, p, f"1UBQ, vacuum MD 300 K, {(k + 1) * 5} ps")
json.dump(manifest, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
