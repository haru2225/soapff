"""The references of the library oplsaa and the test structures, typed by Foyer (OPLS-AA, oplsaa.xml).

    python make_refs.py LIBDIR TESTDIR        (rdkit, openbabel, foyer, openmm; LIBDIR = ../../src/soapff/libraries)

References (LIBDIR/oplsaa/*.extxyz, listed in [references] of LIBDIR/oplsaa.txt), for each of the 11 polymers of
polymers.POLYMERS: a 5-mer as built (Open Babel, an extended chain, then minimised with Foyer's force field), another
5-mer (another tacticity) after 10 ps of MD at 413 K in vacuum, an 8-mer as built, and 5 frames (10 to 50 ps) of MD at 413 K
in vacuum of another 8-mer (the chain folds up: the conformations of a melt; without them some H of PEO in a melt cell were
taken for the H of PVA).
Test structures (TESTDIR/*.extxyz with Foyer's types and charges, and *.json with the molecules), none of them a reference:
a long chain of each polymer (PE C100, the others 20 repeat units, another tacticity), two random copolymers
(ethylene-propylene, MMA-MA), and the last frame of 10 ps of MD at 413 K and at 600 K in vacuum of each long chain.  Also,
for the user library of user_library.py (extra/): a PE-b-PEO block copolymer (its junction is in no reference of oplsaa)
and a short PE-b-PEO, the reference added to the user library.  The amorphous cells are made by
make_cells.py.  The minimisation and the MD use Foyer's force field in OpenMM (polymers.system).  Open Babel's --gen3d is
not deterministic: a new run gives other (equally good) structures."""
import os
import re
import sys

import numpy as np

import polymers as P

LIB, T = (os.path.abspath(x) for x in sys.argv[1:3])
REF = os.path.join(LIB, "oplsaa")
os.makedirs(REF, exist_ok=True)
os.makedirs(T, exist_ok=True)


def report(kind, name, s, box=None):
    q = np.array([a.charge for a in s.atoms])
    print(f"{kind} {name:18s} {len(q):5d} atoms  net {q.sum():+.4f}" + (f"  box {box[0]:.2f} A" if box is not None else ""),
          flush=True)


# ------------------------------------------------------------------ references
refs = []
for p in P.POLYMERS:
    for n, seed, md in ((5, 1, False), (5, 2, True), (8, 3, False)):
        smi = P.chain_smiles(p, n, seed=seed)
        mol, xyz, smi3 = P.embed(smi)
        xyz = P.relax(mol, xyz)
        if md:
            xyz = P.shake_up(mol, xyz, temp=413.0, ps=10.0, seed=seed)
        name = f"{p}_{n}mer" + ("_md413K" if md else "")
        s = P.save(os.path.join(REF, name + ".extxyz"), [mol], xyz, comment=f"{P.POLYMERS[p]}, {n}-mer"
                   + (", 10 ps at 413 K" if md else ""), smiles=[smi3])
        refs.append(f"oplsaa/{name}.extxyz")
        report("ref ", name, s)
    # the conformations of a melt: 5 frames (10 to 50 ps) of MD at 413 K in vacuum of another 8-mer, which folds up
    smi = P.chain_smiles(p, 8, seed=4)
    mol, xyz, smi3 = P.embed(smi)
    for k, x in enumerate(P.frames(mol, P.relax(mol, xyz), temp=413.0, every_ps=10.0, n=5, seed=4)):
        name = f"{p}_8mer_md413K_{k + 1}"
        s = P.save(os.path.join(REF, name + ".extxyz"), [mol], x, comment=f"{P.POLYMERS[p]}, 8-mer, {10 * (k + 1)} ps at 413 K",
                   smiles=[smi3])
        refs.append(f"oplsaa/{name}.extxyz")
        report("ref ", name, s)
txt = os.path.join(LIB, "oplsaa.txt")
lib = open(txt).read()
lib = re.sub(r"\[references\].*", "[references]\n# oligomers typed by Foyer (examples/oplsaa/make_refs.py)\n"
             + "\n".join(refs) + "\n", lib, flags=re.S)
open(txt, "w").write(lib)

# ------------------------------------------------------------------ test structures
chains = {f"chain_{p}": P.chain_smiles(p, 50 if p == "PE" else 20, seed=11) for p in P.POLYMERS}
rng = np.random.default_rng(5)
chains["copoly_EP"] = P.chain_smiles(None, 30, seed=12, units=list(rng.choice(["PE", "PP"], 30)))
chains["copoly_MMA_MA"] = P.chain_smiles(None, 20, seed=13, units=list(rng.choice(["PMMA", "PMA"], 20)))
for k, (name, smi) in enumerate(chains.items()):
    mol, xyz, smi3 = P.embed(smi)
    xyz = P.relax(mol, xyz)
    report("test", name, P.save(os.path.join(T, name + ".extxyz"), [mol], xyz, comment=name, smiles=[smi3]))
    if name.startswith("chain_"):
        for temp, seed in ((413.0, 30 + k), (600.0, 60 + k)):           # MD frames: the melt temperature, and hot
            x = P.shake_up(mol, xyz, temp=temp, ps=10.0, seed=seed)
            md = f"md{temp:.0f}K_" + name[6:]
            report("test", md, P.save(os.path.join(T, md + ".extxyz"), [mol], x,
                                      comment=f"{name} after 10 ps at {temp:.0f} K", smiles=[smi3]))

# the example of a user library (user_library.py): a PE-b-PEO block copolymer, and a short PE-b-PEO with its junction
EXTRA = os.path.join(os.path.dirname(T), "extra")
os.makedirs(EXTRA, exist_ok=True)
mol, xyz, smi3 = P.embed("CC" * 10 + "OCC" * 10 + "OC")
xyz = P.relax(mol, xyz)
report("extr", "block_PE_PEO", P.save(os.path.join(EXTRA, "block_PE_PEO.extxyz"), [mol], xyz, comment="block_PE_PEO",
                                      smiles=[smi3]))
smi = "CC" * 4 + "OCC" * 3 + "OC"
mol, xyz, smi3 = P.embed(smi)
xyz = P.relax(mol, xyz)
report("extr", "PE4_PEO3", P.save(os.path.join(EXTRA, "PE4_PEO3.extxyz"), [mol], xyz, comment="PE-b-PEO junction",
                                  smiles=[smi3]))
