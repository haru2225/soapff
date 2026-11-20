"""The polymer chains of the clay + polymer + water tests: the PEO and PVA 20-mers of ../oplsaa/test (chain_PEO, chain_PVA),
made compact by Langevin MD at 600 K in vacuum inside a sphere (Foyer's OPLS-AA in OpenMM, ../oplsaa/polymers.shake_up),
minimised, and typed by Foyer (the right answer of the tests).  Needs foyer, openmm and parmed (the environment of
../oplsaa).

    python chains.py            -> chains/PEO.extxyz, chains/PVA.extxyz (+ .json, as ../oplsaa/polymers.save)"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "oplsaa"))
import polymers  # noqa: E402

RADIUS = {"PEO": 8.0, "PVA": 7.0}        # A: the sphere of the compact chain


def main():
    os.makedirs(os.path.join(HERE, "chains"), exist_ok=True)
    for p in ("PEO", "PVA"):
        mols, xyz, _, types0, q0 = polymers.load(os.path.join(HERE, "..", "oplsaa", "test", f"chain_{p}.extxyz"))
        mol = mols[0]
        x = polymers.shake_up(mol, xyz, temp=600.0, ps=30.0, seed=7, radius=RADIUS[p])
        x = polymers.relax(mol, x)
        x = x - x.mean(0)
        s = polymers.save(os.path.join(HERE, "chains", f"{p}.extxyz"), [mol], x,
                          comment=f"{p} 20-mer of ../oplsaa/test/chain_{p}, compact (600 K in a sphere of {RADIUS[p]:g} A), "
                                  f"minimised; types and charges by Foyer")
        same = [a.type for a in s.atoms] == list(types0)
        ext = x.max(0) - x.min(0)
        print(f"{p}: {len(x)} atoms, extent {ext.round(1)} A, Foyer types as in chain_{p}: {same}")


if __name__ == "__main__":
    main()
