"""Polymer chains for the library oplsaa and its tests: SMILES of the chains, 3D structures (Open Babel), amorphous cells
(packmol), and the OPLS-AA types and charges that Foyer gives them (the right answer of the tests).
Needs rdkit, openbabel, foyer, parmed, openmm and packmol (conda-forge).

A chain is H-[repeat unit]_n-H written as SMILES (a random tacticity for the vinyl polymers).  A molecule is kept as
(elements, bonds) in the atom order of its structure, and a structure file (extxyz) has a .json next to it with its molecules,
so that the topology for Foyer and OpenMM is rebuilt exactly."""
import json
import logging
import os
import subprocess
import sys
import tempfile
import warnings

import numpy as np

warnings.filterwarnings("ignore")
logging.getLogger("foyer").setLevel(logging.ERROR)

POLYMERS = {
    "PE":   "polyethylene",
    "PP":   "atactic polypropylene",
    "PIB":  "polyisobutylene",
    "PB":   "cis-1,4-polybutadiene",
    "PI":   "cis-1,4-polyisoprene",
    "PEO":  "poly(ethylene oxide), CH3O- ends",
    "PMMA": "atactic poly(methyl methacrylate)",
    "PMA":  "atactic poly(methyl acrylate)",
    "PVA":  "atactic poly(vinyl alcohol)",
    "PAN":  "atactic polyacrylonitrile",
    "PTFE": "polytetrafluoroethylene",
}
VINYL = {"PP": "C", "PMA": "C(=O)OC", "PVA": "O", "PAN": "C#N"}          # -CH2-CH(R)-


def unit(p, rng, last=False):
    """SMILES of one repeat unit of polymer p (the last unit of the chain without a stereo mark)"""
    if p == "PE":
        return "CC"
    if p in VINYL:
        r = VINYL[p]
        return f"CC({r})" if last else f"C[C{'@' if rng.random() < 0.5 else '@@'}H]({r})"
    if p == "PMMA":
        return "CC(C)(C(=O)OC)" if last else f"C[C{'@' if rng.random() < 0.5 else '@@'}](C)(C(=O)OC)"
    if p == "PIB":
        return "CC(C)(C)"
    if p == "PB":
        return "C/C=C\\C"
    if p == "PI":
        return "C/C(C)=C\\C"
    raise KeyError(p)


def chain_smiles(p, n, seed=0, units=None):
    """H-[unit]_n-H of polymer p (PEO: CH3O-[CH2CH2O]_n-CH3, PTFE: F-[CF2CF2]_n-F); units: the polymer of each unit (a
    copolymer of vinyl and olefin units)"""
    rng = np.random.default_rng(seed)
    if p == "PEO":
        return "CO" + "CCO" * n + "C"
    if p == "PTFE":
        return "F" + "C(F)(F)C(F)(F)" * n + "F"
    units = units or [p] * n
    return "".join(unit(u, rng, last=(i == len(units) - 1)) for i, u in enumerate(units))


# ------------------------------------------------------------------ molecules
def obabel():
    """the obabel executable of this Python environment and the environment variables of its plugins"""
    d = os.path.dirname(sys.executable)
    exe = os.path.join(d, "obabel") if os.path.isfile(os.path.join(d, "obabel")) else "obabel"
    env = dict(os.environ)
    root = os.path.dirname(d)
    for var, sub in (("BABEL_LIBDIR", "lib"), ("BABEL_DATADIR", "share")):
        base = os.path.join(root, sub, "openbabel")
        if os.path.isdir(base):
            env[var] = os.path.join(base, sorted(os.listdir(base))[-1])
    return exe, env


def embed(smi):
    """(elements, bonds), positions (A) and the SMILES of the 3D structure (its stereo) of a chain: Open Babel --gen3d
    (rule-based, an extended chain; RDKit's ETKDG is too slow for long chains with nitriles).  The connectivity is checked
    against the SMILES."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    exe, env = obabel()
    r = subprocess.run([exe, f"-:{smi}", "-osdf", "--gen3d"], capture_output=True, text=True, env=env)
    m = Chem.MolFromMolBlock(r.stdout, removeHs=False)
    assert m is not None and m.GetNumConformers(), r.stderr[-500:]
    flat = lambda x: Chem.MolToSmiles(x, isomericSmiles=False)                    # noqa: E731
    assert flat(Chem.RemoveHs(m)) == flat(Chem.MolFromSmiles(smi)), f"connectivity changed: {smi}"
    Chem.AssignStereochemistryFrom3D(m)
    mol = ([a.GetSymbol() for a in m.GetAtoms()], [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()])
    return mol, m.GetConformer().GetPositions(), Chem.MolToSmiles(Chem.RemoveHs(m))


def topology(mols):
    """an OpenMM topology with one residue per molecule; molecules with the same atoms and bonds get the same residue name,
    so that Foyer types them once"""
    import openmm.app as app
    top = app.Topology()
    ch = top.addChain()
    names = {}
    for el, bonds in mols:
        key = (tuple(el), tuple(map(tuple, bonds)))
        res = top.addResidue(names.setdefault(key, f"M{len(names):02d}"), ch)
        at = [top.addAtom(f"{e}{i + 1}", app.Element.getBySymbol(e), res) for i, e in enumerate(el)]
        for i, j in bonds:
            top.addBond(at[i], at[j])
    return top


def foyer_types(mols, xyz, box=None):
    """the parmed structure typed by Foyer's OPLS-AA (types, charges and parameters); xyz in A, box (A) or None"""
    import foyer
    import parmed as pmd
    s = pmd.openmm.load_topology(topology(mols), xyz=np.asarray(xyz))
    if box is not None:
        s.box = [box[0], box[1], box[2], 90, 90, 90]
    ff = foyer.Forcefield(name="oplsaa")
    return ff.apply(s, use_residue_map=True, assert_dihedral_params=False, assert_improper_params=False, verbose=False)


# ------------------------------------------------------------------ files
def write_extxyz(path, sym, xyz, types, charges, box=None, comment="", bonds=None):
    """extxyz with the per-atom columns fftype and ffcharge (what soapff reads in a reference, and the right answer of a
    test) and, with bonds (0-based pairs), the column bonds (the 1-based ids of the bonded atoms) that a reference needs"""
    lat = (f'Lattice="{box[0]:.6f} 0 0 0 {box[1]:.6f} 0 0 0 {box[2]:.6f}" pbc="T T T"' if box is not None
           else 'pbc="F F F"')
    col = [""] * len(sym)
    if bonds is not None:
        nbr = [[] for _ in sym]
        for i, j in bonds:
            nbr[i].append(j + 1); nbr[j].append(i + 1)
        col = ["  " + (",".join(str(x) for x in sorted(b)) or "-") for b in nbr]
    L = [str(len(sym)), f'{lat} Properties=species:S:1:pos:R:3:fftype:S:1:ffcharge:R:1{":bonds:S:1" if bonds is not None else ""} '
                        f'comment="{comment}"']
    L += [f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f} {t:10s} {q:10.6f}{c}" for e, (x, y, z), t, q, c in zip(sym, xyz, types, charges, col)]
    open(path, "w").write("\n".join(L) + "\n")


def save(path, mols, xyz, box=None, comment="", smiles=None):
    """type the structure with Foyer and write path (.extxyz) and its molecules (.json); returns the parmed structure"""
    s = foyer_types(mols, xyz, box)
    sym = [e for el, _ in mols for e in el]
    pairs, off = [], 0
    for el, bonds in mols:
        pairs += [(off + i, off + j) for i, j in bonds]
        off += len(el)
    write_extxyz(path, sym, xyz, [a.type for a in s.atoms], [a.charge for a in s.atoms], box, comment, bonds=pairs)
    json.dump(dict(molecules=[dict(elements=" ".join(el), bonds=[int(x) for b in bonds for x in b]) for el, bonds in mols],
                   smiles=smiles, box=None if box is None else [float(x) for x in box], comment=comment),
              open(path[:-7] + ".json", "w"))
    return s


def load(path):
    """molecules, positions (A), box (A or None), types, charges of a structure written by save"""
    d = json.load(open(path[:-7] + ".json"))
    mols = [(m["elements"].split(), list(zip(m["bonds"][::2], m["bonds"][1::2]))) for m in d["molecules"]]
    L = open(path).read().splitlines()
    t = [l.split() for l in L[2:2 + int(L[0])]]
    return (mols, np.array([[float(v) for v in x[1:4]] for x in t]), None if d["box"] is None else np.array(d["box"]),
            [x[4] for x in t], np.array([float(x[5]) for x in t]))


MASS = dict(H=1.008, C=12.011, N=14.007, O=15.999, F=18.998)


def mass(mols):
    return sum(MASS[e] for el, _ in mols for e in el)


# ------------------------------------------------------------------ MD with Foyer's force field (OpenMM)
def system(mols, xyz, box=None, cutoff=1.2, pme_tol=1e-6, dispersion=False):
    """the OpenMM System of Foyer's OPLS-AA (ParmEd; geometric mixing as a CustomNonbondedForce): no cutoff in vacuum, PME
    and an LJ cutoff (nm) with a box; flexible bonds"""
    import openmm as mm
    from openmm import app, unit
    s = foyer_types(mols, xyz, box)
    if box is None:
        sysm = s.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, removeCMMotion=False)
    else:
        sysm = s.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=cutoff * unit.nanometer, constraints=None,
                              ewaldErrorTolerance=pme_tol, removeCMMotion=False)
        for f in sysm.getForces():
            if isinstance(f, (mm.NonbondedForce, mm.CustomNonbondedForce)):
                if isinstance(f, mm.NonbondedForce):
                    f.setUseDispersionCorrection(dispersion)
                else:
                    f.setUseLongRangeCorrection(dispersion)
                    f.setUseSwitchingFunction(False)
    # the RB torsions with all coefficients zero (X-C-C#N of a nitrile) add nothing, but OpenMM's force of a torsion is NaN
    # when three of its atoms are on a line (Open Babel builds a nitrile exactly linear): they are left out
    for k, f in enumerate(sysm.getForces()):
        if isinstance(f, mm.RBTorsionForce):
            g = mm.RBTorsionForce()
            g.setUsesPeriodicBoundaryConditions(f.usesPeriodicBoundaryConditions())
            for i in range(f.getNumTorsions()):
                p = f.getTorsionParameters(i)
                if any(c.value_in_unit(unit.kilojoule_per_mole) != 0.0 for c in p[4:]):
                    g.addTorsion(*p)
            sysm.removeForce(k)
            sysm.addForce(g)
            break
    return s, sysm


def platform():
    """the OpenMM platform of the melt MD and its properties: $OPENMM_PLATFORM, else OpenCL (single precision; a GPU) when
    there is one, else CPU"""
    import openmm as mm
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    name = os.environ.get("OPENMM_PLATFORM") or ("OpenCL" if "OpenCL" in names else "CPU")
    return mm.Platform.getPlatformByName(name), ({"Precision": "single"} if name in ("OpenCL", "CUDA") else {})


def clash(mol, xyz, limit=1.2):
    """the shortest distance (A) between two atoms more than 2 bonds apart, if below limit (else None)"""
    el, bonds = mol
    n = len(el)
    nb = [set() for _ in range(n)]
    for i, j in bonds:
        nb[i].add(j)
        nb[j].add(i)
    d = np.linalg.norm(np.asarray(xyz)[:, None] - np.asarray(xyz)[None], axis=-1)
    for i in range(n):
        near = nb[i] | {k for j in nb[i] for k in nb[j]} | {i}
        d[i, list(near)] = np.inf
    return float(d.min()) if d.min() < limit else None


def relax(mol, xyz, tol=0.1):
    """the energy minimum (Foyer's force field, OpenMM, vacuum) next to xyz; stops if atoms are left on top of each other"""
    import openmm as mm
    from openmm import unit
    _, sysm = system([mol], xyz)
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(np.asarray(xyz) * 0.1)
    mm.LocalEnergyMinimizer.minimize(ctx, tol, 0)
    x = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    c = clash(mol, x)
    assert c is None, f"atoms {c:.2f} A apart after the minimisation"
    return x


def minimise(mols, xyz, box, tol=10.0):
    """the energy minimum of a periodic cell (Foyer's force field, PME; the box fixed) next to xyz; the chains stay whole"""
    import openmm as mm
    from openmm import unit
    _, sysm = system(mols, xyz, box)
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001), *platform())
    ctx.setPeriodicBoxVectors(*(mm.Vec3(*(np.eye(3)[i] * box[i] * 0.1)) for i in range(3)))
    ctx.setPositions(np.asarray(xyz) * 0.1)
    mm.LocalEnergyMinimizer.minimize(ctx, tol, 0)
    return ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)


def shake_up(mol, xyz, temp=600.0, ps=30.0, seed=0, radius=None):
    """a random conformation of a chain: Langevin MD in vacuum (Foyer's force field) from the extended chain; with radius
    (A), the atoms are kept within that distance of the centre (flat-bottom restraint), so that the chain becomes compact"""
    import openmm as mm
    from openmm import unit
    _, sysm = system([mol], xyz)
    if radius:
        c = np.asarray(xyz).mean(0) * 0.1
        f = mm.CustomExternalForce("k*step(d-R)*(d-R)^2; d=sqrt((x-cx)^2+(y-cy)^2+(z-cz)^2)")
        for k, v in (("k", 500.0), ("R", radius * 0.1), ("cx", c[0]), ("cy", c[1]), ("cz", c[2])):
            f.addGlobalParameter(k, v)
        for i in range(len(xyz)):
            f.addParticle(i, [])
        sysm.addForce(f)
    integ = mm.LangevinMiddleIntegrator(temp * unit.kelvin, 1.0 / unit.picosecond, 1.0 * unit.femtoseconds)
    integ.setRandomNumberSeed(seed + 1)
    ctx = mm.Context(sysm, integ, mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(np.asarray(xyz) * 0.1)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 500)
    ctx.setVelocitiesToTemperature(temp * unit.kelvin, seed + 1)
    integ.step(int(ps * 1000))
    x = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    return x - x.mean(0)


def frames(mol, xyz, temp=413.0, every_ps=10.0, n=5, seed=0):
    """n frames, every every_ps, of Langevin MD in vacuum (Foyer's force field) from a minimised chain: the conformations a
    chain takes at temp (it folds up in vacuum, as in a melt)"""
    import openmm as mm
    from openmm import unit
    _, sysm = system([mol], xyz)
    integ = mm.LangevinMiddleIntegrator(temp * unit.kelvin, 1.0 / unit.picosecond, 1.0 * unit.femtoseconds)
    integ.setRandomNumberSeed(seed + 1)
    ctx = mm.Context(sysm, integ, mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(np.asarray(xyz) * 0.1)
    ctx.setVelocitiesToTemperature(temp * unit.kelvin, seed + 1)
    out = []
    for _ in range(n):
        integ.step(int(every_ps * 1000))
        x = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        out.append(x - x.mean(0))
    return out


def cell(smiles, density, seed=0, tol=2.0, temp=600.0, ps=15.0):
    """an amorphous cell of the chains `smiles` (one copy each) at the density (g/cm3): every chain made compact by shake_up
    (within a sphere of 1.15 x the radius it has at 0.8 g/cm3), then packed by packmol with periodic boundaries, and the
    cell minimised (minimise: the chains keep their compact shapes, the bonds and angles of the 600 K MD relax).  Returns
    mols, xyz (A, every chain whole), box (A), smiles of the 3D structures"""
    mols, confs, smi3 = [], [], []
    for k, smi in enumerate(smiles):
        mol, xyz, s3 = embed(smi)
        r = 1.15 * (3.0 * mass([mol]) / (0.8 * 0.60221408) / (4.0 * np.pi)) ** (1 / 3)
        mols.append(mol)
        confs.append(shake_up(mol, relax(mol, xyz), temp, ps, seed=seed * 100 + k, radius=r))
        smi3.append(s3)
    L = (mass(mols) / density / 0.60221408) ** (1 / 3)
    with tempfile.TemporaryDirectory() as d:
        inp = [f"tolerance {tol}", "filetype xyz", "output out.xyz", f"seed {seed + 1}", f"pbc 0. 0. 0. {L:.4f} {L:.4f} {L:.4f}"]
        for i, ((el, _), x) in enumerate(zip(mols, confs)):
            open(os.path.join(d, f"m{i}.xyz"), "w").write(
                f"{len(el)}\n\n" + "".join(f"{e} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n" for e, p in zip(el, x)))
            inp += [f"structure m{i}.xyz", "  number 1", "end structure"]
        open(os.path.join(d, "in.inp"), "w").write("\n".join(inp) + "\n")
        with open(os.path.join(d, "in.inp")) as f:
            r = subprocess.run(["packmol"], stdin=f, cwd=d, capture_output=True, text=True)
        if not os.path.isfile(os.path.join(d, "out.xyz")):
            raise RuntimeError("packmol failed:\n" + r.stdout[-2000:])
        xyz = np.array([[float(v) for v in l.split()[1:4]] for l in open(os.path.join(d, "out.xyz")).read().splitlines()[2:]])
    box = np.array([L, L, L])
    return mols, minimise(mols, xyz, box), box, smi3      # every chain whole (it may stick out of the box)
