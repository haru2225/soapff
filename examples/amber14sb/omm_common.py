"""OpenMM side of the protein references: PDB -> pdbfixer -> hydrogens at pH 7 -> AMBER ff14SB / TIP3P templates ->
per-atom AMBER atom class and charge (as OpenMM assigns them)."""
import io
import os

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

FF = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
HERE = os.path.dirname(os.path.abspath(__file__))
PDB = os.path.join(HERE, "pdb")


def rename_class(c):
    """OpenMM class names -> the names of the soapff library"""
    if c == "tip3p-O":
        return "OW"
    if c == "tip3p-H":
        return "HW"
    if c.startswith("tip3p_standard-"):
        return c.split("-", 1)[1]
    return c


def pdb_file(pdbid):
    """pdb/<ID>.pdb, downloaded from the RCSB PDB the first time"""
    path = os.path.join(PDB, pdbid + ".pdb")
    if not os.path.isfile(path):
        import urllib.request
        os.makedirs(PDB, exist_ok=True)
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdbid}.pdb", path)
    return path


def atom_lines(pdbid, chain="A"):
    """ATOM records of the first model of a chain, altloc ' ' or 'A' (no heterogens, no water)"""
    out = []
    for l in open(pdb_file(pdbid)):
        if l.startswith("ENDMDL"):
            break
        # the hydrogens of the file (NMR, neutron) are left out: OpenMM adds them with its own names
        if l.startswith("ATOM") and l[21] == chain and l[16] in " A" and l[76:78].strip() not in ("H", "D"):
            out.append(l[:16] + " " + l[17:])
    return out


def residue_numbers(lines):
    seen = []
    for l in lines:
        r = int(l[22:26])
        if not seen or seen[-1][0] != r:
            seen.append((r, l[17:20]))
    return seen


def fragment(lines, first, last, ace=False, nme=False):
    """the residues first..last; ace: residue first-1 kept as an ACE cap (its CA -> CH3, C, O), nme: residue last+1 kept as
    an NME cap (its N, CA -> C)"""
    out = []
    for l in lines:
        r, name = int(l[22:26]), l[12:16].strip()
        if first <= r <= last:
            out.append(l)
        elif ace and r == first - 1 and name in ("CA", "C", "O"):
            out.append(l[:12] + (" CH3" if name == "CA" else l[12:16]) + l[16:17] + "ACE" + l[20:])
        elif nme and r == last + 1 and name in ("N", "CA"):
            out.append(l[:12] + (" C  " if name == "CA" else l[12:16]) + l[16:17] + "NME" + l[20:])
    return "".join(out) + "END\n"


def prepare(pdb_text, variant=None, pH=7.0):
    """pdbfixer (missing heavy atoms, no missing loops) and OpenMM hydrogens; variant: {residue name: variant} e.g.
    {"HIS": "HID"} forced on every residue of that name"""
    fx = PDBFixer(pdbfile=io.StringIO(pdb_text))
    # ATOM records only were given (no heterogens); ACE / NME caps must stay
    fx.findMissingResidues()
    fx.missingResidues = {}
    fx.findNonstandardResidues()
    fx.replaceNonstandardResidues()
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    mod = app.Modeller(fx.topology, fx.positions)
    var = None
    if variant:
        var = [variant.get(r.name) for r in mod.topology.residues()]
        # a free cysteine only (a disulfide stays CYX)
        ss = {a.residue for b in mod.topology.bonds() for a in b if b[0].name == "SG" and b[1].name == "SG"}
        var = [None if (r in ss and v in ("CYM", "CYS")) else v for r, v in zip(mod.topology.residues(), var)]
    mod.addHydrogens(FF, pH=pH, variants=var)
    return mod.topology, mod.positions


def typed(top, pos):
    """(symbols, positions A, class, charge, template name, atom name, residue index) of every atom"""
    data = FF._SystemData(top)
    tmpl = FF._matchAllResiduesToTemplates(data, top, {}, False)
    atoms = list(top.atoms())
    cls = [rename_class(FF._atomTypes[data.atomType[a]].atomClass) for a in atoms]
    q = [float(data.atomParameters[a]["charge"]) for a in atoms]
    res = [tmpl[a.residue].name for a in atoms]
    xyz = np.array(pos.value_in_unit(unit.angstrom))
    return dict(symbols=[a.element.symbol for a in atoms], positions=xyz, fftype=cls, ffcharge=q, template=res,
                atomname=[a.name for a in atoms], resid=[a.residue.index for a in atoms])


def write_extxyz(path, d, comment=""):
    """extxyz without a cell (a molecule): species, positions, fftype, ffcharge, and the template / atom name / residue for
    the reader"""
    n = len(d["symbols"])
    L = [str(n), f'Properties=species:S:1:pos:R:3:fftype:S:1:ffcharge:R:1:template:S:1:atomname:S:1:resid:I:1 pbc="F F F" '
                 f'comment="{comment}"']
    for i in range(n):
        x, y, z = d["positions"][i]
        L.append(f"{d['symbols'][i]:2s} {x:12.6f} {y:12.6f} {z:12.6f} {d['fftype'][i]:4s} {d['ffcharge'][i]:10.6f} "
                 f"{d['template'][i]:5s} {d['atomname'][i]:5s} {d['resid'][i]:5d}")
    open(path, "w").write("\n".join(L) + "\n")


def system_vacuum(top, rigid_water=False):
    """the OpenMM System of the topology, no cutoff, no constraints (for the energy comparison with LAMMPS)"""
    return FF.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=rigid_water,
                           removeCMMotion=False)


def energies(system, pos):
    """energy of each force (kcal/mol) at the positions (A)"""
    out = {}
    for k, f in enumerate(system.getForces()):
        f.setForceGroup(k)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(np.asarray(pos) * 0.1)
    for k, f in enumerate(system.getForces()):
        e = ctx.getState(getEnergy=True, groups={k}).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        out[type(f).__name__] = out.get(type(f).__name__, 0.0) + e
    out["total"] = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    return out
