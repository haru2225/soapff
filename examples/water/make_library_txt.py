"""The library files of soapff for TIP3P (tip3p.txt) and SPC/E (spce.txt) water with their ions, from the OpenMM force-field
files amber14/tip3p.xml and amber14/spce.xml, and their references (a water box with NaCl, and each ion alone).

    python make_library_txt.py LIBDIR      (a Python with OpenMM)  ->  LIBDIR/tip3p.txt, LIBDIR/spce.txt, LIBDIR/tip3p/, LIBDIR/spce/

Units: OpenMM kJ/mol, nm, rad -> kcal/mol, A, deg; OpenMM harmonic E = k/2 (x - x0)^2 -> LAMMPS harmonic E = K (x - x0)^2 with
K = k/2.  The parameter tables are those of the XML files.  The water box of the references is made by OpenMM's Modeller,
which picks the waters it replaces by ions at random, so the coordinates of a new tip3p_NaCl.extxyz differ from the bundled
one (the types and charges do not)."""
import math
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import openmm as mm
from openmm import app, unit

DATA = os.path.join(os.path.dirname(app.__file__), "data", "amber14")
KJ = 4.184
OUT = sys.argv[1]
IONS = ["Li+", "Na+", "K+", "Rb+", "Cs+", "F-", "Cl-", "Br-", "I-", "Mg2+", "Ca2+", "Zn2+", "Sr2+", "Ba2+"]
ION_EL = {"Li+": "Li", "Na+": "Na", "K+": "K", "Rb+": "Rb", "Cs+": "Cs", "F-": "F", "Cl-": "Cl", "Br-": "Br", "I-": "I",
          "Mg2+": "Mg", "Ca2+": "Ca", "Zn2+": "Zn", "Sr2+": "Sr", "Ba2+": "Ba"}
g = lambda x: f"{x:.10g}"


def name_of(t):
    """OpenMM type or class name -> soapff type name (OW, HW, the ion names)"""
    if t.split("-")[0] in ("tip3p", "spce") and t.endswith("-O"):
        return "OW"
    if t.split("-")[0] in ("tip3p", "spce") and t.endswith("-H"):
        return "HW"
    if "_standard-" in t:
        return t.split("-", 1)[1]
    return t


def tables(xml):
    root = ET.parse(os.path.join(DATA, xml)).getroot()
    T = {}
    for t in root.find("AtomTypes"):
        n = name_of(t.get("name"))
        if "_standard-" in t.get("name") and n not in IONS:
            continue
        T[n] = dict(element=t.get("element"), mass=float(t.get("mass")))
    for a in root.find("NonbondedForce").findall("Atom"):
        n = name_of(a.get("type"))
        if n in T:
            T[n]["epsilon"] = float(a.get("epsilon")) / KJ
            T[n]["sigma"] = float(a.get("sigma")) * 10.0
    B, A = [], []
    for b in root.find("HarmonicBondForce") if root.find("HarmonicBondForce") is not None else []:
        B.append(((name_of(b.get("type1")), name_of(b.get("type2"))), (float(b.get("k")) / 2 / KJ / 100.0, float(b.get("length")) * 10)))
    for a in root.find("HarmonicAngleForce") if root.find("HarmonicAngleForce") is not None else []:
        A.append(((name_of(a.get("type1")), name_of(a.get("type2")), name_of(a.get("type3"))),
                  (float(a.get("k")) / 2 / KJ, math.degrees(float(a.get("angle"))))))
    return T, B, A


def write_tables(f, T, B, A):
    f.write("[types]\n# type  element  mass  epsilon(kcal/mol)  sigma(A)\n")
    for n, v in T.items():
        f.write(f"{n:6s} {v['element']:3s} {g(v['mass']):>9s} {g(v['epsilon']):>14s} {g(v['sigma']):>14s}\n")
    f.write("\n[bonds]\n# t1 t2  K(kcal/mol/A^2)  r0(A)       E = K (r - r0)^2\n")
    for k, v in B:
        f.write(f"{k[0]:5s} {k[1]:5s} {g(v[0]):>12s} {g(v[1]):>10s}\n")
    f.write("\n[angles]\n# t1 t2 t3  K(kcal/mol/rad^2)  theta0(deg)       E = K (theta - theta0)^2\n")
    for k, v in A:
        f.write(f"{k[0]:5s} {k[1]:5s} {k[2]:5s} {g(v[0]):>12s} {g(v[1]):>12s}\n")


LAMMPS = """[lammps]
pair_style     = lj/cut/coul/long 10.0
pair_modify    = mix arithmetic
kspace_style   = pppm 1.0e-5
bond_style     = harmonic
angle_style    = harmonic
dihedral_style = charmm
improper_style = cvff
special_bonds  = amber
rigid          = OW HW
timestep       = 1.0
"""
SOAP = """[soap]
r_cut = 5.0
n_max = 4
l_max = 3
sigma = 0.4
depth = 4
prune = 0.999
threshold = 0.9

[bond_lengths]
# E1 E2 r(A).  H-O by hand: in MD at 300 K the H...O of the hydrogen bonds comes to 1.43-1.50 A and the O-H bonds stretch
# to 1.12 A (flexible water); the other pairs from the references (halfway)
H O 1.25
"""


def typed(ff, top, pos):
    """symbols, positions (A), soapff type and charge of every atom, as OpenMM assigns them, with the template, atom name and
    residue index for the reader"""
    data = ff._SystemData(top)
    tmpl = ff._matchAllResiduesToTemplates(data, top, {}, False)
    atoms = list(top.atoms())
    return dict(symbols=[a.element.symbol for a in atoms], positions=np.array(pos.value_in_unit(unit.angstrom)),
                fftype=[name_of(ff._atomTypes[data.atomType[a]].atomClass) for a in atoms],
                ffcharge=[float(data.atomParameters[a]["charge"]) for a in atoms],
                template=[tmpl[a.residue].name for a in atoms], atomname=[a.name for a in atoms],
                resid=[a.residue.index for a in atoms], bonds=[(b[0].index, b[1].index) for b in top.bonds()])


def write_extxyz(path, d, comment=""):
    """extxyz without a cell: species, positions, fftype, ffcharge, the template / atom name / residue for the reader, and
    the bonds (the 1-based ids of the bonded atoms, '-' for none) that soapff needs in a reference"""
    n = len(d["symbols"])
    nbr = [[] for _ in range(n)]
    for i, j in d["bonds"]:
        nbr[i].append(j + 1); nbr[j].append(i + 1)
    L = [str(n), f'Properties=species:S:1:pos:R:3:fftype:S:1:ffcharge:R:1:template:S:1:atomname:S:1:resid:I:1:bonds:S:1 '
                 f'pbc="F F F" comment="{comment}"']
    for i in range(n):
        x, y, z = d["positions"][i]
        L.append(f"{d['symbols'][i]:2s} {x:12.6f} {y:12.6f} {z:12.6f} {d['fftype'][i]:4s} {d['ffcharge'][i]:10.6f} "
                 f"{d['template'][i]:5s} {d['atomname'][i]:5s} {d['resid'][i]:5d}  {','.join(str(x) for x in sorted(nbr[i])) or '-'}")
    open(path, "w").write("\n".join(L) + "\n")


def water_refs(model, ff, folder):
    os.makedirs(folder, exist_ok=True)
    mod = app.Modeller(app.Topology(), [])
    mod.addSolvent(ff, model=model, boxSize=mm.Vec3(2.0, 2.0, 2.0) * unit.nanometer, positiveIon="Na+", negativeIon="Cl-",
                   ionicStrength=0.5 * unit.molar)
    files = [f"{model}_NaCl.extxyz"]
    write_extxyz(os.path.join(folder, files[0]), typed(ff, mod.topology, mod.positions), f"{model} water, NaCl 0.5 M (OpenMM)")
    for ion in IONS:
        el = ION_EL[ion]
        top = app.Topology()
        res = top.addResidue(el.upper(), top.addChain())
        top.addAtom(el.upper(), app.Element.getBySymbol(el), res)
        f = f"ion_{el}.extxyz"
        write_extxyz(os.path.join(folder, f), typed(ff, top, [mm.Vec3(0, 0, 0)] * unit.nanometer), f"{ion} ({model} ion set)")
        files.append(f)
    return files


os.makedirs(OUT, exist_ok=True)
for model, xml, title, refs in (
        ("tip3p", "tip3p.xml", "TIP3P water (Jorgensen et al. 1983, rigid) with the ions of OpenMM amber14/tip3p.xml: Joung-Cheatham "
                               "monovalent ions (2008) and Li-Merz 12-6 divalent ions (2013)",
         "W. L. Jorgensen et al., J. Chem. Phys. 79, 926 (1983); I. S. Joung, T. E. Cheatham, J. Phys. Chem. B 112, 9020 (2008); "
         "P. Li, B. P. Roberts, D. K. Chakravorty, K. M. Merz, J. Chem. Theory Comput. 9, 2733 (2013)"),
        ("spce", "spce.xml", "SPC/E water (Berendsen et al. 1987, rigid) with the ions of OpenMM amber14/spce.xml: Joung-Cheatham "
                             "monovalent ions for SPC/E (2008) and Li-Merz 12-6 divalent ions (2013)",
         "H. J. C. Berendsen, J. R. Grigera, T. P. Straatsma, J. Phys. Chem. 91, 6269 (1987); I. S. Joung, T. E. Cheatham, "
         "J. Phys. Chem. B 112, 9020 (2008); P. Li et al., J. Chem. Theory Comput. 9, 2733 (2013)")):
    ff = app.ForceField(f"amber14/{xml}")
    files = water_refs(model, ff, os.path.join(OUT, model))
    T, B, A = tables(xml)
    with open(os.path.join(OUT, f"{model}.txt"), "w") as f:
        f.write(f"""# soapff library: {title}.
# Parameters converted from OpenMM {mm.__version__} (amber14/{xml}) by examples/water/make_library_txt.py.
# The water is rigid: system.settings ends with fix shake on the O-H bonds and the H-O-H angle (in a minimisation LAMMPS
# turns the constraints into stiff restraints).  The force constants below are used without it (System.commands(rigid=False)).
# Every ion is an atom without bonds: its element gives its type.

[library]
name        = {model}
engine      = molecular
description = {title}
reference   = {refs}
source      = OpenMM {mm.__version__} amber14/{xml}

{SOAP}
{LAMMPS}
""")
        write_tables(f, T, B, A)
        f.write("\n[references]\n" + "".join(f"{model}/{x}\n" for x in files))
    print(model, len(T), "types", len(B), "bonds", len(A), "angles")
