"""The library files of soapff for AMBER ff14SB (amber14sb.txt), TIP3P (tip3p.txt) and SPC/E (spce.txt) water with their ions,
from the OpenMM force-field files (amber14/protein.ff14SB.xml, amber14/tip3p.xml, amber14/spce.xml), and the water / ion
references of tip3p and spce.

Units: OpenMM kJ/mol, nm, rad -> kcal/mol, A, deg; OpenMM harmonic E = k/2 (x - x0)^2 -> LAMMPS harmonic E = K (x - x0)^2 with
K = k/2; periodic torsion k (1 + cos(n phi - phase)) -> K = k (dihedral fourier / charmm), improper -> cvff
K [1 + d cos(n phi)] with d = -1 for phase 180 (+1 for phase 0).

    python make_library_txt.py LIBDIR REFDIR_OF_AMBER14SB      (OpenMM)"""
import math
import os
import shutil
import sys
import xml.etree.ElementTree as ET

import openmm as mm
from openmm import app, unit

from omm_common import typed, write_extxyz

DATA = os.path.join(os.path.dirname(app.__file__), "data", "amber14")
KJ = 4.184
OUT = sys.argv[1]
AMBER_REFS = sys.argv[2]
IONS = ["Li+", "Na+", "K+", "Rb+", "Cs+", "F-", "Cl-", "Br-", "I-", "Mg2+", "Ca2+", "Zn2+", "Sr2+", "Ba2+"]
ION_EL = {"Li+": "Li", "Na+": "Na", "K+": "K", "Rb+": "Rb", "Cs+": "Cs", "F-": "F", "Cl-": "Cl", "Br-": "Br", "I-": "I",
          "Mg2+": "Mg", "Ca2+": "Ca", "Zn2+": "Zn", "Sr2+": "Sr", "Ba2+": "Ba"}
g = lambda x: f"{x:.10g}"


def name_of(t):
    """OpenMM type name -> soapff type name"""
    if t == "":
        return "X"
    if t.endswith(("p-O", "e-O")) and t.split("-")[0] in ("tip3p", "spce"):
        return "OW"
    if t.endswith(("p-H", "e-H")) and t.split("-")[0] in ("tip3p", "spce"):
        return "HW"
    if "_standard-" in t:
        return t.split("-", 1)[1]
    if t.startswith("protein-"):
        return t[len("protein-"):]
    return t


def tables(xml, keep_ions=True):
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
    B, A, D, I = [], [], [], []
    for b in root.find("HarmonicBondForce") if root.find("HarmonicBondForce") is not None else []:
        B.append(((name_of(b.get("type1")), name_of(b.get("type2"))), (float(b.get("k")) / 2 / KJ / 100.0, float(b.get("length")) * 10)))
    for a in root.find("HarmonicAngleForce") if root.find("HarmonicAngleForce") is not None else []:
        A.append(((name_of(a.get("type1")), name_of(a.get("type2")), name_of(a.get("type3"))),
                  (float(a.get("k")) / 2 / KJ, math.degrees(float(a.get("angle"))))))
    tf = root.find("PeriodicTorsionForce")
    if tf is not None:
        for p in tf.findall("Proper"):
            terms, i = [], 1
            while p.get(f"k{i}") is not None:
                terms.append((float(p.get(f"k{i}")) / KJ, int(p.get(f"periodicity{i}")), math.degrees(float(p.get(f"phase{i}")))))
                i += 1
            D.append((tuple(name_of(p.get(f"type{k}")) for k in range(1, 5)), terms))
        for p in tf.findall("Improper"):
            assert p.get("k2") is None
            ph = math.degrees(float(p.get("phase1")))
            assert abs(ph - 180) < 1e-6 or abs(ph) < 1e-6, ph
            c, t2, t3, t4 = (name_of(p.get(f"type{k}")) for k in range(1, 5))
            # OpenMM: type1 = the central atom; the torsion is (atom of type2, atom of type3, central, atom of type4)
            I.append(((t2, t3, c, t4), (float(p.get("k1")) / KJ, -1 if abs(ph - 180) < 1e-6 else 1, int(p.get("periodicity1")))))
    return T, B, A, D, I


def write_tables(f, T, B, A, D, I):
    f.write("[types]\n# type  element  mass  epsilon(kcal/mol)  sigma(A)\n")
    for n, v in T.items():
        f.write(f"{n:6s} {v['element']:3s} {g(v['mass']):>9s} {g(v['epsilon']):>14s} {g(v['sigma']):>14s}\n")
    f.write("\n[bonds]\n# t1 t2  K(kcal/mol/A^2)  r0(A)       E = K (r - r0)^2\n")
    for k, v in B:
        f.write(f"{k[0]:5s} {k[1]:5s} {g(v[0]):>12s} {g(v[1]):>10s}\n")
    f.write("\n[angles]\n# t1 t2 t3  K(kcal/mol/rad^2)  theta0(deg)       E = K (theta - theta0)^2\n")
    for k, v in A:
        f.write(f"{k[0]:5s} {k[1]:5s} {k[2]:5s} {g(v[0]):>12s} {g(v[1]):>12s}\n")
    if D:
        f.write("\n[dihedrals]\n# t1 t2 t3 t4  K1 n1 d1 [K2 n2 d2 ...]   E = sum K [1 + cos(n phi - d)] (kcal/mol, deg); X = any type\n")
        for k, terms in D:
            f.write(" ".join(f"{x:5s}" for x in k) + "  " + "  ".join(f"{g(K)} {n} {g(d)}" for K, n, d in terms) + "\n")
    if I:
        f.write("\n[impropers]\n# t1 t2 tc t4  K d n   tc = the central atom (bonded to the three others); E = K [1 + d cos(n phi)]\n")
        for k, v in I:
            f.write(" ".join(f"{x:5s}" for x in k) + f"  {g(v[0])} {v[1]} {v[2]}\n")


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

[bonding]
factor = 1.15
ions   = Li Na K Rb Cs F Cl Br I Mg Ca Zn Sr Ba
"""


def water_refs(model, ff, folder):
    os.makedirs(folder, exist_ok=True)
    mod = app.Modeller(app.Topology(), [])
    mod.addSolvent(ff, model=model, boxSize=mm.Vec3(2.0, 2.0, 2.0) * unit.nanometer, positiveIon="Na+", negativeIon="Cl-",
                   ionicStrength=0.5 * unit.molar)
    files = [f"{model}_NaCl.extxyz"]
    write_extxyz(os.path.join(folder, files[0]), typed_ff(ff, mod.topology, mod.positions), f"{model} water, NaCl 0.5 M (OpenMM)")
    for ion in IONS:
        el = ION_EL[ion]
        top = app.Topology()
        res = top.addResidue(el.upper(), top.addChain())
        top.addAtom(el.upper(), app.Element.getBySymbol(el), res)
        f = f"ion_{el}.extxyz"
        write_extxyz(os.path.join(folder, f), typed_ff(ff, top, [mm.Vec3(0, 0, 0)] * unit.nanometer), f"{ion} ({model} ion set)")
        files.append(f)
    return files


def typed_ff(ff, top, pos):
    import omm_common
    old = omm_common.FF
    omm_common.FF = ff
    try:
        d = typed(top, pos)
    finally:
        omm_common.FF = old
    d["fftype"] = [name_of(x) if not x.startswith(("tip3p", "spce")) else name_of(x) for x in d["fftype"]]
    return d


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
    T, B, A, D, I = tables(xml)
    with open(os.path.join(OUT, f"{model}.txt"), "w") as f:
        f.write(f"""# soapff library: {title}.
# Parameters converted from OpenMM {mm.__version__} (amber14/{xml}) by examples/amber14sb/make_library_txt.py.
# The water is rigid (SHAKE on the O-H bonds and the H-O-H angle during the MD; flexible with these force constants in the
# minimisation).  Every ion is an atom without bonds: its element gives its type.

[library]
name        = {model}
engine      = molecular
description = {title}
reference   = {refs}
source      = OpenMM {mm.__version__} amber14/{xml}

{SOAP}
{LAMMPS}
""")
        write_tables(f, T, B, A, D, I)
        f.write("\n[references]\n" + "".join(f"{model}/{x}\n" for x in files))
    print(model, len(T), "types", len(B), "bonds", len(A), "angles")

# amber14sb: the protein on top of tip3p (extends)
T, B, A, D, I = tables("protein.ff14SB.xml")
os.makedirs(os.path.join(OUT, "amber14sb"), exist_ok=True)
refs = sorted(x for x in os.listdir(AMBER_REFS) if x.endswith(".extxyz") and not x.startswith(("ion_", "tip3p_")))
for x in refs:
    shutil.copyfile(os.path.join(AMBER_REFS, x), os.path.join(OUT, "amber14sb", x))
with open(os.path.join(OUT, "amber14sb.txt"), "w") as f:
    f.write(f"""# soapff library: AMBER ff14SB proteins (Maier et al. 2015) with TIP3P water and ions (extends tip3p).
# Parameters converted from OpenMM {mm.__version__} (amber14/protein.ff14SB.xml) by examples/amber14sb/make_library_txt.py.
#
# Types: the SOAP of the covalent neighbourhood of every atom (atoms within 4 bonds).  Charges: ff14SB gives every atom of every
# residue its own charge, so the molecule is cut at the peptide bonds (C-N, [graph] cut) and every atom takes the type and the
# charge of the reference atom whose residue (its atoms and their bond distances) is the same; similarity 1 means the residue is
# in the references.  References: hen egg-white lysozyme (1AKI) and T4 lysozyme (2LZM) prepared by pdbfixer / OpenMM at pH 7,
# every amino acid (and HID, HIE, HIP) as an N- and a C-terminus, ACE / NME caps, HID HIE HIP ASH GLH LYN CYS in capped
# fragments, disulfides (CYX).  Not included: CYM, HYP, NHE caps, N-terminal CYX, non-standard residues, ligands.
# 1-4 interactions: special_bonds amber (LJ 1/2, Coulomb 1/1.2); LJ mixing: Lorentz-Berthelot (pair_modify mix arithmetic).

[library]
name        = amber14sb
extends     = tip3p
engine      = molecular
description = AMBER ff14SB proteins (20 amino acids, HID/HIE/HIP, ASH/GLH/LYN, CYS/CYX, termini, ACE/NME) with TIP3P water and ions
reference   = J. A. Maier, C. Martinez, K. Kasavajhala, L. Wickstrom, K. E. Hauser, C. Simmerling, J. Chem. Theory Comput. 11, 3696 (2015); TIP3P and ions: see tip3p
source      = OpenMM {mm.__version__} amber14/protein.ff14SB.xml

[graph]
cut   = C N
depth = 12

""")
    write_tables(f, T, B, A, D, I)
    f.write("\n[references]\n" + "".join(f"amber14sb/{x}\n" for x in refs))
print("amber14sb", len(T), "types", len(B), "bonds", len(A), "angles", len(D), "propers", len(I), "impropers", len(refs), "references")
