"""The parameters of the soapff library oplsaa from Foyer's OPLS-AA (foyer/forcefields/xml/oplsaa.xml), as text.

    python make_library_txt.py ../../src/soapff/libraries/oplsaa.txt          (needs foyer; the references are added by
                                                                              make_refs.py)
- [types]: every type of the file with an element, its mass, the LJ parameters of NonbondedForce (kJ/mol, nm -> kcal/mol,
  A) and its class (the bonds, angles and dihedrals of OPLS-AA are given by class);
- [bonds], [angles]: harmonic, E = K (x - x0)^2 with K = k / 2 (kJ/mol/nm^2 -> kcal/mol/A^2);
- [dihedrals]: the Ryckaert-Bellemans torsions E = sum_n C_n cos^n(psi), psi = phi - 180 deg, written exactly as the Fourier
  series sum_m K_m [1 + cos(m phi - d_m)] (m = 1..5, d_m = 0 or 180; they differ by a constant per torsion, the forces are the
  same; checked at random angles below); a class '' (any) is X.  The file has no impropers.
OPLS-AA: geometric mixing of sigma and epsilon, 1-4 interactions x 0.5 (LJ and Coulomb)."""
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import foyer

KJ = 4.184
XML = os.path.join(os.path.dirname(foyer.__file__), "forcefields", "xml", "oplsaa.xml")
out = sys.argv[1]
root = ET.parse(XML).getroot()
assert root.get("combining_rule") == "geometric"
nbf = root.find("NonbondedForce")
assert (nbf.get("coulomb14scale"), nbf.get("lj14scale")) == ("0.5", "0.5")
lj = {a.get("type"): (float(a.get("sigma")) * 10.0, float(a.get("epsilon")) / KJ) for a in nbf if a.tag == "Atom"}
assert root.find("PeriodicTorsionForce") is None and not root.find("RBTorsionForce").findall("Improper")

types, classes = [], set()
for t in root.find("AtomTypes"):
    el = t.get("element")
    if not el or t.get("name") not in lj:
        continue
    s, e = lj[t.get("name")]
    types.append((t.get("name"), el, float(t.get("mass")), e, s, t.get("class")))
    classes.add(t.get("class"))
assert "X" not in classes, "a class named X would be read as the wildcard"


def c(x):
    return x if x else "X"


def fourier(C):
    """RB coefficients C0..C5 (kcal/mol, psi = phi - 180) -> [(K, m, d)] with sum_m K [1 + cos(m phi - d)] + const"""
    a = [C[n] * (-1) ** n for n in range(6)]                       # in cos(phi)
    b = {1: a[1] + 0.75 * a[3] + 0.625 * a[5], 2: 0.5 * a[2] + 0.5 * a[4], 3: 0.25 * a[3] + 0.3125 * a[5], 4: 0.125 * a[4],
         5: 0.0625 * a[5]}
    return [(abs(v), m, 0.0 if v > 0 else 180.0) for m, v in b.items() if abs(v) > 1e-12]


rng = np.random.default_rng(0)
phi = rng.uniform(-np.pi, np.pi, 50)


def check(C, terms):
    rb = sum(C[n] * np.cos(phi - np.pi) ** n for n in range(6))
    fo = sum(K * (1 + np.cos(m * phi - np.radians(d))) for K, m, d in terms) if terms else 0 * phi
    diff = rb - fo
    assert np.ptp(diff) < 1e-9 * max(1.0, np.abs(rb).max()), (C, terms)


L = ["# soapff library: OPLS-AA (Jorgensen, Maxwell, Tirado-Rives 1996) as distributed with Foyer (oplsaa.xml), for polymers.",
     f"# Parameters converted from foyer {foyer.__version__} forcefields/xml/oplsaa.xml by examples/oplsaa/make_library_txt.py.",
     "#",
     "# Types and charges in two steps: the bonds from the distance and the bond length of each pair of elements ([bond_lengths];",
     "# by default halfway between the longest bond and the shortest distance of two atoms not bonded in the references, which",
     "# carry the bonds of Foyer), and the type and charge from the SOAP of the covalent neighbourhood of every atom (atoms",
     "# within 4 bonds).  OPLS-AA gives every type one charge, and the type of an atom follows from its neighbours (the SMARTS of",
     "# Foyer look at most 3 bonds away), so the nearest reference atom has the type and the charge.  References: oligomers of",
     "# 11 polymers typed by Foyer (examples/oplsaa/make_refs.py).  Bonds, angles and dihedrals are given by class (the 6th",
     "# column of [types]).",
     "# 1-4 interactions: special_bonds lj/coul 0 0 0.5; LJ mixing: geometric (sigma and epsilon).",
     "",
     "[library]",
     "name      = oplsaa",
     "reference = W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996); types by SMARTS: "
     "C. Klein et al., Comput. Mater. Sci. 167, 215 (2019)",
     "",
     "[soap]", "r_cut = 5.0", "n_max = 4", "l_max = 3", "sigma = 0.4", "depth = 4", "prune = 0.999", "threshold = 0.9", "",
     "[bond_lengths]",
     "# E1 E2 r(A).  H-O by hand: in MD at 300 K the H...O of the hydrogen bonds comes to 1.43-1.50 A (in the references",
     "# not nearer than 1.76 A) and the O-H bonds stretch to 1.12 A; the other pairs from the references (halfway)",
     "H O 1.25", "",
     "[lammps]",
     "pair_style     = lj/cut/coul/long 12.0",
     "pair_modify    = mix geometric tail yes",
     "kspace_style   = pppm 1.0e-5",
     "bond_style     = harmonic",
     "angle_style    = harmonic",
     "dihedral_style = fourier",
     "special_bonds  = lj/coul 0.0 0.0 0.5",
     "timestep       = 1.0",
     "",
     "[types]",
     "# type  element  mass  epsilon(kcal/mol)  sigma(A)  class"]
L += [f"{n:10s} {el:3s} {m:9.4f} {e:12.8g} {s:12.8g}  {k}" for n, el, m, e, s, k in types]
L += ["", "[bonds]", "# c1 c2  K(kcal/mol/A^2)  r0(A)       E = K (r - r0)^2"]
nb = na = nd = 0
for b in root.find("HarmonicBondForce"):
    L.append(f"{c(b.get('class1')):5s} {c(b.get('class2')):5s} {float(b.get('k')) / 2 / KJ / 100:12.8g} "
             f"{float(b.get('length')) * 10:10.6g}")
    nb += 1
L += ["", "[angles]", "# c1 c2 c3  K(kcal/mol/rad^2)  theta0(deg)       E = K (theta - theta0)^2"]
for a in root.find("HarmonicAngleForce"):
    L.append(f"{c(a.get('class1')):5s} {c(a.get('class2')):5s} {c(a.get('class3')):5s} "
             f"{float(a.get('k')) / 2 / KJ:12.8g} {np.degrees(float(a.get('angle'))):12.8g}")
    na += 1
L += ["", "[dihedrals]", "# c1 c2 c3 c4  K1 n1 d1 [K2 n2 d2 ...]   E = sum K [1 + cos(n phi - d)] (kcal/mol, deg); from the "
      "Ryckaert-Bellemans C0..C5 of Foyer", "# (all C = 0: one term K = 0, so that the torsion is known and not reported missing)"]
for p in root.find("RBTorsionForce").findall("Proper"):
    C = [float(p.get(f"c{n}")) / KJ for n in range(6)]
    terms = fourier(C)
    check(C, terms)
    t = terms or [(0.0, 1, 0.0)]
    L.append(" ".join(f"{c(p.get(f'class{i}')):5s}" for i in range(1, 5)) + "  "
             + "  ".join(f"{K:.10g} {m} {d:g}" for K, m, d in t))
    nd += 1
# the references written by make_refs.py (kept when the file exists)
old = open(out).read().split("[references]", 1)[1] if os.path.isfile(out) and "[references]" in open(out).read() else \
    "\n# oligomers typed by Foyer (make_refs.py)\n"
open(out, "w").write("\n".join(L) + "\n\n[references]" + old)
print(f"{out}: {len(types)} types, {len(classes)} classes, {nb} bonds, {na} angles, {nd} dihedrals")
