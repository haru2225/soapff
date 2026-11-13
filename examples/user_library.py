"""extending a library, as a user would: a library that extends amber14sb with one reference (methanol, types and charges only;
ff14SB has the parameters of CT, H1, OH, HO), made with the command line and the Python API, then a methanol / water mixture
read from a CIF file.   python user_library.py WORKDIR"""
import os, subprocess, sys
import numpy as np
from ase.build import molecule
from ase import Atoms
import soapff
W = sys.argv[1]
os.environ["SOAPFF_LIBRARY_PATH"] = os.path.join(W, "libs")
os.makedirs(os.environ["SOAPFF_LIBRARY_PATH"], exist_ok=True)
PY = [sys.executable, "-m", "soapff"]
lib = os.path.join(W, "libs", "methanol.txt")
if os.path.exists(lib):
    os.remove(lib)
print(subprocess.run(PY + ["new", "methanol", "--extends", "amber14sb", "--folder", os.path.join(W, "libs"),
                           "--description", "amber14sb + methanol"], capture_output=True, text=True, check=True).stdout.strip())
# the reference: methanol with its types and charges (RESP-like, sum 0)
m = molecule("CH3OH")
sym = m.get_chemical_symbols()
types, q = [], []
for i, s in enumerate(sym):
    if s == "C": types.append("CT"); q.append(0.1166)
    elif s == "O": types.append("OH"); q.append(-0.6497)
    else:
        d = np.linalg.norm(m.positions[i] - m.positions[[k for k, x in enumerate(sym) if x == "O"][0]])
        if d < 1.2: types.append("HO"); q.append(0.4215)
        else: types.append("H1"); q.append(0.0372)
print("reference", list(zip(sym, types, q)), "sum", round(sum(q), 6))
ref = soapff.add_reference("methanol", m, types=types, charges=q, name="methanol.extxyz")
print("added", ref)
print(subprocess.run(PY + ["info", "methanol", "--show", "0"], capture_output=True, text=True).stdout)
# a refusal: adding to a bundled library
try:
    soapff.add_reference("amber14sb", ref)
except PermissionError as e:
    print("bundled library refused:", str(e)[:80], "...")
# a mixture: 4 methanol + 20 water molecules on a grid in a 12 A box
wat = molecule("H2O")
at = Atoms(cell=[14.0, 14.0, 14.0], pbc=True)
rng = np.random.default_rng(1)
k = 0
for ix in range(3):
    for iy in range(3):
        for iz in range(3):
            mol = (m if k < 4 else wat).copy()
            mol.rotate(rng.uniform(0, 360), rng.normal(size=3))
            mol.translate(np.array([ix, iy, iz]) * 4.6 + 2.3 - mol.positions.mean(0))
            at += mol
            k += 1
at.write(os.path.join(W, "mix.xyz"))        # a plain xyz file (no cell): the box comes from --padding
at.write(os.path.join(W, "mix.cif"))
s = soapff.build("methanol", os.path.join(W, "mix.cif"), os.path.join(W, "mix_out"), md=False, log=lambda *a: None)
tt, qq = np.array(s.types), s.charges
print(s)
print("types", sorted(set(zip(at.get_chemical_symbols(), tt, np.round(qq, 4)))))
print("report", {k: s.report[k] for k in ("molecules", "molecule_kinds", "sources", "terms", "warnings")})
# the same with the parent library: methanol is not in amber14sb -> TypingError
try:
    soapff.build("amber14sb", os.path.join(W, "mix.cif"), os.path.join(W, "mix_parent"), md=False, log=lambda *a: None)
    print("parent: no error ?!")
except soapff.TypingError as e:
    print("parent library refuses:", str(e)[:160], "...")
