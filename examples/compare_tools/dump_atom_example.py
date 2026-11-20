"""A clay structure made apart from soapff's references: the 3-water-layer Na-montmorillonite that comes with atom
(structures/minerals/3WNaMMT_4x6x2.pdb: 4,128 atoms, monoclinic, Mg in 1 of 4 octahedral sites) -> soapff's types
(clayff) and the input of run_atom.m, as dump_clay.py does for the test structures.  The atoms are put in the order of a
.gro file of a solvated clay (the mineral, then each water molecule O H H, then the ions), the order run_atom.m needs.
usage (soapff's Python): python dump_atom_example.py /path/to/atom IN       -> IN/ext__mmt3w.{npz,txt}, IN/ext__mmt3w_out/"""
import collections, os, sys
import numpy as np
import ase.io
from ase.neighborlist import neighbor_list
import soapff

ATOM, IN = sys.argv[1:3]
os.makedirs(IN, exist_ok=True)
at = ase.io.read(os.path.join(ATOM, "structures", "minerals", "3WNaMMT_4x6x2.pdb"))
s = soapff.build("clayff", at, os.path.join(IN, "ext__mmt3w_out"))
t = np.array(s.types)
print(len(at), dict(collections.Counter(t)), "min similarity", round(s.report["min_similarity"], 4))
hw = collections.defaultdict(list)
for i, j in zip(*neighbor_list("ij", at, 1.25)):
    if t[i] == "ow" and t[j] == "hw":
        hw[i].append(j)
order = [i for i in range(len(at)) if t[i] not in ("ow", "hw", "Na")]
for o in np.nonzero(t == "ow")[0]:
    assert len(hw[o]) == 2
    order += [int(o)] + hw[o]
order += [int(i) for i in np.nonzero(t == "Na")[0]]
assert sorted(order) == list(range(len(at)))
# save() of dump_clay.py writes IN/NAME.npz and IN/NAME.txt
src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dump_clay.py")
sys.argv = [src, IN]
ns = {"__name__": "dump_clay", "__file__": src}
exec(open(src).read().split("for g, sub in")[0], ns)
ns["save"]("ext", "mmt3w", at[order], t[order])
