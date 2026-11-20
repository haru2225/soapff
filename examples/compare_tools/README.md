# soapff and atom on the clay structures

soapff does with one method what Foyer does for OPLS-AA polymers and what atom (Holmboe, *Clays Clay Miner.* **67**, 419
(2019); https://github.com/mholmboe/atom) does for ClayFF minerals. For the polymers, the types and charges of every test
structure are Foyer's (`../oplsaa/`, and the forces against Foyer's force field in OpenMM: `../oplsaa/results/forces.json`).
This folder is the clay side: `clayff_atom` of atom 3.0.2 (run in GNU Octave 10.3) on the 16 `clayff` test structures in the
reference arrangement (the 10 (001) slabs with water and the 6 water shapes; 28,311 atoms), compared with soapff's types and
charges (which are the right ones of the validation table).

```
python dump_clay.py IN                                   # soapff's Python: IN/NAME.npz (soapff's types), IN/NAME.txt (atom's input)
for f in IN/*.txt; do
  octave-cli run_atom.m /path/to/atom $PWD/$f $PWD/OUT/$(basename $f)    # atom's types and charges
done
python score_atom.py IN OUT                              # -> OUT/score.json (results/score.json)
```

`clayff_atom(atom, Box_dim, 'clayff', 'spc')`: the defaults of atom (no healing, metal-O cutoff 2.25 Å) with SPC water, which
is Cygan's ClayFF. The input is what a .gro file of a solvated clay gives atom: element names for the mineral, which is one
molecule; each water molecule (`SOL`, `Ow Hw Hw`) and each ion (`ION`) its own molecule. atom's names are mapped to Cygan's with
the correspondence in atom's `clayff_param.m` (`MAP` in `score_atom.py`).

Result (`results/score.json`): every one of the 28,311 atoms has atom's charge, water included, and 27,939 also have atom's
type. The other 372 are O bonded to a tetrahedral Al (obts in Cygan's definition, which soapff gives): atom computes the charge
of an O from the cations bonded to it and gives them the obts charge, but the name `Ob` or `Op` (an apical O next to the
tetrahedral Al; atom's neighbour string is made of the element names, so its `Ops` rule for `AlAlAlt` does not match) or
`Oalsi` (a basal O whose tetrahedral Al comes after it in the atom list, so that atom has not yet renamed that Al `Alt`
when it types the O).

The NaCl droplet of Meng et al. (105,977 atoms) is left out: atom's neighbour search keeps N × N distance matrices, which
do not fit in the 16 GB of the machine used.

## A structure made apart from the references

atom comes with a 3-water-layer Na-montmorillonite, `structures/minerals/3WNaMMT_4x6x2.pdb` (4,128 atoms: 1,920 of the
mineral, 720 water molecules, 48 Na⁺; a monoclinic cell, Mg in 1 of 4 octahedral sites, where the one montmorillonite of
`clayff` has 1 in 8). Nothing of it is in soapff's references.

```
python dump_atom_example.py /path/to/atom IN             # soapff's types (clayff) and atom's input: IN/ext__mmt3w.{npz,txt}
octave-cli run_atom.m /path/to/atom $PWD/IN/ext__mmt3w.txt $PWD/OUT/ext__mmt3w.txt
python score_atom.py IN OUT                              # -> OUT/score.json (results/score_atom_example.json)
```

soapff types it without stopping (smallest similarity 0.951, smallest margin 0.010): 48 mgo with 192 obos and 96 ohs next
to them, and every one of the 4,128 atoms has atom's type and charge.
