# soapff

soapff is a Python package that turns a structure file into a force field for LAMMPS. It assigns a force-field type and a charge to every atom. It does this by comparing each atom's local environment (SOAP descriptor) with reference structures whose types and charges are known. It then writes the LAMMPS data and settings files.

soapff stops there. The simulation itself (minimisation, NVT, NPT, output, analysis) is written in LAMMPS, for example with the LAMMPS Python module:

```python
import soapff
from lammps import lammps

s = soapff.build("oplsaa", "pe_melt.extxyz", "pe_out")  # library, structure file (elements, coordinates, cell)
lmp = lammps()
lmp.commands_list(s.commands())                          # the force field
lmp.commands_string("""
    minimize 1.0e-4 1.0e-6 1000 10000
    velocity all create 413 1 dist gaussian
    fix 1 all npt temp 413 413 100 iso 1 1 1000
    run 100000
""")
lmp.get_thermo("density")
```

- **Structure files:** any format ASE reads (cif, xyz, extxyz, pdb, POSCAR, ...), a LAMMPS data file, or an `ase.Atoms`.
- **Libraries:** text files. Each one names the reference structures and gives the force-field parameters. Users can extend a bundled library or write their own.
- **Bundled libraries:** ClayFF (clay minerals), OPLS-AA (polymers: 11 polymers and their copolymers, types and charges as Foyer gives them), TIP3P and SPC/E (water and ions).
- **Checks:** an atom whose environment is not in the library (similarity below the threshold) stops the build with an error that names the atom. So does a molecule whose charge is not a whole number (molecular engine): when the charges belong to the types, as in OPLS-AA and the water models, a mistyped atom of another charge shows there even when its similarity is above the threshold (a junction of a block copolymer, a conformation the references do not cover). The error lists the atoms of that molecule with the smallest margins first. The similarity of every atom, and its margin over the nearest reference with another type or charge, go to `charges.csv`; the atoms with the smallest margins are listed in `report.json`. In the tests of `examples/oplsaa`, soapff either wrote a force field with every atom right or stopped. What these checks cannot see is a confusion of two types with the same charge, and, in the clayff engine, of two types with the same bonded elements and another charge (the octahedral Mg mgo and mgh). The net charge of the structure is written to `report.json` (`charge_sum`) and the log, not checked: a clay layer or particle may carry a charge.

## Installation

```bash
git clone https://github.com/haru2225/soapff.git
pip install -e soapff          # Python >= 3.9; installs numpy, scipy, ase, dscribe
```

soapff itself does not run LAMMPS and does not import it. To run the simulations from Python, install the LAMMPS Python module, e.g. `pip install "lammps[mpi]"` (LAMMPS 22 Jul 2025 with KSPACE, MOLECULE, EXTRA-MOLECULE, RIGID and MPI; this wheel has no OpenMP threads, so a run in parallel is `mpirun -np 4 python script.py`, where the script reads a folder built beforehand with `soapff.System("out").commands()`), or use any LAMMPS build with its Python module (conda, your own build with GPU or MPI). The files soapff writes also work with the `lmp` executable: `include system.settings` in a LAMMPS input.

The clay examples place water around a mineral with [packmol](https://m3g.github.io/packmol/) (`soapff_cui.water`).

Optional extras:
- `pip install -e "soapff[examples]"`: what the examples need to run their MD and plots (the LAMMPS Python module, matplotlib, RDKit).
- `pip install -e "soapff[surfaces]"`: (001) slabs of the minerals (pymatgen).

## Usage

```python
import soapff

soapff.libraries()                                     # the libraries that can be given by name
s = soapff.build("clayff", "kaolinite.cif")            # writes kaolinite_soapff/ in the current folder
s = soapff.build("oplsaa", "melt_PE.extxyz", "out")     # a polymer melt, into out/
s.files                  # system.data, system.settings, charges.csv, report.json in s.out_dir
s.types, s.charges       # the type and charge of every atom
s.commands()             # the force field as LAMMPS commands (system.settings, with the path of system.data)
s.type_names             # the type names in the order of the LAMMPS atom types (e.g. for dump_modify element)
s.report                 # report.json: similarities, margins, terms, warnings, recommended time step
```

The force field contains everything that belongs to the model: units, styles, coefficients, `special_bonds`, mixing, long-range Coulomb, and SHAKE for rigid water (`fix rigid`, at the end of `system.settings`). `s.commands(rigid=False)` leaves SHAKE out, e.g. for a single-point force of the flexible molecule. Everything else is written in LAMMPS. `examples/quickstart.py` runs an NPT this way, `examples/polymer_smiles.py` makes a PMMA chain from its SMILES (RDKit), types it and runs an NVT, and every example in `examples/` runs its MD with the LAMMPS Python module.

Options of `soapff.build`:

- `cutoff`, `kspace` (e.g. `"ewald 1.0e-5"`), `boundary` (e.g. `"p p f"`: MSM or the slab correction is chosen for the non-periodic faces, and atoms are reflected at fixed faces).
- `strict=False`: write the force field even when some atoms are not in the library.
- `padding`: a structure without a cell, such as an xyz file or a PDB file without CRYST1, is put in a box this much larger (Å).

The command `soapff` does the same:

```bash
soapff list
soapff info oplsaa
soapff build oplsaa cell_PE.extxyz -o out      # out/system.data, out/system.settings: include it in a LAMMPS input
soapff new mine --extends oplsaa
soapff add-reference mine PE4_PEO3.extxyz
```

## Libraries

| Library | Engine | References | Force field |
|---|---|---|---|
| `clayff` | clayff | Bulk crystals of 11 layered minerals and their (001) slabs | ClayFF (Cygan et al. 2004), flexible SPC water |
| `clayff_edge` | clayff | Montmorillonite particle with edges | Underwood & Bourg (2020) |
| `oplsaa` | molecular | 88 oligomer structures of 11 polymers (5-mers, 8-mers and MD frames at 413 K), typed by Foyer | OPLS-AA (Jorgensen et al. 1996; from Foyer 1.3.0 `oplsaa.xml`; the bonded terms by class, the Ryckaert–Bellemans torsions as `dihedral_style fourier`) |
| `tip3p` | molecular | Water, NaCl and 14 ions | TIP3P (rigid), Joung–Cheatham 2008 monovalent and Li–Merz 2013 divalent ions |
| `spce` | molecular | The same, for SPC/E | SPC/E (rigid) and its ions |

**Library files.** Each library is a text file (`src/soapff/libraries/*.txt`) with these sections:

- `[library]`: name, engine, `extends`, description.
- `[soap]`, `[lammps]`: the descriptor (and the similarity threshold) and the LAMMPS styles.
- `[bond_lengths]`: the bond length of each pair of elements (`E1 E2 r`); each must lie between the longest bond and the shortest distance not bonded of the pair in the references (checked when the table is built), and a pair not written gets the halfway value.
- `[types]`, `[bonds]`, `[angles]`, `[dihedrals]`, `[impropers]`: the parameters (molecular engine).
- `[references]`: the reference structures. These are extxyz files with the per-atom columns `fftype`, `ffcharge` and (molecular engine) `bonds` (the 1-based ids of the bonded atoms), or LAMMPS data files with the type names in the Masses section and a Bonds section.

The docstring of `soapff.library` gives the full format. Libraries are looked up by name in `$SOAPFF_LIBRARY_PATH`, then `~/.soapff/libraries`, then the bundled folder (`soapff list --paths`).

**Your own library.** First, `soapff.new_library("mine", extends="oplsaa")` writes `~/.soapff/libraries/mine.txt`. Then `soapff.add_reference("mine", molecule, types=[...], charges=[...], bonds=[...])` adds a reference, and any parameter rows that are needed go into the text file. The table of reference environments is computed once and cached in `~/.cache/soapff`. It is computed again when the references or the descriptor settings change. See `examples/oplsaa/user_library.py` (the junction of a block copolymer, extending `oplsaa`).

**Engines.**
Both engines type in the same two steps:
1. **Bonds from the bond lengths.** Two atoms are bonded when they are nearer than the bond length of their pair of elements. The bond lengths are written in the library file (`[bond_lengths]`). In the bundled libraries they are halfway between the longest bond and the shortest distance of two atoms that are not bonded, both measured in the references (clayff: Si–O 2.41, Al–O 2.66, Mg–O 2.70 Å; oplsaa: C–C 1.99, C–H 1.51, C–O 1.86 Å), and H–O is 1.25 Å, set by hand: in MD at 300 K the H…O of a hydrogen bond comes to 1.43 Å, nearer than in any reference, and an O–H bond stretches to 1.12 Å. When the table is built, each written length is checked to lie between those two distances in the references, so a reference added that a length does not fit stops the build. A pair bonded in the references but not written gets the halfway value, which the log gives as a row to copy into the file. Elements that have no bond in the references (the ions) get no bonds. If a pair's bonds and non-bonded distances overlap in the references, the table is not built.
   The water molecules (an O bonded to two H only, each H to this O only) are then taken apart from the rest: their bonds to other atoms are taken out, so the mineral or the polymer is typed without the water and each water molecule on its own. A water O can come nearer than the Al–O or Mg–O bond length to an edge Al or Mg (1.92 Å in the minimised nano-gouge, `examples/nanogouge`), and ClayFF bonds it to nothing. This only decides which bonds are taken out; the water's type still comes from step 2. It is used when the only O with two H in the references is that of water (all bundled libraries); a library with another such O (an edge M–OH₂) does not take water apart. `report.json` lists the bonds taken out (`water_apart`).
2. **Type and charge from the SOAP of the covalent neighbourhood.** The SOAP of each atom is computed on the atoms within 4 bonds only, so other chains, the next layer, the water and the ions are left out. The nearest reference row of the same element gives the type and the charge.

The build stops when an atom is below the similarity threshold, when an atom's bonds are not those of its type in the references, when a molecule's charge is not a whole number (molecular engine), or when an element is not in the library.
- **clayff:** the ClayFF writer (`soapff_cui`). The mineral, the water and the ions are typed in one build (the water taken apart in step 1); the O–H bonds and the M–O–H angles of ClayFF come from the bonds of step 1.
- **molecular:** any force field with bonds, angles, dihedrals and impropers, given by type or by class. The charge comes with the type from the nearest row, so it is right only when the charge of a type is decided within 4 bonds (as in OPLS-AA, the water models and ClayFF); a force field that gives each residue its own charges (protein force fields) is outside what soapff types.

## Test cases

The scripts and their results are in `examples/`.

| Test | Script | Result |
|---|---|---|
| Clay minerals (`clayff`, `clayff_edge`) | `examples/<folder>/` (`examples/README.md`) | No test structure is a reference. Water on the (001) slabs of 10 minerals and 6 water shapes on a slab (packmol; 692–4,249 atoms), the last frames of 19 MD runs with water at 300 K, a NaCl droplet on montmorillonite (105,977 atoms), the particle with edges, alone, twice and with a water film, and a nano-gouge of 3 such particles, dry and in water (66,900 atoms): 0 wrong types and 0 wrong charges, except 12 edge Mg of the minimised nano-gouge typed `mgb` for `mga` (same charge, LJ and angle constants). The written force fields minimise in LAMMPS. |
| Polymers (`oplsaa`) | `examples/oplsaa/` | 40 structures that are not references (long chains of 11 polymers, MD frames at 413 K and 600 K, random copolymers, 5 amorphous cells; 18,230 atoms): 0 wrong types and 0 wrong charges against Foyer. Forces of the written files in LAMMPS against Foyer's OPLS-AA in OpenMM: within 1.7e−7 kcal/mol/Å in vacuum; RMS 2.2e−6 kcal/mol/Å at most in the periodic cells (Ewald / PME 1e−8). Melt densities of PE, PP and PIB at 413 K in LAMMPS: within 1.5% of OpenMM after 500 ps, and within 0.7% of OpenMM's 2 ns run. The 153 frames of these LAMMPS melts (424,422 atoms), typed again from their coordinates and box: 0 wrong. A user library with one more reference (a short PE-b-PEO typed by Foyer) types a PE-b-PEO block copolymer, whose junction is in no reference of `oplsaa`: 0 wrong. |
| Water | `examples/md_check.py` | NPT with the LAMMPS Python module. TIP3P: 299.5 K, 0.965(3) g/cm³. SPC/E: 299.3 K, 0.985(5) g/cm³ (10 Å cutoff, no tail correction). |
| Swelling curve | `examples/montmorillonite_swelling/` | d001 of Na montmorillonite against its water content, NPT with the `clayff` force field, compared with Cygan et al. (2004). |

Besides soapff, the clay examples need packmol (water) and pymatgen (surfaces). The polymer examples need Foyer, OpenMM, ParmEd, RDKit, Open Babel and packmol to remake the library, the test structures and the OpenMM side. The water libraries are remade from OpenMM's files by `examples/water/make_library_txt.py`.

**How the MD of the tests was run.** The examples run LAMMPS through its Python module (`pip install "lammps[mpi]"`, LAMMPS 22 Jul 2025). Some results kept in `examples/` were run before that change, with the same LAMMPS commands written to an input file and the `lmp` executable (the clay MD on another computer over SSH; see each folder's README). Those runs used Ewald where the Homebrew LAMMPS crashed when PPPM was set up a second time.

## Layout

```
src/soapff/              the package: api.py (build, System), library.py (library files), molecular.py, clayff.py,
                         lammps.py (system.data and system.settings of the molecular engine)
src/soapff/libraries/    the bundled libraries (*.txt), their reference structures and the prebuilt ClayFF tables
src/soapff_cui/          the ClayFF writer used by the clayff engine (SOAP table, ClayFF types and charges, water, surfaces)
examples/                quick start, your own library, short MD checks, and all the test cases (examples/README.md); the MD
                         of the clay examples (clay_md.py)
```

## License

MIT (`LICENSE`).
