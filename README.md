# soapff

soapff is a Python package that turns a structure file into a LAMMPS simulation. It assigns a force-field type and a charge to every atom. It does this by comparing each atom's local environment (SOAP descriptor) with reference structures whose types and charges are known. It then writes the LAMMPS data, settings and input files.

```python
import soapff
s = soapff.build("amber14sb", "protein_in_water.pdb")   # library, structure file
s.run(np=4)                                              # LAMMPS
```

- **Structure files:** any format ASE reads (cif, xyz, extxyz, pdb, POSCAR, ...), a LAMMPS data file, or an `ase.Atoms`.
- **Libraries:** text files. Each one names the reference structures and gives the force-field parameters. Users can extend a bundled library or write their own.
- **Bundled libraries:** ClayFF (clay minerals), AMBER ff14SB (proteins), TIP3P and SPC/E (water and ions).
- **Checks:** an atom whose environment is not in the library stops the build with an error that names the atom. soapff never silently guesses a type.

## Installation

```bash
git clone https://github.com/haru2225/soapff.git
pip install -e soapff          # Python >= 3.9; installs numpy, scipy, ase, dscribe
```

Running the simulations needs LAMMPS. soapff looks for `lmp`, `lmp_mpi` or `lmp_serial` on the PATH; set `$SOAPFF_LMP` (and `$SOAPFF_MPIRUN`) to choose another executable.

The ClayFF engine adds water with [packmol](https://m3g.github.io/packmol/) when asked to.

Optional extras:
- `pip install -e "soapff[analysis]"`: analysis of finished runs (matplotlib, sovapy).
- `pip install -e "soapff[surfaces]"`: (001) slabs of the minerals (pymatgen).

## Usage

```python
import soapff

soapff.libraries()                                     # the libraries that can be given by name
s = soapff.build("clayff", "kaolinite.cif")            # writes kaolinite_soapff/ in the current folder
s = soapff.build("amber14sb", "protein_in_water.pdb", "out", nvt_ps=100, npt_ps=100)
s.files                  # system.data, system.settings, in.lammps, charges.csv, report.json in s.out_dir
s.types, s.charges       # the type and charge of every atom
s.run(np=4)              # LAMMPS in s.out_dir
s.thermo()["PotEng"]     # the thermo output of the last run
lmp = s.lammps()         # or the LAMMPS Python module, with in.lammps run
```

Options of `soapff.build`:

- `md=False`: `in.lammps` stops after the energy minimisation.
- MD stages: `temp`, `relax_ps`, `nvt_ps`, `npt_ps`, `press`, `baro`, `dump_ps`, `timestep`.
- `cutoff`, `kspace` (e.g. `"ewald 1.0e-5"`), `boundary`.
- `strict=False`: write the input even when some atoms are not in the library.
- `padding`: a structure without a cell, such as an xyz file or a PDB file without CRYST1, is put in a box this much larger (Å).

The command `soapff` does the same:

```bash
soapff list
soapff info amber14sb
soapff build amber14sb protein.pdb -o out --run --np 4
soapff run out --np 4
soapff new mine --extends amber14sb
soapff add-reference mine ligand.extxyz
```

## Libraries

| Library | Engine | References | Force field |
|---|---|---|---|
| `clayff` | clayff | Bulk crystals of 11 layered minerals and their (001) slabs | ClayFF (Cygan et al. 2004), flexible SPC water |
| `clayff_edge` | clayff | Montmorillonite particle with edges | Underwood & Bourg (2020) |
| `amber14sb` | molecular | Lysozymes (1AKI, 2LZM) with protonation variants, N- and C-terminal fragments of every amino acid, ACE/NME caps; plus `tip3p` | AMBER ff14SB (from OpenMM `amber14/protein.ff14SB.xml`) |
| `tip3p` | molecular | Water, NaCl and 14 ions | TIP3P (rigid), Joung–Cheatham 2008 monovalent and Li–Merz 2013 divalent ions |
| `spce` | molecular | The same, for SPC/E | SPC/E (rigid) and its ions |

**Library files.** Each library is a text file (`src/soapff/libraries/*.txt`) with these sections:

- `[library]`: name, engine, `extends`, description.
- `[soap]`, `[graph]`, `[bonding]`, `[lammps]`: the descriptor, the bond detection and the LAMMPS styles.
- `[types]`, `[bonds]`, `[angles]`, `[dihedrals]`, `[impropers]`: the parameters (molecular engine).
- `[references]`: the reference structures. These are extxyz files with the per-atom columns `fftype` and `ffcharge`, or LAMMPS data files with the type names in the Masses section.

The docstring of `soapff.library` gives the full format. Libraries are looked up by name in `$SOAPFF_LIBRARY_PATH`, then `~/.soapff/libraries`, then the bundled folder (`soapff list --paths`).

**Your own library.** First, `soapff.new_library("mine", extends="amber14sb")` writes `~/.soapff/libraries/mine.txt`. Then `soapff.add_reference("mine", molecule, types=[...], charges=[...])` adds a reference, and any parameter rows that are needed go into the text file. The table of reference environments is computed once and cached in `~/.cache/soapff`. It is computed again when the references or the descriptor settings change. See `examples/user_library.py`.

**Engines.**
- **clayff:** the ClayFF writer (`soapff_cui`). The mineral is typed first; the water and the ions are typed in a second step with a solution library.
- **molecular:** any force field with bonds, angles, dihedrals and impropers. It works in three steps:
  1. Bonds come from the distances (covalent radii × 1.15).
  2. The type comes from the SOAP vector of each atom's covalent neighbourhood.
  3. For libraries whose charges depend on more than the neighbourhood (`[graph]`, such as the ff14SB residues), the charge comes from an exact match of the bonded fragment.

## Test cases

The scripts and their results are in `examples/`.

| Test | Script | Result |
|---|---|---|
| Clay minerals (`clayff`) | `examples/clayff/test_clayff.py` | The 11 bulk crystals and 11 (001) slabs as 2×2×1 supercells, also with random displacements (652–3,264 atoms), plus kaolinite read from a CIF file: 0 wrong types and 0 wrong charges in 45 structures. Kaolinite with water (1,659 atoms): every atom typed. |
| Proteins (`amber14sb`) | `examples/amber14sb/` | 12 structures that are not references (ubiquitin, BPTI, crambin, Trp-cage, fragments, MD frames): 0 wrong types and 0 wrong charges against OpenMM. Vacuum energies against OpenMM: bond, angle and nonbonded terms within 1e−4 kcal/mol, torsions within 0.1 kcal/mol. |
| Protein in water | `examples/amber14sb/solvated_lammps.py` | Ubiquitin in TIP3P with 0.15 M NaCl (17,018 atoms) from a PDB file: 0 wrong types and 0 wrong charges. Total energy differs from OpenMM PME by 0.21 kcal/mol (of −36,286). The written `in.lammps` runs. |
| Water | `examples/md_check.py` | TIP3P: 300.5 K, 0.977 g/cm³. SPC/E: 298.9 K, 0.988 g/cm³ (NPT, 10 Å cutoff, no tail correction). |
| Your own library | `examples/user_library.py` | Methanol added to a library that extends `amber14sb`. A methanol/water mixture read from a CIF file is typed correctly; the parent library refuses it. |
| More clay cases | `examples/<folder>/` | Larger crystals, MD frames, basal slabs, interfaces with water, water shapes, a 105,977-atom droplet, the particle with edges, the swelling curve of Na montmorillonite, diffusion, mechanical properties, and a supercomputer workflow. See `examples/README.md`. |

Besides soapff, the clay examples need packmol (water) and pymatgen (surfaces). The protein examples need OpenMM and pdbfixer to remake the library and the reference energies.

**Note on PPPM.** The tests were run with Ewald, because the Homebrew LAMMPS used for them crashes when PPPM is set up a second time (FFTW3 with threads). Other builds run the `in.lammps` soapff writes as it is.

## Layout

```
src/soapff/              the package: api.py (build, System), library.py (library files), molecular.py, clayff.py, lammps.py,
                         remote.py (LAMMPS on another computer over SSH)
src/soapff/libraries/    the bundled libraries (*.txt), their reference structures and the prebuilt ClayFF tables
src/soapff_cui/          the ClayFF writer used by the clayff engine (SOAP table, ClayFF types and charges, water, surfaces, analysis)
examples/                quick start, your own library, short MD checks, and all the test cases (examples/README.md)
```
