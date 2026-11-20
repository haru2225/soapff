# Polymers with OPLS-AA (`oplsaa`)

The library `oplsaa` assigns OPLS-AA types and charges (Jorgensen, Maxwell, Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)) from the SOAP of each atom alone. It uses no SMARTS and no bond orders. Its parameters come from Foyer's `oplsaa.xml` (Foyer 1.3.0). Its 88 reference structures are oligomers of 11 polymers, typed by Foyer.

This folder makes the library and tests it in four ways:
- the types and charges against Foyer;
- the forces against Foyer's force field in OpenMM;
- the melt MD against Foyer's force field in OpenMM;
- the melt densities against the measurements of Fetters, Lohse, Richter, Witten and Zirkel (Macromolecules 27, 4639 (1994)).

The polymers are listed in `polymers.POLYMERS`:
- PE;
- atactic PP;
- PIB;
- cis-1,4-PB and cis-1,4-PI;
- PEO, with CH3O ends;
- atactic PMMA, PMA, PVA and PAN;
- PTFE.

## Files

| File | What it does |
|---|---|
| `polymers.py` | Code shared by the scripts below. It covers the SMILES of the chains, the 3D structures (Open Babel), and Foyer's types and charges. It also builds the OpenMM System of Foyer's force field and does the minimisation, the MD in vacuum and the amorphous cells (packmol). |
| `make_library_txt.py` | Writes the parameters of `src/soapff/libraries/oplsaa.txt` from Foyer's `oplsaa.xml`. The Ryckaert–Bellemans torsions become `dihedral_style fourier`, and the bonded terms are given by class. |
| `make_refs.py` | Makes the 88 references (`src/soapff/libraries/oplsaa/`) and the test structures (`test/`), all typed by Foyer. It also writes the two structures of `user_library.py` (`extra/`). |
| `make_cells.py` | Makes the 5 amorphous cells of `test/`. Each cell is made in four steps: the chains are built, made compact by MD at 600 K in a sphere, packed by packmol at 0.25 g/cm³, and minimised. |
| `test_typing.py` | Runs soapff on every test structure, given only the elements and coordinates, and compares the result with Foyer. Writes `results/typing.json`. |
| `depth.py` | Types the test structures with `oplsaa` (SOAP of the atoms within 4 bonds) and with a library that keeps the atoms within 3 bonds only, and counts the types exchanged. Writes `results/depth.json`. |
| `test_frames.py` | Types every frame of the LAMMPS melt MD (`md_lammps.py`) from its elements, coordinates and box, and compares with Foyer's types and bonds of the start. Writes `results/frames.json`. |
| `user_library.py` | Builds a user library that extends `oplsaa` with one reference, a short PE-b-PEO (`extra/PE4_PEO3.extxyz`), and types a longer PE-b-PEO (`extra/block_PE_PEO.extxyz`) with it. Writes `results/user_library.json`. |
| `smiles_check.py` | The check behind `../polymer_smiles.py`: the PMMA 10-mer that RDKit makes from its SMILES, typed by Foyer (`extra/PMMA_10mer_rdkit.extxyz`, written by `smiles_check.py foyer`), against soapff given its elements and positions only. Writes `results/smiles_check.json`: all 152 atoms get Foyer's type and charge. |
| `forces_openmm.py`, `forces_lammps.py` | Compute the force on every atom twice: from Foyer's OPLS-AA in OpenMM, and from soapff's force field in LAMMPS (the LAMMPS Python module, `run 0`). Write `results/forces.json`. |
| `md_prepare.py`, `md_lammps.py`, `md_openmm.py`, `md_density.py` | Run NPT MD of the PE, PP and PIB melts at 413 K and 1 bar, from the same start, both with soapff's force field in LAMMPS and with Foyer's force field in OpenMM. `md_lammps.py` shows how an NPT is written with the LAMMPS Python module after `soapff.build`. Write `md/` and `results/md_density.json`. |

## Requirements

- **Tests:** soapff (this repository) and the LAMMPS Python module (`pip install "lammps[mpi]"`).
- **Library, references, test structures and OpenMM side:** Foyer 1.3.0, OpenMM, ParmEd, RDKit, Open Babel (`obabel`) and packmol. The versions used were OpenMM 8.6.1, ParmEd 4.3.1, RDKit 2026.09.2, Open Babel 3.2.1 and LAMMPS 22 Jul 2025 (the Python module from PyPI, `lammps` 2025.7.22).
- **Open Babel's `--gen3d` is not deterministic.** Running `make_refs.py` again gives other structures, which are just as valid, so the numbers below change a little. The structures that were used are in the repository.

## Commands

```bash
cd examples/oplsaa
# the library (Foyer environment)
python make_library_txt.py ../../src/soapff/libraries/oplsaa.txt
python make_refs.py ../../src/soapff/libraries test
python make_cells.py test
python smiles_check.py foyer                   # the RDKit PMMA of ../polymer_smiles.py, typed by Foyer
# typing (soapff)
python test_typing.py test WORK results/typing.json
python user_library.py extra/block_PE_PEO.extxyz extra/PE4_PEO3.extxyz WORK results/user_library.json
python smiles_check.py WORK
python depth.py test WORK results/depth.json
# forces
python forces_openmm.py test WORK/openmm.npz                          # Foyer environment
python forces_lammps.py test WORK/openmm.npz WORK results/forces.json
# melt MD
python md_prepare.py test md PE PP PIB                                # Foyer environment (a GPU if there is one)
mpirun -np 4 python md_lammps.py md 500 PE PP PIB                     # LAMMPS Python module on 4 MPI processes
python md_openmm.py md 500 PE PP PIB                                  # Foyer environment
python md_density.py md 100 results/md_density.json PE PP PIB
python test_frames.py md results/frames.json PE PP PIB                 # the frames of the LAMMPS melt MD
```

## References of the library

There are 8 structures for each polymer (88 structures, 6,264 atoms; 4,612 rows of the SOAP table after pruning):
- a 5-mer and an 8-mer as built: an extended chain from Open Babel, minimised with Foyer's force field;
- a 5-mer of another tacticity after 10 ps of MD at 413 K in vacuum;
- 5 frames (10, 20, 30, 40 and 50 ps) of MD at 413 K in vacuum of another 8-mer, which folds up.

The MD frames of the folding 8-mers put the conformations of chains in a melt among the references (in the PEO frames, 75–100% of the O–C–C–O dihedrals are gauche).

Each reference carries Foyer's bonds as well as its types and charges (the per-atom column `bonds` of the extxyz). From them the library learns two things:
- **The bond length of each pair of elements.** By default it is halfway between the longest bond and the shortest distance of two atoms not bonded (within 4 Å) in the references. `[bond_lengths]` in `oplsaa.txt` sets H–O by hand (see below).
- **The elements and numbers of the bonded atoms of each type**, used to check the bonds (opls_185, the H of a C bonded to an ether O, has one C).

| Pair | Longest bond (Å) | Shortest not bonded (Å) | Bond length (Å) |
|---|---|---|---|
| C–C | 1.716 | 2.260 | 1.988 |
| C–H | 1.217 | 1.797 | 1.507 |
| C–O | 1.545 | 2.179 | 1.862 |
| C–N | 1.210 | 2.562 | 1.886 |
| C–F | 1.425 | 2.110 | 1.767 |
| H–O | 1.021 | 1.761 | 1.25, by hand (halfway: 1.39) |

The longest bonds are those of the 413 K MD frames of the references. H–O is set by hand because of hydrogen bonds. In MD at 300 K, an H comes to 1.43–1.50 Å from an O it is not bonded to: 1.50 Å for PVA in TIP3P water, 1.43 Å for clay with flexible SPC water. Meanwhile O–H bonds stretch to 1.12 Å. 1.25 Å lies between the two. The tip3p and spce libraries use the same value.

Types and charges then come in two steps:
1. The bonds: two atoms closer than the bond length of their pair of elements are bonded. A pair never bonded in the references (H–H, O–O) is never bonded.
2. The type and charge: from the SOAP of the atoms within 4 bonds, the nearest row of the same element.

## Results

### Types and charges (`results/typing.json`)

None of these structures is a reference.

| Structures | Atoms | Wrong types | Wrong charges | Lowest similarity |
|---|---|---|---|---|
| 11 long chains: PE C100, the others 20 units, with another tacticity than the references | 122–302 | 0 | 0 | 0.984–1.000 |
| the same after 10 ps at 413 K in vacuum | 122–302 | 0 | 0 | 0.974–0.995 |
| the same after 10 ps at 600 K in vacuum | 122–302 | 0 | 0 | 0.983–0.992 |
| random copolymers: ethylene–propylene and MMA–MA | 218, 269 | 0 | 0 | 0.987, 0.978 |
| 5 amorphous cells: PE, PP, PIB, PEO and PMMA, periodic | 1,192–2,892 | 0 | 0 | 0.980–0.995 |

In total, the 40 structures hold 18,230 atoms, and every one gets Foyer's type and charge. The bonds made from the bond lengths are Foyer's bonds in all 40 structures.

Foyer's `oplsaa.xml` has no parameters for the PB and PI backbone dihedrals CM–CT–CT–CM and HC–CM–CT–CT. Like Foyer, soapff leaves them out and reports them in `report.json`.

### The checks

soapff stops (`TypingError`) in four cases:
- an atom is below the similarity threshold;
- an element is not in the library;
- an atom's bonds are not those of its type (the elements and numbers that type has in the references);
- a molecule's charge is not a whole number (more than `[soap] charge_tolerance`, 0.001 e, from one).

None of the tests here stopped. The charge and bond checks do not see a confusion of two types with the same charge and the same bonded atoms. opls_185 and opls_469 are both +0.03 e and bonded to one C, with ε 0.03 and 0.015 kcal/mol. None happened with the settings of `oplsaa`; with SOAP depth 3, 87 H were typed so (see Depth below).

### User library (`results/user_library.json`)

The junction of a PE-b-PEO block copolymer is in no reference of `oplsaa`. `user_library.py` makes a library that extends `oplsaa` with one reference, a short PE-b-PEO typed by Foyer (`extra/PE4_PEO3.extxyz`: C8 and 3 EO units, 51 atoms), and types a longer one with it (`extra/block_PE_PEO.extxyz`: C20 and 10 EO units, 136 atoms): every atom gets Foyer's type and charge, with similarity 0.999 at least. No parameter line was added: `oplsaa` has the types and terms it needs.

### MD frames of the melts (`results/frames.json`)

`test_frames.py` types every frame of the 500 ps NPT MD of the PE, PP and PIB melts at 413 K in LAMMPS (`md_lammps.py`; a frame every 10 ps, 51 per polymer) from its elements, positions and box alone, and compares with Foyer's types and charges of the start.

| Melt | Frames | Atoms typed | Wrong | Frames that stop | Lowest similarity | Lowest margin |
|---|---|---|---|---|---|---|
| PE (C100 × 9) | 51 | 138,618 | 0 | 0 | 0.984 | 0.0038 |
| PP (50 units × 6) | 51 | 138,312 | 0 | 0 | 0.978 | 0.0020 |
| PIB (40 units × 6) | 51 | 147,492 | 0 | 0 | 0.987 | 0.0091 |

The bonds come from the bond lengths (C–C 1.99 Å, C–H 1.51 Å). In every frame they are Foyer's bonds of the start (`frames_with_bonds_differing`: 0). Some C–H bonds stretch to 1.26 Å at 413 K, still well below 1.51 Å.

### Depth (`results/depth.json`)

The SOAP of step 2 is computed on the atoms within 4 bonds. With 3 bonds (the other settings of `oplsaa`), 87 H of the PEO, PMA and PMMA structures and of the MMA–MA copolymer are mistyped: the H of the ester O–CH3 (opls_469) and the H on a C bonded to an ether O (opls_185) are exchanged. Both have +0.03 e, so the charges stay right and soapff does not stop. What tells them apart is the carbonyl O, four bonds from the H. With 4 bonds every atom of the 40 structures is right (`depth.py`; the 87 H are in 11 of the 40 structures, and none stops).
