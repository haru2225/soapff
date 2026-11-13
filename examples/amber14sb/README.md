# The amber14sb library of soapff: how it was made and checked

The libraries `amber14sb`, `tip3p` and `spce` bundled with soapff (`src/soapff/libraries/`) were written from the OpenMM
force-field files (`amber14/protein.ff14SB.xml`, `amber14/tip3p.xml`, `amber14/spce.xml`). Their reference structures are
proteins typed by OpenMM. This folder holds the scripts that made them, a test set that is not in the references, and the
checks against OpenMM.

## Scripts

| Script | Python | What it does |
|---|---|---|
| `omm_common.py` | OpenMM, pdbfixer | PDB → pdbfixer → hydrogens at pH 7 → the ff14SB/TIP3P templates OpenMM matches → the AMBER atom class and charge of each atom. PDB files are downloaded to `pdb/` the first time. |
| `make_refs.py OUT` | OpenMM, pdbfixer | The references (`OUT/refs`) and the test set (`OUT/test`) as extxyz files with the per-atom columns `fftype` and `ffcharge`. |
| `make_library_txt.py LIBDIR OUT/refs` | OpenMM | `amber14sb.txt`, `tip3p.txt`, `spce.txt` and the water/ion references, converted from the OpenMM XML files. |
| `test_typing.py test WORK results/typing.json` | soapff | soapff on the test set, compared with OpenMM's types and charges. |
| `energy_openmm.py test results/energies_openmm.json` | OpenMM | OpenMM energies of the test set: vacuum, no cutoff, each force separately. |
| `energy_lammps.py test results/energies_openmm.json WORK results/energies_lammps.json` | soapff, LAMMPS | The same energies from the files soapff writes, using `lj/cut/coul/cut` with a 999 Å cutoff and boundary `f f f`. |
| `solvated.py` | OpenMM, pdbfixer | Ubiquitin in TIP3P with 0.15 M NaCl. It writes `ubq_water.pdb` and `ubq_water.extxyz`, and the OpenMM PME energy of these files to `results/energies_solvated.json`. |
| `solvated_lammps.py . WORK results/solvated.json` | soapff, LAMMPS | soapff on `ubq_water.pdb`: the types, the energy against OpenMM, and the `in.lammps` soapff writes (minimisation, then NVT with SHAKE water). |

## The references and the test set

- **References** (64 proteins and fragments, plus the water and ion references of `tip3p`):
  - **Whole proteins:** hen egg-white lysozyme (1AKI) and T4 lysozyme (2LZM), with every histidine as HID, HIE and HIP. T4 lysozyme also appears with ASH, GLH and LYN.
  - **Terminal fragments:** for each of the 20 amino acids (and HID, HIE, HIP), one fragment where it is the N terminus and one where it is the C terminus.
  - **Capped fragments:** fragments capped with ACE and NME.
  - Fragment provenance is listed in `results/references_manifest.json`.
- **Test set** (`test/`, not in the references):
  - ubiquitin (1UBQ), BPTI (5PTI), crambin (1CRN) and Trp-cage (1L2Y, NMR);
  - fragments of ubiquitin, free and ACE/NME-capped;
  - three frames (5, 10 and 15 ps) of a 300 K vacuum MD of ubiquitin with flexible bonds.

## How the molecular engine types a protein

1. **Bonds** come from the distances (ASE covalent radii × 1.15). Ions are never bonded.
2. **Type:** the SOAP vector of the covalent neighbourhood (the atoms within 4 bonds) is compared with the reference rows, and the nearest row gives the type.
3. **Charge:** ff14SB charges depend on the residue, so the molecule is cut at the peptide C–N bonds (`[graph] cut = C N`). The charge comes from a fragment signature: the types of the atoms within 12 bonds, by bond distance. This signature must match a reference exactly (similarity ≥ 0.9999). Otherwise soapff stops (`soapff.TypingError`) and names the atoms.

## Results

**Typing** (`results/typing.json`): **0 wrong types and 0 wrong charges** on all 12 test structures (92–1231 atoms). This includes the NMR structure and the MD frames. The lowest SOAP similarity is 0.966; every fragment similarity is 1.0000. Typing takes 0.1–1.4 s per structure. Computing the library table takes about 10 s, once (cached in `~/.cache/soapff`).

**Energies in vacuum** (`results/energies_lammps.json`, LAMMPS − OpenMM, kcal/mol):

| Structure | Atoms | Bond | Angle | Torsion | Nonbonded | Total |
|---|---|---|---|---|---|---|
| 1UBQ | 1231 | −1e−7 | +4e−7 | −0.0019 | +2e−5 | −0.0019 |
| 1UBQ, MD 15 ps | 1231 | +4e−7 | +4e−7 | −0.0974 | +3e−5 | −0.0973 |
| 5PTI | 892 | −1e−7 | −3e−7 | −0.0042 | +2e−5 | −0.0042 |
| 1CRN | 642 | −1e−7 | +5e−7 | −0.0420 | +1e−5 | −0.0420 |
| 1L2Y | 304 | +8e−8 | +2e−7 | −0.0006 | +5e−6 | −0.0006 |
| fragments (4) | 92–119 | < 5e−7 | < 5e−7 | 0.0000 | < 4e−6 | 0.0000 |

- The proper torsions are identical.
- The torsion difference of the whole proteins comes from 3 of the 216 impropers of ubiquitin. Both programs order the outer atoms of an improper by type (OpenMM's `amber` ordering). When two outer atoms have the same type, OpenMM breaks the tie by residue and then by the atom's place in the residue template; soapff has no templates and uses the atom index. In these three impropers the two orders differ, so the dihedral angle, and with it the energy, differs.

**Ubiquitin in water** (`results/solvated.json`; 17,018 atoms, 5,253 TIP3P waters, 14 Na⁺ and 14 Cl⁻):

- **Typing:** the PDB file written by OpenMM gives **0 wrong types and 0 wrong charges**, in 2.7 s.
- **Energy:** LAMMPS (Ewald 1e−6) against OpenMM (PME 1e−6), no dispersion correction, flexible water:

  | Term | LAMMPS − OpenMM (kcal/mol) |
  |---|---|
  | Bond | +0.000 |
  | Angle | +0.000 |
  | Torsion | −0.002 |
  | Nonbonded | −0.21 of −37,631 |
  | Total | −0.21 of −36,286 (6e−6) |

  With LAMMPS's `g_ewald` set to OpenMM's α, the total difference is −0.10.
- **MD:** the `in.lammps` soapff writes (Ewald 1e−5, serial LAMMPS): the minimisation stops on the energy tolerance after 274 iterations; then 1 ps NVT with SHAKE on the 5,253 waters (`b 56 a 119`), which starts at 300 K. The temperature is 279 K after 1 ps, still rising back from the minimised structure, and `final.data` is written. This took 18 min.
- **Use the box of the file.** The Modeller box is 56.7416 Å; CRYST1 rounds it to 56.742 Å. The water box of `addSolvent` is not relaxed at the boundary, so this rounding alone changes the nonbonded energy by about 20 kcal/mol. The OpenMM energy is therefore computed from the files.

## Limits

- **Residues:**
  - Covered: the 20 amino acids, HID/HIE/HIP, ASH/GLH/LYN, CYS/CYX, N- and C-termini, and ACE/NME caps.
  - Not covered: CYM, HYP, NHE, and ligands, nucleic acids, sugars and lipids. To add them, put references with their types and charges into your own library (`soapff.new_library("mine", extends="amber14sb")`, then `soapff.add_reference`), plus the parameter rows the new types need.
- **Equivalent atoms with different charges:** in a few ff14SB terminal residues, atoms that are equivalent by symmetry carry different charges (NLEU CD1/CD2, NPHE CD1/CD2 and CE1/CE2). Their environments cannot tell them apart, so soapff gives both the mean. The molecule's total charge is unchanged.
- **Torsions:** dihedrals are written with `dihedral_style charmm` (weight 0: 1-4 pairs come from `special_bonds amber`). Terms with K = 0 are left out, as OpenMM does.
- **PPPM crash in the Homebrew LAMMPS:** the Homebrew LAMMPS used here (22 Jul 2025, FFTW3 with threads) crashes in `FFT3d::~FFT3d` → `fftw_plan_awake` whenever PPPM is set up a second time, for example in the run after a minimisation. The checks above therefore use Ewald (`soapff.build(..., kspace="ewald 1.0e-5")`). This is a problem of that build, not of the input.
