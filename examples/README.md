# Examples and test cases of soapff

Run these with a Python where `pip install -e ".[examples]"` was run: soapff and the LAMMPS Python module (`pip install "lammps[mpi]"`). soapff writes the force field only; every example writes its own minimisation and MD in LAMMPS, through the LAMMPS Python module (`from lammps import lammps`); the minimisations of the clay examples also run with the `lmp` executable given by `LMP=/path/to/lmp` (`clay_md.run`). The scripts of the first table are run from the root of the repository; the clay examples are run inside their own folder (`python run.py`). No test structure is a reference of the library it is typed with, except the particle with edges and the particle twice of `montmorillonite_edge/` (the references of `clayff_edge` themselves).

## The package and its libraries

| File | What it shows |
|---|---|
| `quickstart.py [--run]` | Lists the libraries. Builds the force field of ClayFF kaolinite/water from an extxyz file, of a PMMA chain with oplsaa from a plain xyz file (no types, no cell, so it is boxed), and of a PE melt with oplsaa. Shows `charges.csv` and the force field as LAMMPS commands. `--run` then runs 10 ps of NPT of the PE melt at 413 K with the LAMMPS Python module. |
| `polymer_smiles.py [ps]` | A polymer from its SMILES: a PMMA 10-mer made 3D with H by RDKit (only the elements and positions go to soapff), its OPLS-AA force field (`oplsaa`; every type and charge is the one Foyer gives the same coordinates: `oplsaa/smiles_check.py`), then minimisation and NVT at 413 K (default 100 ps; 24 s on one core) with the LAMMPS Python module. Needs `pip install rdkit`. |
| `md_check.py WORK` | Short MD with soapff's force fields. TIP3P and SPC/E boxes of 216 waters (SHAKE from the force field; minimisation, 2 ps NVT relaxation, 5 ps NVT, 20 ps NPT). Kaolinite with water (ClayFF, `clay_md.py`: minimisation, 0.5 + 2 ps NVT). Results in `results/md.json`. |
| `clay_md.py` | The MD protocol of the clay examples as LAMMPS commands (two-stage minimisation, NVT relaxation, NVT with trajectory and MSD, NPT), for `lmp.commands_list`, or written to an input file for the `lmp` executable. |
| `clayff/water_md_frames.py` | `clayff` on the last frames of 19 MD runs with water at 300 K (`library_water/refs/*/final.data` and `montmorillonite_2W/diffusion/start/final.data`), against the types of their MD input; the longest O-H bond and the shortest H…O of the hydrogen bonds, which decide the H-O bond length (1.25 Å). Results in `results/clayff_water_md_frames.json`. |
| `clayff/water_apart_check.py` | Step 1 on the 37 clay structures with water: the bonds between water and the rest that taking the water molecules apart takes out (none; the nano-gouge, `nanogouge/`, is the case where there are). Results in `results/clayff_water_apart_check.json`. |
| `clayff/obss_check.py WORK` | A type with no reference: one Si of the dry Na montmorillonite reference changed to Al, so that an O touches both a tetrahedral and an octahedral substitution (ClayFF obss, in no `clayff` reference). The nearest rows are obts (similarity 0.987), but the bonded elements of the O match no obts reference atom, and the build stops with a `TypingError` that names the atoms. Results in `results/clayff_obss_check.json`. |
| `oplsaa/` | Polymers with OPLS-AA: how the library was made from Foyer's `oplsaa.xml` and oligomers typed by Foyer, and the checks against Foyer (types and charges of chains, MD frames, copolymers and amorphous cells), against Foyer's force field in OpenMM (forces, melt MD) and against the measured melt densities (Fetters et al. 1994); a user library for a block-copolymer junction. See `oplsaa/README.md`. |
| `compare_tools/` | atom 3.0.2 (Holmboe 2019; ClayFF by coordination rules, run in GNU Octave) on the 16 `clayff` test structures in the reference arrangement: soapff gives every atom atom's charge, and atom's type except some O next to a tetrahedral Al. See `compare_tools/README.md`. |
| `paper_figures/` | Figs. 3 and 7 of the paper: the SOAP vectors of the O references on their first two principal components (`pca.py`), and the typing of the test structures with the SOAP settings changed one at a time (`soap_params.py`, its summary in `figdata/soap_params/results.json`). |
| `water/make_library_txt.py LIBDIR` | How the `tip3p` and `spce` libraries were made from OpenMM's `amber14/tip3p.xml` and `amber14/spce.xml` (needs OpenMM). |

## Clay minerals (`clayff`, `clayff_edge`)

The references of the two libraries are in `src/soapff/libraries/clayff/` (`bulk/`, `basal/`, `water/`, `edge/`). Each folder here has a `README.md`, its inputs, the script that redoes the calculation, and its results (`results.json`).

| Test case | Folder | What it holds |
|---|---|---|
| Water and ion references, MD frames | `library_water/` | How the water and ion references of the libraries were made (packmol liquids, films and a droplet in vacuum), and the 300 K MD of the liquids, of the 11 (001) slabs with water and of the 2-layer hydrate, whose last frames are the test frames of `clayff/water_md_frames.py`. |
| Interfaces with water | `interfaces_001_water/` | Water on the (001) slabs of 10 minerals (packmol): force field and minimisation. |
| Water shapes | `interface_shapes/` | 6 initial structures made with packmol's constraints (sphere, hemisphere, cylinder, cavity, xygauss, NaCl film): force field and minimisation. |
| Large droplet | `meng2022_droplet/` | The NaCl droplet of Meng et al. 2022 on montmorillonite (001), 105,977 atoms. |
| Particle with edges | `montmorillonite_edge/` | The `clayff_edge` library on the particle, the particle twice and the particle with a water film and NaCl: types, charges, coefficients, single point against the original files, minimisations. |
| Nano-gouge | `nanogouge/` | 3 particles with edges in a periodic box, dry and in water (66,900 atoms), typed with `clayff_edge` against the types of the builder (`typing_test.py`). |
| 2-layer hydrate | `montmorillonite_2W/` | Na montmorillonite with two water layers (1,554 atoms): its force field and a 300 K NPT, the last frame of which is one of the test frames of `clayff/water_md_frames.py`. |
| Swelling curve | `montmorillonite_swelling/` | The Na montmorillonite of Cygan et al. 2004 and its swelling curve (10 water contents; `swelling.py`). |

Notes on these examples:
- **Force field:** the charges of Table 1 of Cygan et al. 2004 as they are, the LJ parameters unrounded, and the flexible SPC water. The edge examples use the particle's own force field.
- **Where the MD ran:** the MD of `montmorillonite_2W/`, `montmorillonite_swelling/` and the systems of `library_water/refs/` was run on another computer (Linux, Open MPI, PPPM 1e-5) with LAMMPS 7 Feb 2024 (Update 1); the water, NaCl and KCl systems of `library_water/refs/` with LAMMPS 29 Aug 2024 (the first line of each `log.lammps`). These runs used the `in.lammps` files written by an earlier soapff, whose commands are those the scripts now write with `clay_md.py`; a rerun with another LAMMPS version may relax to slightly different minima.
- **Minimisations:** those of the interface examples and the particle with edges were last run with the `lmp` executable of LAMMPS 22 Jul 2025 (`clay_md.run`): the force field, `kspace_style ewald 1.0e-4` (PPPM stopped at the start of a minimisation with the Homebrew LAMMPS used for the stored results), then the two-stage minimisation of `clay_md.py`. The same commands are written to `in.min`. The `in.lammps` files in some folders were written by the earlier soapff and are kept as records.
- **Not kept:** trajectories (`*.lammpstrj`) and packmol's xyz files of the large structures.

## Results of `md_check.py` (`results/md.json`; LAMMPS 22 Jul 2025 Python module, serial, PPPM 1e-5)

| System | Atoms | Run | Result |
|---|---|---|---|
| TIP3P, 216 waters | 648 | 27 ps | mean temperature of the NPT 299.5 K; density over the second half of the NPT 0.965(3) g/cm³ (10 Å cutoff, no tail correction; TIP3P is about 0.98 with it) |
| SPC/E, 216 waters | 648 | 27 ps | 299.3 K; 0.985(5) g/cm³ (SPC/E is about 0.99 with the tail correction) |
| kaolinite / water (ClayFF) | 1659 | minimisation, 0.5 ps relaxation, 2 ps NVT (0.5 fs) | mean temperature over the second half of the NVT 299.2 K; `final.data` written |

The numbers in brackets are standard errors from 5 blocks. SHAKE is part of the water force field, so it is on during the minimisation too (LAMMPS turns the constraints into stiff restraints there and warns about it). An earlier run of the same protocol with the `lmp` executable, Ewald 1e-5 and SHAKE switched on after the minimisation gave 0.977 and 0.988 g/cm³. A box of 216 waters is small: the density of single frames varies by ±0.014 g/cm³ (standard deviation), so 10 ps averages of two runs can differ by about 0.01.

For a density closer to the published values in NPT, add the long-range dispersion correction: `pair_modify mix arithmetic tail yes` in the `[lammps]` section of your own library, or in `system.settings`.
