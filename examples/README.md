# Examples and test cases of soapff

Run these from the root of the repository, with a Python where `pip install -e .` was run.

| File | What it shows |
|---|---|
| `quickstart.py [--run]` | Lists the libraries. Builds ClayFF kaolinite/water from an extxyz file, then Trp-cage with amber14sb from a plain xyz file (no types, no cell, so it is boxed). Shows `charges.csv`. `--run` runs LAMMPS. |
| `user_library.py WORK` | Your own library. `soapff new methanol --extends amber14sb` (command line), then `soapff.add_reference` with a methanol molecule, its types and its charges. ff14SB already has the CT/H1/OH/HO parameters. A methanol/water mixture read from a CIF file is then typed. The parent library refuses it, and adding to a bundled library is refused. |
| `clayff/test_clayff.py [WORK]` | Test case of the clay minerals. The 11 bulk crystals and the 11 (001) slabs of the `clayff` references are built as 2×2×1 supercells, and again with random displacements (0.05 Å; each H moves with its O, plus 0.02 Å). Kaolinite is read from a CIF file, and kaolinite with water (`clayff/kaolinite_water.extxyz`) is typed. Every atom is compared with the type and charge of the reference. Results in `results/clayff.json`. |
| `md_check.py WORK` | Short LAMMPS runs of what soapff writes. TIP3P and SPC/E boxes of 216 waters (2 ps NVT relaxation, 5 ps NVT, 20 ps NPT, SHAKE). Kaolinite with water (ClayFF, 0.5 + 2 ps NVT). Results in `results/md.json`. |
| `amber14sb/` | How the protein library was made from OpenMM's ff14SB, and the checks against OpenMM (types, charges, energies in vacuum and in water). See `amber14sb/README.md`. |

## Results of `clayff/test_clayff.py` (`results/clayff.json`)

| Structures | Atoms | Wrong types | Wrong charges (> 0.01 e) | Lowest similarity |
|---|---|---|---|---|
| 11 bulk crystals, 2×2×1 | 652–1,632 | 0 | 0 | 1.0000 |
| the same, displaced | 652–1,632 | 0 | 0 | 0.9934–0.9963 |
| 11 (001) slabs, 2×2×1 | 1,280–3,264 | 0 | 0 | 0.9990–0.9996 |
| the same, displaced | 1,280–3,264 | 0 | 0 | 0.9902–0.9961 |
| kaolinite from a CIF file | 408 | 0 | 0 | 1.0000 |

Kaolinite with water (1,659 atoms): all 816 mineral atoms and 843 water atoms are typed by SOAP; no atom falls back to a rule; net charge 0.

## Results of `md_check.py` (`results/md.json`; LAMMPS 22 Jul 2025, serial, Ewald 1e-5)

| System | Atoms | Run | Result |
|---|---|---|---|
| TIP3P, 216 waters | 648 | 27 ps | mean temperature of the NPT 300.5 K; density over the second half of the NPT 0.977 g/cm³ (10 Å cutoff, no tail correction; TIP3P is about 0.98–1.00) |
| SPC/E, 216 waters | 648 | 27 ps | 298.9 K; 0.988 g/cm³ (SPC/E is about 0.99–1.00) |
| kaolinite / water (ClayFF) | 1659 | minimisation, 0.5 ps relaxation, 2 ps NVT (0.5 fs) | mean temperature over the second half of the NVT 295.7 K; `final.data` written |

These checks use Ewald, not the PPPM the libraries choose. The local Homebrew LAMMPS (FFTW3 with threads) crashes when PPPM is set up a second time, as in the run after a minimisation. Other builds run the `in.lammps` soapff writes as it is.

For a density closer to the published values in NPT, add the long-range dispersion correction: `pair_modify mix arithmetic tail yes` in the `[lammps]` section of your own library, or in `system.settings`.
