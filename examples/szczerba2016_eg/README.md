# Example: ethylene glycol in a smectite (Szczerba & Kalinichev 2016)

The SWy-1 montmorillonite model of Szczerba & Kalinichev, *Clays Clay Miner.* **64**, 488 (2016), with Ca²⁺, ethylene glycol (EG) and water in the interlayers, typed with `clayff_oplsaa` (ClayFF and OPLS-AA as one force field, Lorentz–Berthelot for every pair, as their CLAYFF + OPLS-aa); and the density of liquid EG with the OPLS-AA soapff writes, against their Table 2.

EG is not among the references of `oplsaa`, but its atoms are: the HO–CH₂– at the ends of the PVA references. The bundled `oplsaa` and `clayff_oplsaa` give every EG atom Foyer's type and charge without more references. The EG molecule typed by Foyer (`eg_refs.py`) is only the molecule packed into the models and the liquid, and their right answer.

| Step | Script | Environment | Output |
|---|---|---|---|
| 1. one EG molecule, typed by Foyer | `eg_refs.py` | foyer, openmm, parmed (as `../oplsaa`) | `refs/EG.extxyz` |
| 2. SWy-1 + EG + water models, right answer | `build.py` | soapff (+ packmol) | `structures/SWy-1_EG{1.7,1.5}_start.extxyz` |
| 3. typing of the models with `clayff_oplsaa` | `test_typing.py` | soapff | `results/typing_start.json` |
| 4. liquid EG: 230 molecules, typed with `oplsaa`, LAMMPS input | `eg_liquid.py build` | soapff (+ packmol) | `liquid/` |
| 5. liquid EG: 0.1 ns NVT, 0.4 ns NPT at 293 K and 1 bar | `cd liquid && mpirun -np 4 lmp -in in.md` | LAMMPS | `liquid/density.dat` |
| 6. the same with their cutoff, 8.5 Å, no tail correction | `eg_liquid.py build85`, then `cd liquid_8.5 && mpirun -np 4 lmp -in in.md` | soapff, LAMMPS | `liquid_8.5/density.dat` |
| 6b. the same box typed with `clayff_oplsaa`, the library of the models (Lorentz–Berthelot, no tail correction), with their cutoff, 8.5 Å | `eg_liquid.py build_clay`, then `cd liquid_clayff_oplsaa && mpirun -np 6 lmp -in in.md` | soapff, LAMMPS | `liquid_clayff_oplsaa/density.dat` |
| 7. density, last 0.3 ns | `eg_liquid.py result` | numpy | `results/eg_liquid.json` |

The models (as theirs): the pyrophyllite crystal of `clayff` repeated to 8 × 4 × 2 unit cells (2560 atoms); in each layer 18 of the octahedral Al replaced by Mg (0.28 per half unit cell, phuc), as far apart as can be (smallest Mg–Mg distance 7.9 Å); per interlayer 9 Ca²⁺, 51 water (0.8 phuc) and 109 or 96 EG (1.7 or 1.5 phuc). The right answer: the mineral, the ClayFF type of the bond-graph rule (`soapff_cui.core.rule_types`); water ow/hw; Ca²⁺; EG, Foyer's.

## Results

Typing (`results/typing_start.json`):

| Model | Atoms | Wrong | Lowest similarity (ClayFF atoms / OPLS-AA atoms) |
|---|---|---|---|
| 1.7 EG phuc | 5,064 | 0 | 0.963 / 0.971 |
| 1.5 EG phuc | 4,804 | 0 | 0.963 / 0.971 |

The net charge of each model is +0.0072 e: ClayFF Table 1 gives −0.9998 e per octahedral Mg (36 Mg, 18 Ca²⁺).

Density of liquid EG at 293 K and 1 bar (`results/eg_liquid.json`; the mean of the last 0.3 ns of the NPT, standard error of 5 blocks). Every atom of the box gets Foyer's type and charge with both libraries:

| Library, interactions | Density (g/cm³) |
|---|---|
| `clayff_oplsaa` with the cutoff 8.5 Å: Lorentz–Berthelot, no tail correction | 1.040 ± 0.002 |
| `oplsaa` with the cutoff 8.5 Å, no tail correction: geometric mixing | 1.053 ± 0.001 |
| `oplsaa`: geometric mixing, cutoff 12 Å, tail correction | 1.080 ± 0.002 |
| Szczerba & Kalinichev, Table 2, OPLS-aa | 1.060 |
| experiment (their Table 2) | 1.115 |

The MD used LAMMPS 22 Jul 2025 (PyPI `lammps[mpi]`), 4 MPI processes (6 for `liquid_clayff_oplsaa`). With the force field of the smectite models the density is 2 % below theirs; most of the difference is the mixing rule (1.053 g/cm³ with geometric mixing and the same cutoff).
