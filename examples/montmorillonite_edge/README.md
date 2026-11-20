# Example: the library `clayff_edge` (Montmorillonite edge, Underwood & Bourg 2020) on the particle with edges

The montmorillonite particle with edges on every side (`src/soapff/libraries/clayff/edge/`: `clay_montmorillonite.data` and `.settings`, a single layer of 2,292 atoms; `clay_montmorillonite_2.data`, the particle twice along c, 4,584 atoms) typed with the `clayff_edge` library (`soapff.build("clayff_edge", ...)`), whose force field is taken whole from those files. The particle files were converted from the GROMACS files of Underwood and Bourg (J. Phys. Chem. C 124, 3702 (2020); `new.gro` and `v07.itp`, the first line of the data file); the LJ parameters of the settings file are the ClayFF values put on each type (the edge types take those of their parent type).

```
pip install "lammps[mpi]"     # the LAMMPS Python module (or LMP=/path/to/lmp: the lmp executable)
python run.py                 # packmol, when particle_water.extxyz is missing
```

LAMMPS runs with its Python module, or with the `lmp` executable of `$LMP` (the single point and the minimisations below).

1. `clay_montmorillonite.extxyz`, `clay_montmorillonite_2.extxyz`: the particle files reduced to elements, coordinates and box. soapff's files for them (`out_clay_montmorillonite*/`: `structure.extxyz`, `system.data`, `system.settings`, `charges.csv`, `report.json`) are compared with the particle files: types, charges, every pair_coeff, the bond and angle coefficients, the numbers of bonds and angles.
2. `single_point/`: a LAMMPS single point (Ewald 1e-5, cutoff 10 Å, special_bonds 0 0 0) of the original files (`original.data`, a copy with the blank line after "Angle Coeffs" that LAMMPS needs; `original.settings`) and of the written `system.data` / `system.settings` of `out_clay_montmorillonite/` (with `kspace_style ewald 1.0e-5` instead of PPPM); logs `log.original`, `log.written` (and their inputs `log.original.in`, `log.written.in`, when run with the `lmp` executable).
3. Minimisations (Ewald 1e-4) of the particle, the particle twice, and the particle with a water film 2–12 Å above it and 4 Na⁺ + 4 Cl⁻ (`particle_water.extxyz`, packmol seed 1; `out_particle_water/`). The particle, the water and the ions are typed in one build with `clayff_edge` (its references include water and ions), the water molecules taken apart from the particle in step 1, and compared with the right answer (the particle file, then ow, hw, Na and Cl for the added atoms). The water and Cl⁻ take ClayFF (flexible SPC), the particle and Na⁺ the edge model. Each minimisation is the force field (`system.settings`), then `kspace_style ewald 1.0e-4`, then the two-stage minimisation of `../clay_md.py` (`clay_md.minimize`). Its commands are written to `in.min` in its folder (after `include system.settings`; `lmp -in in.min` runs it with the `lmp` executable), its log is `log.min`.

## Results (`results.json`)

| | Atoms | Wrong type | Wrong charge | Smallest similarity | pair_coeff differing | Bonds, angles | Minimisation stopped by |
|---|---|---|---|---|---|---|---|
| particle | 2,292 | 0 | 0 | 1.000 | none | 288, 624 (as in the file) | energy tolerance |
| particle twice | 4,584 | 0 | 0 | 1.000 | none | 576, 1,248 | energy tolerance |

Written bonded coefficients: O-H 553.935 kcal/mol/Å², 1.0 Å; M-O-H 30.0191 kcal/mol/rad², 109.47°; edge metal-O-metal 800 kcal/mol/rad², 100° (all as in the files).

Single point (kcal/mol), the original files and the written ones: pe -509528.707277, bonds 7.090368, angles 5917.865731, LJ 38869.102638, real-space Coulomb -335112.681592, k-space -219210.084421, the same in every term to 1e-6.

The particle with water: 2,560 water molecules, 4 Na⁺ and 4 Cl⁻ (9,980 atoms: 2,268 mineral atoms, 7,680 water atoms, 32 ions with the 24 Na⁺ of the particle). The 2,292 atoms of the particle get the types and charges of the particle files, and the added atoms 2,560 ow, 5,120 hw, 4 Na and 4 Cl. The smallest similarity is 1.000, no bond between the water and the particle was taken out (`water_apart`), the net charge is 0, and the minimisation stopped by the energy tolerance.

All the files here (`results.json`, `run.log`, the `out_*` folders, `single_point/`, `particle_water.extxyz`) were written by `run.py` with the `lmp` executable of LAMMPS 22 Jul 2025 (the `particle_water.extxyz` of the earlier packmol run used again).
