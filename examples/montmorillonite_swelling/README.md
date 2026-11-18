# Example: swelling curve of Na montmorillonite (same calculation as Cygan et al. 2004)

For the Na montmorillonite Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆·nH₂O of Cygan et al. (2004, J. Phys. Chem. B 108, 1255) (4 unit cells, one layer), NPT was run while water was added to the interlayer, and the relation between the basal spacing d001 and the amount of water was obtained.

## 1. Structure (`build_structure.py`)

A pyrophyllite layer (the library reference) is replicated 2x2x1, 1 tetrahedral Si is changed to Al, 2 octahedral Al atoms are changed to Mg, and 3 Na⁺ ions are placed in the interlayer. The Mg positions are the ones with the lowest static ClayFF energy (with the box also relaxed) among all 120 ways of choosing 2 of the 16 octahedral sites (4 placements of Na⁺ for each; the 29 arrangements in which an O bonds to both a tetrahedral Al and a Mg (the ClayFF obss type) are excluded because the library does not have them). All results are written to `screening.csv`, and the chosen structure (163 atoms, dry) to `dry.extxyz`.

`python build_structure.py [nproc]` does the relaxations with the LAMMPS Python module (`pip install "lammps[mpi]"`), one serial LAMMPS per process. The stored `screening.csv` and `dry.extxyz` come from an earlier run of the same commands with the `lmp` executable and an earlier LAMMPS. With LAMMPS 22 Jul 2025 the relaxed energies of some pairs come out different (another local minimum after the box relaxation: up to 13 kcal/mol for the 6 pairs checked; the chosen pair Mg 2, 42 gives -39219.01 kcal/mol instead of -39218.87), so a rerun may rank the pairs differently.

Na⁺ in a dry interlayer is an environment that no other reference has, so this structure (with ClayFF types and standard charges, `dry_reference.data`) has been added to the library.

## 2. Calculation (`make_run.py`)

`python make_run.py` writes `run/` (the folder is replaced; `--out DIR` writes elsewhere):

1. The force field of `dry.extxyz` with the `clayff` library (all atoms above the similarity threshold, 0 mismatches in the type cross-check).
2. The structure has no water type, so a water type is added to `system.data` and `system.settings`, and `water.mol` is written.
3. `in.lammps`: the minimization of `../clay_md.py`, then the swelling curve of `swelling.py`, with the numbers of water molecules "4 8 12 16 20 25 32 40 48 59". Each state is equilibrated with NPT for 30 ps, and the following 20 ps are recorded.

`python make_run.py --run` also runs it here, with the LAMMPS Python module (`pip install "lammps[mpi]"`). The commands are the same as in `in.lammps`. In the input file the loops are `label` / `jump SELF` (they need a file); with the Python module they are Python loops (`swelling.run`). Both give the same `swelling.dat` (checked on a short run). `in.lammps` can also be run with the `lmp` executable: `cd run && mpirun -np 4 lmp -in in.lammps -log log.lammps`.

`swelling.py` holds the swelling-curve code (`swelling_schedule`, `swelling_lines` for an input file, `run` for the Python module). It was part of soapff (`soapff_cui.core`) before soapff was reduced to writing the force field.

The stored run in `run/` was made with the `in.lammps` of the earlier soapff and the `lmp` executable (`soapff run --np 4, an option of the earlier soapff`: 4 MPI processes on another computer, LAMMPS 7 Feb 2024, about 5 minutes). Its commands are the same as those `make_run.py` writes now (the swelling block is identical; the earlier file also set `neighbor 2.0 bin` and `neigh_modify every 1 delay 0 check yes`, the LAMMPS defaults).

For each state, water molecules are inserted into the interlayer (the region between the clay and the box edge, on the side across the periodic boundary) with `create_atoms ... mol` (at least 2 Å from other atoms; if they do not fit, the box is stretched by 0.5 Å in z), the clay is moved back to the center of the box, and then a minimization, a 2 ps relaxation with damping, and NPT (300 K, 1 atm) are run. The NPT moves x, y, z, and the xy tilt, and keeps the tilts across the layer, xz and yz, fixed (when water enters, the layer slides and the tilt becomes large, which makes LAMMPS stop; d001 = lz does not depend on this).

The charges are the values of Table 1 of Cygan et al. (2004) as they are, and the Lennard-Jones parameters are written without rounding. Because the O charges of the table are rounded to 4 decimals, the net charge of this system is +0.0002 e (LAMMPS neutralizes it with a uniform background). Water is the flexible ClayFF SPC water (timestep 0.5 fs, PPPM 1e-5).

`run/` contains the inputs (`in.lammps`, `system.data`, `system.settings`, `water.mol`) and the outputs (`log.lammps`, `swelling.dat`, and the last structure of each state `swelling_states.lammpstrj` (large, so it is not in git)).

## 3. Analysis (`analyze_swelling.py`)

From `swelling.dat`, the mean and standard deviation of lz (= d001) over the last 20 ps of each state, and the water content (g H₂O / g clay; the clay is all atoms other than water, 2,945 g/mol), are obtained (`run/swelling.csv`, `swelling.png`). `python analyze_swelling.py [folder]` (default `run`).

| Water (molecules) | g/g | d001 (Å) |
|---|---|---|
| 4 | 0.024 | 11.56 |
| 8 | 0.049 | 12.03 |
| 12 | 0.073 | 12.19 |
| 16 | 0.098 | 12.29 |
| 20 | 0.122 | 12.47 |
| 25 | 0.153 | 13.66 |
| 32 | 0.196 | 14.89 |
| 40 | 0.245 | 15.69 |
| 48 | 0.294 | 17.27 |
| 59 | 0.361 | 19.05 |

The 1-layer hydrate (about 12.0-12.5 Å at 0.05-0.12 g/g), the expansion to the 2-layer hydrate (about 13.7-14.9 Å at 0.15-0.2 g/g), and the continued swelling with more water were the same as in the MD of Fig. 8 of Cygan et al.
