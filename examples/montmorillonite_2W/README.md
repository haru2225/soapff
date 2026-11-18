# Example: Na montmorillonite with a 2-layer hydrate (1,554 atoms)

`initial.xyz` is a 3x2x1 supercell (31.12 x 35.95 x 14.76 Å) of the 32-water state of `../montmorillonite_swelling/` (the last structure of the NPT run with 32 water molecules, `run/swelling_states.lammpstrj`), reduced to elements, coordinates, and box only: one 2:1 layer with the composition of Cygan et al. (2004), Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆ per unit (layer charge 0.75 e per O₂₀(OH)₄, 3 Na⁺ per unit), and 8.0 water molecules per O₂₀(OH)₄ (2-layer hydrate, d₀₀₁ ≈ 14.8 Å). It has 960 clay atoms, 18 Na⁺ ions, and 192 water molecules (576 atoms); it contains neither charges nor types.

The last frame of its MD (`diffusion/start/final.data`) is one of the 19 MD frames with water on which `../clayff/water_md_frames.py` tests `clayff`; its force field (`diffusion/start/system.data`) is the start of the `hydrate_2W` MD of `../library_water/`.

## Procedure

`python diffusion/run_md.py` builds the force field of `initial.xyz` with soapff (the bundled `clayff` library) and writes the MD in LAMMPS (`diffusion/start/`): the minimization, a relaxation NVT of 10 ps, NVT 10 ps, and NPT (perpendicular to the plane only, 1 atm) 100 ps at 300 K, the protocol of `../clay_md.py`, with a trajectory every 10 ps and `compute msd` of the water O and Na⁺ every 0.1 ps. The time step is 0.5 fs, and the force field is ClayFF (Tables 1 and 2 of Cygan et al. 2004: the charges of Table 1 as they are, unrounded LJ parameters) with the flexible SPC water. It runs with the LAMMPS Python module (`pip install "lammps[mpi]"`; serial, in the Python process), or, with `--lmp "mpirun -np 3 lmp"`, with an `lmp` executable; `start/in.lammps` holds the same commands (`include system.settings`).

The stored run was made with the force field and the `in.lammps` of an earlier soapff, on another computer (Linux, Open MPI, LAMMPS 7 Feb 2024, 3 MPI processes). `run_md.py` writes the same commands (the earlier input also had `neighbor 2.0 bin` and `neigh_modify every 1 delay 0 check yes`, the LAMMPS defaults) the same `system.settings`, and a `system.data` with the same atoms, types, charges, bonds and angles (the bonds and angles listed in another order). The trajectory (`traj.lammpstrj`) is large and is not included; the MSD files (`msd_Ow.dat`, `msd_Na.dat`), the inputs, the log, and the last structure (`final.data`) are.

## Force field (`diffusion/start/report.json`)

| Item | Value |
|---|---|
| Atoms | 1,554 (clay 960, water 576, ions 18) |
| Types and charges | st +2.1 (186), at +1.575 (6), ao +1.575 (84), mgo +1.36 (12), ob -1.05 (408), obts -1.1688 (24), obos -1.1808 (48), oh -0.95 (72), ohs -1.0808 (24), ho +0.425 (96), Na +1.0 (18), ow -0.82 (192), hw +0.41 (384) |
| Bonds / angles | 480 / 384 (water O-H 384 and H-O-H 192; clay O-H 96 and M-O-H 192) |
| Total charge | +0.0012 e (the rounding of Table 1: +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al, 6 units; not shifted; LAMMPS adds a uniform neutralizing background) |

## MD (`diffusion/start/log.lammps`, 300 K)

- In the relaxation NVT the temperature is 264-316 K (second-half average 298 K), and in the NVT 291-310 K (second-half average 299 K).
- In NPT the interlayer distance d₀₀₁ (the box height) stays at 14.7-15.1 Å (average 14.87 Å over the last 50 ps). The x and y box lengths are fixed, so the pressure of the whole system is not controlled in-plane (average -150 atm in the second half).
