# Example: the water and ion references of the clay libraries, and the MD frames with water

The libraries `clayff` and `clayff_edge` type the mineral, the water and the ions in one build. Their references include water and the ions of ClayFF Table 1 (`src/soapff/libraries/clayff/water/`): those files are the `system.data` of the packmol structures made here (initial structures, not relaxed). The 300 K MD of the same liquids, of the 11 (001) slabs with water and of the 2-layer hydrate gives the test frames of `../clayff/water_md_frames.py`; none of these frames is a reference.

| System | Composition | Atoms | Use |
|---|---|---|---|
| `water` | 512 H₂O in a 24.83 Å cube | 1,536 | reference (`water.data`); MD frame |
| `NaCl`, `KCl`, `CsCl` | 496 H₂O + 8 M⁺ + 8 Cl⁻ (about 0.9 mol/L) | 1,504 | references; MD frames |
| `CaCl2`, `BaCl2` | 494 H₂O + 6 M²⁺ + 12 Cl⁻ (about 0.7 mol/L) | 1,500 | references; MD frames |
| `vac_water_film`, `vac_water_film_thin` | water film 20 and 10 Å thick, 24.83 Å periodic in x and y, 50 Å cell | 1,236, 618 | references (a water at a free surface) |
| `vac_water_droplet` | water sphere of radius 13 Å in a 50 Å cube | 924 | reference |
| `vac_NaCl_film`, `vac_KCl_film`, `vac_CsCl_film`, `vac_CaCl2_film`, `vac_BaCl2_film` | the 20 Å film with 7 + 7 ions (5 + 10 for Ca, Ba) | 1,250-1,251 | references |
| `if_<mineral>` | the (001) slab of each of the 11 minerals (`src/soapff/libraries/clayff/basal/`) with the vacuum above it filled with water (1 g/cm³, packmol seed 2; the test structures of `../interfaces_001_water/` use seed 1) | 692-1,659 | MD frames |
| `hydrate_2W` | the 2-layer Na-montmorillonite hydrate of `../montmorillonite_2W/` (`hydrate_2W_start.data`) | 1,554 | MD frame |

Water: flexible SPC; ions: ClayFF Table 1. In the libraries the rows of these systems are pruned at similarity 0.99. A water molecule is typed on its own bonds (taken apart from the mineral in step 1) and an ion on itself, so the liquids give the water and the ions in water, and the films and the droplet in vacuum give a water or an ion at a free surface, as at the surface of a water film or a droplet on a mineral.

## Files

- `build_refs_vacuum.py` -> `refs/vac_*/`: the water and the solutions in vacuum (packmol, initial structures; the force field files only, no MD).
- `build_refs.py` -> `refs/<system>/` for the liquids, `if_<mineral>` and `hydrate_2W`: the structure, its force field (`soapff_cui.core.make` with the bundled `clayff` table) and a 300 K MD with the LAMMPS Python module (`pip install "lammps[mpi]"`), with the commands of `../clay_md.py`: minimisation, 5 ps NVT relaxation, then the liquids and the hydrate 5 ps NVT + 20 ps NPT at 1 atm (liquids isotropic, hydrate full cell), the interfaces 25 ps NVT; 0.5 fs, frames every 1 ps. `python build_refs.py [--jobs 3] [system ...]`. The same commands are written to `refs/<system>/in.lammps`; `python build_refs.py --inputs-only [system ...]` writes only the inputs (`system.data`, `system.settings`, `in.lammps`) for the `lmp` executable.
- `refs/<system>/`: the inputs and the MD output (`log.lammps`, `final.data`, `msd_*.dat`). The stored systems were typed by an earlier soapff, which typed the water and the ions by rule; the bundled `clayff` library gives them the same types and charges. The stored MD runs were made with the `in.lammps` that soapff wrote before, which has the same commands as `build_refs.py` now: `water`, `NaCl` and `KCl` with the `lmp` executable here (LAMMPS 29 Aug 2024), the others on another computer (LAMMPS 7 Feb 2024). The trajectories (`traj.lammpstrj`) are large and are not included.
- `hydrate_2W_start.data`: the `system.data` of `../montmorillonite_2W/diffusion/start/`, the start of the `hydrate_2W` MD.
