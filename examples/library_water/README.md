# Example: the libraries with water and aqueous ions (solution libraries) and the two-step typing

The default solution libraries are built from **initial structures only** (what packmol builds: not relaxed, no MD frame):

| File (also in `src/soapff/libraries/clayff/prebuilt/`) | For | Rows | Contents |
|---|---|---|---|
| `library_clayff_water.npz` | ClayFF library | 4,124 | the 3,265 rows of the water-free mineral library (11 bulk layered minerals and their basal (0 0 1) slabs) + 320 rows of the packmol liquids + 539 rows of the packmol systems in vacuum (below) |
| `library_clayff_edge_water.npz` | edge library (Underwood & Bourg 2020) | 3,151 | the two particle files of the edge library (2,292 rows) + the same liquids and vacuum systems (859 rows) |

| System | Composition | Atoms |
|---|---|---|
| `water` | 512 H2O, 24.83 Å cube | 1,536 |
| `NaCl`, `KCl`, `CsCl` | 496 H2O + 8 M+ + 8 Cl- | 1,504 |
| `CaCl2`, `BaCl2` | 494 H2O + 6 M2+ + 12 Cl- | 1,500 |
| `vac_water_film`, `vac_water_film_thin` | water film 20 and 10 Å thick, 24.83 Å periodic in x and y, 50 Å cell (`build_refs_vacuum.py`, seed 2) | 1,236, 618 |
| `vac_water_droplet` | water sphere of radius 13 Å in a 50 Å cube | 924 |
| `vac_NaCl_film`, `vac_KCl_film`, `vac_CsCl_film`, `vac_CaCl2_film`, `vac_BaCl2_film` | the 20 Å film with 7 + 7 ions (5 + 10 for Ca, Ba) | 1,250-1,251 |

Water: flexible SPC; ions: ClayFF Table 1 (also in the edge library: the model of the particle has neither Cl- nor water). Rows of the systems with water are pruned at similarity 0.99, slab rows at 0.999. The 2-layer hydrate (its structure comes from an NPT run) is not in these libraries.

Typing (`soapff_cui.core.make(..., solution_library=...)`, `soapff.clayff.prepare` chooses the solution library of the library): step 1 types the mineral atoms with the water-free library on the structure without the water molecules (and without the atoms marked "added"); step 2 types the water, the added ions and the step-1 atoms below the threshold with the solution library: the water and the added ions on the structure of the water and the ions alone (the mineral taken away for the SOAP, the atoms kept), a mineral atom on the whole structure (water and ions may only match water or ion types; a mineral atom may match any type except water). The libraries have no mineral-water interface: the mineral is typed first, on the dry structure, and the water on its own, so neither a water-covered mineral nor a water at a mineral surface is needed; the water films and the droplet in vacuum give the environment of a water at a free surface. (The first version of the solution libraries had 411 interface rows (ClayFF: the 11 `if_<mineral>` packmol structures) and 534 rows of a particle with a water film (edge); they are needed only to type a structure with water by its whole-structure SOAP alone, see Results.)

The two libraries are not copied into this directory (they are in `src/soapff/libraries/clayff/prebuilt/`; `build_library_initial.py` and `build_library_clayff_edge_water.py` write them here again).

The libraries built with the frames of a 300 K MD (the previous default, 7,435 rows; `build_library.py`, `build_refs.py`, `test.py`, `soap_only.py`, `mineral_in_water.py`) are kept as `library_clayff_water_md.npz`; its MD frames are not used by the new libraries.

## Files

- `build_refs_vacuum.py` -> `refs/vac_*/`: the water and the solutions in vacuum (packmol, initial structures).
- `build_library_initial.py` -> `library_clayff_water.npz`, `build_library_initial.json`; `build_library_clayff_edge_water.py` -> `library_clayff_edge_water.npz`, `build_library_clayff_edge_water.json` (rows per type, one-to-one check).
- `build_refs.py`: the reference systems (inputs; the MD of the previous library); `run_remote.py`: the same MD on another computer over SSH (`soapff.remote`; the password is read from the environment variable `SOAPFF_PW` only). `refs/<system>/`: inputs (the packmol structure of each system is its `system.data`; the MD output is not used by the new libraries).
- `water_common.py`: the nearest row of a library for every atom of a structure (whole-structure SOAP), the true types of a structure with water.
- `initial_only.py` -> `initial_only.json`: the 10 packmol interfaces (`../interfaces_001_water/`) and the 6 shapes (`../interface_shapes/`), the effect of the water on the SOAP of the mineral atoms (with the MD library; the water-free library on the whole structure and on the structure without the water).
- `test_two_step.py` -> `test_two_step.json`: the same 16 initial structures (28,311 atoms) with `soapff.clayff.prepare` (the two steps) and with the whole-structure SOAP alone, against `library_clayff_water.npz`.
- `test_edge_water.py` -> `test_edge_water.json`: the edge particle and the particle twice with a water film (seeds 2 and 3; no such structure in the library), the two steps and the whole-structure SOAP, against `library_clayff_edge_water.npz`.

## Results

- One-to-one check (nearest other row of the same element): `library_clayff_water.npz` 4,124 / 4,124; `library_clayff_edge_water.npz` 3,151 / 3,151 rows with the same charge and type.
- `test_two_step.json`, 16 initial structures (28,311 atoms; 4,037 water O, 422 ions), the two steps (similarity threshold 0.9): 0 wrong types, 0 fallback, 0 atoms by the rule. The water and the ions alone against the library: 0 wrong, 0 below 0.9 (minimum 0.954 for the 10 interfaces; 0.915 for the 6 shapes: a Cl- at the top of the water film of `ions` with another Cl- 2.2 Å away). The whole-structure SOAP alone (no first step, no removal of the mineral): 6 wrong types, 8 atoms below 0.9 (none of them a mineral atom). Without the water films and the droplet in vacuum, 3 atoms of the shapes would be below 0.9 in the second step (minimum 0.880; 24 below 0.95).
- Mineral atoms and ions of the 10 interfaces (4,328 atoms) with the whole-structure SOAP: water-free library 2 below 0.9 (minimum 0.885) and 1 wrong type; with the water removed (step 1) 0 and 0 (minimum 0.999) (`initial_only.json`; 23 atoms were below 0.95: the hydroxyl H of kaolinite (11) and brucite (7) that point to the water and 5 interlayer K+ and Cs+).
- `test_edge_water.json`: the two steps give 0 wrong types, 0 fallback and 0 atoms by the rule on both structures (9,980 and 12,280 atoms); the water and the ions alone against the library 0 wrong, 0 below 0.9 (minimum 0.973); the whole-structure SOAP alone: 0 and 2 wrong, 0 and 1 below 0.9; against the water-free edge library no atom of the particle falls below 0.9 (minimum 0.948 and 0.923; with 0.95: 1 and 2 atoms).
