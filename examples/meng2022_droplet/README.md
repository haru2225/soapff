# Example: the NaCl droplet of Meng et al. 2022 on montmorillonite (0 0 1)

The initial structure of Meng et al., Appl. Clay Sci. 224, 106513 (2022), their Droplet-1: a droplet of 9,171 water molecules with 112 Na⁺ and 112 Cl⁻ (3.83 wt%), about 80 Å across and about 5 Å above the (0 0 1) surface of montmorillonite, in a 209.2 x 217.4 x 201 Å box.

`paper/`: the (0 0 1) slab of dry Na montmorillonite (the bulk reference of the clayff library, `src/soapff/libraries/clayff/bulk/montmorillonite_dry_Na.data`; pymatgen, the vacuum above it) with c set to 201 Å, repeated 20 x 12 (206.9 x 215.5 Å, 78,240 atoms); the water, Na⁺ and Cl⁻ inside one sphere (packmol `inside sphere`) that holds 9,171 molecules at 1 g/cm³ (radius 40.3 Å), its bottom 5 Å above the surface; the force field (`soapff.build("clayff", ...)`), compared with the right answer (the ClayFF types of the bond-graph rule and the charges of Table 1, `soapff_cui.core.clayff_answer`).

```
python run.py     # packmol when paper/structure.extxyz is missing
```

`paper/` keeps the structure `structure.extxyz` (105,977 atoms), the files soapff writes for it (`system.data`, `system.settings`, `charges.csv`, `report.json`) and packmol's input and log (`packmol/`); packmol's large xyz files (`solid.xyz`, `packed.xyz`) are written by `run.py` but not kept in the repository. The structure is large, so no LAMMPS calculation was run for this test.

## Results (`results.json`)

| | Cell (Å) | Atoms | Placed | Nearest atom to the droplet (Å) | Droplet bottom above the surface (Å) | Wrong type or charge | Net charge (e) |
|---|---|---|---|---|---|---|---|
| `paper/` | 206.93 x 215.53 x 201.00 | 105,977 | 9,171 water, 112 Na⁺, 112 Cl⁻ | 6.38 | 5.28 | 0 | +0.096 |

packmol converged (37 s). Every atom has similarity 1.000 in step 2. The net charge is the rounding of ClayFF Table 1 over the 480 units of the slab (+0.0002 e each).
