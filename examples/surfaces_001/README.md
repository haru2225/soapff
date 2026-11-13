# Example: basal (0 0 1) slabs of the 11 layered minerals

`run.py` does two tests with the bulk references (`../library/refs_bulk/`):

1. **Which faces can be made.** Every face with |h|, |k|, |l| ≤ 1 (13 faces) of each mineral is cut with pymatgen (`soapff_cui.surfaces.make_slab`). A face is made only when it cuts no cation–O bond. Of the 143 faces, only the 11 (0 0 1) faces are made; the 132 others are refused.
2. **Slabs in a new cell.** Each (0 0 1) slab is cut again from the bulk doubled along a, with 30 Å of vacuum (the default). It is written as elements, coordinates and box only (`slabs/<mineral>_001.extxyz`) and typed with the `clayff` library. The truth is the type each atom had in the bulk and its ClayFF Table 1 charge.

```
python run.py        # pymatgen must be installed  ->  slabs/, out/<mineral>/, results.json
```

`out/<mineral>/` holds the slab (`structure.extxyz`) and the files soapff writes for it (`system.data`, `system.settings`, `in.lammps`, `charges.csv`, `report.json`).

## Results (`results.json`)

| Test | Result |
|---|---|
| Faces | 143 tried, 11 made (all (0 0 1)), 132 refused |
| (0 0 1) slabs, 11 minerals, 640–1,632 atoms | 0 wrong types, 0 fallback, lowest SOAP similarity 0.999 |
