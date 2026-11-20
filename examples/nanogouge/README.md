# Nano-gouge: many edge-bearing montmorillonite particles in water (paper: Table "clay test cases", Fig. "nanogouge")

The initial state of a nano-gouge after Underwood & Bourg, J. Phys. Chem. C 124, 3702 (2020): N montmorillonite particles with edges on every side (the particle of the optional library `clayff_edge`, `../montmorillonite_edge/out_clay_montmorillonite/`, 2,292 atoms, ~6 nm disk, ClayFF + Lammers-type edges, 24 Na⁺) at random positions and orientations in a periodic box, in flexible SPC water. soapff types it from its elements, coordinates and box alone.

- `build_gouge.py`: the builder (numpy and scipy only). The particles are added one by one with a random position and orientation (no two atoms of two particles nearer than 3.5 Å); of 300 random configurations the most random one (the lowest order parameter S) is kept, and the box is filled with water at a water/clay mass ratio m_w/m_c = 3. Its data file carries the types of the particle (those of the particle's files) and ow, hw. `python build_gouge.py -n 3 --seed 1 -o gouge_N3`.
- `typing_test.py`: builds the gouge, types it with the edge library (`soapff.build("clayff_edge", ...)`, strict) and compares the types and charges with those of the builder's data file; for the smallest N also a LAMMPS single point (Ewald 1e-5, cutoff 10 Å, special_bonds 0 0 0) of soapff's force field against the builder's files. Needs soapff and the LAMMPS Python module (`pip install "lammps[mpi]"`).

```
python typing_test.py WORK 3          # -> typing_results.json; WORK/gouge_N3.*, WORK/out_N3_{wet,wet_min,dry}/
```

`minimised/`: the minimised structure of N = 3 (`wet_min`) and soapff's result on it, so that the typing can be checked without LAMMPS (the minimisation depends on the LAMMPS version): `builder_min.data.gz` (the builder's LAMMPS data file after the minimisation, LAMMPS 22 Jul 2025, with the builder's types and charges), `charges.csv.gz` and `report.json` (soapff's types, charges, similarities and margins; `soapff build clayff_edge` on the elements, coordinates and box of that file). They were written again with the current soapff, which gives the same types and charges as the earlier one (`typing_results.json`); with the table cached, typing these 66,900 atoms took 8 s on an Apple M5 (22 s in `typing_results.json`, with the earlier soapff).

## Results (`typing_results.json`)

N = 3 (`--seed 1`, S = 0.128; box 73.9 x 73.9 x 123.1 Å; 6,876 clay atoms, 20,008 water molecules, 66,900 atoms):

| Structure | Result |
|---|---|
| `dry`: the particles without the water (edges next to the edges and faces of the other particles, in vacuum) | all 6,876 atoms right (similarity 1.000); single point equal to the builder's in every term |
| `wet`: as built (water on a jittered grid with random orientations) | stops: in 16 places a water H is 1.11-1.25 Å from the O of another water molecule (overlapping molecules of the builder), 18 H below similarity 0.9 |
| `wet_min`: after a minimisation with the builder's force field (`../clay_md.py`, PPPM 1e-4) | all types and charges right but 12 edge Mg typed `mgb` for `mga` (same charge, LJ and angle constants in this force field); smallest similarity 0.912 (particles), 0.997 (water); single point equal to the builder's (total -1,851,521.5420 against -1,851,521.5417 kcal/mol: the positions written to extxyz at 1e-8 Å) |

In the minimised structure water O come to 1.92 Å of edge Al and 2.19 Å of edge Mg, nearer than the Al-O and Mg-O bond lengths (2.65 Å). soapff takes the water molecules apart from the particles in step 1 (an O bonded to two H only; in the references only the O of water has two H), so these 23 water-metal pairs (17 Al-O, 6 Mg-O) are not bonded and the particles and the water are typed each on its own bonds (`water_apart` in `typing_results.json`).
