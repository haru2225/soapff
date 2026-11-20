# Example: clay + polymer + water (the library `clayff_oplsaa`)

A smectite layer with water, Na⁺ and two polymer chains (PEO or PVA) in its interlayer, typed from the elements, positions and cell alone with `clayff_oplsaa` (ClayFF and OPLS-AA as one force field), and compared with the right answer. The test also checks that the force field written is ClayFF for the mineral, water and ions and OPLS-AA for the chains, with the Lorentz–Berthelot rules for every pair.

`clayff_oplsaa` extends `oplsaa` with the ClayFF types, the O–H bonds and the H–O–H and M–O–H angles of ClayFF, and the references of `clayff`; it is typed by the molecular engine. The metals and bridging O of ClayFF are named in no row of `[bonds]`, so their bonds of step 1 only give the SOAP its neighbourhood (the force field has no metal–O bond), and the mineral layer they hold is not checked for a whole-number charge (`report.json`: `bonds_not_in_force_field`, `frameworks`).

| Step | Script | Environment | Output |
|---|---|---|---|
| 1. compact 20-mer chains (600 K Langevin MD in a sphere), typed by Foyer | `chains.py` | foyer, openmm, parmed (as `../oplsaa`) | `chains/{PEO,PVA}.extxyz` |
| 2. layer, substitutions, packmol, right answer | `build.py` | soapff (+ packmol) | `structures/{PEO,PVA}_start.extxyz` |
| 3. typing of the start structures | `test_typing.py` | soapff | `results/typing.json` |
| 4. the force field of `clayff_oplsaa` vs those of `clayff` and `oplsaa` | `check_ff.py` | soapff | `results/check_ff.json` |

The layer is the 2:1 layer of the pyrophyllite crystal of `clayff`, repeated to 31.0 × 35.9 Å, with substitutions at random apart from each other: PEO, 16 octahedral Mg (a montmorillonite); PVA, 18 tetrahedral Al (a beidellite). The interlayer is 22 Å between the O planes. The right answer: the mineral, the ClayFF type of the bond-graph rule (`soapff_cui.core.rule_types`) with the charges of ClayFF Table 1; water ow/hw; Na⁺; the chains, Foyer's.

## Results

Typing (`results/typing.json`):

| System | Atoms | Wrong | Lowest similarity (ClayFF atoms / OPLS-AA atoms) |
|---|---|---|---|
| PEO | 3,104 | 0 | 0.963 / 0.988 |
| PVA | 3,098 | 0 | 0.992 / 0.986 |

Force field (`results/check_ff.json`): the force field of `clayff_oplsaa` against those `clayff` writes for the mineral, water and ions and `oplsaa` writes for the chains (each in the cell of the structure), atom by atom and term by term:
- the type and charge of every atom and the LJ of every type: the same;
- the bonds, angles and dihedrals (atoms and coefficients): those of the two together (PEO 1,612 bonds, 1,348 angles, 612 dihedrals; PVA 1,602, 1,324, 822);
- `pair_modify mix arithmetic` (Lorentz–Berthelot for every pair; `oplsaa` alone mixes by geometric means) and `special_bonds lj/coul 0.0 0.0 0.5` (ClayFF has no 1-4 pair).
