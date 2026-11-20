"""ClayFF helpers of the examples (soapff itself types and writes every library the same way, soapff.build):

  core      reading structures (LAMMPS data, extxyz, cif, POSCAR), the ClayFF types of the bond-graph rule (the right
            answer of the clay tests), the boundary helpers of soapff.lammps and a ClayFF writer (the slabs of surfaces)
  clayff    the ClayFF parameters (Tables 1 and 2) and the flexible SPC water
  _clayff   the ClayFF types from the bond graph
  water     water and ions added to a structure with packmol (the test structures of the examples)
  surfaces  slabs of a crystal along a face (pymatgen; the basal (0 0 1) references)
"""

# part of soapff
__version__ = "0.1.0"
