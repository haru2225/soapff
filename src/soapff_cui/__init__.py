"""The engine "clayff" of soapff (ClayFF, Cygan et al. 2004, with its flexible SPC water), used through soapff.build:

  core      reading structures (LAMMPS data, extxyz, cif, POSCAR), the two steps with a table (make_steps) and the LAMMPS
            writer (system.data, system.settings, charges.csv, report.json)
  steps     the table of the two steps (bond lengths and SOAP rows of the references) and the typing (type_atoms)
  clayff    the ClayFF parameters (Tables 1 and 2) and the flexible SPC water
  _clayff   the ClayFF types from the bond graph (for the cross-check of report.json and the reference structures)
  water     water and ions added to a structure with packmol (the test structures of the examples)
  surfaces  slabs of a crystal along a face (pymatgen; the basal (0 0 1) references)
"""

# part of soapff
__version__ = "0.1.0"
