"""soapff: the force field of any structure for LAMMPS, from a library of reference structures (SOAP).

soapff types every atom and writes the force field (system.data, system.settings); the rest (minimisation, MD, analysis) is
written in LAMMPS, here with its Python module (pip install "lammps[mpi]"):

    import soapff
    from lammps import lammps

    s = soapff.build("oplsaa", "cell_PE.extxyz", "pe_out")   # library, structure file (any format ASE reads, or LAMMPS data)
    s.types[:3], s.charges[:3]                               # the type and charge of every atom (charges.csv)
    lmp = lammps()
    lmp.commands_list(s.commands())                          # the force field (system.settings, read_data system.data)
    lmp.commands_string('''
    minimize 1.0e-4 1.0e-6 1000 10000
    velocity all create 413 1 dist gaussian
    fix 1 all npt temp 413 413 100 iso 1 1 1000
    run 100000
    ''')

    soapff.libraries()         # the libraries that can be given by name
    soapff.new_library("mine", extends="clayff")         # ~/.soapff/libraries/mine.txt
    soapff.add_reference("mine", "annite.data")           # a reference structure with its types and charges

The command line does the same: soapff list | info LIB | build LIB STRUCTURE -o OUT | new NAME --extends LIB |
add-reference LIB FILE   (python -m soapff --help); a LAMMPS input then includes OUT/system.settings.

Each library is a text file (see soapff.library): the reference structures and the parameters of the force field (ClayFF,
OPLS-AA, TIP3P / SPC/E water and ions, ...).  Every library is typed and written the same way."""
from .api import System, add_reference, build, new_library, read_structure
from .library import Library, LibraryError
from .library import available as _available
from .library import load as load_library
from .molecular import TypingError, write_reference

__version__ = "1.0.1"


def libraries():
    """[{name, reference, path}] of the libraries that can be given by name"""
    return _available()


__all__ = ["build", "System", "libraries", "load_library", "Library", "LibraryError", "TypingError", "read_structure",
           "new_library", "add_reference", "write_reference"]
