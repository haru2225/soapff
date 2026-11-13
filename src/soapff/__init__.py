"""soapff: the force field of any structure for LAMMPS, from a library of reference structures (SOAP).

    import soapff
    s = soapff.build("clayff", "kaolinite.cif")          # library, structure file (any format ASE reads, or LAMMPS data)
    s = soapff.build("amber14sb", "protein.pdb", md=False)
    s.files                    # system.data, system.settings, in.lammps, charges.csv, report.json in s.out_dir
    s.run(np=4)                # LAMMPS (lmp executable, mpirun) in s.out_dir
    s.thermo()["PotEng"]       # the thermo output of the last run
    lmp = s.lammps()           # or the LAMMPS Python module, with in.lammps run

    soapff.libraries()         # the libraries that can be given by name
    soapff.new_library("mine", extends="clayff")         # ~/.soapff/libraries/mine.txt
    soapff.add_reference("mine", "annite.data")           # a reference structure with its types and charges

The command line does the same: soapff list | info LIB | build LIB STRUCTURE -o OUT [--run --np 4] | run OUT | new NAME
--extends LIB | add-reference LIB FILE   (python -m soapff --help).

Each library is a text file (see soapff.library): the reference structures and, for the molecular engine, the force-field
parameters.  Engines: "clayff" (soapff-cui of soapff-gui, unchanged: ClayFF, Cygan et al. 2004) and "molecular" (any force
field written in the library: AMBER ff14SB proteins, TIP3P / SPC/E water and ions, ...)."""
from .api import System, add_reference, build, new_library, read_structure, read_thermo
from .library import Library, LibraryError
from .library import available as _available
from .library import load as load_library
from .molecular import TypingError, write_reference

__version__ = "0.1.0"


def libraries():
    """[{name, engine, description, path}] of the libraries that can be given by name"""
    return _available()


def run(out_dir, np=1, lmp=None, mpirun=None, input="in.lammps", echo=True):
    """run LAMMPS in an output folder of soapff (or any folder with an input); returns the path of log.lammps"""
    from .lammps import run as _run
    return _run(out_dir.out_dir if isinstance(out_dir, System) else out_dir, np_=np, lmp=lmp, mpirun=mpirun, input=input,
                echo=echo)


__all__ = ["build", "System", "libraries", "load_library", "Library", "LibraryError", "TypingError", "read_structure",
           "new_library", "add_reference", "write_reference", "run", "read_thermo"]
