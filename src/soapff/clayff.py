"""engine "clayff": the ClayFF writer soapff_cui (soapff_cui.core: SOAP table, ClayFF of Cygan et al. 2004 with its flexible
SPC water, the mineral typed first and the water and ions in a second step with the solution library, and its checks).  The library file gives the references and the build options of the SOAP
table; the force field itself is ClayFF (soapff_cui/clayff.py), or the one of the references when [library] forcefield names
a LAMMPS data and settings file (the montmorillonite edge model)."""
import os

from soapff_cui import core

from .library import load

# the bundled tables of the clayff libraries (npz files of soapff_cui.core), for code that calls soapff_cui.core directly;
# soapff.build("clayff", ...) finds them through the library files (prebuilt = ...)
PREBUILT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "libraries", "clayff", "prebuilt")
LIBRARY = os.path.join(PREBUILT, "library_clayff_bulk.npz")                  # clayff: 11 minerals and their (001) slabs
LIBRARY_EDGE = os.path.join(PREBUILT, "library_clayff_edge.npz")             # clayff_edge: montmorillonite particle with edges
LIBRARY_WATER = os.path.join(PREBUILT, "library_clayff_water.npz")           # their solution libraries: the water and the ions,
LIBRARY_EDGE_WATER = os.path.join(PREBUILT, "library_clayff_edge_water.npz")  # typed after the mineral
SOLUTION_LIBRARY = {LIBRARY: LIBRARY_WATER, LIBRARY_EDGE: LIBRARY_EDGE_WATER}


def prepare(structure, workdir, library=None, log=print, **opts):
    """soapff_cui.core.make with a table (npz; default LIBRARY) and its solution table (SOLUTION_LIBRARY), stopping when an
    atom falls back to a rule (as soapff.build with strict=True): structure -> workdir/{system.data, system.settings,
    in.lammps, charges.csv, report.json}; returns the report (dict)"""
    opts.setdefault("stop_on_fallback", True)
    library = library or LIBRARY
    sol = SOLUTION_LIBRARY.get(os.path.abspath(library))
    if sol and os.path.isfile(sol):
        opts.setdefault("solution_library", sol)
    return core.make(structure, workdir, library, log=log, **opts)


def _forcefield(lib):
    spec = lib.library.get("forcefield")
    if not spec:
        return None
    tok = spec.split()
    if len(tok) != 2:
        raise ValueError(f"{lib.path}: forcefield = DATA SETTINGS (a LAMMPS data file with type names and its settings file)")
    data, settings = tok
    return core.read_forcefield(data, settings, lib.library.get("forcefield_name", lib.name),
                                lib.library.get("forcefield_source", lib.library.get("reference", "")))


def build_table(lib, out, log=print):
    """the SOAP table of the references (soapff_cui.core.build_library), with the prune threshold of each reference"""
    soap = {k: v for k, v in lib.soap_params().items()}
    wis = lib.soap.get("water_in_soap", "no").lower() in ("yes", "true", "1")
    core.build_library([r["path"] for r in lib.references], [], out=out, soap=soap,
                       prune=[r.get("prune") for r in lib.references], water_in_soap=wis, protect=0,
                       forcefield=_forcefield(lib), log=log)
    return out


def make(lib, structure, out_dir, log=print, strict=True, **opts):
    """soapff_cui.core.make with the table of the library and of its solution library ([library] solution)"""
    table = lib.table_path(log=log)
    sol = lib.library.get("solution")
    sol_table = None
    if sol:
        sol_table = load(sol).table_path(log=log)
    opts.setdefault("threshold", float(lib.soap.get("threshold", core.THRESHOLD)))
    opts.setdefault("stop_on_fallback", strict)
    return core.make(structure, out_dir, table, solution_library=sol_table, log=log, **opts)
