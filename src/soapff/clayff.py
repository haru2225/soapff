"""engine "clayff": the ClayFF writer soapff_cui (soapff_cui.core: ClayFF of Cygan et al. 2004 with its flexible SPC water).
The types come from the two steps of the molecular engine (soapff_cui.steps): the bonds from the bond lengths of the library
([bond_lengths], checked against its references), and the SOAP of the atoms within 4 bonds gives the type and charge; the mineral, the water and the
ions are typed together.  The library file gives the references
and the build options of the table; the force field itself is ClayFF (soapff_cui/clayff.py), or the one of the references
when [library] forcefield names a LAMMPS data and settings file (the montmorillonite edge model)."""
import os

from soapff_cui import core, steps

# the bundled tables of the clayff libraries (npz files of soapff_cui.steps), for code that calls soapff_cui.core directly;
# soapff.build("clayff", ...) finds them through the library files (prebuilt = ...)
PREBUILT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "libraries", "clayff", "prebuilt")
LIBRARY = os.path.join(PREBUILT, "library_clayff_bulk.npz")                  # clayff: 11 minerals, (001) slabs, water, ions
LIBRARY_EDGE = os.path.join(PREBUILT, "library_clayff_edge.npz")             # clayff_edge: montmorillonite particle with edges


def prepare(structure, workdir, library=None, log=print, **opts):
    """soapff_cui.core.make with a table (npz; default LIBRARY), stopping when an atom is below the threshold (as
    soapff.build with strict=True): structure -> workdir/{system.data, system.settings, charges.csv, report.json}; returns
    the report (dict)"""
    opts.setdefault("strict", True)
    return core.make(structure, workdir, library or LIBRARY, log=log, **opts)


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
    """the table of the references for the two steps (soapff_cui.steps.build_table: the bond lengths and the SOAP rows), with
    the prune threshold of each reference and the bond lengths of [bond_lengths]"""
    steps.build_table([r["path"] for r in lib.references], out, lib.soap_params(),
                      prune=[r.get("prune") for r in lib.references], forcefield=_forcefield(lib),
                      depth=int(lib.soap.get("depth", steps.DEPTH)), bond_lengths=lib.bond_lengths, log=log)
    return out


def make(lib, structure, out_dir, log=print, strict=True, **opts):
    """soapff_cui.core.make with the table of the library"""
    table = lib.table_path(log=log)
    opts.setdefault("threshold", float(lib.soap.get("threshold", core.THRESHOLD)))
    opts.setdefault("strict", strict)
    return core.make(structure, out_dir, table, log=log, **opts)
