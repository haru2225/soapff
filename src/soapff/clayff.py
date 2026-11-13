"""engine "clayff": soapff-cui, the force-field writer of soapff-gui, unchanged (soapff_cui.core: SOAP library, ClayFF of
Cygan et al. 2004 with its flexible SPC water, the mineral typed first and the water and ions in a second step with the
solution library, all the checks of soapff-gui).  The library file gives the references and the build options of the SOAP
table; the force field itself is ClayFF (soapff_cui/clayff.py), or the one of the references when [library] forcefield names
a LAMMPS data and settings file (the montmorillonite edge model)."""
import os

from soapff_cui import core

from .library import load


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
