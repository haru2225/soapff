"""engine "clayff": the ClayFF writer soapff_cui (soapff_cui.core: ClayFF of Cygan et al. 2004 with its flexible SPC water).
The types come from the two steps of the molecular engine (soapff_cui.steps): the bonds from the bond lengths of the library
([bond_lengths], checked against its references), and the SOAP of the atoms within 4 bonds gives the type and charge; the mineral, the water and the
ions are typed together.  The library file gives the references
and the build options of the table, and the parameters of the force field in the tables of the molecular engine: [types]
(mass, LJ), [bonds] (the O-H bonds), [angles] (H-O-H of water, metal-O-H, the edge metal-O-metal), [lammps]; the charges are
those of the reference atoms.  Which bonds and angles there are is ClayFF's rule (soapff_cui.core.topology_from_bonds)."""
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
    """the force field of the library file for soapff_cui.core: the mass, element and LJ of every type, the rows of [bonds]
    and [angles] and the [lammps] settings (the charges are those of the reference atoms); None for a library without
    [types] (ClayFF of soapff_cui/clayff.py)"""
    if not lib.types:
        return None
    return dict(name=lib.library.get("forcefield_name", lib.name),
                source=lib.library.get("forcefield_source", lib.library.get("reference", "")),
                mass={t: v["mass"] for t, v in lib.types.items()}, element={t: v["element"] for t, v in lib.types.items()},
                lj={t: [v["epsilon"], v["sigma"]] for t, v in lib.types.items()},
                bonds=[[*k, *v] for k, v in lib.bonds], angles=[[*k, *v] for k, v in lib.angles], lammps=dict(lib.lammps))


def build_table(lib, out, log=print):
    """the table of the references for the two steps (soapff_cui.steps.build_table: the bond lengths and the SOAP rows), with
    the prune threshold of each reference and the bond lengths of [bond_lengths]"""
    steps.build_table([r["path"] for r in lib.references], out, lib.soap_params(),
                      prune=[r.get("prune") for r in lib.references], forcefield=_forcefield(lib),
                      depth=int(lib.soap.get("depth", steps.DEPTH)), bond_lengths=lib.bond_lengths, log=log)
    return out


def make(lib, structure, out_dir, log=print, strict=True, **opts):
    """soapff_cui.core.make with the table of the library and its force field (the cutoff and kspace of [lammps] by default)"""
    table = lib.table_path(log=log)
    opts.setdefault("threshold", float(lib.soap.get("threshold", core.THRESHOLD)))
    opts.setdefault("strict", strict)
    ps = lib.lammps.get("pair_style", "").split()
    if len(ps) > 1:
        opts.setdefault("cutoff", float(ps[1]))
    if lib.lammps.get("kspace_style"):
        opts.setdefault("kspace", lib.lammps["kspace_style"])
    return core.make(structure, out_dir, table, log=log, forcefield=_forcefield(lib), **opts)
