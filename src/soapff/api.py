"""soapff.build(library, structure): a structure file of any format in, the force field of the library for LAMMPS out
(system.data, system.settings).  What is done with it (minimisation, MD, ...) is written in LAMMPS, e.g. with its Python
module:  lmp = lammps(); lmp.commands_list(system.commands()); lmp.command("fix 1 all npt ...")."""
import collections
import csv
import json
import os

import numpy as np

from . import library as _library
from . import molecular as M


# ------------------------------------------------------------------------------------------------- structures
def read_structure(path, **ase_kwargs):
    """any structure file -> ase.Atoms: every format ASE reads (cif, xyz, extxyz, pdb, POSCAR/CONTCAR, gro, mol2 (ASE), sdf,
    xsf, ...; ase_kwargs go to ase.io.read, e.g. format= or index=) and LAMMPS data files (.data, .lmp, .lammps; elements from
    the type names in the comments of Masses, else from the masses)"""
    if not isinstance(path, (str, os.PathLike)):
        raise TypeError("give the path of a structure file (or pass an ase.Atoms to soapff.build)")
    path = os.fspath(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    if path.endswith((".data", ".lmp", ".lammps")) or ".data" in os.path.basename(path) and not ase_kwargs:
        from soapff_cui.core import read_lammps_data
        return read_lammps_data(path)[0]
    import ase.io
    at = ase.io.read(path, **ase_kwargs)
    if isinstance(at, list):
        at = at[-1]
    return at


def _has_cell(at):
    """a usable cell: three vectors, none shorter than 2 A (PDB files of NMR structures carry a 1 A placeholder cell)"""
    c = np.array(at.cell)
    return bool(at.cell.rank == 3 and np.linalg.norm(c, axis=1).min() > 2.0)


def boxed(at, padding=10.0):
    """an Atoms without a usable cell put in an orthogonal box `padding` A larger than the atoms on every side (periodic)"""
    at = at.copy()
    lo, hi = at.positions.min(0), at.positions.max(0)
    at.positions = at.positions - lo + padding
    at.cell = np.diag(hi - lo + 2.0 * padding)
    at.pbc = True
    return at


# ------------------------------------------------------------------------------------------------- the result
class System:
    """an output folder of soapff: system.data, system.settings (the force field for LAMMPS), charges.csv, report.json"""

    def __init__(self, out_dir, report=None, atoms=None):
        self.out_dir = os.path.abspath(out_dir)
        self.report = report if report is not None else json.load(open(os.path.join(self.out_dir, "report.json")))
        self.atoms = atoms

    def __repr__(self):
        r = self.report
        return (f"<soapff System {r.get('library_name', os.path.basename(str(r.get('library', ''))))}: {r.get('atoms')} atoms, "
                f"net charge {r.get('charge_sum', 0.0):+.4f} e, {self.out_dir}>")

    @property
    def files(self):
        return {k: os.path.join(self.out_dir, f) for k, f in (("data", "system.data"), ("settings", "system.settings"),
                                                              ("charges", "charges.csv"), ("report", "report.json"))}

    def commands(self, rigid=True):
        """the LAMMPS commands of the force field (system.settings, with the path of system.data made absolute), for the
        LAMMPS Python module:  lmp.commands_list(s.commands()).  rigid=False leaves out the SHAKE of the rigid molecules
        (fix rigid; e.g. a single point of the flexible molecule)."""
        L = []
        for line in open(self.files["settings"]):
            c = line.split("#")[0].strip()
            if not c:
                continue
            if c.split()[0] == "read_data":
                c = f'read_data "{self.files["data"]}"' + (" " + " ".join(c.split()[2:]) if len(c.split()) > 2 else "")
            if not rigid and c.startswith("fix rigid "):
                continue
            L.append(c)
        return L

    @property
    def type_names(self):
        """the type names in the order of the LAMMPS atom types (type 1 first)"""
        tid = self.report["type_ids"]
        return sorted(tid, key=tid.get)

    @property
    def type_ids(self):
        """{type name: LAMMPS atom type}"""
        return dict(self.report["type_ids"])

    @property
    def elements(self):
        """the element of each LAMMPS atom type (type 1 first), e.g. for dump_modify element"""
        return list(self.report["type_elements"])

    def table(self):
        """the rows of charges.csv (one dict per atom: id, element, type, charge, ...)"""
        return list(csv.DictReader(open(self.files["charges"])))

    @property
    def types(self):
        rows = self.table()
        key = "type" if rows and "type" in rows[0] else "clayff_type"
        return [r[key] for r in rows]

    @property
    def charges(self):
        return np.array([float(r["charge"]) for r in self.table()])


# ------------------------------------------------------------------------------------------------- build
def build(library, structure, out_dir=None, *, cutoff=None, kspace=None, boundary="p p p", strict=True, padding=10.0,
          log=print, **options):
    """type every atom of `structure` with `library` and write its force field for LAMMPS.

    library    a library name (soapff.libraries() lists them: clayff, clayff_edge, oplsaa, tip3p, spce, ...) or
               the path of a library .txt file
    structure  a structure file of any format ASE reads (cif, xyz, extxyz, pdb, POSCAR, ...) or a LAMMPS data file, or an
               ase.Atoms
    out_dir    the output folder (default: <structure name>_soapff in the current folder): system.data, system.settings,
               charges.csv, report.json
    cutoff     the cutoff of LJ and of the real-space Coulomb (A; default: the library's)
    kspace     the long-range Coulomb, e.g. "pppm 1.0e-5" (default: the library's; MSM or slab are chosen for non-periodic
               boundaries)
    boundary   the LAMMPS boundary, e.g. "p p f"
    strict     stop (soapff.TypingError) when an atom is not in the library (its environment has no counterpart in the
               references) instead of writing the force field with the nearest reference
    padding    molecular engine: a structure without a cell (xyz, pdb without CRYST1) is put in a box this much larger (A)
    options    clayff engine: the other options of soapff_cui.core.make (e.g. reserve_water, water_rule)
    Returns a soapff.System: its commands() are the force field for the LAMMPS Python module (lmp.commands_list), and its
    folder holds the same as files (include system.settings in a LAMMPS input).  Minimisation, MD and analysis are not
    soapff's: they are written in LAMMPS."""
    lib = _library.load(library)
    if out_dir is None:
        stem = os.path.splitext(os.path.basename(os.fspath(structure)))[0] if isinstance(structure, (str, os.PathLike)) else "structure"
        out_dir = f"{stem}_soapff"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    if lib.engine == "clayff":
        from . import clayff
        kw = dict(boundary=boundary)
        if cutoff is not None:
            kw["cutoff"] = cutoff
        if kspace is not None:
            kw["kspace"] = kspace
        kw.update(options)
        src = os.fspath(structure) if isinstance(structure, (str, os.PathLike)) else structure
        if isinstance(src, str):
            from soapff_cui.core import read_structure as _read
            try:
                _read(src)                         # a file soapff-cui reads itself (with a cell): keep its path in report.json
            except Exception:
                src = read_structure(src)
        rep = clayff.make(lib, src, out_dir, log=log, strict=strict, **kw)
        rep["library_name"] = lib.name
        rep["library_file"] = lib.path
        rep["software"] = _library.software_versions()
        json.dump(rep, open(os.path.join(out_dir, "report.json"), "w"), indent=1)
        return System(out_dir, rep)
    if options:
        raise TypeError(f"unknown options for the molecular engine: {sorted(options)}")
    return _build_molecular(lib, structure, out_dir, cutoff=cutoff, kspace=kspace, boundary=boundary, strict=strict,
                            padding=padding, log=log)


def _build_molecular(lib, structure, out_dir, cutoff, kspace, boundary, strict, padding, log):
    from . import lammps as W
    src = os.fspath(structure) if isinstance(structure, (str, os.PathLike)) else "Atoms"
    at = read_structure(structure) if isinstance(structure, (str, os.PathLike)) else structure.copy()
    warn = []
    if not _has_cell(at):
        at = boxed(at, padding)
        warn.append(f"the structure has no cell: put in an orthogonal box {padding:g} A larger than the atoms on every side "
                    f"({' x '.join(f'{x:.1f}' for x in at.cell.lengths())} A, periodic)")
        log("NOTE " + warn[-1])
    else:
        at.pbc = True
    table_file = lib.table_path(log=log)
    table = M.load_table(table_file)
    res = M.type_atoms(lib, table, at, log=log)
    nb = res["nb"]
    if res["water"] and res["water"]["bonds_taken_out"]:
        log("NOTE the water typed apart from the rest: bonds taken out between water molecules and the rest (number, "
            "shortest A): " + ", ".join(f"{k} {n} ({d:.2f})" for k, (n, d) in res["water"]["bonds_taken_out"].items()))
    el = np.array(at.get_chemical_symbols())
    types, q = res["type"], res["q"]
    thr = float(lib.soap.get("threshold", 0.9))
    known = np.array([t != "" for t in types])
    low1 = known & (res["source"] == "soap") & (res["sim1"] < thr)
    odd = np.zeros(len(at), bool)
    odd[res["odd"]] = True
    bad = ~known | low1 | odd
    if bad.any():
        idx = np.nonzero(bad)[0]
        species = set(str(x) for x in table["allspecies"])

        def why(i):
            if not known[i]:
                if el[i] not in species:
                    return "element not in the library"
                return f"bonded to {M.coordination(nb, el, i)}: no atom of the library has these bonds"
            if odd[i]:
                return (f"bonded to {M.coordination(nb, el, i)}, but its type {types[i]} has "
                        f"{' or '.join(res['want'][types[i]])} (similarity {res['sim1'][i]:.3f})")
            return f"SOAP similarity {res['sim1'][i]:.3f} < {thr} (nearest {res['label'][i]})"
        detail = "; ".join(f"{el[i]}{i + 1}: {why(i)}" for i in idx[:8]) + (" ..." if len(idx) > 8 else "")
        msg = (f"{len(idx)} atoms are not in the library {lib.name}: {detail}.  Add a reference structure with these "
               f"environments to the library (soapff.add_reference), or check the structure (missing hydrogens, overlapping "
               f"molecules, close contacts)")
        if strict:
            raise M.TypingError(msg)
        warn.append(msg)
        log("WARNING " + msg)
        if not known.all():
            raise M.TypingError("atoms without a type cannot be written: " + msg)
    tlist = [str(t) for t in types]
    topo = M.topology(nb, tlist, el, lib)
    miss = {k: v for k, v in topo["missing"].items() if not k.startswith("dihedral")}
    if miss:
        raise M.TypingError(f"no parameters in {lib.path} for: " + ", ".join(f"{k} (x{v})" for k, v in sorted(miss.items())[:12]))
    nodih = {k: v for k, v in topo["missing"].items() if k.startswith("dihedral")}
    if nodih:
        warn.append(f"{sum(nodih.values())} dihedrals have no parameters and are left out (as AMBER / OpenMM do): "
                    + ", ".join(sorted(nodih)[:6]))
    # margin: how much nearer the reference row that gave the type and charge was than the nearest row with another type or
    # charge (atoms typed by SOAP only); a small margin marks an atom to check (e.g. an environment in no reference)
    margin = np.where(res["source"] == "soap", res["margin"], np.nan)
    # molecules and their charges.  The charges of a force field whose charges belong to the types (OPLS-AA, water models,
    # ions) add up to a whole number in every molecule; a molecule whose charge does not has a mistyped atom, even when every
    # similarity is above the threshold (e.g. a junction of two blocks in no reference): stop before writing anything.
    # A molecule held by bonds of step 1 that are not bonds of the force field (a mineral layer of ClayFF, whose charges of
    # Table 1 do not add up to a whole number) is not checked; its charge is in the report
    mol = M.molecules(nb)
    qm = np.bincount(mol, weights=q)
    frac = np.abs(qm - np.round(qm))
    framework = np.unique([mol[i] for i, _ in topo["not_ff"]]).astype(int)
    frac[framework] = 0.0
    tol = float(lib.soap.get("charge_tolerance", 1e-3))
    odd = np.nonzero(frac > tol)[0]
    if len(odd):
        def suspects(k):
            idx = np.nonzero(mol == k)[0]
            idx = idx[np.argsort(np.where(np.isnan(margin[idx]), np.inf, margin[idx]))][:4]
            return ", ".join(f"{el[i]}{i + 1} {tlist[i]} (margin {margin[i]:.4f}, nearest {res['label'][i]})" for i in idx)
        msg = (f"{len(odd)} molecules have a charge that is not a whole number, so atoms in them are mistyped: "
               + "; ".join(f"molecule {k + 1} ({int((mol == k).sum())} atoms) charge {qm[k]:+.4f} e, smallest margins: "
                           f"{suspects(k)}" for k in odd[:3]) + (" ..." if len(odd) > 3 else "")
               + ".  Add a reference structure with these environments to the library (soapff.add_reference)")
        if strict:
            raise M.TypingError(msg)
        warn.append(msg)
        log("WARNING " + msg)
    info = W.write(out_dir, at, nb, tlist, q, topo, lib, cutoff=cutoff, kspace=kspace, boundary=boundary)
    total = float(q.sum())
    if abs(total) > 1e-6:
        warn.append(f"net charge {total:+.4f} e (LAMMPS adds a uniform neutralising background for the long-range Coulomb)")
    with open(os.path.join(out_dir, "charges.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "element", "type", "charge", "molecule", "source", "similarity_soap", "nearest_reference",
                    "margin_soap", "bonded_to"])
        for i in range(len(at)):
            w.writerow([i + 1, el[i], tlist[i], f"{q[i]:+.6f}", mol[i] + 1, res["source"][i], f"{res['sim1'][i]:.5f}",
                        res["label"][i],
                        "" if np.isnan(margin[i]) else f"{margin[i]:.5f}",
                        " ".join(str(b + 1) for b, _ in sorted(nb[i], key=lambda x: x[0]))])
    small = [int(i) for i in np.argsort(np.where(np.isnan(margin), np.inf, margin))[:5] if np.isfinite(margin[i])]
    per = {}
    for t in sorted(set(tlist)):
        m = np.array(tlist) == t
        per[t] = dict(atoms=int(m.sum()), charges=sorted({round(float(x), 6) for x in q[m]}),
                      min_similarity_soap=float(res["sim1"][m].min()))
    kinds = {}
    for k in np.unique(mol):
        key = f"{int((mol == k).sum())} atoms, charge {qm[k]:+.4f}"
        kinds[key] = kinds.get(key, 0) + 1
    rep = dict(structure=src, library=lib.path, library_name=lib.name, engine="molecular", table=table_file,
               software=_library.software_versions(), atoms=len(at),
               charge_sum=total, molecules=int(mol.max() + 1), molecule_kinds=dict(sorted(kinds.items(), key=lambda kv: -kv[1])[:20]),
               min_similarity_soap=float(res["sim1"][known].min()) if known.any() else None,
               min_margin_soap=float(margin[small[0]]) if small else None,
               smallest_margins=[dict(id=i + 1, element=str(el[i]), type=tlist[i], margin=round(float(margin[i]), 5),
                                      similarity_soap=round(float(res["sim1"][i]), 5), nearest_reference=str(res["label"][i]))
                                 for i in small],
               sources={s: int((res["source"] == s).sum()) for s in sorted(set(res["source"]))},
               bond_lengths={k: v["r"] for k, v in table["bond_lengths"].items()}, water_apart=res["water"],
               terms=dict(bonds=info["nbonds"], angles=info["nangles"], dihedrals=info["ndihedrals"], impropers=info["nimpropers"]),
               dihedrals_without_parameters=nodih, box=info["box"], boundary=boundary, kspace=info["kspace"],
               per_type=per, warnings=warn, type_ids=info["type_ids"], type_elements=info["elements"],
               rigid=dict(bonds=info["shake_bonds"], angles=info["shake_angles"]),
               timestep_fs=float(lib.lammps.get("timestep", W.TIMESTEP_FS)))
    if topo["not_ff"]:
        # the bonds of step 1 that are not bonds of the force field (a metal of ClayFF and its O), and the charges of the
        # molecules they hold (not checked)
        rep["bonds_not_in_force_field"] = dict(collections.Counter(
            M.pair_key(tlist[i], tlist[j]) for i, j in topo["not_ff"]).most_common())
        rep["frameworks"] = [dict(molecule=int(k) + 1, atoms=int((mol == k).sum()), charge=round(float(qm[k]), 6))
                             for k in framework]
    json.dump(rep, open(os.path.join(out_dir, "report.json"), "w"), indent=1)
    log(f"{len(at)} atoms, {rep['molecules']} molecules; {rep['terms']}; net charge {total:+.4f} e; wrote {out_dir}/system.data, "
        f"system.settings, charges.csv, report.json")
    typed = at.copy()
    typed.set_array("fftype", np.array(tlist))
    typed.set_array("ffcharge", q)
    return System(out_dir, rep, typed)


# ------------------------------------------------------------------------------------------------- user libraries
def user_library_dir():
    d = os.path.join(os.path.expanduser("~"), ".soapff", "libraries")
    os.makedirs(d, exist_ok=True)
    return d


def new_library(name, extends=None, engine="molecular", folder=None, description=""):
    """write a library file NAME.txt (in ~/.soapff/libraries unless `folder`): extends=BASE starts from a library (its
    references and parameters), to which references and parameter rows are added; returns its path"""
    folder = folder or user_library_dir()
    path = os.path.join(folder, name + ".txt")
    if os.path.exists(path):
        raise FileExistsError(path)
    if extends:
        engine = _library.load(extends).engine
    L = [f"# soapff library {name}", "[library]", f"name        = {name}"]
    L += [f"extends     = {extends}"] if extends else [f"engine      = {engine}"]
    L += [f"description = {description or ('references added to ' + extends if extends else '')}", ""]
    if not extends and engine == "molecular":
        L += ["[soap]", "r_cut = 5.0", "n_max = 4", "l_max = 3", "sigma = 0.4", "depth = 4", "prune = 0.999", "threshold = 0.9", "",
              "[bond_lengths]", "# E1 E2 r(A): two atoms nearer than r are bonded.  r must lie between the longest bond and the "
              "shortest distance not bonded of the pair in the references; a pair not written gets the halfway value (the "
              "log of the first build gives it as a row to copy here)",
              "# H...O of hydrogen bonds comes to 1.43 A in MD, O-H bonds stretch to 1.12 A",
              "H O 1.25", "",
              "[lammps]", "pair_style     = lj/cut/coul/long 10.0", "pair_modify    = mix arithmetic",
              "kspace_style   = pppm 1.0e-5", "bond_style     = harmonic", "angle_style    = harmonic",
              "dihedral_style = fourier", "improper_style = cvff", "special_bonds  = lj/coul 0.0 0.0 0.5", "", "[types]",
              "# type  element  mass  epsilon(kcal/mol)  sigma(A)", "", "[bonds]", "# t1 t2  K(kcal/mol/A^2)  r0(A)", "",
              "[angles]", "# t1 t2 t3  K(kcal/mol/rad^2)  theta0(deg)", "", "[dihedrals]",
              "# t1 t2 t3 t4  K1 n1 d1 [K2 n2 d2 ...]   (E = sum K [1 + cos(n phi - d)])", "", "[impropers]",
              "# t1 t2 tc t4  K d n   (tc central; E = K [1 + d cos(n phi)])", ""]
    L += ["[references]", "# file [prune=0.999]   (relative to this file; extxyz with fftype, ffcharge and bonds, or LAMMPS data)",
          ""]
    open(path, "w").write("\n".join(L))
    return path


def add_reference(library, structure, types=None, charges=None, bonds=None, name=None, prune=None):
    """add a reference structure to a library file (never to a bundled one: make your own with
    soapff.new_library(name, extends=...)).  structure: a file that already carries the types, charges and bonds (extxyz with
    fftype / ffcharge / bonds, or a LAMMPS data file with type names and Bonds), or an ase.Atoms / any structure file with
    `types`, `charges` and `bonds` (the bonded pairs, 0-based atom indices) given.
    The file is copied (or written) next to the library, in <library>_refs/, and listed in [references].  Returns its path."""
    path = _library.find(os.fspath(library) if not isinstance(library, _library.Library) else library.path)
    if os.path.commonpath([path, _library.BUNDLED]) == _library.BUNDLED:
        raise PermissionError(f"{path} is bundled with soapff: make your own library with "
                              f"soapff.new_library('mine', extends='{os.path.basename(path)[:-4]}') and add the reference to it")
    lib_dir = os.path.join(os.path.dirname(path), os.path.basename(path)[:-4] + "_refs")
    os.makedirs(lib_dir, exist_ok=True)
    if types is None:
        src = os.fspath(structure)
        at0 = M.read_reference(src)[0]                    # it must carry the types, the charges
        M.reference_bonds(at0, src)                       # and the bonds
        name = name or os.path.basename(src)
        dest = os.path.join(lib_dir, name)
        if os.path.abspath(src) != os.path.abspath(dest):
            import shutil
            shutil.copyfile(src, dest)
    else:
        at = read_structure(structure) if isinstance(structure, (str, os.PathLike)) else structure
        if charges is None or len(types) != len(at) or len(charges) != len(at):
            raise ValueError("give one type and one charge per atom")
        if bonds is None and "bonds" not in at.arrays:
            raise ValueError("give the bonds of the reference: bonds=[(i, j), ...] (0-based atom indices)")
        name = name or ((os.path.splitext(os.path.basename(os.fspath(structure)))[0] if isinstance(structure, (str, os.PathLike))
                         else f"ref{len(os.listdir(lib_dir)) + 1}") + ".extxyz")
        dest = os.path.join(lib_dir, name)
        M.write_reference(dest, at, types, charges, bonds=bonds)
    rel = os.path.relpath(dest, os.path.dirname(path))
    with open(path, "a") as f:
        txt = open(path).read()
        if "[references]" not in txt:
            f.write("\n[references]\n")
        elif not txt.endswith("\n"):
            f.write("\n")
        f.write(rel + (f"  prune={prune}" if prune else "") + "\n")
    _library.load(path)                                    # the library still reads
    return dest
