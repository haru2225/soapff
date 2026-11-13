"""soapff.build(library, structure): a structure file of any format in, the LAMMPS input of the force field of the library out."""
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
    """an output folder of soapff: system.data, system.settings, in.lammps, charges.csv, report.json"""

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
                                                              ("input", "in.lammps"), ("charges", "charges.csv"),
                                                              ("report", "report.json"))}

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

    @property
    def input(self):
        return open(self.files["input"]).read()

    def write_input(self, text, name="in.lammps"):
        """replace (or add) a LAMMPS input in the folder, e.g. one written by hand after include system.settings"""
        open(os.path.join(self.out_dir, name), "w").write(text)
        return os.path.join(self.out_dir, name)

    def run(self, np=1, lmp=None, mpirun=None, input="in.lammps", log="log.lammps", echo=True):
        """run LAMMPS here (lmp executable; np > 1 with mpirun) and return the path of log.lammps"""
        from .lammps import run
        return run(self.out_dir, np_=np, lmp=lmp, mpirun=mpirun, input=input, log=log, echo=echo)

    def lammps(self, input="in.lammps", run_input=True, cmdargs=("-screen", "none")):
        """the LAMMPS Python module started in this folder (with `input` run first when run_input)"""
        from .lammps import session
        return session(self.out_dir, input=input, run_input=run_input, cmdargs=cmdargs)

    def thermo(self, log="log.lammps"):
        """the thermo output of the last run in log.lammps: {column: numpy array}"""
        return read_thermo(os.path.join(self.out_dir, log))

    def analyze(self, **kw):
        """the analyses of soapff-gui (soapff_cui.analysis: diffusion from compute msd, SOVA g(r) and S(Q), XRD)"""
        from soapff_cui.analysis import analyze
        return analyze(self.out_dir, **kw)


def read_thermo(path, which=-1):
    """the thermo table of a LAMMPS log (`which`: the run, -1 = the last) -> {column: array}"""
    blocks, cur, head = [], None, None
    for line in open(path):
        p = line.split()
        if p[:1] == ["Step"]:
            head, cur = p, []
            blocks.append((head, cur))
            continue
        if cur is not None:
            if line.startswith("Loop time") or not p:
                cur = None
                continue
            try:
                cur.append([float(x) for x in p])
            except ValueError:
                cur = None
    if not blocks:
        return {}
    head, rows = blocks[which]
    a = np.array([r for r in rows if len(r) == len(head)])
    return {h: a[:, k] for k, h in enumerate(head)} if len(a) else {h: np.zeros(0) for h in head}


# ------------------------------------------------------------------------------------------------- build
def build(library, structure, out_dir=None, *, md=True, temp=300.0, relax_ps=10.0, nvt_ps=50.0, npt_ps=0.0, press=1.0,
          baro="aniso", dump_ps=1.0, timestep=None, cutoff=None, kspace=None, boundary="p p p", strict=True, padding=10.0,
          log=print, **options):
    """type every atom of `structure` with `library` and write the LAMMPS input.

    library    a library name (soapff.libraries() lists them: clayff, clayff_edge, amber14sb, tip3p, spce, ...) or the path of
               a library .txt file
    structure  a structure file of any format ASE reads (cif, xyz, extxyz, pdb, POSCAR, ...) or a LAMMPS data file, or an
               ase.Atoms
    out_dir    the output folder (default: <structure name>_soapff in the current folder)
    md         False: in.lammps stops after the energy minimisation; True: then NVT (relax_ps not recorded, nvt_ps recorded) and
               NPT (npt_ps > 0) at temp (K), press (atm), barostat coupling baro (iso, aniso, tri, z)
    strict     stop (soapff.TypingError) when an atom is not in the library (its environment has no counterpart in the
               references) instead of writing the input with the nearest reference
    padding    molecular engine: a structure without a cell (xyz, pdb without CRYST1) is put in a box this much larger (A)
    options    clayff engine: the other options of soapff_cui.core.make (e.g. reserve_water, water_rule)
    Returns a soapff.System (out_dir, report, files, run(), lammps())."""
    lib = _library.load(library)
    if out_dir is None:
        stem = os.path.splitext(os.path.basename(os.fspath(structure)))[0] if isinstance(structure, (str, os.PathLike)) else "structure"
        out_dir = f"{stem}_soapff"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    if lib.engine == "clayff":
        from . import clayff
        kw = dict(md=md, temp=temp, relax_ps=relax_ps, nvt_ps=nvt_ps, npt_ps=npt_ps, press=press, baro=baro, dump_ps=dump_ps,
                  boundary=boundary)
        if timestep is not None:
            kw["timestep_fs"] = timestep
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
        try:
            rep = clayff.make(lib, src, out_dir, log=log, strict=strict, **kw)
        except Exception as e:
            from soapff_cui.core import FallbackError
            if isinstance(e, FallbackError):
                raise M.TypingError(str(e)) from e
            raise
        rep["library_name"] = lib.name
        json.dump(rep, open(os.path.join(out_dir, "report.json"), "w"), indent=1)
        return System(out_dir, rep)
    if options:
        raise TypeError(f"unknown options for the molecular engine: {sorted(options)}")
    return _build_molecular(lib, structure, out_dir, md=md, temp=temp, relax_ps=relax_ps, nvt_ps=nvt_ps, npt_ps=npt_ps,
                            press=press, baro=baro, dump_ps=dump_ps, timestep=timestep, cutoff=cutoff, kspace=kspace,
                            boundary=boundary, strict=strict, padding=padding, log=log)


def _build_molecular(lib, structure, out_dir, md, temp, relax_ps, nvt_ps, npt_ps, press, baro, dump_ps, timestep, cutoff,
                     kspace, boundary, strict, padding, log):
    from . import lammps as W
    from soapff_cui.core import parse_boundary
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
    _, periodic = parse_boundary(boundary)
    if npt_ps > 0 and not all(periodic):
        raise ValueError(f"NPT needs a periodic box (boundary {boundary})")
    table = M.load_table(lib.table_path(log=log))
    factor, ions = M.bonding_options(lib)
    nb = M.bond_graph(at, factor, ions)
    el = np.array(at.get_chemical_symbols())
    res = M.assign(lib, table, at, nb, log=log)
    types, q = res["type"], res["q"]
    thr = float(lib.soap.get("threshold", 0.9))
    graph = "G" in table
    known = np.array([t != "" for t in types])
    low1 = known & (res["source"] == "soap") & (res["sim1"] < thr)
    low2 = known & (res["source"] == "graph") & (res["sim2"] < M.GRAPH_EXACT)
    bad = ~known | low1 | low2
    if bad.any():
        idx = np.nonzero(bad)[0]
        why = lambda i: ("element not in the library" if not known[i] else
                         f"fragment not in the library (similarity {res['sim2'][i]:.4f}, nearest {res['label'][i]})" if low2[i] else
                         f"SOAP similarity {res['sim1'][i]:.3f} < {thr} (nearest {res['label'][i]})")
        detail = "; ".join(f"{el[i]}{i + 1}: {why(i)}" for i in idx[:8]) + (" ..." if len(idx) > 8 else "")
        # an H bonded to two atoms: two molecules overlap (e.g. water placed across the periodic boundary)
        over = [i for i in range(len(at)) if el[i] == "H" and len(nb[i]) > 1]
        clash = ""
        if over:
            clash = (f"  The structure has overlapping atoms: {len(over)} H atoms are within bonding distance of two atoms ("
                     + "; ".join(f"H{i + 1}: " + ", ".join(f"{el[j]}{j + 1} {np.linalg.norm(d):.2f} A" for j, d in nb[i])
                                 for i in over[:3]) + (" ..." if len(over) > 3 else "") + "): remove or move those molecules.")
        msg = (f"{len(idx)} atoms are not in the library {lib.name}: {detail}.{clash}  Otherwise add a reference structure "
               f"with these environments to the library (soapff.add_reference), or check the structure (missing hydrogens, "
               f"close contacts)")
        if strict:
            raise M.TypingError(msg)
        warn.append(msg)
        log("WARNING " + msg)
        if not known.all():
            raise M.TypingError("atoms of elements the library does not have cannot be written: " + msg)
    tlist = [str(t) for t in types]
    topo = M.topology(nb, tlist, el, lib)
    miss = {k: v for k, v in topo["missing"].items() if not k.startswith("dihedral")}
    if miss:
        raise M.TypingError(f"no parameters in {lib.path} for: " + ", ".join(f"{k} (x{v})" for k, v in sorted(miss.items())[:12]))
    nodih = {k: v for k, v in topo["missing"].items() if k.startswith("dihedral")}
    if nodih:
        warn.append(f"{sum(nodih.values())} dihedrals have no parameters and are left out (as AMBER / OpenMM do): "
                    + ", ".join(sorted(nodih)[:6]))
    info = W.write(out_dir, at, nb, tlist, q, topo, lib, cutoff=cutoff, kspace=kspace, boundary=boundary)
    open(os.path.join(out_dir, "in.lammps"), "w").write(
        W.input_file(lib, info, md=md, temp=temp, relax_ps=relax_ps, nvt_ps=nvt_ps, npt_ps=npt_ps, press=press, baro=baro,
                     dump_ps=dump_ps, timestep_fs=timestep))
    # molecules and their charges
    mol = M.molecules(nb)
    qm = np.bincount(mol, weights=q)
    frac = np.abs(qm - np.round(qm))
    odd = np.nonzero(frac > 1e-3)[0]
    if len(odd):
        warn.append(f"{len(odd)} molecules have a charge that is not a whole number (largest deviation {frac.max():.4f} e)")
        log("WARNING " + warn[-1])
    total = float(q.sum())
    if abs(total) > 1e-6:
        warn.append(f"net charge {total:+.4f} e (LAMMPS adds a uniform neutralising background for the long-range Coulomb)")
    with open(os.path.join(out_dir, "charges.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "element", "type", "charge", "molecule", "source", "similarity_soap", "similarity_fragment",
                    "nearest_reference"])
        for i in range(len(at)):
            w.writerow([i + 1, el[i], tlist[i], f"{q[i]:+.6f}", mol[i] + 1, res["source"][i], f"{res['sim1'][i]:.5f}",
                        "" if np.isnan(res["sim2"][i]) else f"{res['sim2'][i]:.5f}", res["label"][i]])
    per = {}
    for t in sorted(set(tlist)):
        m = np.array(tlist) == t
        per[t] = dict(atoms=int(m.sum()), charges=sorted({round(float(x), 6) for x in q[m]}),
                      min_similarity_soap=float(res["sim1"][m].min()))
    kinds = {}
    for k in np.unique(mol):
        key = f"{int((mol == k).sum())} atoms, charge {qm[k]:+.4f}"
        kinds[key] = kinds.get(key, 0) + 1
    rep = dict(structure=src, library=lib.path, library_name=lib.name, engine="molecular", atoms=len(at),
               charge_sum=total, molecules=int(mol.max() + 1), molecule_kinds=dict(sorted(kinds.items(), key=lambda kv: -kv[1])[:20]),
               min_similarity_soap=float(res["sim1"][known].min()) if known.any() else None,
               min_similarity_fragment=(float(np.nanmin(res["sim2"])) if graph and np.isfinite(res["sim2"]).any() else None),
               sources={s: int((res["source"] == s).sum()) for s in sorted(set(res["source"]))},
               terms=dict(bonds=info["nbonds"], angles=info["nangles"], dihedrals=info["ndihedrals"], impropers=info["nimpropers"]),
               dihedrals_without_parameters=nodih, box=info["box"], boundary=boundary, kspace=info["kspace"],
               per_type=per, warnings=warn,
               md=dict(md=md, temp=temp, relax_ps=relax_ps, nvt_ps=nvt_ps, npt_ps=npt_ps, press=press, baro=baro,
                       dump_ps=dump_ps, timestep_fs=float(timestep or lib.lammps.get("timestep", W.TIMESTEP_FS))))
    json.dump(rep, open(os.path.join(out_dir, "report.json"), "w"), indent=1)
    log(f"{len(at)} atoms, {rep['molecules']} molecules; {rep['terms']}; net charge {total:+.4f} e; wrote {out_dir}/system.data, "
        f"system.settings, in.lammps, charges.csv, report.json")
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
              "[bonding]", "factor = 1.15", "", "[lammps]", "pair_style     = lj/cut/coul/long 10.0", "pair_modify    = mix arithmetic",
              "kspace_style   = pppm 1.0e-5", "bond_style     = harmonic", "angle_style    = harmonic",
              "dihedral_style = fourier", "improper_style = cvff", "special_bonds  = lj/coul 0.0 0.0 0.5", "", "[types]",
              "# type  element  mass  epsilon(kcal/mol)  sigma(A)", "", "[bonds]", "# t1 t2  K(kcal/mol/A^2)  r0(A)", "",
              "[angles]", "# t1 t2 t3  K(kcal/mol/rad^2)  theta0(deg)", "", "[dihedrals]",
              "# t1 t2 t3 t4  K1 n1 d1 [K2 n2 d2 ...]   (E = sum K [1 + cos(n phi - d)])", "", "[impropers]",
              "# t1 t2 tc t4  K d n   (tc central; E = K [1 + d cos(n phi)])", ""]
    L += ["[references]", "# file [prune=0.999]   (relative to this file; extxyz with fftype and ffcharge, or LAMMPS data)", ""]
    open(path, "w").write("\n".join(L))
    return path


def add_reference(library, structure, types=None, charges=None, name=None, prune=None):
    """add a reference structure to a library file (never to a bundled one: make your own with
    soapff.new_library(name, extends=...)).  structure: a file that already carries the types and charges (extxyz with fftype /
    ffcharge, or a LAMMPS data file with type names), or an ase.Atoms / any structure file with `types` and `charges` given.
    The file is copied (or written) next to the library, in <library>_refs/, and listed in [references].  Returns its path."""
    path = _library.find(os.fspath(library) if not isinstance(library, _library.Library) else library.path)
    if os.path.commonpath([path, _library.BUNDLED]) == _library.BUNDLED:
        raise PermissionError(f"{path} is bundled with soapff: make your own library with "
                              f"soapff.new_library('mine', extends='{os.path.basename(path)[:-4]}') and add the reference to it")
    lib_dir = os.path.join(os.path.dirname(path), os.path.basename(path)[:-4] + "_refs")
    os.makedirs(lib_dir, exist_ok=True)
    if types is None:
        src = os.fspath(structure)
        M.read_reference(src)                             # it must carry the types and charges
        name = name or os.path.basename(src)
        dest = os.path.join(lib_dir, name)
        if os.path.abspath(src) != os.path.abspath(dest):
            import shutil
            shutil.copyfile(src, dest)
    else:
        at = read_structure(structure) if isinstance(structure, (str, os.PathLike)) else structure
        if charges is None or len(types) != len(at) or len(charges) != len(at):
            raise ValueError("give one type and one charge per atom")
        name = name or ((os.path.splitext(os.path.basename(os.fspath(structure)))[0] if isinstance(structure, (str, os.PathLike))
                         else f"ref{len(os.listdir(lib_dir)) + 1}") + ".extxyz")
        dest = os.path.join(lib_dir, name)
        M.write_reference(dest, at, types, charges)
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
