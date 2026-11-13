"""soapff from the command line (python -m soapff, or the soapff command):

    soapff list                                        the libraries that can be given by name
    soapff info amber14sb [--table]                    a library: file, engine, references, types (--table: compute its table)
    soapff build amber14sb protein.pdb -o out [--run --np 4]
    soapff run out [--np 4] [--lmp lmp] [--mpirun mpirun]
    soapff new mine --extends clayff                   ~/.soapff/libraries/mine.txt
    soapff add-reference mine annite.data              a reference with its types and charges (extxyz fftype/ffcharge, LAMMPS data)
"""
import argparse
import json
import os
import sys


def _floats(s):
    return float(s) if s is not None else None


def cmd_list(a):
    from . import libraries
    rows = libraries()
    if not rows:
        print("no libraries")
        return 0
    w = max(len(r["name"]) for r in rows)
    for r in rows:
        print(f"{r['name']:<{w}}  {r['engine']:<9}  {r['description']}")
        if a.paths:
            print(f"{'':<{w}}  {'':<9}  {r['path']}")
    return 0


def cmd_info(a):
    from .library import load
    lib = load(a.library)
    chain, p = [], lib.parent
    while p is not None:
        chain.append(p.name)
        p = p.parent
    print(f"name         {lib.name}")
    print(f"file         {lib.path}")
    print(f"engine       {lib.engine}")
    if chain:
        print(f"extends      {' -> '.join(chain)}")
    print(f"description  {lib.description}")
    if lib.library.get("reference"):
        print(f"reference    {lib.library['reference']}")
    print(f"references   {len(lib.references)} structures")
    for r in lib.references[:a.show]:
        extra = " ".join(f"{k}={v}" for k, v in r.items() if k != "path")
        print(f"             {os.path.relpath(r['path'], os.path.dirname(lib.path))}  {extra}".rstrip())
    if len(lib.references) > a.show:
        print(f"             ... ({len(lib.references) - a.show} more; --show N)")
    if lib.engine == "molecular":
        print(f"types        {len(lib.types)}: {' '.join(sorted(lib.types))}")
        print(f"terms        {len(lib.bonds)} bonds, {len(lib.angles)} angles, {len(lib.dihedrals)} dihedrals, "
              f"{len(lib.impropers)} impropers")
    for sec in ("soap", "graph", "bonding", "lammps"):
        d = getattr(lib, sec)
        if d:
            print(f"[{sec}]".ljust(13) + "  ".join(f"{k}={v}" for k, v in d.items()))
    path = lib.table_path(build=a.table)
    print(f"table        {path or 'not computed yet (computed by the first build, or soapff info --table)'}")
    return 0


def cmd_build(a):
    from . import build
    kw = dict(md=not a.no_md, temp=a.temp, relax_ps=a.relax_ps, nvt_ps=a.nvt_ps, npt_ps=a.npt_ps, press=a.press, baro=a.baro,
              dump_ps=a.dump_ps, timestep=a.timestep, cutoff=a.cutoff, kspace=a.kspace, boundary=a.boundary,
              strict=not a.no_strict, padding=a.padding)
    for o in a.option or []:
        if "=" not in o:
            raise SystemExit(f"--option {o!r}: give KEY=VALUE")
        k, v = o.split("=", 1)
        try:
            kw[k] = json.loads(v)
        except ValueError:
            kw[k] = v
    s = build(a.library, a.structure, a.out, **kw)
    print(s)
    if a.run:
        log = s.run(np=a.np, lmp=a.lmp, mpirun=a.mpirun)
        print("LAMMPS log:", log)
    return 0


def cmd_run(a):
    from .lammps import run
    print("LAMMPS log:", run(a.folder, np_=a.np, lmp=a.lmp, mpirun=a.mpirun, input=a.input, echo=not a.quiet))
    return 0


def cmd_new(a):
    from . import new_library
    path = new_library(a.name, extends=a.extends, engine=a.engine, folder=a.folder, description=a.description or "")
    print(path)
    return 0


def cmd_add_reference(a):
    from . import add_reference
    print(add_reference(a.library, a.structure, name=a.name, prune=a.prune))
    return 0


def parser():
    p = argparse.ArgumentParser(prog="soapff", description="the force field of any structure for LAMMPS, from a library of "
                                "reference structures (SOAP)")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("list", help="the libraries that can be given by name")
    s.add_argument("--paths", action="store_true", help="also the file of each library")
    s.set_defaults(func=cmd_list)

    s = sub.add_parser("info", help="a library: file, engine, references, types, parameters")
    s.add_argument("library", help="name or path of a library .txt file")
    s.add_argument("--table", action="store_true", help="compute the SOAP table of the library now if it is not cached")
    s.add_argument("--show", type=int, default=10, help="references listed (default 10)")
    s.set_defaults(func=cmd_info)

    s = sub.add_parser("build", help="type a structure with a library and write the LAMMPS input")
    s.add_argument("library", help="name (soapff list) or path of a library .txt file")
    s.add_argument("structure", help="any structure file ASE reads (cif, xyz, extxyz, pdb, POSCAR, ...) or a LAMMPS data file")
    s.add_argument("-o", "--out", help="output folder (default: <structure>_soapff)")
    s.add_argument("--no-md", action="store_true", help="in.lammps stops after the minimisation")
    s.add_argument("--temp", type=float, default=300.0, help="K (default 300)")
    s.add_argument("--relax-ps", type=float, default=10.0, help="NVT relaxation, not recorded (default 10 ps)")
    s.add_argument("--nvt-ps", type=float, default=50.0, help="NVT, recorded (default 50 ps)")
    s.add_argument("--npt-ps", type=float, default=0.0, help="NPT after the NVT (default 0: none)")
    s.add_argument("--press", type=float, default=1.0, help="atm (default 1)")
    s.add_argument("--baro", default="aniso", choices=["iso", "aniso", "tri", "z"], help="NPT coupling (default aniso)")
    s.add_argument("--dump-ps", type=float, default=1.0, help="dump interval (default 1 ps)")
    s.add_argument("--timestep", type=float, help="fs (default: that of the library)")
    s.add_argument("--cutoff", type=float, help="pair cutoff, A (default: that of the library)")
    s.add_argument("--kspace", help="e.g. 'pppm 1.0e-5', 'ewald 1.0e-6', none (default: that of the library)")
    s.add_argument("--boundary", default="p p p", help="LAMMPS boundary (default 'p p p')")
    s.add_argument("--no-strict", action="store_true", help="write the input even when atoms are not in the library "
                   "(they take the nearest reference; see report.json)")
    s.add_argument("--padding", type=float, default=10.0, help="box padding (A) of a structure without a cell (default 10)")
    s.add_argument("--option", action="append", metavar="KEY=VALUE", help="clayff engine: another option of soapff_cui "
                   "make (repeatable), e.g. --option reserve_water=true")
    s.add_argument("--run", action="store_true", help="then run LAMMPS in the output folder")
    s.add_argument("--np", type=int, default=1, help="MPI processes of --run (default 1)")
    s.add_argument("--lmp", help="LAMMPS executable (default: $SOAPFF_LMP, else found on PATH)")
    s.add_argument("--mpirun", help="MPI launcher (default: $SOAPFF_MPIRUN, else mpirun)")
    s.set_defaults(func=cmd_build)

    s = sub.add_parser("run", help="run LAMMPS in an output folder")
    s.add_argument("folder")
    s.add_argument("--np", type=int, default=1)
    s.add_argument("--lmp")
    s.add_argument("--mpirun")
    s.add_argument("--input", default="in.lammps")
    s.add_argument("--quiet", action="store_true", help="do not echo the LAMMPS output")
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("new", help="write a new library file (in ~/.soapff/libraries unless --folder)")
    s.add_argument("name")
    s.add_argument("--extends", help="start from this library (its references and parameters)")
    s.add_argument("--engine", default="molecular", choices=["molecular", "clayff"])
    s.add_argument("--folder")
    s.add_argument("--description")
    s.set_defaults(func=cmd_new)

    s = sub.add_parser("add-reference", help="add a reference structure (with its types and charges) to a library")
    s.add_argument("library", help="name or path of a library of yours (not a bundled one)")
    s.add_argument("structure", help="extxyz with the per-atom columns fftype and ffcharge, or a LAMMPS data file with "
                   "the type names in the comments of Masses")
    s.add_argument("--name", help="file name in <library>_refs/")
    s.add_argument("--prune", type=float, help="similarity above which a row of this reference adds nothing")
    s.set_defaults(func=cmd_add_reference)
    return p


def main(argv=None):
    a = parser().parse_args(argv)
    from .library import LibraryError
    from .molecular import TypingError
    try:
        return a.func(a)
    except (LibraryError, TypingError, FileNotFoundError, FileExistsError, PermissionError) as e:
        print(f"soapff: {type(e).__name__}: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
