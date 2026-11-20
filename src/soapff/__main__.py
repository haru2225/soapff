"""soapff from the command line (python -m soapff, or the soapff command):

    soapff list                                        the libraries that can be given by name
    soapff info oplsaa [--table]                       a library: file, references, types (--table: compute its table)
    soapff build oplsaa cell_PE.extxyz -o out          the force field: out/system.data, out/system.settings (include it in
                                                       a LAMMPS input, or lmp.commands_list(System.commands()) in Python)
    soapff new mine --extends clayff                   ~/.soapff/libraries/mine.txt
    soapff add-reference mine annite.data              a reference with its types and charges (extxyz fftype/ffcharge, LAMMPS data)
"""
import argparse
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
        print(f"{r['name']:<{w}}  {r['reference']}")
        if a.paths:
            print(f"{'':<{w}}  {r['path']}")
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
    if chain:
        print(f"extends      {' -> '.join(chain)}")
    if lib.reference:
        print(f"reference    {lib.reference}")
    print(f"references   {len(lib.references)} structures")
    for r in lib.references[:a.show]:
        extra = " ".join(f"{k}={v}" for k, v in r.items() if k != "path")
        print(f"             {os.path.relpath(r['path'], os.path.dirname(lib.path))}  {extra}".rstrip())
    if len(lib.references) > a.show:
        print(f"             ... ({len(lib.references) - a.show} more; --show N)")
    if lib.types:
        print(f"types        {len(lib.types)}: {' '.join(sorted(lib.types))}")
        print(f"terms        {len(lib.bonds)} bonds, {len(lib.angles)} angles, {len(lib.dihedrals)} dihedrals, "
              f"{len(lib.impropers)} impropers")
    for sec in ("soap", "lammps"):
        d = getattr(lib, sec)
        if d:
            print(f"[{sec}]".ljust(13) + "  ".join(f"{k}={v}" for k, v in d.items()))
    path = lib.table_path(build=a.table)
    print(f"table        {path or 'not computed yet (computed by the first build, or soapff info --table)'}")
    return 0


def cmd_build(a):
    from . import build
    kw = dict(cutoff=a.cutoff, kspace=a.kspace, boundary=a.boundary, strict=not a.no_strict, padding=a.padding,
              reserve_water=a.reserve_water)
    s = build(a.library, a.structure, a.out, **kw)
    print(s)
    print(f"in a LAMMPS input:  include {s.files['settings']}   (run LAMMPS in {s.out_dir}: read_data system.data)")
    return 0


def cmd_new(a):
    from . import new_library
    path = new_library(a.name, extends=a.extends, folder=a.folder, reference=a.reference or "")
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

    s = sub.add_parser("info", help="a library: file, references, types, parameters")
    s.add_argument("library", help="name or path of a library .txt file")
    s.add_argument("--table", action="store_true", help="compute the SOAP table of the library now if it is not cached")
    s.add_argument("--show", type=int, default=10, help="references listed (default 10)")
    s.set_defaults(func=cmd_info)

    s = sub.add_parser("build", help="type a structure with a library and write its force field for LAMMPS (system.data, system.settings)")
    s.add_argument("library", help="name (soapff list) or path of a library .txt file")
    s.add_argument("structure", help="any structure file ASE reads (cif, xyz, extxyz, pdb, POSCAR, ...) or a LAMMPS data file")
    s.add_argument("-o", "--out", help="output folder (default: <structure>_soapff)")
    s.add_argument("--cutoff", type=float, help="pair cutoff, A (default: that of the library)")
    s.add_argument("--kspace", help="e.g. 'pppm 1.0e-5', 'ewald 1.0e-6', none (default: that of the library)")
    s.add_argument("--boundary", default="p p p", help="LAMMPS boundary (default 'p p p')")
    s.add_argument("--no-strict", action="store_true", help="write the force field even when atoms are not in the library "
                   "(they take the nearest reference; see report.json)")
    s.add_argument("--padding", type=float, default=10.0, help="box padding (A) of a structure without a cell (default 10)")
    s.add_argument("--reserve-water", action="store_true", help="also write the water types of the library and room for "
                   "their bonds, so that water can be added in LAMMPS")
    s.set_defaults(func=cmd_build)

    s = sub.add_parser("new", help="write a new library file (in ~/.soapff/libraries unless --folder)")
    s.add_argument("name")
    s.add_argument("--extends", help="start from this library (its references and parameters)")
    s.add_argument("--folder")
    s.add_argument("--reference", help="the papers of the force field")
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
