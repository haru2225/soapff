"""soapff libraries: a text file (.txt) that names the force field, the reference structures (atoms with a known type and
charge) and the parameters, so that a library can be read, copied and extended by hand.

Format (one file per library; '#' starts a comment):

    [library]                    key = value
    name        = oplsaa
    engine      = molecular      # molecular (any force field written in the tables below) or clayff (soapff_cui, ClayFF)
    extends     = NAME           # optional: start from another library (its references and parameters come first)
    description = ...
    reference   = ...            # the papers of the force field
    [soap]                       key = value: r_cut, n_max, l_max, sigma (dscribe SOAP), depth (molecular: atoms within
                                 # `depth` bonds), prune (similarity above which a reference row adds nothing), threshold,
                                 # charge_tolerance (molecular: a molecule whose charge is further than this from a whole
                                 # number stops the build; default 0.001 e)
    [lammps]                     key = value: pair_style, pair_modify, kspace_style, bond_style, angle_style, dihedral_style,
                                 # improper_style, special_bonds, rigid (types made rigid by SHAKE), timestep (fs)
    [types]        type element mass epsilon sigma [class]       (kcal/mol, A; LJ 12-6, sigma not R_min; class: optional,
                                                                  a name the rows below may use in place of the type, as
                                                                  the classes of OPLS-AA; default the type itself)
    [bonds]        t1 t2 K r0                                    (E = K (r - r0)^2; kcal/mol/A^2, A)
    [angles]       t1 t2 t3 K theta0                             (E = K (theta - theta0)^2; kcal/mol/rad^2, deg)
    [dihedrals]    t1 t2 t3 t4 K1 n1 d1 [K2 n2 d2 ...]           (E = sum K [1 + cos(n phi - d)]; kcal/mol, deg; X = any)
    [impropers]    t1 t2 tc t4 K d n                             (tc: the central atom, bonded to the three others;
                                                                  E = K [1 + d cos(n phi)], LAMMPS cvff; X = any)
    [bond_lengths] E1 E2 r                                       (A: two atoms of these elements nearer than r are bonded
                                                                  (step 1).  r must lie between the longest bond of the pair
                                                                  and the shortest distance of two of its atoms that are not
                                                                  bonded, in the references (else the table is not built); a
                                                                  pair bonded in the references and not written gets the
                                                                  halfway value, and the log gives it as a row to copy here)
    [references]   FILE [prune=S]                                (paths relative to this file: extxyz with the per-atom
                                                                  columns fftype, ffcharge and (molecular) bonds, the
                                                                  1-based ids of the bonded atoms; or a LAMMPS data file
                                                                  with the type names in the comments of Masses and Bonds)

Libraries are found by name in $SOAPFF_LIBRARY_PATH, ~/.soapff/libraries and the libraries bundled with soapff (in this
order), or given as a path.  The table of SOAP vectors of the references is computed once and kept in ~/.cache/soapff
($SOAPFF_CACHE); any change of the references or of [soap] makes a new one."""
import hashlib
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLED = os.path.join(HERE, "libraries")
# changes when the way the table of a library is computed changes (per engine: the prebuilt ClayFF tables stay valid)
DESCRIPTOR_VERSION = {"clayff": "3", "molecular": "5"}  # molecular 5: the prune= of each reference
KEYVALUE = ("library", "soap", "lammps")
TABLES = ("types", "bonds", "angles", "dihedrals", "impropers", "bond_lengths", "references")


class LibraryError(ValueError):
    """a library file that cannot be read or a library that cannot be found"""


def search_path():
    """the folders where libraries are looked for by name"""
    env = [p for p in os.environ.get("SOAPFF_LIBRARY_PATH", "").split(os.pathsep) if p]
    return env + [os.path.join(os.path.expanduser("~"), ".soapff", "libraries"), BUNDLED]


def cache_dir():
    d = os.environ.get("SOAPFF_CACHE") or os.path.join(os.path.expanduser("~"), ".cache", "soapff")
    os.makedirs(d, exist_ok=True)
    return d


def find(name):
    """path of a library given by path or by name"""
    if os.path.isfile(name):
        return os.path.abspath(name)
    base = name[:-4] if name.endswith(".txt") else name
    for d in search_path():
        p = os.path.join(d, base + ".txt")
        if os.path.isfile(p):
            return os.path.abspath(p)
    raise LibraryError(f"library {name!r} not found (looked in {', '.join(search_path())}); available: "
                       f"{', '.join(l['name'] for l in available()) or 'none'}")


def available():
    """[{name, path, engine, description}] of every library that can be given by name (the first of a name wins)"""
    out, seen = [], set()
    for d in search_path():
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if not f.endswith(".txt") or f[:-4] in seen:
                continue
            try:
                raw = parse(os.path.join(d, f))
            except LibraryError:
                continue
            lib = raw["library"]
            if lib.get("hidden", "no") == "yes":
                continue
            seen.add(f[:-4])
            out.append(dict(name=f[:-4], path=os.path.join(d, f), engine=lib.get("engine", "molecular"),
                            description=lib.get("description", "")))
    return out


def parse(path):
    """one library file -> {section: {key: value}} for the key = value sections and {section: [(tokens, file, line)]} for
    the tables (no `extends`)"""
    out = {s: {} for s in KEYVALUE}
    out.update({s: [] for s in TABLES})
    sec = None
    for k, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        m = re.fullmatch(r"\[\s*([A-Za-z_]+)\s*\]", line)
        if m:
            sec = m.group(1).lower()
            if sec == "bonding":
                raise LibraryError(f"{path}:{k}: [bonding] is no longer used: the bond lengths are rows "
                                   f"[bond_lengths] E1 E2 r")
            if sec not in out:
                raise LibraryError(f"{path}:{k}: unknown section [{sec}] (sections: {', '.join(KEYVALUE + TABLES)})")
            continue
        if sec is None:
            raise LibraryError(f"{path}:{k}: a line before the first [section]")
        if sec in KEYVALUE:
            if "=" not in line:
                raise LibraryError(f"{path}:{k}: [{sec}] takes key = value lines, got {line!r}")
            key, val = (x.strip() for x in line.split("=", 1))
            out[sec][key.lower()] = val
        else:
            out[sec].append((line.split(), path, k))
    return out


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


class Library:
    """a library read from its text file (with the libraries it extends)"""

    def __init__(self, path, _chain=()):
        self.path = os.path.abspath(path)
        if self.path in _chain:
            raise LibraryError(f"{path}: 'extends' makes a loop ({' -> '.join(_chain + (self.path,))})")
        raw = parse(self.path)
        here = os.path.dirname(self.path)
        kv = raw["library"]
        # the keys that name files are made absolute here (relative to the file that declares them)
        if "forcefield" in kv:
            kv["forcefield"] = " ".join(_resolve(x, here) for x in kv["forcefield"].split())
        if "prebuilt" in kv:
            kv["prebuilt"] = _resolve(kv["prebuilt"], here)
        parent = raw["library"].get("extends")
        if parent:
            base = Library(find(parent) if not os.path.isfile(os.path.join(os.path.dirname(self.path), parent))
                           else os.path.join(os.path.dirname(self.path), parent), _chain + (self.path,))
            merged = {s: dict(getattr(base, s)) for s in KEYVALUE}
            for s in KEYVALUE:
                merged[s].update(raw[s])
            merged["library"].pop("extends", None)
            merged["library"].pop("prebuilt", None)            # a prebuilt table is that of the parent only
            merged["library"].pop("prebuilt_sha", None)
            merged["library"].pop("prebuilt_dscribe", None)
            merged["library"].pop("hidden", None)
            # rows of this file first (they win over the parent's for the same types), then the parent's
            tables = {s: raw[s] + base._rows[s] for s in TABLES if s != "references"}
            tables["references"] = base._rows["references"] + raw["references"]
            self.parent = base
        else:
            merged = {s: raw[s] for s in KEYVALUE}
            tables = {s: raw[s] for s in TABLES}
            self.parent = None
        self._rows = tables
        for s in KEYVALUE:
            setattr(self, s, merged[s])
        self.name = self.library.get("name") or os.path.basename(self.path)[:-4]
        self.engine = self.library.get("engine", "molecular")
        if self.engine not in ("molecular", "clayff"):
            raise LibraryError(f"{self.path}: engine must be molecular or clayff, got {self.engine!r}")
        self.description = self.library.get("description", "")
        self.references = self._references()
        self.bond_lengths = self._bond_lengths()
        if self.engine == "molecular":
            self.types = self._types()
            self.bonds = self._table("bonds", 2, 2)
            self.angles = self._table("angles", 3, 2)
            self.dihedrals = self._dihedrals()
            self.impropers = self._table("impropers", 4, 3)

    def __repr__(self):
        return f"<soapff Library {self.name} ({self.engine}, {len(self.references)} references) {self.path}>"

    # ------------------------------------------------------------------------------------------- tables
    def _references(self):
        out = []
        for tok, src, k in self._rows["references"]:
            path = tok[0] if os.path.isabs(tok[0]) else os.path.join(os.path.dirname(src), tok[0])
            if not os.path.isfile(path):
                raise LibraryError(f"{src}:{k}: reference {tok[0]} not found ({path})")
            opts = {}
            for o in tok[1:]:
                if "=" not in o:
                    raise LibraryError(f"{src}:{k}: reference option {o!r} is not key=value")
                key, val = o.split("=", 1)
                opts[key] = float(val) if key == "prune" else val
            out.append(dict(path=os.path.abspath(path), **opts))
        return out

    def _bond_lengths(self):
        """[bond_lengths] E1 E2 r (A): the bond lengths of the library (the first row of a pair wins: this file before the
        parent)"""
        out = {}
        for tok, src, k in self._rows["bond_lengths"]:
            if len(tok) != 3:
                raise LibraryError(f"{src}:{k}: [bond_lengths] rows are: element element length(A)")
            out.setdefault("-".join(sorted(tok[:2])), float(tok[2]))
        return out

    def _types(self):
        out = {}
        for tok, src, k in reversed(self._rows["types"]):          # the first row of a type wins (this file before the parent)
            if len(tok) not in (5, 6):
                raise LibraryError(f"{src}:{k}: [types] rows are: type element mass epsilon sigma [class]")
            out[tok[0]] = dict(element=tok[1], mass=float(tok[2]), epsilon=float(tok[3]), sigma=float(tok[4]),
                               cls=tok[5] if len(tok) == 6 else tok[0])
        return out

    def _table(self, sec, n, nval):
        out = []
        seen = set()
        for tok, src, k in self._rows[sec]:
            if len(tok) != n + nval:
                raise LibraryError(f"{src}:{k}: [{sec}] rows are {n} types and {nval} numbers, got {' '.join(tok)!r}")
            key = tuple(tok[:n])
            if key in seen:
                continue
            seen.add(key)
            out.append((key, tuple(float(x) for x in tok[n:])))
        return out

    def _dihedrals(self):
        out, seen = [], set()
        for tok, src, k in self._rows["dihedrals"]:
            if len(tok) < 7 or (len(tok) - 4) % 3:
                raise LibraryError(f"{src}:{k}: [dihedrals] rows are t1 t2 t3 t4 K1 n1 d1 [K2 n2 d2 ...]")
            key = tuple(tok[:4])
            if key in seen:
                continue
            seen.add(key)
            v = [float(x) for x in tok[4:]]
            out.append((key, [(v[i], int(round(v[i + 1])), v[i + 2]) for i in range(0, len(v), 3)]))
        return out

    # ------------------------------------------------------------------------------------------- the SOAP table
    def soap_params(self):
        s = self.soap
        return dict(r_cut=float(s.get("r_cut", 5.0)), n_max=int(s.get("n_max", 4)), l_max=int(s.get("l_max", 3)),
                    sigma=float(s.get("sigma", 0.4)))

    def recipe(self, software=True):
        """everything the table of the library depends on (its hash names the cached table), with the version of dscribe
        (software=False: without it, for the prebuilt tables, whose dscribe version is checked on its own)"""
        keep = {k: v for k, v in self.soap.items() if k not in ("threshold", "charge_tolerance")}
        ff = [_sha(p) for p in self.library.get("forcefield", "").split()]
        return dict(engine=self.engine, version=DESCRIPTOR_VERSION[self.engine], soap=keep,
                    forcefield=ff,
                    references=[dict(sha=_sha(r["path"]), **{k: v for k, v in r.items() if k != "path"})
                                for r in self.references],
                    types=sorted(self.types) if self.engine == "molecular" else None,
                    bond_lengths=sorted(self.bond_lengths.items()),
                    **(dict(dscribe=dscribe_version()) if software else {}))

    def recipe_sha(self, software=True):
        return hashlib.sha256(json.dumps(self.recipe(software), sort_keys=True).encode()).hexdigest()

    def table_path(self, build=True, log=print):
        """the npz table of the library (the bond lengths and the SOAP rows): the prebuilt one when the
        recipe is unchanged, else the cached one, built when missing"""
        pre = self.library.get("prebuilt")
        if pre and self.library.get("prebuilt_sha", "") == self.recipe_sha(software=False)[:16] and os.path.isfile(pre):
            v = self.library.get("prebuilt_dscribe")
            if v and v != dscribe_version() and self.path not in _WARNED:
                _WARNED.add(self.path)
                log(f"WARNING library {self.name}: the prebuilt table was checked with dscribe {v}, this is dscribe "
                    f"{dscribe_version()}; if its SOAP vectors differ, the similarities differ (delete prebuilt_sha in "
                    f"{self.path} to compute the table again)")
            return pre
        sha = self.recipe_sha()
        out = os.path.join(cache_dir(), f"{self.name}-{sha[:16]}.npz")
        if not os.path.isfile(out):
            if not build:
                return None
            log(f"soapff: library {self.name}: computing the SOAP table of {len(self.references)} references "
                f"(once; kept in {out})")
            tmp = out[:-4] + f".tmp{os.getpid()}.npz"
            if self.engine == "clayff":
                from .clayff import build_table
            else:
                from .molecular import build_table
            build_table(self, tmp, log=log)
            os.replace(tmp, out)
        return out


_WARNED = set()


def dscribe_version():
    """the version of dscribe (the SOAP of the tables depends on it)"""
    try:
        from importlib.metadata import version
        return version("dscribe")
    except Exception:
        return None


def software_versions():
    """the versions of soapff and dscribe, for report.json"""
    from . import __version__
    return dict(soapff=__version__, dscribe=dscribe_version())


def _resolve(spec, base):
    """a path relative to the library file, or pkg:PACKAGE/rel/path inside an installed package"""
    if spec.startswith("pkg:"):
        pkg, rel = spec[4:].split("/", 1)
        import importlib.util
        mod = importlib.util.find_spec(pkg)
        if mod is None or not mod.submodule_search_locations:
            return ""
        return os.path.join(list(mod.submodule_search_locations)[0], rel)
    return spec if os.path.isabs(spec) else os.path.join(base, spec)


def load(name_or_path):
    """a library by name (oplsaa, clayff, ...) or by the path of its .txt file"""
    if isinstance(name_or_path, Library):
        return name_or_path
    return Library(find(str(name_or_path)))
