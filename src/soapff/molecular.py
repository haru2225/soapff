"""The typing of soapff, for any force field written in the tables of a library (types, bonds, angles, dihedrals, impropers).

Every atom gets its type and charge from the most similar atom of the references (atoms with a known type, charge and bonds),
in two steps:
  1. bonds from the elements and the distances: two atoms are bonded when they are nearer than the bond length of their two
     elements.  The bond lengths are written in the library ([bond_lengths]); each must lie between the longest bond and
     the shortest distance of two atoms of these elements that are not bonded, in the references (else the table is not
     built).  A pair bonded in the references but not written gets the halfway value.  Elements never bonded in the
     references (ions) get no bond.
     When the only O with two H in the references is that of water (water_rule), the water molecules (an O bonded to two H,
     each H to this O only) are then taken apart from the rest: their bonds to other atoms are taken out (water_apart), so
     that the rest is typed without the water and each water molecule on its own bonds.
  2. SOAP of the covalent neighbourhood of the atom: the atoms within `depth` bonds of it ([soap] depth), the rest of the
     structure left out, so that the type does not depend on the conformation or on the molecules around it.  The type and
     charge are those of the most similar row of the same element.  An atom with no bond (an ion) takes the type of the
     reference atom of its element that has no bond either.
Every atom must then have bonds that an atom of its type has in the references (else the build stops).  The angles,
dihedrals and impropers come from the bonds, their parameters from the tables (X = any type; the most specific row wins).
The bonds of step 1 of a type that no row of [bonds] names (the metals and bridging O of ClayFF, which has no metal-O
bond) are not bonds of the force field: they only give the SOAP its neighbourhood (topology)."""
import collections
import itertools
import json
import os

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list


class TypingError(ValueError):
    """atoms whose environment is not in the library, or terms without parameters"""


# ------------------------------------------------------------------------------------------------- references
def read_reference(path):
    """a reference structure -> (Atoms, types, charges, labels).  extxyz (or any file ASE reads with the per-atom arrays
    fftype and ffcharge) or a LAMMPS data file with the type names in the comments of Masses."""
    name = os.path.basename(path)
    if path.endswith((".data", ".lmp", ".lammps")):
        from soapff_cui.core import read_lammps_data
        at, q, types = read_lammps_data(path)
        if types is None:
            raise ValueError(f"{path}: the Masses section must name every type (e.g. '1 12.011  # CT')")
        at.pbc = True
        return at, np.array(types), np.asarray(q, float), np.array([f"{name}:{i + 1}" for i in range(len(at))])
    import ase.io
    at = ase.io.read(path)
    a = at.arrays
    t = a.get("fftype", a.get("type"))
    q = a.get("ffcharge", a.get("q"))
    if t is None or q is None:
        raise ValueError(f"{path}: a reference needs the per-atom columns fftype and ffcharge")
    if "template" in a and "atomname" in a:
        lab = [f"{name}:{x}:{y}" for x, y in zip(a["template"], a["atomname"])]
    else:
        lab = [f"{name}:{i + 1}" for i in range(len(at))]
    return at, np.array([str(x) for x in t]), np.asarray(q, float), np.array(lab)


def write_reference(path, atoms, types, charges, bonds=None, extra=None):
    """write a reference structure (extxyz with the columns fftype, ffcharge and bonds) that a library can list in
    [references].  bonds: the bonded pairs (0-based atom indices); None keeps the column bonds of `atoms`"""
    import ase.io
    at = atoms.copy()
    at.calc = None
    col = bonds_column(len(at), bonds) if bonds is not None else at.arrays.get("bonds")
    if col is None:
        raise ValueError("a reference needs its bonds: give bonds=[(i, j), ...] (0-based atom indices)")
    for k in list(at.arrays):
        if k not in ("numbers", "positions"):
            del at.arrays[k]
    at.set_array("fftype", np.array([str(t) for t in types]))
    at.set_array("ffcharge", np.asarray(charges, float))
    at.set_array("bonds", np.asarray(col))
    for k, v in (extra or {}).items():
        at.set_array(k, np.asarray(v))
    ase.io.write(path, at, format="extxyz")
    return path


# ------------------------------------------------------------------------------------------------- bonds
def bonds_column(n, pairs):
    """the per-atom column bonds of a reference: the 1-based ids of the bonded atoms, comma separated ('-' for none)"""
    nbr = [[] for _ in range(n)]
    for i, j in pairs:
        nbr[int(i)].append(int(j) + 1)
        nbr[int(j)].append(int(i) + 1)
    return np.array([",".join(str(x) for x in sorted(b)) or "-" for b in nbr])


def neighbours(at, pairs):
    """neighbours of every atom, [(j, vector i->j)], for the bonded pairs (the nearest image in a periodic cell)"""
    nb = [[] for _ in range(len(at))]
    pairs = sorted({(min(int(i), int(j)), max(int(i), int(j))) for i, j in pairs})
    if not pairs:
        return nb
    P = np.array(pairs)
    D = at.positions[P[:, 1]] - at.positions[P[:, 0]]
    if at.pbc.any():
        from ase.geometry import find_mic
        D = find_mic(D, at.cell, at.pbc)[0]
    for (i, j), d in zip(P, D):
        nb[i].append((int(j), d))
        nb[j].append((int(i), -d))
    return nb


def _lammps_bonds(path):
    """the bonded pairs (0-based, atoms in the order of their ids) of the Bonds section of a LAMMPS data file"""
    L = open(path).read().splitlines()
    sec = lambda h: next((k for k, l in enumerate(L) if l.split("#")[0].strip() == h), None)      # noqa: E731

    def rows(k):
        i = k + 1
        while i < len(L) and not L[i].strip():
            i += 1
        out = []
        while i < len(L) and L[i].strip():
            out.append(L[i].split("#")[0].split())
            i += 1
        return out
    ids = sorted(int(p[0]) for p in rows(sec("Atoms")))
    rank = {a: k for k, a in enumerate(ids)}
    kb = sec("Bonds")
    return [] if kb is None else [(rank[int(p[2])], rank[int(p[3])]) for p in rows(kb)]


FRAMEWORK = ("Si", "Al", "Mg", "Fe", "Ti")      # the cations of a mineral framework (ClayFF: no metal-O bond in its files)
RMO = 2.5                                       # A: a framework cation of a reference is bonded to the O within RMO


def reference_bonds(at, path):
    """the bonds of a reference: the column bonds of an extxyz, or the Bonds of a LAMMPS data file and, as ClayFF has no
    metal-O bond, each framework cation (Si Al Mg Fe Ti) bonded to the O within RMO of it (in the minerals of clayff every
    such cation has 4, 5 or 6 O within 2.19 A and the next O beyond 3.14 A) -> neighbours"""
    if path.endswith((".data", ".lmp", ".lammps")):
        pairs = list(_lammps_bonds(path))
        el = np.array(at.get_chemical_symbols())
        if np.isin(el, FRAMEWORK).any():
            I, J = neighbor_list("ij", at, RMO)
            pairs += [(int(i), int(j)) for i, j in zip(I, J) if el[i] in FRAMEWORK and el[j] == "O"]
        return neighbours(at, pairs)
    col = at.arrays.get("bonds")
    if col is None:
        raise ValueError(f"{path}: a reference needs its bonds: the per-atom column bonds (the 1-based ids of the bonded "
                         f"atoms, comma separated, '-' for none; soapff.write_reference(..., bonds=...) writes it)")
    pairs = [(i, int(x) - 1) for i, c in enumerate(col) for x in str(c).split(",") if x not in ("-", "")]
    bad = [(i, j) for i, j in pairs if not 0 <= j < len(at) or (j, i) not in set(pairs)]
    if bad:
        raise ValueError(f"{path}: the column bonds is not symmetric or names atoms that do not exist: atoms "
                         + ", ".join(f"{i + 1}-{j + 1}" for i, j in bad[:5]))
    return neighbours(at, pairs)


def coordination(nb, el, i):
    """the bonded atoms of atom i, counted by element ('C2 H2')"""
    c = collections.Counter(el[b] for b, _ in nb[i])
    return " ".join(f"{e}{c[e]}" for e in sorted(c)) or "none"


def pair_key(a, b):
    return "-".join(sorted((str(a), str(b))))


def bond_lengths(refs, override=None, rmax=4.0):
    """the bond length of every pair of elements bonded in the references: the one written in the library ([bond_lengths],
    override: {"C-H": r}), else halfway between the longest bond and the shortest distance between two atoms of these
    elements that are not bonded (within rmax), over all the references.
    -> {"C-H": dict(r, longest_bond, shortest_nonbonded, rule)}; ValueError when a pair has a bond longer than a distance
    between two atoms of the pair that are not bonded (the distance cannot tell them apart), or when a written length does
    not lie between the two (it would not give the bonds of the references: a reference added that the length does not fit)"""
    longest, shortest = {}, {}
    for at, nb in refs:
        el = np.array(at.get_chemical_symbols())
        for i, x in enumerate(nb):
            for j, v in x:
                k = pair_key(el[i], el[j])
                longest[k] = max(longest.get(k, 0.0), float(np.linalg.norm(v)))
    for at, nb in refs:
        el = np.array(at.get_chemical_symbols())
        bonded = {(i, j) for i, x in enumerate(nb) for j, _ in x}
        I, J, d = neighbor_list("ijd", at, rmax, self_interaction=False)
        for i, j, r in zip(I, J, d):
            if i < j and (i, j) not in bonded:
                k = pair_key(el[i], el[j])
                if k in longest and r < shortest.get(k, np.inf):
                    shortest[k] = float(r)
    override = {pair_key(*k.split("-")): float(v) for k, v in (override or {}).items()}
    out, bad, wrong = {}, [], []
    for k in sorted(set(longest) | set(override)):
        lo, hi = longest.get(k), shortest.get(k)
        if lo is not None and hi is not None and lo >= hi:
            bad.append(f"{k}: bond up to {lo:.3f} A, unbonded from {hi:.3f} A")
            continue
        if k in override:
            r, how = override[k], "library file"
            if (lo is not None and r <= lo) or (hi is not None and r >= hi):
                wrong.append(f"{k} {r:g} A (the references: bonds up to {lo:.3f} A"
                             + ("" if hi is None else f", two atoms not bonded from {hi:.3f} A, halfway {0.5 * (lo + hi):.2f} A")
                             + ")")
        elif hi is None:
            r, how = 1.2 * lo, "1.2 x longest bond (no unbonded pair within %g A)" % rmax
        else:
            r, how = 0.5 * (lo + hi), "halfway"
        out[k] = dict(r=round(r, 4), longest_bond=None if lo is None else round(lo, 4),
                      shortest_nonbonded=None if hi is None else round(hi, 4), rule=how)
    if bad:
        raise ValueError("the references have bonds longer than distances between atoms that are not bonded, so the bonds "
                         "cannot come from the distance: " + "; ".join(bad) + ".  Check the bonds of these references")
    if wrong:
        raise ValueError("[bond_lengths] gives lengths that do not separate the bonds of the references from the atoms that "
                         "are not bonded: " + "; ".join(wrong) + ".  Write a length between the two")
    return out


def bond_lengths_log(bl):
    """the bond lengths for the log; the pairs not written in [bond_lengths] also as rows to copy into it"""
    s = ", ".join(f"{k} {v['r']:.3f}" for k, v in bl.items())
    miss = [f"{' '.join(k.split('-'))} {v['r']:.2f}" for k, v in bl.items() if v["rule"] != "library file"]
    return s + (f"; not in [bond_lengths], from the references: {'; '.join(miss)}" if miss else "")


def bonds_by_distance(at, lengths):
    """step 1: two atoms are bonded when they are nearer than the bond length of their elements (lengths: {"C-H": r} or the
    dict of bond_lengths); the nearest image in a periodic cell.  -> neighbours [(j, vector i->j)]"""
    n = len(at)
    nb = [[] for _ in range(n)]
    L = {k: (v["r"] if isinstance(v, dict) else float(v)) for k, v in lengths.items()}
    if not n or not L:
        return nb
    I, J, d, D = neighbor_list("ijdD", at, max(L.values()), self_interaction=False)
    el = np.array(at.get_chemical_symbols())
    sym, ei = np.unique(el, return_inverse=True)            # the bond length of a pair: a table of the elements (memory)
    R = np.array([[L.get(pair_key(a, b), -1.0) for b in sym] for a in sym])
    r = R[ei[I], ei[J]]
    best = {}
    for i, j, dd, v in zip(I[d < r], J[d < r], d[d < r], D[d < r]):
        if (i, j) not in best or dd < best[(i, j)][0]:
            best[(i, j)] = (dd, v)
    for (i, j), (_, v) in sorted(best.items()):
        nb[i].append((int(j), v))
    return nb


def _n_h(bonds):
    """the number of H in a coordination ('Al1 H2' -> 2)"""
    return sum(int(x[1:]) for x in str(bonds).split() if x[:1] == "H" and x[1:].isdigit())


def water_rule(rows):
    """True when the water of a structure is told by its bonds alone: the references have a water molecule (an O bonded to
    two H and to nothing else) and no other O bonded to two H (as in the clay libraries and oplsaa_water: an O of a mineral or
    a polymer has at most one H).  rows: (element, bonds) of the reference atoms, e.g. ('O', 'H2'), ('O', 'Al2 H1')"""
    o2h = {str(b) for e, b in rows if str(e) == "O" and _n_h(b) == 2}
    return o2h == {"H2"}


def water_apart(at, nb):
    """step 1 for the water: the water molecules (an O bonded to exactly two H, each H bonded to this O only) are taken apart
    from the rest of the structure, so that the mineral or the polymer is typed without the water and each water molecule on
    its own bonds.  The bonds between a water molecule and the rest are taken out: an O of water nearer than the bond length
    to a metal (a water on an edge Mg or Al, which ClayFF bonds to nothing).  -> (nb, the O of the water molecules,
    {pair: [number of bonds taken out, shortest]})"""
    el = np.array(at.get_chemical_symbols())
    mol = -np.ones(len(nb), int)
    ow = []
    for o in np.nonzero(el == "O")[0]:
        hs = [b for b, _ in nb[o] if el[b] == "H"]
        if len(hs) == 2 and all(len(nb[h]) == 1 for h in hs):
            mol[[o] + hs] = len(ow)
            ow.append(int(o))
    cut = {}
    out = []
    for i, x in enumerate(nb):
        keep = []
        for j, v in x:
            if mol[i] == mol[j]:
                keep.append((j, v))
            elif i < j:
                k = pair_key(el[i], el[j])
                d = float(np.linalg.norm(v))
                c = cut.setdefault(k, [0, d])
                c[0] += 1
                c[1] = min(c[1], d)
        out.append(keep)
    return out, np.array(ow, int), {k: [n, round(d, 3)] for k, (n, d) in sorted(cut.items())}


def molecules(nb):
    """molecule index of every atom (connected by bonds)"""
    n = len(nb)
    mol = -np.ones(n, int)
    k = 0
    for s in range(n):
        if mol[s] >= 0:
            continue
        mol[s] = k
        stack = [s]
        while stack:
            a = stack.pop()
            for b, _ in nb[a]:
                if mol[b] < 0:
                    mol[b] = k
                    stack.append(b)
        k += 1
    return mol


def unwrapped(at, nb):
    """positions with every molecule whole (each atom moved next to the atom it is bonded to), for the image flags"""
    pos = at.positions.copy()
    seen = np.zeros(len(at), bool)
    for s in range(len(at)):
        if seen[s]:
            continue
        seen[s] = True
        stack = [s]
        while stack:
            a = stack.pop()
            for b, d in nb[a]:
                if not seen[b]:
                    seen[b] = True
                    pos[b] = pos[a] + d
                    stack.append(b)
    return pos


# ------------------------------------------------------------------------------------------------------------- SOAP
def space_rows(at, species, params):
    """not used to type (for comparisons): normalised SOAP vectors of every atom on all the atoms within r_cut (periodic when the structure is); the
    atoms of elements not in `species` are left out (their rows are zero)"""
    from dscribe.descriptors import SOAP
    el = np.array(at.get_chemical_symbols())
    ok = np.isin(el, list(species))
    sub = at if ok.all() else at[ok]
    soap = SOAP(species=list(species), periodic=bool(sub.pbc.any()), sparse=False, r_cut=params["r_cut"],
                n_max=params["n_max"], l_max=params["l_max"], sigma=params["sigma"])
    X = np.zeros((len(at), soap.get_number_of_features()), np.float32)
    if ok.any():
        x = np.atleast_2d(soap.create(sub))
        X[ok] = x / np.linalg.norm(x, axis=1, keepdims=True)
    return X


def _soap(species, p):
    from dscribe.descriptors import SOAP
    return SOAP(species=list(species), periodic=False, sparse=False, r_cut=p["r_cut"], n_max=p["n_max"], l_max=p["l_max"],
                sigma=p["sigma"])


def soap_rows(at, nb, species, params, depth, centers, chunk=2000):
    """normalised SOAP vectors of the atoms `centers`, each computed on its covalent neighbourhood (the atoms within `depth`
    bonds, positions relative to it through the bonds; atoms of elements not in `species` left out).  The neighbourhoods of a
    block of atoms are put side by side, far apart, in one structure, so that dscribe computes them at once.  An atom reached
    through two periodic images (a crystal in a small cell) is in the neighbourhood once per image."""
    soap = _soap(species, params)
    el = np.array(at.get_chemical_symbols())
    ok = np.isin(el, list(species))
    xyz = at.positions
    inv = np.linalg.inv(at.cell.array) if at.pbc.any() and at.cell.rank == 3 else None
    X = np.zeros((len(centers), soap.get_number_of_features()), np.float32)
    for k0 in range(0, len(centers), chunk):
        part = centers[k0:k0 + chunk]
        envs = []
        for c in part:
            c = int(c)
            rel = {(c, (0, 0, 0)): np.zeros(3)}         # (atom, periodic image) -> position relative to the centre
            front = list(rel)
            for _ in range(depth):
                nxt = []
                for ka in front:
                    for b, d in nb[ka[0]]:
                        r = rel[ka] + d
                        img = (0, 0, 0) if inv is None else tuple(np.rint((r - xyz[b] + xyz[c]) @ inv).astype(int))
                        if (b, img) not in rel:
                            rel[(b, img)] = r
                            nxt.append((b, img))
                front = nxt
            keys = [k for k in rel if ok[k[0]]]
            envs.append(([k[0] for k in keys], np.array([rel[k] for k in keys])))
        ext = max(float(np.abs(p).max()) for _, p in envs) + params["r_cut"] + 1.0
        g = int(np.ceil(len(part) ** (1 / 3)))
        syms, pos, cen = [], [], []
        for m, (idx, rel) in enumerate(envs):
            off = 2.0 * ext * np.array([m % g, (m // g) % g, m // (g * g)], float)
            cen.append(len(syms))                     # the centre is the first atom of its neighbourhood
            syms += list(el[idx])
            pos.append(rel + off)
        x = np.atleast_2d(soap.create(Atoms(syms, positions=np.vstack(pos)), centers=cen))
        X[k0:k0 + len(part)] = x / np.linalg.norm(x, axis=1, keepdims=True)
    return X


# ------------------------------------------------------------------------------------------------- the table of a library
def _prune(X, keys, thr):
    """indices of the rows to keep: a row is dropped when a kept row with the same key is at least thr similar
    (thr: one value, or one per row: the prune= of the row's reference, else [soap] prune)"""
    thr = np.broadcast_to(np.asarray(thr, float), (len(X),))
    keep = np.zeros(len(X), bool)
    for key in sorted(set(keys)):
        idx = np.nonzero(keys == key)[0]
        kept = []
        for i in idx:
            if kept:
                s = X[kept] @ X[i]
                if s.max() >= thr[i]:
                    continue
            kept.append(i)
        keep[kept] = True
    return np.nonzero(keep)[0]


def build_table(lib, out, log=print):
    """the bond lengths (step 1), the bonds of every type, the SOAP rows of the covalent neighbourhood (step 2), the rows of
    the atoms with no bond, of every reference atom"""
    depth = int(lib.soap.get("depth", 4))
    params = lib.soap_params()
    refs = []
    for r in lib.references:
        at, t, q, lab = read_reference(r["path"])
        bad = sorted(set(t) - set(lib.types))
        if bad:
            raise ValueError(f"{r['path']}: types {bad} are not in [types] of {lib.path}")
        refs.append((r, at, t, q, lab, reference_bonds(at, r["path"])))
    # the bonds of every type in the references (an atom must have bonds of its type)
    coord = collections.defaultdict(set)
    for r, at, t, q, lab, nb in refs:
        el = at.get_chemical_symbols()
        for i in range(len(at)):
            coord[str(t[i])].add(coordination(nb, el, i))
    coord = {k: sorted(v) for k, v in coord.items()}
    thr = float(lib.soap.get("prune", 0.999))
    # step 1: the bond lengths
    bl = bond_lengths([(at, nb) for _, at, _, _, _, nb in refs], lib.bond_lengths)
    log("  bond lengths: " + bond_lengths_log(bl))
    # step 2: SOAP of the covalent neighbourhood
    species = sorted({e for _, at, _, _, _, nb in refs for i, e in enumerate(at.get_chemical_symbols()) if nb[i]})
    X, T, Q, E, L, P = [], [], [], [], [], []
    iso = {}
    for r, at, t, q, lab, nb in refs:
        el = np.array(at.get_chemical_symbols())
        bonded = np.array([len(x) > 0 for x in nb])
        for i in np.nonzero(~bonded)[0]:
            iso.setdefault(el[i], (t[i], q[i], lab[i]))
        cen = np.nonzero(bonded)[0]
        if len(cen):
            X.append(soap_rows(at, nb, species, params, depth, cen))
            T += list(t[cen]); Q += list(q[cen]); E += list(el[cen]); L += list(lab[cen])
            P += [float(r.get("prune", thr))] * len(cen)             # prune= of this reference, else [soap] prune
        log(f"  {os.path.basename(r['path'])}: {len(at)} atoms")
    X = np.vstack(X) if X else np.zeros((0, 1), np.float32)
    T, Q, E, L = np.array(T), np.array(Q), np.array(E), np.array(L)
    keys = np.array([f"{e}|{t}|{q:.6f}" for e, t, q in zip(E, T, Q)])
    keep = _prune(X, keys, np.array(P)) if len(X) else np.zeros(0, int)
    log(f"  SOAP rows: {len(X)} -> {len(keep)} (pruned at similarity {thr}, or the prune= of a reference, within the same "
        "element, type and charge)")
    data = dict(X=X[keep], type=T[keep], q=Q[keep], element=E[keep], label=L[keep], species=np.array(species),
                soap=json.dumps(dict(params, depth=depth)), iso_element=np.array(list(iso)),
                iso_type=np.array([v[0] for v in iso.values()]), iso_q=np.array([v[1] for v in iso.values()]),
                iso_label=np.array([v[2] for v in iso.values()]), library=lib.path,
                allspecies=np.array(sorted({e for _, at, *_ in refs for e in at.get_chemical_symbols()})),
                coord=json.dumps(coord), bond_lengths=json.dumps(bl))
    np.savez_compressed(out, **data)
    log(f"library {lib.name}: {len(bl)} bond lengths, {len(keep)} SOAP rows, {len(iso)} elements without bonds; species {species}")
    return out


def load_table(path):
    d = np.load(path, allow_pickle=False)
    t = {k: d[k] for k in d.files}
    t["soap"] = json.loads(str(t["soap"]))
    t["coord"] = json.loads(str(t["coord"]))
    t["bond_lengths"] = json.loads(str(t["bond_lengths"]))
    return t


def _nearest(F, el, X, E, block=4000, keys=None):
    """index (into X) of the most similar row of the same element for every row of F, and the similarity; with keys (one
    per row of X), also the similarity of the most similar row with another key (NaN if there is none)"""
    j = -np.ones(len(F), int)
    s = np.zeros(len(F))
    s2 = np.full(len(F), np.nan)
    for e in set(el):
        m = np.nonzero(el == e)[0]
        r = np.nonzero(E == e)[0]
        if not len(r):
            continue
        for k in range(0, len(m), block):
            b = m[k:k + block]
            S = F[b] @ X[r].T
            a = S.argmax(1)
            j[b] = r[a]
            s[b] = S[np.arange(len(b)), a]
            if keys is not None:
                kr = keys[r]
                other = np.where(kr[None, :] == kr[a][:, None], -np.inf, S).max(1)
                s2[b] = np.where(np.isfinite(other), other, np.nan)
    return (j, s, s2) if keys is not None else (j, s)


# ------------------------------------------------------------------------------------------------- typing
def type_atoms(lib, table, at, log=print):
    """the two steps (module docstring): the dict of assign with nb (the bonds of step 1) and odd (the atoms whose bonds are
    not bonds that an atom of their type has in the references), want (the bonds of every type)"""
    nb = bonds_by_distance(at, table["bond_lengths"])
    want = {t: c for t, c in table["coord"].items()}
    elt = dict(zip(table["type"], table["element"])) | dict(zip(table["iso_type"], table["iso_element"]))
    water = None
    if water_rule((elt.get(t, ""), c) for t, cs in want.items() for c in cs):
        nb, ow, cut = water_apart(at, nb)
        water = dict(molecules=len(ow), bonds_taken_out=cut)
    res = assign(lib, table, at, nb, log=log)
    el = at.get_chemical_symbols()
    odd = [i for i in range(len(at)) if res["type"][i] and coordination(nb, el, i) not in want.get(res["type"][i], [])]
    res.update(nb=nb, odd=np.array(odd, int), want=want, water=water)
    return res


def assign(lib, table, at, nb, log=print):
    """type and charge of every atom.  Returns dict(type, q, sim1, margin, label, source).
    margin: sim1 minus the similarity of the nearest SOAP row with another type or charge (NaN if the library has none for
    that element); a small margin means that another answer was almost as near"""
    el = np.array(at.get_chemical_symbols())
    n = len(at)
    types = np.array([""] * n, dtype=object)
    q = np.full(n, np.nan)
    sim1 = np.zeros(n)
    margin = np.full(n, np.nan)
    label = np.array([""] * n, dtype=object)
    source = np.array(["none"] * n, dtype=object)
    bonded = np.array([len(x) > 0 for x in nb])
    # atoms without a bond: the reference atom of the same element without a bond
    iso = {e: (t, x, l) for e, t, x, l in zip(table["iso_element"], table["iso_type"], table["iso_q"], table["iso_label"])}
    for i in np.nonzero(~bonded)[0]:
        if el[i] in iso:
            types[i], q[i], label[i] = iso[el[i]]
            sim1[i] = 1.0
            source[i] = "ion"
    # step 2: SOAP of the covalent neighbourhood
    species = list(table["species"])
    cen = np.nonzero(bonded & np.isin(el, species))[0]
    if len(cen):
        p = table["soap"]
        F = soap_rows(at, nb, species, p, int(p["depth"]), cen)
        keys = np.unique([f"{t}|{x:.6f}" for t, x in zip(table["type"], table["q"])], return_inverse=True)[1]
        j, s, s2 = _nearest(F, el[cen], table["X"], table["element"], keys=keys.ravel())
        ok = j >= 0
        c = cen[ok]
        types[c] = table["type"][j[ok]]
        q[c] = table["q"][j[ok]]
        sim1[c] = s[ok]
        margin[c] = s[ok] - s2[ok]
        label[c] = table["label"][j[ok]]
        source[c] = "soap"
    return dict(type=types, q=q, sim1=sim1, margin=margin, label=label, source=source)


# ------------------------------------------------------------------------------------------------- topology and parameters
def _match(key, rows, wild="X", ckey=None):
    """the row of `rows` [(types, values)] for the types `key` (or the reverse): exact first, else the one with the fewest
    wildcards, the first in the list among equals.  A row names each atom by its type or by its class (`ckey`: the classes
    of the types of `key`, [types] class)."""
    best, nbest = None, 99
    ckey = key if ckey is None else ckey
    both = ((key, ckey), (key[::-1], ckey[::-1]))
    for tk, v in rows:
        for k, c in both:
            if all(a == b or a == cc or a == wild for a, b, cc in zip(tk, k, c)):
                nw = sum(a == wild for a in tk)
                if nw < nbest:
                    best, nbest = (tk, v), nw
                break
        if nbest == 0:
            break
    return best


def water_types(lib, table):
    """the water of the library, for the water types reserved in system.data (soapff.build reserve_water): the O type of the
    water molecules of the references (an O bonded to two H and nothing else), the H type of its H-O-H row of [angles], and
    the rows of the O-H bond and the H-O-H angle -> dict(o, h, bond, angle)"""
    el = {t: v["element"] for t, v in lib.types.items()}
    cls = lambda key: tuple(lib.types[t].get("cls", t) for t in key)         # noqa: E731
    found = []
    for o, cs in table["coord"].items():
        if el.get(o) != "O" or list(cs) != ["H2"]:
            continue
        for h, hc in table["coord"].items():
            if el.get(h) == "H" and list(hc) == ["O1"]:
                b = _match((o, h), lib.bonds, ckey=cls((o, h)))
                a = _match((h, o, h), lib.angles, ckey=cls((h, o, h)))
                if b and a:
                    found.append(dict(o=o, h=h, bond=b[0], angle=a[0]))
    if len(found) != 1:
        raise ValueError(f"reserve_water: the library {lib.name} has {len(found)} water models in its references "
                         f"({', '.join(f['o'] + '/' + f['h'] for f in found) or 'none'}); one is needed")
    return found[0]


def bonding_types(lib, types):
    """the types (of `types`) that have bonds in the force field: named, by type or class, in a row of [bonds].  The others
    (ions, and the metals and bridging O of ClayFF, which has no metal-O bond) have none: their bonds of step 1 only give
    the SOAP its neighbourhood"""
    named = {t for key, _ in lib.bonds for t in key}
    return {t for t in set(types) if t in named or (t in lib.types and lib.types[t].get("cls", t) in named)}


def topology(nb, types, el, lib):
    """bonds, angles, proper dihedrals and impropers of the bond graph with their parameter rows.  A bond of step 1 is a
    bond of the force field when both its types have bonds in the force field (bonding_types; it then needs a row of
    [bonds]); the other bonds of step 1 (a metal of ClayFF and its O) are not written.  An angle on two bonds of the force
    field needs a row of [angles]; an angle on a bond that is not (a metal - O - H of ClayFF) is written when a row gives it.
    Dihedrals and impropers are made of bonds of the force field.  Returns dict of lists of (atom indices, parameter key),
    the missing terms and the bonds of step 1 that are not bonds of the force field (not_ff, (i, j) pairs)."""
    n = len(nb)
    nbr = [sorted(b for b, _ in nb[i]) for i in range(n)]
    has = bonding_types(lib, types)
    ffb = [[b for b in nbr[i] if types[i] in has and types[b] in has] for i in range(n)]
    out = dict(bonds=[], angles=[], dihedrals=[], impropers=[], not_ff=[])
    missing = collections.Counter()
    bcache, acache, dcache = {}, {}, {}
    cls = {t: lib.types[t].get("cls", t) if t in lib.types else t for t in set(types)}
    ck = lambda key: tuple(cls[t] for t in key)                                   # noqa: E731
    for i in range(n):
        for j in nbr[i]:
            if i < j:
                if j not in ffb[i]:
                    out["not_ff"].append((i, j))
                    continue
                k = (types[i], types[j])
                if k not in bcache:
                    bcache[k] = _match(k, lib.bonds, ckey=ck(k))
                if bcache[k] is None:
                    missing["bond " + "-".join(k)] += 1
                else:
                    out["bonds"].append(((i, j), bcache[k][0]))
    for j in range(n):
        for i, k in itertools.combinations(nbr[j], 2):
            key = (types[i], types[j], types[k])
            if key not in acache:
                acache[key] = _match(key, lib.angles, ckey=ck(key))
            if acache[key] is None:
                if i in ffb[j] and k in ffb[j]:
                    missing["angle " + "-".join(key)] += 1
            else:
                out["angles"].append(((i, j, k), acache[key][0]))
    nbr = ffb                                   # dihedrals and impropers: the bonds of the force field
    for j in range(n):
        for k in nbr[j]:
            if j > k:
                continue
            for i in nbr[j]:
                if i == k:
                    continue
                for l in nbr[k]:
                    if l == j or l == i:
                        continue
                    key = (types[i], types[j], types[k], types[l])
                    if key not in dcache:
                        dcache[key] = _match(key, lib.dihedrals, ckey=ck(key))
                    if dcache[key] is None:
                        missing["dihedral " + "-".join(key)] += 1      # reported, not an error (no term, as in AMBER)
                    else:
                        out["dihedrals"].append(((i, j, k, l), dcache[key][0]))
    # impropers: the order of OpenMM's 'amber' ordering, with the atom index in place of the residue and template index
    for c in range(n):
        if len(nbr[c]) < 3:
            continue
        for trio in itertools.combinations(nbr[c], 3):
            m = _improper(c, trio, types, el, lib.impropers, cls=cls)
            if m:
                out["impropers"].append(m)
    out["missing"] = dict(missing)
    return out


def _improper(c, trio, types, el, rows, wild="X", cls=None):
    match = None
    cls = cls or {}
    ok = lambda x, a: x == wild or x == types[a] or x == cls.get(types[a])        # noqa: E731
    for tk, v in rows:                  # tk = (t1, t2, tc, t4): tc the central atom
        has_w = wild in tk
        if match is not None and has_w:
            continue                   # a definition without wildcards wins over those with
        if not ok(tk[2], c):
            continue
        for a2, a3, a4 in itertools.permutations(trio):
            if all(ok(x, a) for x, a in zip((tk[0], tk[1], tk[3]), (a2, a3, a4))):
                if not has_w:
                    if types[a2] == types[a4] and a2 > a4:
                        a2, a4 = a4, a2
                    if types[a3] == types[a4] and a3 > a4:
                        a3, a4 = a4, a3
                    if types[a2] == types[a3] and a2 > a3:
                        a2, a3 = a3, a2
                else:
                    if el[a2] == el[a4] and a2 > a4:
                        a2, a4 = a4, a2
                    if el[a3] == el[a4] and a3 > a4:
                        a3, a4 = a4, a3
                    if a2 > a3:
                        a2, a3 = a3, a2
                match = ((a2, a3, c, a4), tk)
                break
    return match
