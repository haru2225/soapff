"""engine "molecular": any force field written in the tables of a library (types, bonds, angles, dihedrals, impropers).

Every atom gets the type and the charge of the most similar atom of the references, in two passes:
  1. SOAP (dscribe) of the covalent neighbourhood of the atom: the atoms within `depth` bonds of it ([soap] depth), the rest of
     the structure left out, so that the type does not depend on the conformation or on the molecules around it.  An atom
     with no bond (an ion) takes the type of the reference atom of its element that has no bond either.
  2. (when the library has [graph]) the molecule is cut into fragments at the bonds between the types of [graph] cut (both
     sides keeping at least two heavy atoms, not in a ring of up to 9 atoms; for a protein: the peptide bonds, so a fragment
     is a residue), and the atom is described by the number of atoms of each type at each bond distance inside its fragment
     (the atoms across a cut bond are counted apart).  This does not depend on the conformation at all and gives the charges
     of force fields whose charges belong to a residue (AMBER): the reference atom with the same fragment has similarity 1.
Bonds come from the covalent radii (x [bonding] factor; the elements of [bonding] ions are never bonded), the angles,
dihedrals and impropers from the bonds, their parameters from the tables (X = any type; the most specific row wins)."""
import collections
import itertools
import json
import os

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.neighborlist import neighbor_list

ION_DEFAULT = ("Li", "Na", "K", "Rb", "Cs", "F", "Cl", "Br", "I", "Mg", "Ca", "Zn", "Sr", "Ba", "Mn", "Fe", "Co", "Ni", "Cu", "Cd")
GRAPH_EXACT = 0.9999                 # pass 2: the fragment of the atom is in the references (same counts)
RING = 8                             # a bond in a ring of up to RING + 1 atoms is never cut


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


def write_reference(path, atoms, types, charges, extra=None):
    """write a reference structure (extxyz with the columns fftype and ffcharge) that a library can list in [references]"""
    import ase.io
    at = atoms.copy()
    at.calc = None
    for k in list(at.arrays):
        if k not in ("numbers", "positions"):
            del at.arrays[k]
    at.set_array("fftype", np.array([str(t) for t in types]))
    at.set_array("ffcharge", np.asarray(charges, float))
    for k, v in (extra or {}).items():
        at.set_array(k, np.asarray(v))
    ase.io.write(path, at, format="extxyz")
    return path


# ------------------------------------------------------------------------------------------------- bonds
def bonding_options(lib):
    b = lib.bonding
    return float(b.get("factor", 1.15)), tuple(b.get("ions", " ".join(ION_DEFAULT)).split())


def bond_graph(at, factor=1.15, ions=ION_DEFAULT):
    """neighbours of every atom: list of (j, vector i->j) for the covalent bonds (distance < factor x (r_i + r_j), the
    covalent radii of ASE); the elements `ions` have no bond"""
    el = at.get_chemical_symbols()
    cut = np.array([0.0 if e in ions else covalent_radii[n] * factor for e, n in zip(el, at.numbers)])
    nb = [[] for _ in range(len(at))]
    if not len(at) or cut.max() == 0:
        return nb
    I, J, D = neighbor_list("ijD", at, cut, self_interaction=False)
    for i, j, d in zip(I, J, D):
        if cut[i] > 0 and cut[j] > 0:
            nb[i].append((int(j), d))
    return nb


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


# ------------------------------------------------------------------------------------------------- pass 1: SOAP
def _soap(species, p):
    from dscribe.descriptors import SOAP
    return SOAP(species=list(species), periodic=False, sparse=False, r_cut=p["r_cut"], n_max=p["n_max"], l_max=p["l_max"],
                sigma=p["sigma"])


def soap_rows(at, nb, species, params, depth, centers, chunk=2000):
    """normalised SOAP vectors of the atoms `centers`, each computed on its covalent neighbourhood (the atoms within `depth`
    bonds, positions relative to it through the bonds; atoms of elements not in `species` left out).  The neighbourhoods of a
    block of atoms are put side by side, far apart, in one structure, so that dscribe computes them at once."""
    soap = _soap(species, params)
    el = np.array(at.get_chemical_symbols())
    ok = np.isin(el, list(species))
    X = np.zeros((len(centers), soap.get_number_of_features()), np.float32)
    for k0 in range(0, len(centers), chunk):
        part = centers[k0:k0 + chunk]
        envs = []
        for c in part:
            rel = {int(c): np.zeros(3)}
            front = [int(c)]
            for _ in range(depth):
                nxt = []
                for a in front:
                    for b, d in nb[a]:
                        if b not in rel:
                            rel[b] = rel[a] + d
                            nxt.append(b)
                front = nxt
            idx = [i for i in rel if ok[i]]
            envs.append((idx, np.array([rel[i] for i in idx])))
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


# ------------------------------------------------------------------------------------------------- pass 2: fragments
def cut_pairs(lib):
    tok = lib.graph.get("cut", "").split()
    if len(tok) % 2:
        raise ValueError(f"{lib.path}: [graph] cut takes pairs of types, got {lib.graph.get('cut')!r}")
    return {frozenset((tok[i], tok[i + 1])) if tok[i] != tok[i + 1] else frozenset((tok[i],)) for i in range(0, len(tok), 2)}


def _side(nb, el, a, b):
    """heavy atoms within RING bonds of a on its side of the bond a-b; None when the bond is in a ring of up to RING + 1 atoms"""
    seen = {a}
    front = [a]
    heavy = int(el[a] != "H")
    for _ in range(RING):
        nxt = []
        for x in front:
            for y, _ in nb[x]:
                if x == a and y == b:
                    continue
                if y == b:
                    return None
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
                    heavy += el[y] != "H"
        front = nxt
    return heavy


def cut_bonds(nb, el, types, pairs):
    """the bonds (both directions) between the types of a cut pair, both sides keeping >= 2 heavy atoms, not in a small ring"""
    cut = set()
    if not pairs:
        return cut
    for a in range(len(nb)):
        for b, _ in nb[a]:
            if a < b and frozenset((types[a], types[b])) in pairs:
                sa, sb = _side(nb, el, a, b), _side(nb, el, b, a)
                if sa is not None and sb is not None and sa >= 2 and sb >= 2:
                    cut.add((a, b))
                    cut.add((b, a))
    return cut


def graph_rows(nb, types, cut, tindex, depth, centers):
    """counts of (type, bond distance <= depth, across a cut bond or not) inside the fragment of each centre (sparse: rows,
    cols); the atoms across a cut bond are counted but not walked through"""
    rows, cols = [], []
    w = depth + 1
    for r, c in enumerate(centers):
        c = int(c)
        dist = {c: 0}
        rows.append(r); cols.append((2 * tindex[types[c]]) * w)
        front = [c]
        for d in range(1, depth + 1):
            nxt = []
            for a in front:
                for b, _ in nb[a]:
                    if b in dist:
                        continue
                    dist[b] = d
                    stub = (a, b) in cut
                    rows.append(r); cols.append((2 * tindex[types[b]] + stub) * w + d)
                    if not stub:
                        nxt.append(b)
            front = nxt
    return np.array(rows, int), np.array(cols, int)


def _dense(rows, cols, n, ncol):
    X = np.zeros((n, ncol), np.float32)
    np.add.at(X, (rows, cols), 1.0)
    return X / np.linalg.norm(X, axis=1, keepdims=True)


# ------------------------------------------------------------------------------------------------- the table of a library
def _prune(X, keys, thr):
    """indices of the rows to keep: a row is dropped when a kept row with the same key is at least thr similar"""
    keep = np.zeros(len(X), bool)
    for key in sorted(set(keys)):
        idx = np.nonzero(keys == key)[0]
        kept = []
        for i in idx:
            if kept:
                s = X[kept] @ X[i]
                if s.max() >= thr:
                    continue
            kept.append(i)
        keep[kept] = True
    return np.nonzero(keep)[0]


def build_table(lib, out, log=print):
    """the SOAP rows (pass 1), the rows of the atoms with no bond, and the fragment rows (pass 2) of every reference atom"""
    factor, ions = bonding_options(lib)
    depth = int(lib.soap.get("depth", 4))
    params = lib.soap_params()
    refs = []
    for r in lib.references:
        at, t, q, lab = read_reference(r["path"])
        bad = sorted(set(t) - set(lib.types))
        if bad:
            raise ValueError(f"{r['path']}: types {bad} are not in [types] of {lib.path}")
        refs.append((r, at, t, q, lab, bond_graph(at, factor, ions)))
    species = sorted({e for _, at, _, _, _, nb in refs for i, e in enumerate(at.get_chemical_symbols()) if nb[i]})
    X, T, Q, E, L = [], [], [], [], []
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
        log(f"  {os.path.basename(r['path'])}: {len(at)} atoms")
    X = np.vstack(X) if X else np.zeros((0, 1), np.float32)
    T, Q, E, L = np.array(T), np.array(Q), np.array(E), np.array(L)
    graph = bool(lib.graph)
    thr = float(lib.soap.get("prune", 0.999))
    keys = np.array([f"{e}|{t}" if graph else f"{e}|{t}|{q:.6f}" for e, t, q in zip(E, T, Q)])
    keep = _prune(X, keys, thr) if len(X) else np.zeros(0, int)
    log(f"  SOAP rows: {len(X)} -> {len(keep)} (pruned at similarity {thr} within the same element, type"
        + ("" if graph else " and charge") + ")")
    data = dict(X=X[keep], type=T[keep], q=Q[keep], element=E[keep], label=L[keep], species=np.array(species),
                soap=json.dumps(dict(params, depth=depth)), iso_element=np.array(list(iso)),
                iso_type=np.array([v[0] for v in iso.values()]), iso_q=np.array([v[1] for v in iso.values()]),
                iso_label=np.array([v[2] for v in iso.values()]), library=lib.path)
    if graph:
        gdepth = int(lib.graph.get("depth", 12))
        pairs = cut_pairs(lib)
        tindex = {t: k for k, t in enumerate(sorted(lib.types))}
        ncol = 2 * len(tindex) * (gdepth + 1)
        G, GT, GQ, GE, GL = [], [], [], [], []
        for r, at, t, q, lab, nb in refs:
            el = np.array(at.get_chemical_symbols())
            cen = np.nonzero([len(x) > 0 for x in nb])[0]
            if not len(cen):
                continue
            rr, cc = graph_rows(nb, list(t), cut_bonds(nb, el, list(t), pairs), tindex, gdepth, cen)
            G.append(_dense(rr, cc, len(cen), ncol))
            GT += list(t[cen]); GQ += list(q[cen]); GE += list(el[cen]); GL += list(lab[cen])
        G = np.vstack(G)
        GT, GQ, GE, GL = np.array(GT), np.array(GQ), np.array(GE), np.array(GL)
        # identical fragment rows with the same element, type and charge: one is enough
        seen, first, owner = set(), [], {}
        clash = set()
        for i in range(len(G)):
            sig = np.round(G[i], 5).tobytes()
            k = (sig, GE[i], GT[i], round(float(GQ[i]), 6))
            if k in seen:
                continue
            seen.add(k)
            first.append(i)
            if (sig, GE[i]) in owner:
                clash.add((sig, GE[i]))        # the same fragment environment with another type or charge
            owner.setdefault((sig, GE[i]), i)
        first = np.array(first)
        # the same fragment environment with two charges: atoms that the bonds cannot tell apart (e.g. CD1 / CD2 of the
        # N-terminal LEU of ff14SB, -0.4106 / -0.4104).  With the same type they get the mean charge (the residue keeps its
        # total charge); with two types the first is kept and reported.
        GQ = GQ.astype(float).copy()
        two_types = []
        for key in clash:
            rows = [i for i in first if (np.round(G[i], 5).tobytes(), GE[i]) == key]
            if len({GT[i] for i in rows}) == 1:
                GQ[rows] = np.mean(GQ[rows])
                log(f"  equivalent atoms, mean charge {np.mean(GQ[rows]):+.6f}: " + ", ".join(f"{GL[i]}" for i in rows))
            else:
                two_types.append(key)
                log("  WARNING the same fragment environment has two types (the first is used): "
                    + ", ".join(f"{GL[i]} {GT[i]}" for i in rows))
        unique = [i for i in first if not any((np.round(G[i], 5).tobytes(), GE[i]) == k and i != owner[k] for k in clash
                                              if k not in two_types)]
        first = np.array(unique)
        log(f"  fragment rows: {len(G)} -> {len(first)}")
        data.update(G=G[first], gtype=GT[first], gq=GQ[first], gelement=GE[first], glabel=GL[first],
                    gindex=json.dumps(tindex), gdepth=gdepth, gcut=json.dumps(sorted(sorted(p) for p in pairs)))
    np.savez_compressed(out, **data)
    log(f"library {lib.name}: {len(keep)} SOAP rows" + (f", {len(data['G'])} fragment rows" if graph else "")
        + f", {len(iso)} elements without bonds; species {species}")
    return out


def load_table(path):
    d = np.load(path, allow_pickle=False)
    t = {k: d[k] for k in d.files}
    t["soap"] = json.loads(str(t["soap"]))
    if "gindex" in t:
        t["gindex"] = json.loads(str(t["gindex"]))
        t["gcut"] = {frozenset(p) for p in json.loads(str(t["gcut"]))}
        t["gdepth"] = int(t["gdepth"])
    return t


def _nearest(F, el, X, E, block=4000):
    """index (into X) of the most similar row of the same element for every row of F, and the similarity"""
    j = -np.ones(len(F), int)
    s = np.zeros(len(F))
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
    return j, s


# ------------------------------------------------------------------------------------------------- typing
def assign(lib, table, at, nb, log=print):
    """type and charge of every atom.  Returns dict(type, q, sim1, sim2, label, source) (sim2 NaN without [graph])"""
    el = np.array(at.get_chemical_symbols())
    n = len(at)
    types = np.array([""] * n, dtype=object)
    q = np.full(n, np.nan)
    sim1 = np.zeros(n)
    sim2 = np.full(n, np.nan)
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
    # pass 1: SOAP of the covalent neighbourhood
    species = list(table["species"])
    cen = np.nonzero(bonded & np.isin(el, species))[0]
    if len(cen):
        p = table["soap"]
        F = soap_rows(at, nb, species, p, int(p["depth"]), cen)
        j, s = _nearest(F, el[cen], table["X"], table["element"])
        ok = j >= 0
        c = cen[ok]
        types[c] = table["type"][j[ok]]
        q[c] = table["q"][j[ok]]
        sim1[c] = s[ok]
        label[c] = table["label"][j[ok]]
        source[c] = "soap"
    # pass 2: the fragment of the atom
    if "G" in table and len(cen):
        tindex = table["gindex"]
        ncol = 2 * len(tindex) * (table["gdepth"] + 1)
        for it in range(2):            # the types of pass 1, then those of the first round (pass 1 may confuse a few types)
            t = [x if x in tindex else None for x in types]
            ok = np.array([x is not None for x in t])
            cc = cen[ok[cen]]
            if not len(cc):
                break
            tt = [x if x is not None else next(iter(tindex)) for x in t]
            rr, kk = graph_rows(nb, tt, cut_bonds(nb, el, tt, table["gcut"]), tindex, table["gdepth"], cc)
            Gx = _dense(rr, kk, len(cc), ncol)
            j, s = _nearest(Gx, el[cc], table["G"], table["gelement"])
            good = j >= 0
            c = cc[good]
            types[c] = table["gtype"][j[good]]
            q[c] = table["gq"][j[good]]
            sim2[c] = s[good]
            label[c] = table["glabel"][j[good]]
            source[c] = "graph"
    return dict(type=types, q=q, sim1=sim1, sim2=sim2, label=label, source=source)


# ------------------------------------------------------------------------------------------------- topology and parameters
def _match(key, rows, wild="X"):
    """the row of `rows` [(types, values)] for the types `key` (or the reverse): exact first, else the one with the fewest
    wildcards, the first in the list among equals"""
    best, nbest = None, 99
    rk = key[::-1]
    for tk, v in rows:
        for k in (key, rk):
            if all(a == b or a == wild for a, b in zip(tk, k)):
                nw = sum(a == wild for a in tk)
                if nw < nbest:
                    best, nbest = (tk, v), nw
                break
        if nbest == 0:
            break
    return best


def topology(nb, types, el, lib):
    """bonds, angles, proper dihedrals and impropers of the bond graph with their parameter rows.  Returns dict of lists of
    (atom indices, parameter key) and the missing terms."""
    n = len(nb)
    nbr = [sorted(b for b, _ in nb[i]) for i in range(n)]
    out = dict(bonds=[], angles=[], dihedrals=[], impropers=[])
    missing = collections.Counter()
    bcache, acache, dcache = {}, {}, {}
    for i in range(n):
        for j in nbr[i]:
            if i < j:
                k = (types[i], types[j])
                if k not in bcache:
                    bcache[k] = _match(k, lib.bonds)
                if bcache[k] is None:
                    missing["bond " + "-".join(k)] += 1
                else:
                    out["bonds"].append(((i, j), bcache[k][0]))
    for j in range(n):
        for i, k in itertools.combinations(nbr[j], 2):
            key = (types[i], types[j], types[k])
            if key not in acache:
                acache[key] = _match(key, lib.angles)
            if acache[key] is None:
                missing["angle " + "-".join(key)] += 1
            else:
                out["angles"].append(((i, j, k), acache[key][0]))
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
                        dcache[key] = _match(key, lib.dihedrals)
                    if dcache[key] is None:
                        missing["dihedral " + "-".join(key)] += 1      # reported, not an error (no term, as in AMBER)
                    else:
                        out["dihedrals"].append(((i, j, k, l), dcache[key][0]))
    # impropers: the order of OpenMM's 'amber' ordering, with the atom index in place of the residue and template index
    for c in range(n):
        if len(nbr[c]) < 3:
            continue
        for trio in itertools.combinations(nbr[c], 3):
            m = _improper(c, trio, types, el, lib.impropers)
            if m:
                out["impropers"].append(m)
    out["missing"] = dict(missing)
    return out


def _improper(c, trio, types, el, rows, wild="X"):
    match = None
    for tk, v in rows:                  # tk = (t1, t2, tc, t4): tc the central atom
        has_w = wild in tk
        if match is not None and has_w:
            continue                   # a definition without wildcards wins over those with
        if not (tk[2] == wild or tk[2] == types[c]):
            continue
        for a2, a3, a4 in itertools.permutations(trio):
            if all(x == wild or x == types[a] for x, a in zip((tk[0], tk[1], tk[3]), (a2, a3, a4))):
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
