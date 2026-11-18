"""The two steps of the molecular engine (soapff.molecular) for clay minerals, water and ions in one pass:
  1. bonds from the elements and the distances (soapff.molecular.bonds_by_distance): two atoms are bonded when they are nearer
     than the bond length of their elements, written in the library ([bond_lengths]: halfway between the longest bond and
     the shortest distance of two atoms of these elements that are not bonded in the references, H-O by hand; checked against
     the references by soapff.molecular.bond_lengths).  Ions have no bond.
     The water molecules (an O bonded to two H only) are then taken apart from the mineral (soapff.molecular.water_apart):
     a water O nearer than the bond length to an edge Al or Mg is not bonded to it.
  2. SOAP of the atoms within `depth` bonds only (the atoms beyond, the water and ions near a surface and the next layer
     included, are left out): the row of the nearest reference atom gives the ClayFF type and charge.
The bonds of the references: the O-H bonds of their LAMMPS data files and, as ClayFF has no metal-O bond, each framework
cation (Si Al Mg Fe Ti) bonded to the O within RMO of it (in the 11 crystals and the edge particle every such cation has 4,
5 or 6 O within 2.19 A and the next O beyond 3.14 A); ions have no bond.  Only the references use RMO."""
import collections
import json
import os

import numpy as np
from ase.neighborlist import neighbor_list

from soapff import molecular as M

FRAMEWORK = ("Si", "Al", "Mg", "Fe", "Ti")
RMO = 2.5                   # A: framework cation - O bonds of the references
DEPTH = 4                   # step 2: the atoms within DEPTH bonds
CHUNK = 20000               # step 2: atoms whose SOAP is computed and compared at a time (memory)
VERSION = 2


def coordination(el, nb, i):
    """the bonded atoms of atom i, counted by element ('H1 Mg3', 'none')"""
    return M.coordination(nb, el, i)


def _counts(s):
    """'H1 Mg3' -> {'H': 1, 'Mg': 3}"""
    out = {}
    for x in s.split():
        if x != "none":
            e = x.rstrip("0123456789")
            out[e] = int(x[len(e):])
    return out


def reference_bonds(at, path):
    """the bonds of a reference (LAMMPS data file): its O-H bonds and the framework cation - O pairs within RMO -> neighbours"""
    el = np.array(at.get_chemical_symbols())
    pairs = list(M._lammps_bonds(path))
    I, J = neighbor_list("ij", at, RMO)
    pairs += [(int(i), int(j)) for i, j in zip(I, J) if el[i] in FRAMEWORK and el[j] == "O"]
    return M.neighbours(at, pairs)


def _prune_new(x, el, tn, q, kept, threshold):
    """the rows of one new reference to keep: those that no kept row (dict (element, type|bonds, charge) -> rows) of the same
    element, type, bonds and charge, and no row kept before them in this reference, matches with similarity >= threshold"""
    sel = np.zeros(len(x), bool)
    qr = np.round(q, 6)
    for key in set(zip(el, tn, qr)):
        idx = np.nonzero((el == key[0]) & (tn == key[1]) & (qr == key[2]))[0]
        cand = idx
        if key in kept:
            cand = idx[(x[idx] @ kept[key].T).max(1) < threshold]
        chosen = []
        for i in cand:
            if not chosen or (x[chosen] @ x[i]).max() < threshold:
                chosen.append(i)
        sel[chosen] = True
    return sel


def build_table(refs, out, soap, prune=None, forcefield=None, depth=DEPTH, read=None, bond_lengths=None, log=print):
    """the table of the references (LAMMPS data files with the ClayFF type names in the Masses comments): the bond lengths
    (step 1; bond_lengths: {"Mg-O": r}, [bond_lengths] of the library) and s3_* (step 2: SOAP on the atoms within `depth` bonds, with the bonds
    of each row's atom).  prune: one threshold per reference (None: kept whole): a row is left out when a row kept before it
    (of this or an earlier reference) with the same element, type, charge and bonds is at least that similar."""
    from .core import read_lammps_data
    read = read or read_lammps_data
    thr = list(prune) if isinstance(prune, (list, tuple)) else [prune] * len(refs)
    data = [read(p) for p in refs]
    species = sorted({e for at, _, _ in data for e in at.get_chemical_symbols()})
    soap = {k: soap[k] for k in ("r_cut", "n_max", "l_max", "sigma")}
    D = collections.defaultdict(list)
    kept = {}
    allnb = []
    for path, t_k, (at, q, names) in zip(refs, thr, data):
        el = np.array(at.get_chemical_symbols())
        names = np.array(names)
        q = np.asarray(q, float)
        nb = reference_bonds(at, path)
        allnb.append((at, nb))
        co = np.array([coordination(el, nb, i) for i in range(len(at))])
        key = np.array([f"{t}|{c}" for t, c in zip(names, co)])
        F = M.soap_rows(at, nb, species, soap, depth, np.arange(len(at)))
        sel = np.ones(len(at), bool) if t_k is None else _prune_new(F, el, key, q, kept, t_k)
        for k in set(zip(el[sel], key[sel], np.round(q[sel], 6))):
            m = sel & (el == k[0]) & (key == k[1]) & (np.round(q, 6) == k[2])
            kept[k] = np.vstack([kept[k], F[m]]) if k in kept else F[m]
        D["X"].append(F[sel]); D["type"] += list(names[sel]); D["q"] += list(q[sel]); D["element"] += list(el[sel])
        D["bonds"] += list(co[sel]); D["origin"] += [f"{os.path.basename(path)}:{i + 1}" for i in np.nonzero(sel)[0]]
        log(f"  {os.path.basename(path)}: {len(at)} atoms, {int(sel.sum())} rows" + ("" if t_k is None else f" (pruned at {t_k})"))
    bl = M.bond_lengths(allnb, bond_lengths)
    log("  bond lengths: " + M.bond_lengths_log(bl))
    arrays = {f"s3_{k}": (np.vstack(v) if k == "X" else np.array(v)) for k, v in D.items()}
    np.savez_compressed(out, steps=VERSION, depth=depth, rmo=RMO, species=np.array(species), bond_lengths=json.dumps(bl),
                        source=np.array([os.path.abspath(r) for r in refs]), **soap, **arrays,
                        **({"forcefield": json.dumps(forcefield)} if forcefield else {}))
    log(f"table: {len(bl)} bond lengths, {len(arrays['s3_q'])} rows from {len(refs)} references, species {species} -> {out}")
    return out


def load_table(path):
    d = np.load(path, allow_pickle=False)
    T = {k: d[k] for k in d.files}
    T["soap"] = {k: (float(d[k]) if k in ("r_cut", "sigma") else int(d[k])) for k in ("r_cut", "n_max", "l_max", "sigma")}
    T["species"] = [str(x) for x in d["species"]]
    T["forcefield"] = json.loads(str(d["forcefield"])) if "forcefield" in d.files else None
    T["depth"] = int(d["depth"])
    if "bond_lengths" not in d.files:
        raise ValueError(f"{path}: a table of the former three steps (no bond lengths): build it again")
    T["bond_lengths"] = json.loads(str(d["bond_lengths"]))
    return T


def type_atoms(at, T):
    """the two steps for the structure `at` -> dict of arrays: type, q (of the row), sim (step 2), margin (to the most similar
    row of another type or charge), bonds (the bonds found, by element), want (the bonds of the row's reference atom), origin
    (reference atom of the row), known (element in the table); nb (the bonds) and water (the water molecules taken apart
    and the bonds between them and the rest taken out; None when the table has no water, or an O with two H that is not)"""
    el = np.array(at.get_chemical_symbols())
    sp = T["species"]
    known = np.isin(el, sp)
    nb = M.bonds_by_distance(at, T["bond_lengths"])
    water = None
    if M.water_rule(zip(T["s3_element"], T["s3_bonds"])):      # the water apart: the mineral without it, each molecule alone
        nb, ow, cut = M.water_apart(at, nb)
        water = dict(molecules=len(ow), bonds_taken_out=cut)
    keys = np.array([f"{t}|{q:.6f}" for t, q in zip(T["s3_type"], T["s3_q"])])
    # the SOAP of an atom depends only on its own neighbourhood: computed and compared CHUNK atoms at a time, so that the
    # memory does not grow with the structure (a few 100,000 atoms)
    n = len(at)
    j3, s3, s2 = np.zeros(n, int), np.zeros(n), np.full(n, np.nan)
    for k0 in range(0, n, CHUNK):
        c = np.arange(k0, min(n, k0 + CHUNK))
        F3 = M.soap_rows(at, nb, sp, T["soap"], T["depth"], c)
        j3[c], s3[c], s2[c] = M._nearest(F3, el[c], T["s3_X"], T["s3_element"], keys=keys)
    j3c = np.maximum(j3, 0)
    bonds = np.array([coordination(el, nb, i) for i in range(len(at))])
    pick = lambda a, j, fill="": np.where(known, a[j], fill)     # noqa: E731
    return dict(type=pick(T["s3_type"], j3c), q=np.where(known, T["s3_q"][j3c], np.nan), sim=np.where(known, s3, 0.0),
                margin=np.where(known, s3 - np.nan_to_num(s2, nan=-1.0), np.nan), bonds=bonds,
                want=pick(T["s3_bonds"], j3c), origin=pick(T["s3_origin"], j3c), known=known, nb=nb, water=water)
