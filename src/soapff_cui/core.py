"""The engine "clayff" of soapff: reading structures, the two steps with a table of steps.build_table (make_steps), and the
LAMMPS writer of ClayFF (Cygan et al. 2004) with its flexible SPC water."""
import collections
import csv
import itertools
import json
import os

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list

from ._clayff.chem import MASS
from ._clayff.ff import ELEMENT, ff_types, lammps_frame

from .clayff import (ANGLE_HOH, ANGLE_MOH, ANGLE_MOM_EDGE, BOND_OH, CHARGE_SCHEMES, EDGE_HYDROXYL_O, EDGE_METALS, EDGE_PARENT,
                     LJ, WATER_LABELS, WATER_MODELS, WATER_Q,
                     type_charges)
from .clayff import Q as CLAYFF_Q

THRESHOLD = 0.9
TYPE_ALIAS = {"OW": "ow", "HW": "hw", "OH": "oh", "HO": "ho"}   # library type names -> ClayFF type names (clayff.py)
TIMESTEP_FS = 0.5                                                  # recommended time step (report.json): flexible SPC water, 0.5 fs (as in ClayFF studies)
IONS = ("Li", "Na", "K", "Rb", "Cs", "Ca", "Sr", "Ba", "Cl", "Br", "F", "I")
FORMAL_ION = {"Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Ca": 2, "Sr": 2, "Ba": 2, "Cl": -1, "Br": -1, "F": -1, "I": -1}


# ------------------------------------------------------------------------------------------------- reading structures
def read_lammps_data(path):
    """atom_style full (or charge) data file -> (Atoms, charges, type names or None).  Elements from the type names, else the
    masses."""
    L = open(path).read().splitlines()
    box = {}
    tilt = (0.0, 0.0, 0.0)
    vec, origin = {}, np.zeros(3)
    for l in L[:60]:
        p = l.split()
        if len(p) >= 4 and p[2] in ("xlo", "ylo", "zlo"):
            box[p[2][0]] = (float(p[0]), float(p[1]))
        if len(p) >= 6 and p[3] == "xy":
            tilt = tuple(map(float, p[:3]))
        if len(p) >= 4 and p[3] in ("avec", "bvec", "cvec"):              # the general triclinic box (LAMMPS, OVITO)
            vec[p[3][0]] = np.array(p[:3], float)
        if len(p) >= 5 and p[3:5] == ["abc", "origin"]:
            origin = np.array(p[:3], float)
    i = L.index(next(l for l in L if l.strip().startswith("Masses"))) + 2
    mass, name = {}, {}
    while L[i].strip():
        p = L[i].split()
        mass[int(p[0])] = float(p[1])
        name[int(p[0])] = L[i].split("#")[1].strip() if "#" in L[i] else None
        i += 1
    el_of = {}
    for t, m in mass.items():
        if name[t] in ELEMENT:            # a known type name fixes the element
            el_of[t] = ELEMENT[name[t]]
            continue
        el_of[t] = min(MASS, key=lambda e: abs(MASS[e] - m)) if min(abs(MASS[e] - m) for e in MASS) < 0.1 else _mass_to_el(m)
    head = [k for k, l in enumerate(L) if l.strip().startswith("Atoms")][0]
    i = head + 2
    style = L[head].split("#")[1].strip() if "#" in L[head] else ""
    rows = []
    while i < len(L) and L[i].strip():
        p = L[i].split("#")[0].split()
        if style == "charge" or (not style and len(p) in (6, 9)):        # atom_style charge: id type q x y z [ix iy iz]
            rows.append((int(p[0]), int(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[5])))
        else:                                                           # atom_style full: id mol type q x y z [ix iy iz]
            rows.append((int(p[0]), int(p[2]), float(p[3]), float(p[4]), float(p[5]), float(p[6])))
        i += 1
    rows.sort()
    if len(vec) == 3:
        cell = [vec["a"], vec["b"], vec["c"]]
        lo = origin
    else:
        lx = box["x"][1] - box["x"][0]; ly = box["y"][1] - box["y"][0]; lz = box["z"][1] - box["z"][0]
        cell = [[lx, 0, 0], [tilt[0], ly, 0], [tilt[1], tilt[2], lz]]
        lo = np.array([box["x"][0], box["y"][0], box["z"][0]])
    at = Atoms([el_of[r[1]] for r in rows], positions=np.array([r[3:6] for r in rows]) - lo, cell=cell, pbc=True)
    names = [name[r[1]] for r in rows] if all(name.values()) else None
    return at, np.array([r[2] for r in rows]), names


def _section(L, head, raw=False):
    """the lines of a LAMMPS data file section (after its header and any blank lines, up to the next blank line), split
    without their comments (raw: as they are)"""
    i = next(k for k, l in enumerate(L) if l.strip().startswith(head)) + 1
    while i < len(L) and not L[i].strip():
        i += 1
    out = []
    while i < len(L) and L[i].strip():
        out.append(L[i] if raw else L[i].split("#")[0].split())
        i += 1
    return out


def read_forcefield(data, settings, name, source):
    """the force field of a reference given as LAMMPS files, to be used whole instead of ClayFF for its types (build_library
    forcefield=): data with type names in the Masses comments, charges in Atoms, Bond Coeffs (one bond type: the hydroxyl
    O-H), Angle Coeffs and Angles; settings with pair_coeff i i epsilon sigma.  Returns {name, source, charge {type: q},
    lj {type: [epsilon, sigma]}, bond_oh [k, r0], angle_moh [k, theta0], angle_mom_edge [k, theta0]}: the angle type with
    an H at an end is M-O-H, the other one the edge metal-O-metal.  Every atom of a type must carry the same charge."""
    L = open(data).read().splitlines()
    tname = {int(l.split()[0]): l.split("#")[1].strip() for l in _section(L, "Masses", raw=True)}
    at, q, names = read_lammps_data(data)
    charge = {}
    for t, x in zip(names, q):
        if t in charge and abs(charge[t] - x) > 1e-9:
            raise ValueError(f"{data}: type {t} carries the charges {charge[t]} and {x}")
        charge[t] = float(x)
    head = next(l for l in L if l.strip().startswith("Atoms"))
    full = "full" in head or "charge" not in head
    atype = {int(p[0]): tname[int(p[2] if full else p[1])] for p in _section(L, "Atoms")}
    bc = _section(L, "Bond Coeffs")
    if len(bc) != 1:
        raise ValueError(f"{data}: {len(bc)} bond types; only the hydroxyl O-H is expected")
    acoef = {int(p[0]): [float(p[1]), float(p[2])] for p in _section(L, "Angle Coeffs")}
    kind = {}
    for p in _section(L, "Angles"):
        ends = (atype[int(p[2])], atype[int(p[4])])
        kind.setdefault(int(p[1]), set()).add("moh" if any(ELEMENT.get(e) == "H" for e in ends) else "mom_edge")
    if any(len(v) != 1 for v in kind.values()) or sorted(x for v in kind.values() for x in v) != ["moh", "mom_edge"]:
        raise ValueError(f"{data}: the angle types are not one M-O-H and one metal-O-metal type ({kind})")
    ang = {next(iter(v)): acoef[k] for k, v in kind.items()}
    lj = {}
    for l in open(settings).read().splitlines():
        p = l.split("#")[0].split()
        if len(p) >= 5 and p[0] == "pair_coeff" and p[1] == p[2]:
            lj[tname[int(p[1])]] = [float(p[3]), float(p[4])]
    missing = sorted(set(tname.values()) - set(lj))
    if missing:
        raise ValueError(f"{settings}: no pair_coeff for {missing}")
    return dict(name=name, source=source, charge=charge, lj=lj, bond_oh=[float(bc[0][1]), float(bc[0][2])],
                angle_moh=ang["moh"], angle_mom_edge=ang["mom_edge"])


def _mass_to_el(m):
    from ase.data import atomic_masses, chemical_symbols
    k = int(np.argmin(np.abs(atomic_masses[1:] - m))) + 1
    return chemical_symbols[k]


def read_structure(path):
    """any structure with a cell -> Atoms (charges/types from a LAMMPS data file are ignored here)."""
    if path.endswith((".data", ".lammps", ".lmp")) or "data" in os.path.basename(path).split(".")[-1]:
        at, _, _ = read_lammps_data(path)
        return at
    import ase.io
    at = ase.io.read(path)
    if not at.cell.any():
        raise ValueError(f"{path}: no cell; give a periodic structure (extxyz with Lattice, cif, POSCAR or LAMMPS data)")
    at.pbc = True
    return at


# ------------------------------------------------------------------------------------------------- roles and checks
def roles(at):
    """water (an O with two H at 0.9-1.1 A making an H-O-H angle of 95-125 deg, the shape of a (flexible) water molecule, and those
    H; also when it overlaps the mineral in a built structure, or sits on a metal as an OH2, which the edge model of
    Underwood & Bourg treats as water), ion (isolated ions), mineral (everything else).  Also returns warnings for atoms
    closer than physically possible."""
    el = np.array(at.get_chemical_symbols())
    n = len(at)
    I, J, D = neighbor_list("ijd", at, 2.4)
    nbr = [[] for _ in range(n)]
    for i, j, d in zip(I, J, D):
        nbr[i].append((j, d))
    role = np.array(["mineral"] * n, dtype="U12")
    warn = []
    for h in np.nonzero(el == "H")[0]:
        os_ = [(j, d) for j, d in nbr[h] if el[j] == "O" and d < 1.2]
        if len(os_) > 1:
            warn.append(f"H{h + 1} is within 1.2 A of {len(os_)} O atoms (" + ", ".join(f"O{j + 1} {d:.2f} A" for j, d in os_)
                        + "): overlapping atoms, minimise the structure first")
    for o in np.nonzero(el == "O")[0]:
        hs = [j for j, d in nbr[o] if el[j] == "H" and 0.9 <= d <= 1.1 and role[j] != "water"]
        if len(hs) < 2:
            continue
        pair = None
        for a, b in itertools.combinations(hs, 2):
            if 95.0 <= at.get_angle(a, o, b, mic=True) <= 125.0:
                pair = (a, b)
                break
        if pair and len(hs) == 2:
            role[o] = "water"; role[list(pair)] = "water"
    for i in np.nonzero(np.isin(el, IONS))[0]:
        role[i] = "ion"
    # overlapping atoms: closer than 0.75 A (an O-H bond 0.6-1.2 A excepted), or two atoms that are not H closer than 1.2 A
    close = [(i, j, d) for i, j, d in zip(I, J, D) if i < j and (
        (d < 0.75 and not ({el[i], el[j]} == {"O", "H"} and d > 0.6)) or (d < 1.2 and "H" not in (el[i], el[j])))]
    for i, j, d in close[:20]:
        warn.append(f"atoms {i + 1} ({el[i]}) and {j + 1} ({el[j]}) are {d:.2f} A apart: overlapping atoms in the structure")
    return role, warn


def rule_types(at, role):
    """the ClayFF types (and Table 1 charges) that the bond-graph rule of ClayFF gives (soapff_cui._clayff.ff.ff_types; water
    flexible SPC, ions formal): only compared with the types of the two steps (report.json type_check, charges.csv bond_type)"""
    at = at.copy()
    at.set_array("role", np.array(["water" if r == "water" else "ion" if r == "ion" else "bulk" for r in role], dtype="U12"))
    types = ff_types(at)
    q = np.array([CLAYFF_Q.get(t, 0.0) for t in types])
    qO, qH = WATER_Q
    el = at.get_chemical_symbols()
    for i, r in enumerate(role):
        if r == "water":
            q[i] = qO if el[i] == "O" else qH
        elif r == "ion":
            q[i] = FORMAL_ION.get(el[i], q[i])
    return q, types


# ------------------------------------------------------------------------------------------------- type masks
FRAMEWORK_TYPES = ("ao", "at", "mgo", "mgh", "st", "cao", "cah", "feo", "lio") + EDGE_METALS


def hydroxyl_o(t, extra=()):
    """mask of the hydroxyl O types (oh, ohs, and the edge O with one H of the optional edge library) among the type array t"""
    return np.array([x in ("oh", "ohs", *EDGE_HYDROXYL_O, *extra) for x in t], bool)


def framework_cation(t, si=True):
    """mask of the framework cations; si=False leaves Si out"""
    base = FRAMEWORK_TYPES if si else tuple(x for x in FRAMEWORK_TYPES if x != "st")
    return np.array([x in base for x in t], bool)


# ------------------------------------------------------------------------------------------------- boundary conditions
def parse_boundary(boundary):
    """'p p f' -> (['p', 'p', 'f'], [True, True, False]).  Each axis: p, or one or two of f / s / m (lo and hi face),
    as in the LAMMPS boundary command."""
    b = boundary.split()
    if len(b) != 3 or any(not (x == "p" or (1 <= len(x) <= 2 and set(x) <= set("fsm"))) for x in b):
        raise ValueError(f"boundary must be three of p / f / s / m per axis (e.g. 'p p p', 'p p f'), got {boundary!r}")
    return b, [x == "p" for x in b]


def center_nonperiodic(at, periodic):
    """along every non-periodic axis, move the widest empty gap to the box faces so the material sits in the middle
    (a slab given across the periodic boundary would otherwise be cut).  Returns (Atoms, {axis: gap width in A})."""
    s = at.get_scaled_positions(wrap=True)
    gaps = {}
    for k in range(3):
        if periodic[k]:
            continue
        x = np.sort(s[:, k])
        d = np.diff(np.r_[x, x[0] + 1.0])
        j = int(np.argmax(d))
        s[:, k] = np.mod(s[:, k] - (x[j] + d[j] / 2), 1.0)        # gap centre -> face
        h = at.cell.volume / np.linalg.norm(np.cross(*[at.cell[i] for i in range(3) if i != k]))   # box height along k
        gaps["xyz"[k]] = float(d[j] * h)
    at = at.copy()
    at.set_scaled_positions(s)
    return at, gaps


def kspace_setup(kspace, boundary, periodic, tilted):
    """long-range Coulomb for the boundary: (pair_style coul suffix, settings lines, note or None).
    p p p: kspace as given (PPPM).  Only z non-periodic ('f') in an orthogonal box: PPPM/Ewald with the slab correction.
    Anything else (a tilted box, x or y non-periodic, shrink-wrapped faces): MSM, the LAMMPS solver for non-periodic
    boundaries (PPPM/Ewald slab corrections do not support triclinic boxes).
    MSM is run at 1e-4 with the pair cutoff kept (cutoff/adjust no), whatever the kspace accuracy: on a 320-atom
    pyrophyllite slab its energy is as close to PPPM-slab 1e-6 as PPPM-slab 1e-5 is (2.7 vs 4.1 kcal/mol), at ~18x the
    PPPM cost; MSM 1e-5 lets LAMMPS grow the Coulomb cutoff to ~30 A (neighbour list overflow) or, with the cutoff kept,
    costs ~200x."""
    style = kspace.split()[0]
    if all(periodic):
        return "long", [f"kspace_style {kspace}"], None
    if periodic[0] and periodic[1] and boundary[2] == "f" and not tilted and style in ("pppm", "ewald"):
        return "long", [f"kspace_style {kspace}", "kspace_modify slab 3.0"], f"{style} with the slab correction (slab 3.0)"
    why = "tilted (triclinic) box" if tilted else f"boundary {' '.join(boundary)}"
    return "msm", ["kspace_style msm 1.0e-4", "kspace_modify cutoff/adjust no"], \
        f"msm 1.0e-4 ({why}: PPPM/Ewald cannot handle it; ~10-20x slower than PPPM)"


# ------------------------------------------------------------------------------------------------- topology and writer
def topology_from_types(at, types):
    """bonds and angles from the types (not from distances alone, so that overlapping atoms do not create false bonds):
    each water O takes its two nearest water H, each hydroxyl H its nearest hydroxyl O; M-O-H angles for the framework
    cations within 2.5 A (Al) / 2.6 A (Mg) of a hydroxyl O; H-O-H for water."""
    t = np.array(types)
    pos, cell = at.positions, np.array(at.cell)
    inv = np.linalg.inv(cell)

    def d(i, J):
        v = pos[J] - pos[i]
        f = v @ inv
        v = (f - np.round(f)) @ cell
        return np.linalg.norm(v, axis=-1)
    bonds, angles = [], []
    ow, hw = np.nonzero(t == "ow")[0], np.nonzero(t == "hw")[0]
    free = set(hw.tolist())
    pairs = sorted(((float(dd), o, h) for o in ow for h, dd in zip(hw, d(o, hw)) if dd < 1.5))
    nh = {o: [] for o in ow}
    for dd, o, h in pairs:
        if h in free and len(nh[o]) < 2:
            nh[o].append(h); free.discard(h)
    for o, hs in nh.items():
        for h in hs:
            bonds.append((o, h, "OHw"))
        if len(hs) == 2:
            angles.append((hs[0], o, hs[1], ("HOHw",)))
    ohO = np.nonzero(hydroxyl_o(t, extra=("oaq",)))[0]
    fw = np.nonzero(framework_cation(t))[0]
    for h in np.nonzero(t == "ho")[0]:
        if not len(ohO):
            continue
        dd = d(h, ohO); o = ohO[int(np.argmin(dd))]
        bonds.append((o, h, "OH"))
        dm = d(o, fw) if len(fw) else np.array([])
        cat = [m for m, x in zip(fw, dm) if x < {"Si": 2.0, "Al": 2.5, "Mg": 2.6}.get(ELEMENT[t[m]], 2.6)]
        site = "bulk" if len(cat) >= 2 else "surface"
        for m in cat:
            angles.append((m, o, h, ("MOH", ELEMENT[t[m]], site)))
    # edge model (optional library): each edge octahedral metal - O - octahedral metal angle (O within 2.5 A of Al, 2.6 A of Mg)
    em = np.nonzero(np.isin(t, EDGE_METALS))[0]
    if len(em):
        octa = np.nonzero(framework_cation(t, si=False))[0]
        oxy = np.nonzero(np.array([ELEMENT.get(x) == "O" for x in t]))[0]
        rmax = lambda m: {"Al": 2.5, "Mg": 2.6}[ELEMENT[t[m]]]
        seen = set()
        for m in em:
            for o in oxy[d(m, oxy) < rmax(m)]:
                for c in octa:
                    if c != m and d(o, [c])[0] < rmax(c) and (min(m, c), o, max(m, c)) not in seen:
                        seen.add((min(m, c), o, max(m, c)))
                        angles.append((m, o, c, ("MOMe",)))
    unbonded = [int(h) for h in free] + [int(o) for o, hs in nh.items() if len(hs) != 2]
    return bonds, angles, unbonded


def topology_from_bonds(at, types, nb):
    """bonds and angles from the bonds of step 1 (steps.type_atoms: the bond lengths of the library): the O-H bonds (of a water O: OHw),
    the H-O-H angle of each water, an M-O-H angle for each O-H bond and each framework cation bonded to that O, and (edge
    model) an edge metal - O - octahedral metal angle for each O bonded to both"""
    t = np.array(types)
    el = np.array(at.get_chemical_symbols())
    fw = framework_cation(t)
    octa = framework_cation(t, si=False)
    bonds, angles, unbonded = [], [], []
    for o in np.nonzero(el == "O")[0]:
        hs = sorted(int(b) for b, _ in nb[o] if el[b] == "H")
        if t[o] == "ow":
            bonds += [(int(o), h, "OHw") for h in hs]
            if len(hs) == 2:
                angles.append((hs[0], int(o), hs[1], ("HOHw",)))
            else:
                unbonded += [int(o)] + hs
            continue
        cat = sorted(int(b) for b, _ in nb[o] if fw[b])
        for h in hs:
            bonds.append((int(o), h, "OH"))
            site = "bulk" if len(cat) >= 2 else "surface"
            angles += [(m, int(o), h, ("MOH", ELEMENT[t[m]], site)) for m in cat]
    seen = set()
    for m in np.nonzero(np.isin(t, EDGE_METALS))[0]:
        for o, _ in nb[m]:
            for c, _ in nb[o]:
                key = (min(int(m), int(c)), int(o), max(int(m), int(c)))
                if c != m and octa[c] and key not in seen:
                    seen.add(key)
                    angles.append((int(m), int(o), int(c), ("MOMe",)))
    unbonded += [int(h) for h in np.nonzero(t == "hw")[0] if not any(t[b] == "ow" for b, _ in nb[h])]
    return bonds, angles, sorted(set(unbonded))


def _row_coefficients(ff, types, bonds, angles, reserve_water):
    """the coefficients of every bond and angle type of write_system from the rows of the library ([bonds] t1 t2 K r0,
    [angles] t1 t2 t3 K theta0; X = any type, the most specific row wins): {key: (K, r0 or theta0)}.  The terms of one type
    must have one row (else ValueError)"""
    from soapff.molecular import _match
    brows = [(tuple(r[:2]), tuple(float(x) for x in r[2:])) for r in ff["bonds"]]
    arows = [(tuple(r[:3]), tuple(float(x) for x in r[3:])) for r in ff["angles"]]
    t = list(types)
    want_b = [((t[o], t[h]), k) for o, h, k in bonds] + ([(("ow", "hw"), "OHw")] if reserve_water else [])
    want_a = [((t[a], t[b], t[c]), k) for a, b, c, k in angles] + ([(("hw", "ow", "hw"), ("HOHw",))] if reserve_water else [])
    out = []
    for want, rows_, what in ((want_b, brows, "bonds"), (want_a, arows, "angles")):
        co, cache = {}, {}
        for key, k in want:
            if key not in cache:
                cache[key] = _match(key, rows_)
            if cache[key] is None:
                raise ValueError(f"no row of [{what}] in the library for {'-'.join(key)}")
            co.setdefault(k, set()).add(cache[key][1])
        bad = {k: v for k, v in co.items() if len(v) > 1}
        if bad:
            raise ValueError(f"the {what} of one ClayFF kind have different rows in the library: {bad}")
        out.append({k: next(iter(v)) for k, v in co.items()})
    return out


def write_system(at, q, types, out_dir, cutoff=10.0, kspace="pppm 1.0e-5", boundary="p p p", reserve_water=False,
                 charges="clayff2004", water="spc", ff=None, nb=None):
    """reserve_water: the water types (ow, hw), the water O-H bond and H-O-H angle types and room for their bonds are
    written even when the structure has no water, so that water can be added in LAMMPS (swelling curve).
    ff: the force field of the library file (soapff.clayff: the mass, element and LJ of every type, the rows of [bonds] and
    [angles], [lammps]): every coefficient is taken from it, the bonded ones from the row of the types of each bond and
    angle; or (read_forcefield, the tables of former versions) a force field whose LJ, hydroxyl O-H bond, M-O-H and edge
    M-O-M angles replace ClayFF's for its types; None: ClayFF (soapff_cui/clayff.py).
    nb: the bonds of step 1 (steps.type_atoms): the bonds and angles come from them (topology_from_bonds); without
    nb, from the types and the distances (topology_from_types)."""
    bnd, periodic = parse_boundary(boundary)
    bonds, angles, unbonded = topology_from_bonds(at, types, nb) if nb is not None else topology_from_types(at, types)
    rows = bool(ff) and "bonds" in ff                  # the force field of a library file
    elem = (lambda x: ff["element"].get(x) or ELEMENT[x]) if rows else (lambda x: ELEMENT[x])   # noqa: E731
    tnames = sorted(set(types) | ({"ow", "hw"} if reserve_water else set()), key=lambda x: (elem(x), x))
    tid = {x: k + 1 for k, x in enumerate(tnames)}
    bkeys = sorted(set(b[2] for b in bonds) | ({"OHw"} if reserve_water else set())); bid = {k: n + 1 for n, k in enumerate(bkeys)}
    akeys = sorted(set(a[3] for a in angles) | ({("HOHw",)} if reserve_water else set())); aid = {k: n + 1 for n, k in enumerate(akeys)}
    mol = np.ones(len(at), int); nxt = 2
    for o, h, k in bonds:
        if k == "OHw":
            if mol[o] == 1:
                mol[o] = nxt; nxt += 1
            mol[h] = mol[o]
    for i, x in enumerate(types):
        if x in ("Na", "K", "Ca", "Cs", "Cl", "Ba"):
            mol[i] = nxt; nxt += 1
    p, xl = lammps_frame(at, wrap=True)
    lx, ly, lz, xy, xz, yz = p.get_lammps_prism()
    L = ["soapff-cui: ClayFF types/LJ/bonded terms, charges from SOAP", "", f"{len(at)} atoms", f"{len(tnames)} atom types",
         f"{len(bonds)} bonds", f"{len(bkeys)} bond types", f"{len(angles)} angles", f"{len(akeys)} angle types"]
    if reserve_water:
        L += ["2 extra bond per atom", "1 extra angle per atom", "2 extra special per atom"]
    L += ["",
         f"0.0 {lx:.8f} xlo xhi", f"0.0 {ly:.8f} ylo yhi", f"0.0 {lz:.8f} zlo zhi"]
    if any(abs(v) > 1e-8 for v in (xy, xz, yz)):
        L.append(f"{xy:.8f} {xz:.8f} {yz:.8f} xy xz yz")
    L += ["", "Masses", ""] + [f"{tid[x]} {(ff['mass'][x] if rows else MASS[ELEMENT[x]]):.5f}  # {x}" for x in tnames]
    # image flags: an H bonded across the cell boundary gets the image that puts it next to its O
    img = np.zeros((len(at), 3), int)
    frac = at.get_scaled_positions(wrap=True)
    for o, h, k in bonds:
        img[h] = -np.round(frac[h] - frac[o]).astype(int)
    cut = [(o, h) for o, h, k in bonds if any(img[h, a] and not periodic[a] for a in range(3))]
    if cut:
        o, h = cut[0]
        raise ValueError(f"{len(cut)} O-H bonds cross a non-periodic boundary (e.g. O{o + 1}-H{h + 1}): the structure has no "
                         f"empty gap along that axis; add vacuum or use a periodic boundary")
    coul, kspace_lines, kspace_note = kspace_setup(kspace, bnd, periodic, any(abs(v) > 1e-8 for v in (xy, xz, yz)))
    L += ["", "Atoms  # full", ""] + [f"{i + 1} {mol[i]} {tid[types[i]]} {q[i]:.8f} {xl[i, 0]:.8f} {xl[i, 1]:.8f} {xl[i, 2]:.8f}"
                                      f" {img[i, 0]} {img[i, 1]} {img[i, 2]}" for i in range(len(at))]
    if bonds:
        L += ["", "Bonds", ""] + [f"{k + 1} {bid[b[2]]} {b[0] + 1} {b[1] + 1}" for k, b in enumerate(bonds)]
    if angles:
        L += ["", "Angles", ""] + [f"{k + 1} {aid[a[3]]} {a[0] + 1} {a[1] + 1} {a[2] + 1}" for k, a in enumerate(angles)]
    open(os.path.join(out_dir, "system.data"), "w").write("\n".join(L) + "\n")
    if rows:
        lm = ff.get("lammps", {})
        S = [f"# soapff-cui settings: {ff['name']} ({ff['source']}), the parameters of the library file; types and charges in "
             f"system.data from SOAP (the reference atoms)", "units real", "atom_style full", f"boundary {' '.join(bnd)}",
             f"bond_style {lm.get('bond_style', 'harmonic')}", f"angle_style {lm.get('angle_style', 'harmonic')}",
             f"pair_style lj/cut/coul/{coul} {cutoff}", f"pair_modify {lm.get('pair_modify', 'mix arithmetic tail no')}",
             f"special_bonds {lm.get('special_bonds', 'lj/coul 0.0 0.0 0.0')}", "read_data system.data"]
        bco, aco = _row_coefficients(ff, types, bonds, angles, reserve_water)
        for x in tnames:
            e, sg = ff["lj"][x]
            S.append(f"pair_coeff {tid[x]} {tid[x]} {e!r} {sg!r}  # {x}")
        for k in bkeys:
            K, r0 = bco[k]
            S.append(f"bond_coeff {bid[k]} {K!r} {r0!r}  # {'O-H (water)' if k == 'OHw' else 'O-H'}")
        for k in akeys:
            K, th = aco[k]
            S.append(f"angle_coeff {aid[k]} {K!r} {th!r}  # {'-'.join(map(str, k))}"
                     + ("  (edge metal-O-metal)" if k[0] == "MOMe" else ""))
        S += kspace_lines
        open(os.path.join(out_dir, "system.settings"), "w").write("\n".join(S) + "\n")
        return dict(type_ids=tid, water_bond=bid.get("OHw"), water_angle=aid.get(("HOHw",)), nbonds=len(bonds),
                    nangles=len(angles), unbonded=unbonded, nmol=int(nxt - 1), kspace=kspace_lines, kspace_note=kspace_note,
                    tilted=any(abs(v) > 1e-8 for v in (xy, xz, yz)))
    S = [(f"# soapff-cui settings: {ff['name']} ({ff['source']}) for its types; the other types and the {WATER_LABELS[water]} "
          f"water: ClayFF (Cygan et al. 2004, Tables 1 and 2); types in system.data from SOAP" if ff else
          f"# soapff-cui settings: ClayFF (Cygan et al. 2004, Tables 1 and 2) with {WATER_LABELS[water]} water; types in system.data "
          f"from SOAP, charges: {charges}"),
         "units real", "atom_style full", f"boundary {' '.join(bnd)}",
         "bond_style harmonic", "angle_style harmonic", f"pair_style lj/cut/coul/{coul} {cutoff}", "pair_modify mix arithmetic tail no",
         "special_bonds lj/coul 0.0 0.0 0.0", "read_data system.data"]
    own = ff["lj"] if ff else {}
    for x in tnames:
        # Table 1 as published: epsilon = D0, sigma = R0 / 2^(1/6) written in full double precision (not rounded); the types of
        # the library's force field: its own values
        e, sg = own[x] if x in own else LJ[x]
        S.append(f"pair_coeff {tid[x]} {tid[x]} {e!r} {sg!r}  # {x}" + (f" ({ff['name']})" if x in own else ""))
    for k in bkeys:
        K, r0 = ff["bond_oh"] if (ff and k != "OHw") else BOND_OH
        S.append(f"bond_coeff {bid[k]} {K} {r0}  # {'O-H (water)' if k == 'OHw' else 'O-H'}"
                 + (f" ({ff['name']})" if ff and k != "OHw" else ""))
    for k in akeys:
        if k[0] in ("HOH", "HOHw"):
            K, th = ANGLE_HOH
        elif k[0] == "MOMe":
            K, th = ff["angle_mom_edge"] if ff else ANGLE_MOM_EDGE
        else:
            K, th = ff["angle_moh"] if ff else ANGLE_MOH
        S.append(f"angle_coeff {aid[k]} {K} {th}  # {'-'.join(map(str, k))}"
                 + ("  (edge metal-O-metal)" if k[0] == "MOMe" else "") + (f" ({ff['name']})" if ff and k[0] not in ("HOH", "HOHw") else ""))
    S += kspace_lines
    open(os.path.join(out_dir, "system.settings"), "w").write("\n".join(S) + "\n")
    return dict(type_ids=tid, water_bond=bid.get("OHw"), water_angle=aid.get(("HOHw",)), nbonds=len(bonds),
                nangles=len(angles), unbonded=unbonded, nmol=int(nxt - 1), kspace=kspace_lines, kspace_note=kspace_note,
                tilted=any(abs(v) > 1e-8 for v in (xy, xz, yz)))


# ------------------------------------------------------------------------------------------------- water added in LAMMPS
# (the swelling curve: examples/montmorillonite_swelling): the water types reserved, and the molecule file of one water
SPC_OH, SPC_HOH = 1.0, 109.47                               # flexible SPC water of ClayFF: O-H (A), H-O-H (deg) at rest


def reserve_water_types(out_dir, log=print):
    """rewrite system.data / system.settings of an output folder with the water types reserved (write_system
    reserve_water), from the structure of report.json and the charges and types of charges.csv.
    Returns the write_system info."""
    rep = json.load(open(os.path.join(out_dir, "report.json")))
    at = read_structure(rep["structure"])
    rows = list(csv.DictReader(open(os.path.join(out_dir, "charges.csv"))))
    if len(rows) != len(at):
        raise ValueError(f"charges.csv has {len(rows)} atoms, the structure {len(at)}")
    q = np.array([float(r["charge"]) for r in rows])
    types = [r["clayff_type"] for r in rows]
    kspace = " ".join(rep.get("kspace", [""])[0].split()[1:]) if rep.get("kspace") else "pppm 1.0e-5"
    info = write_system(at, q, types, out_dir, kspace=kspace or "pppm 1.0e-5", boundary=rep.get("boundary", "p p p"),
                        reserve_water=True, charges=rep.get("charges", "clayff2004"), water=rep.get("water_model", "spc"),
                        ff=rep.get("forcefield"))
    log(f"system.data / system.settings rewritten with the water types reserved (ow {info['type_ids']['ow']}, "
        f"hw {info['type_ids']['hw']}, O-H bond {info['water_bond']}, H-O-H angle {info['water_angle']})")
    return info


def water_molecule(type_ids, bond_type, angle_type, water="spc"):
    """LAMMPS molecule file of one flexible water (charges of the water model `water`, spc) for create_atoms ... mol."""
    h = np.radians(SPC_HOH / 2)
    xs = [(0.0, 0.0, 0.0), (SPC_OH * np.sin(h), SPC_OH * np.cos(h), 0.0), (-SPC_OH * np.sin(h), SPC_OH * np.cos(h), 0.0)]
    t = [type_ids["ow"], type_ids["hw"], type_ids["hw"]]
    qq = [WATER_MODELS[water][0], WATER_MODELS[water][1], WATER_MODELS[water][1]]
    L = [f"# {WATER_LABELS[water]} water, written by soapff-cui", "", "3 atoms", "2 bonds", "1 angles", "", "Coords", ""]
    L += [f"{k + 1} {x:.6f} {y:.6f} {z:.6f}" for k, (x, y, z) in enumerate(xs)]
    L += ["", "Types", ""] + [f"{k + 1} {v}" for k, v in enumerate(t)]
    L += ["", "Charges", ""] + [f"{k + 1} {v:.4f}" for k, v in enumerate(qq)]
    L += ["", "Bonds", "", f"1 {bond_type} 1 2", f"2 {bond_type} 1 3", "", "Angles", "", f"1 {angle_type} 2 1 3"]
    return "\n".join(L) + "\n"


def make(structure, out_dir, library, **opts):
    """the force field of a structure with a table of the two steps (steps.build_table; the prebuilt tables of the clayff
    libraries, or library.table_path): make_steps"""
    return make_steps(structure, out_dir, library, **opts)



def _write(at, q, types, out_dir, bnd, periodic, cutoff, kspace, reserve_water, charges, water, ff, warn, log, nb=None):
    """the non-periodic faces moved to the widest gap, system.data and system.settings (write_system), the walls of the
    fixed faces and the neighbour list of the shrink-wrapped ones -> (the structure as written, write_system's info, gaps)"""
    gaps = {}
    if not all(periodic):
        at, gaps = center_nonperiodic(at, periodic)
        for ax, g in gaps.items():
            log(f"boundary {' '.join(bnd)}: widest empty gap along {ax} = {g:.2f} A, moved to the box faces")
            if g < cutoff:
                warn.append(f"the empty gap along the non-periodic {ax} axis is only {g:.2f} A (< cutoff {cutoff} A): atoms on "
                            f"the two faces interacted through the periodic boundary of the given structure, now they do not")
                log("WARNING " + warn[-1])
    info = write_system(at, q, types, out_dir, cutoff=cutoff, kspace=kspace, boundary=" ".join(bnd), reserve_water=reserve_water,
                        charges=charges, water=water, ff=ff, nb=nb)
    if info["kspace_note"]:
        log(f"boundary {' '.join(bnd)}: long-range Coulomb by {info['kspace_note']}")
    # atoms reaching a fixed ('f') face would be lost: reflect them there
    walls = " ".join(f"{ax}{side} EDGE" for ax, b in zip("xyz", bnd) if b != "p"
                     for side, c in zip(("lo", "hi"), b * 2 if len(b) == 1 else b) if c == "f")
    extra = ["# the fixed faces of the boundary: atoms are reflected there", f"fix wall all wall/reflect {walls}"] if walls else []
    # a shrink-wrapped ('s') face sits on the outermost atom; between rebuilds atoms move past it and MSM stops
    # ("Out of range atoms").  Rebuilding every step re-shrinks the box before each force call (on a 320-atom slab no
    # slower than 'm' or 'f': MSM dominates the cost).
    if any("s" in b for b in bnd):
        extra += ["# shrink-wrapped faces with MSM: the neighbour list is rebuilt every step", "neigh_modify every 1 delay 0 check no"]
        log(f"boundary {' '.join(bnd)}: 's' faces with MSM -> neighbour list rebuilt every step (neigh_modify check no)")
    if extra:
        with open(os.path.join(out_dir, "system.settings"), "a") as f:
            f.write("\n".join(extra) + "\n")
    if info["unbonded"]:
        warn.append(f"{len(info['unbonded'])} water O/H could not be paired into molecules (atoms {[i + 1 for i in info['unbonded'][:10]]})")
        log("WARNING " + warn[-1])
    return at, info, gaps


ION_TYPES = tuple(FORMAL_ION) + ("Na", "K", "Ca", "Cs", "Cl", "Ba")


def make_steps(structure, out_dir, library, threshold=THRESHOLD, cutoff=10.0, kspace="pppm 1.0e-5", boundary="p p p",
               check=True, reserve_water=False, strict=False, charges="clayff2004", water="spc", forcefield=None, log=print):
    """make with a table of the two steps (steps.build_table): the mineral, the water and the ions in one pass; the water
    molecules (an O bonded to two H only, soapff.molecular.water_apart) are taken apart from the mineral, so that the mineral
    is typed without them and each molecule on its own bonds.  Stops when an element is not in the table; when an atom's
    similarity (step 2) is below `threshold` (strict; else the atom keeps its row and is reported: charges.csv
    source soap-low); and (check) when an atom does not have the bonds of the reference atom of its row (the bonds of step 1
    disagree with the type of step 2).
    The bonds and angles of the force field come from the bonds of step 1 (topology_from_bonds).
    forcefield: the force field of the library file (soapff.clayff._forcefield), with the charges of the reference atoms;
    default: the one of the table, else ClayFF (soapff_cui/clayff.py) with the charges of Table 1."""
    from soapff.molecular import TypingError

    from . import steps
    bnd, periodic = parse_boundary(boundary)
    if charges not in CHARGE_SCHEMES or water not in WATER_MODELS:
        raise ValueError(f"charges must be one of {list(CHARGE_SCHEMES)}, water one of {list(WATER_MODELS)}")
    os.makedirs(out_dir, exist_ok=True)
    at = read_structure(structure) if isinstance(structure, str) else structure
    T = steps.load_table(library)
    _, warn = roles(at)                                  # only its warnings: overlapping atoms
    for w_ in warn:
        log("WARNING " + w_)
    el = np.array(at.get_chemical_symbols())
    missing = sorted(set(el) - set(T["species"]))
    if missing:
        raise TypingError(f"the elements {missing} are not in the library {os.path.abspath(library)} (it has {T['species']}): "
                            f"add a reference structure with them")
    r = steps.type_atoms(at, T)
    nb = r["nb"]
    if r["water"] and r["water"]["bonds_taken_out"]:
        log("NOTE the water typed apart from the mineral: bonds taken out between water molecules and the rest (number, "
            "shortest A): " + ", ".join(f"{k} {n} ({d:.2f})" for k, (n, d) in r["water"]["bonds_taken_out"].items()))
    types = [TYPE_ALIAS.get(str(x), str(x)) for x in r["type"]]
    tarr = np.array(types)
    sim = r["sim"]
    low = sim < threshold
    if low.any():
        idx = np.nonzero(low)[0]
        detail = (", ".join(f"{el[i]}{i + 1} ({sim[i]:.3f}, nearest {types[i]} {r['origin'][i]})" for i in idx[:10])
                  + (" ..." if len(idx) > 10 else ""))
        if strict:
            raise TypingError(f"{len(idx)} atoms are below similarity {threshold} (step 2): {detail}; library "
                              f"{os.path.abspath(library)}. The library lacks these environments: add a reference structure "
                              f"with them to the library")
        warn.append(f"{len(idx)} atoms below similarity {threshold} kept the nearest row (min {sim[low].min():.3f}): {detail}; "
                    f"see charges.csv (source soap-low)")
        log("WARNING " + warn[-1])
    odd = np.nonzero(r["bonds"] != r["want"])[0]
    problems = []
    if len(odd):
        problems.append(f"{len(odd)} atoms do not have the bonds of the reference atom of their type: " + ", ".join(
            f"{el[i]}{i + 1} {types[i]} ({sim[i]:.3f}) bonded to {r['bonds'][i]}, its reference atom {r['origin'][i]} to "
            f"{r['want'][i]}" for i in odd[:8])
            + (" ..." if len(odd) > 8 else ""))
        for p_ in problems:
            log("CHECK " + p_)
        if check:
            json.dump(dict(structure=structure if isinstance(structure, str) else "Atoms", consistency=problems),
                      open(os.path.join(out_dir, "report.json"), "w"), indent=1)
            raise TypingError("the bonds found do not fit the types (see report.json): " + "; ".join(problems)
                              + ".  The library lacks these environments, or the structure has overlapping atoms")
        warn += [f"consistency: {p_}" for p_ in problems]
    ff = forcefield if forcefield is not None else T["forcefield"]
    if ff and "bonds" in ff:                   # the force field of a library file: the charges of the reference atoms
        q = np.asarray(r["q"], float).copy()
    else:
        tq = type_charges(charges, water, types=types)
        if ff:
            tq.update(ff["charge"])
        q = np.array([tq.get(t, r["q"][i]) for i, t in enumerate(types)], float)
    dq = float(np.nanmax(np.abs(q - r["q"])))
    if dq > 1e-4:
        warn.append(f"the charges of the types differ from those of the reference atoms by up to {dq:.4f} e")
        log("WARNING " + warn[-1])
    role = np.where(np.isin(tarr, ["ow", "hw"]), "water", np.where(np.isin(tarr, ION_TYPES), "ion", "mineral"))
    # for comparison only (not used): the type the bond-graph rule of ClayFF gives (charges.csv column bond_type)
    _, types_geo = rule_types(at, role)
    geo = np.array([TYPE_ALIAS.get(str(t), str(t)) for t in types_geo])
    cls = lambda t: ("bridge" if t in ("ob", "obts", "obos", "obss") else "oh" if t in ("oh", "ohs") else  # noqa: E731
                     "water" if t in ("ow", "oaq", "hw") else "h" if t == "ho" else "cation")
    same_class = np.array([cls(a) == cls(b) and cls(a) in ("bridge", "oh", "cation") and a not in EDGE_PARENT
                           for a, b in zip(tarr, geo)])
    differ = same_class & (tarr != geo)
    type_check = dict(compared=int(same_class.sum()), differ=int(differ.sum()),
                      kinds=dict(collections.Counter(f"{g}->{t}" for g, t in zip(geo[differ], tarr[differ]))))
    # the net charge is written (report.json charge_sum, the log), not checked: a clay layer or particle may carry a charge
    # (LAMMPS adds a uniform neutralising background for the long-range Coulomb)
    total = float(q.sum())
    at, info, gaps = _write(at, q, types, out_dir, bnd, periodic, cutoff, kspace, reserve_water, charges, water, ff, warn, log,
                            nb=nb)
    with open(os.path.join(out_dir, "charges.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "element", "role", "clayff_type", "charge", "source", "matched_library_type", "similarity", "bond_type",
                    "bonds", "margin", "reference_atom"])
        for i in range(len(at)):
            w.writerow([i + 1, el[i], role[i], types[i], f"{q[i]:+.6f}", "soap-low" if low[i] else "soap", r["type"][i],
                        f"{sim[i]:.5f}", geo[i], r["bonds"][i], f"{r['margin'][i]:.5f}", r["origin"][i]])
    per = {}
    for t in sorted(set(types)):
        m = tarr == t
        per[t] = dict(atoms=int(m.sum()), charges=sorted({round(float(x), 4) for x in q[m]}), min_similarity=float(sim[m].min()), bonds=sorted(set(r["bonds"][m])))
    rep = dict(structure=structure if isinstance(structure, str) else "Atoms", atoms=len(at), library=os.path.abspath(library),
               soap=T["soap"], depth=T["depth"], threshold=threshold, boundary=" ".join(bnd), kspace=info["kspace"],
               vacuum_gap=gaps, charges=charges, water_model=water, forcefield=ff,
               charge_sum=float(q.sum()), below_threshold=int(low.sum()),
               min_similarity=float(sim.min()), min_margin=float(np.nanmin(r["margin"])) if len(at) else None,
               bond_lengths={k: v["r"] for k, v in T["bond_lengths"].items()}, bonds=info["nbonds"], angles=info["nangles"],
               type_check=type_check, water_apart=r["water"], warnings=warn,
               roles={k: int((role == k).sum()) for k in ("mineral", "water", "ion")}, per_type=per, consistency=problems,
               type_ids=info["type_ids"], type_elements=[ELEMENT[t] for t, _ in sorted(info["type_ids"].items(), key=lambda kv: kv[1])],
               rigid=dict(bonds=[], angles=[]), timestep_fs=TIMESTEP_FS)
    json.dump(rep, open(os.path.join(out_dir, "report.json"), "w"), indent=1)
    log(f"{len(at)} atoms ({rep['roles']}); similarity >= {rep['min_similarity']:.3f}; {info['nbonds']} bonds, "
        f"{info['nangles']} angles; "
        f"net charge {total:+.4f} e; wrote {out_dir}/system.data, system.settings, charges.csv, report.json")
    return rep
