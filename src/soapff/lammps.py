"""The LAMMPS files of soapff: system.data (atoms, charges, bonds, angles, dihedrals, impropers) and system.settings (units,
styles, coefficients, special_bonds, kspace, the walls of fixed faces, and SHAKE for the rigid molecules of the library).
These are the force field and nothing else: minimisation, velocities, thermostats, barostats, runs and output are written
by the user, in LAMMPS (e.g. its Python module: lmp.commands_list(system.commands()))."""
import os

import numpy as np

from soapff_cui._clayff.ff import lammps_frame
from soapff_cui.core import center_nonperiodic, kspace_setup, parse_boundary

TIMESTEP_FS = 1.0          # recommended time step (fs) when the library gives none (report.json)


def _num(x):
    return f"{x:.10g}"


def write(out_dir, at, nb, types, q, topo, lib, cutoff=None, kspace=None, boundary="p p p", pair_style=None,
          dihedral_style=None, reserve=None, log=print):
    """system.data and system.settings of a typed structure; returns the type ids, the SHAKE bonds and angles and the counts.
    Along a non-periodic axis the widest empty gap is moved to the box faces (the material in the middle).
    reserve: the water of the library (molecular.water_types): its types, its bond and angle types and room for their bonds
    are written even when the structure has no water, so that water can be added in LAMMPS."""
    from .molecular import unwrapped
    lm = lib.lammps
    bnd, periodic = parse_boundary(boundary)
    warnings, gaps = [], {}
    if not all(periodic):
        at, gaps = center_nonperiodic(at, periodic)
        for ax, g in gaps.items():
            log(f"boundary {' '.join(bnd)}: widest empty gap along {ax} = {g:.2f} A, moved to the box faces")
            rc = float((pair_style or lm.get("pair_style", "lj/cut/coul/long 10.0")).split()[1]) if cutoff is None else cutoff
            if g < rc:
                warnings.append(f"the empty gap along the non-periodic {ax} axis is only {g:.2f} A (< cutoff {rc:g} A): atoms "
                                f"on the two faces interacted through the periodic boundary of the given structure, now they "
                                f"do not")
                log("WARNING " + warnings[-1])
    rtypes = {reserve["o"], reserve["h"]} if reserve else set()
    tnames = sorted(set(types) | rtypes, key=lambda t: (lib.types[t]["element"], t))
    tid = {t: k + 1 for k, t in enumerate(tnames)}
    bkeys = sorted({k for _, k in topo["bonds"]} | ({reserve["bond"]} if reserve else set()))
    bid = {k: i + 1 for i, k in enumerate(bkeys)}
    akeys = sorted({k for _, k in topo["angles"]} | ({reserve["angle"]} if reserve else set()))
    aid = {k: i + 1 for i, k in enumerate(akeys)}
    bpar = dict(lib.bonds)
    apar = dict(lib.angles)
    dpar = dict(lib.dihedrals)
    ipar = dict(lib.impropers)
    dstyle = dihedral_style or lm.get("dihedral_style", "fourier")
    # dihedral types: one per row (fourier: all its terms in one type) or one per term (charmm, weight 0: the 1-4 pairs come
    # from special_bonds)
    dihedrals = []
    dtypes = {}
    for idx, key in topo["dihedrals"]:
        terms = [t for t in dpar[key] if t[0] != 0.0]     # terms with K = 0 add nothing (OpenMM leaves them out too)
        if not terms:
            continue
        if dstyle.split()[0] == "charmm":
            for K, n, d in terms:
                if abs(d - round(d)) > 1e-6:
                    raise ValueError(f"dihedral {'-'.join(key)}: phase {d} is not a whole degree; use dihedral_style fourier")
                tk = (key, K, n, int(round(d)))
                dtypes.setdefault(tk, len(dtypes) + 1)
                dihedrals.append((idx, dtypes[tk]))
        else:
            dtypes.setdefault(key, len(dtypes) + 1)
            dihedrals.append((idx, dtypes[key]))
    ikeys = sorted({k for _, k in topo["impropers"]})
    iid = {k: i + 1 for i, k in enumerate(ikeys)}
    # molecules (of the bonds of step 1: a mineral layer is one) and image flags (each molecule of the bonds of the force
    # field whole; the bonds of step 1 that are not, a metal of ClayFF and its O, leave the atoms in the cell)
    from .molecular import molecules
    mol = molecules(nb) + 1
    cell = np.array(at.cell)
    off = {(min(i, j), max(i, j)) for i, j in topo.get("not_ff", [])}
    nb_ff = [[(j, d) for j, d in x if (min(i, j), max(i, j)) not in off] for i, x in enumerate(nb)] if off else nb
    pos_u = unwrapped(at, nb_ff)
    frac = pos_u @ np.linalg.inv(cell)
    img = np.where(periodic, np.floor(frac), 0).astype(int)
    cross = [k for k in range(3) if not periodic[k] and ((frac[:, k] < -1e-9).any() or (frac[:, k] >= 1.0).any())]
    if cross:
        raise ValueError(f"molecules cross the non-periodic boundary along {', '.join('xyz'[k] for k in cross)}: the structure has "
                         f"no empty gap along that axis wide enough; add vacuum or use a periodic boundary")
    wrapped = at.copy()
    wrapped.positions = (frac - img) @ cell
    p, xl = lammps_frame(wrapped, wrap=False)
    lx, ly, lz, xy, xz, yz = p.get_lammps_prism()
    tilted = any(abs(v) > 1e-8 for v in (xy, xz, yz))
    L = [f"soapff: library {lib.name}, types and charges from the references (SOAP)", "",
         f"{len(at)} atoms", f"{len(tnames)} atom types", f"{len(topo['bonds'])} bonds", f"{len(bkeys)} bond types",
         f"{len(topo['angles'])} angles", f"{len(akeys)} angle types", f"{len(dihedrals)} dihedrals",
         f"{len(dtypes)} dihedral types", f"{len(topo['impropers'])} impropers", f"{len(ikeys)} improper types"]
    if reserve:
        L += ["2 extra bond per atom", "1 extra angle per atom", "2 extra special per atom"]
    L += ["",
         f"0.0 {lx:.8f} xlo xhi", f"0.0 {ly:.8f} ylo yhi", f"0.0 {lz:.8f} zlo zhi"]
    if tilted:
        L.append(f"{xy:.8f} {xz:.8f} {yz:.8f} xy xz yz")
    L += ["", "Masses", ""] + [f"{tid[t]} {lib.types[t]['mass']:.5f}  # {t}" for t in tnames]
    L += ["", "Atoms  # full", ""] + [f"{i + 1} {mol[i]} {tid[types[i]]} {q[i]:.6f} {xl[i, 0]:.8f} {xl[i, 1]:.8f} {xl[i, 2]:.8f}"
                                      f" {img[i, 0]} {img[i, 1]} {img[i, 2]}" for i in range(len(at))]
    if topo["bonds"]:
        L += ["", "Bonds", ""] + [f"{k + 1} {bid[key]} {a + 1} {b + 1}" for k, ((a, b), key) in enumerate(topo["bonds"])]
    if topo["angles"]:
        L += ["", "Angles", ""] + [f"{k + 1} {aid[key]} {a + 1} {b + 1} {c + 1}"
                                   for k, ((a, b, c), key) in enumerate(topo["angles"])]
    if dihedrals:
        L += ["", "Dihedrals", ""] + [f"{k + 1} {t} {a + 1} {b + 1} {c + 1} {d + 1}" for k, ((a, b, c, d), t) in enumerate(dihedrals)]
    if topo["impropers"]:
        L += ["", "Impropers", ""] + [f"{k + 1} {iid[key]} {a + 1} {b + 1} {c + 1} {d + 1}"
                                      for k, ((a, b, c, d), key) in enumerate(topo["impropers"])]
    open(os.path.join(out_dir, "system.data"), "w").write("\n".join(L) + "\n")

    # settings
    ps = (pair_style or lm.get("pair_style", "lj/cut/coul/long 10.0")).split()
    if cutoff is not None:
        ps = ps[:1] + [str(cutoff)]
    long_range = ps[0].endswith("coul/long")
    ks_lines, note = [], None
    if long_range:
        coul, ks_lines, note = kspace_setup(kspace or lm.get("kspace_style", "pppm 1.0e-5"), bnd, periodic, tilted)
        if coul == "msm":
            ps[0] = ps[0].replace("coul/long", "coul/msm")
    S = [f"# soapff settings: library {lib.name} ({lib.reference}); types and charges in "
         f"system.data from the references of the library", "units real", "atom_style full", f"boundary {' '.join(bnd)}",
         f"pair_style {' '.join(ps)}"]
    if lm.get("pair_modify"):
        S.append(f"pair_modify {lm['pair_modify']}")
    S += [f"bond_style {lm.get('bond_style', 'harmonic')}", f"angle_style {lm.get('angle_style', 'harmonic')}"]
    if dtypes:
        S.append(f"dihedral_style {dstyle}")
    if ikeys:
        S.append(f"improper_style {lm.get('improper_style', 'cvff')}")
    S += [f"special_bonds {lm.get('special_bonds', 'lj/coul 0.0 0.0 0.0')}", "read_data system.data"]
    for t in tnames:
        x = lib.types[t]
        S.append(f"pair_coeff {tid[t]} {tid[t]} {_num(x['epsilon'])} {_num(x['sigma'])}  # {t}")
    for k in bkeys:
        K, r0 = bpar[k]
        S.append(f"bond_coeff {bid[k]} {_num(K)} {_num(r0)}  # {'-'.join(k)}")
    for k in akeys:
        K, th = apar[k]
        S.append(f"angle_coeff {aid[k]} {_num(K)} {_num(th)}  # {'-'.join(k)}")
    for k, i in sorted(dtypes.items(), key=lambda kv: kv[1]):
        if dstyle.split()[0] == "charmm":
            key, K, n, d = k
            S.append(f"dihedral_coeff {i} {_num(K)} {n} {d} 0.0  # {'-'.join(key)}")
        else:
            terms = [t for t in dpar[k] if t[0] != 0.0]
            S.append(f"dihedral_coeff {i} {len(terms)} " + " ".join(f"{_num(K)} {n} {_num(d)}" for K, n, d in terms)
                     + f"  # {'-'.join(k)}")
    for k in ikeys:
        K, d, n = ipar[k]
        S.append(f"improper_coeff {iid[k]} {_num(K)} {int(round(d))} {int(round(n))}  # {'-'.join(k)}")
    S += ks_lines
    # a fixed ('f') face would lose the atoms that reach it: they are reflected there.  A shrink-wrapped ('s') face sits on
    # the outermost atom; between rebuilds atoms move past it and MSM stops ("Out of range atoms"): the neighbour list is
    # rebuilt every step (on a 320-atom slab no slower than 'm' or 'f': MSM dominates the cost)
    walls = " ".join(f"{ax}{side} EDGE" for ax, b in zip("xyz", bnd) if b != "p"
                     for side, c in zip(("lo", "hi"), b * 2 if len(b) == 1 else b) if c == "f")
    if walls:
        S += ["# the fixed faces of the boundary: atoms are reflected there", f"fix wall all wall/reflect {walls}"]
    if any("s" in b for b in bnd):
        S += ["# shrink-wrapped faces with MSM: the neighbour list is rebuilt every step", "neigh_modify every 1 delay 0 check no"]
        log(f"boundary {' '.join(bnd)}: 's' faces with MSM -> neighbour list rebuilt every step (neigh_modify check no)")
    if note:
        log(f"boundary {' '.join(bnd)}: long-range Coulomb by {note}")
    # SHAKE: the bonds and angles whose atoms are all of the rigid types (rigid water: part of the model).  LAMMPS keeps the
    # constraints as stiff restraints during a minimisation.  A single point of the flexible molecule leaves it out
    # (System.commands(rigid=False), or unfix rigid).
    rigid = set(lm.get("rigid", "").split())
    shake_b = [bid[k] for k in bkeys if rigid and set(k) <= rigid]
    shake_a = [aid[k] for k in akeys if rigid and set(k) <= rigid]
    if shake_b or shake_a:
        S += [f"# the rigid molecules of the model (library: rigid = {lm.get('rigid', '')}): SHAKE",
              "fix rigid all shake 1.0e-4 20 0" + (" b " + " ".join(map(str, shake_b)) if shake_b else "")
              + (" a " + " ".join(map(str, shake_a)) if shake_a else "")]
    open(os.path.join(out_dir, "system.settings"), "w").write("\n".join(S) + "\n")
    return dict(type_ids=tid, elements=[lib.types[t]["element"] for t in tnames], shake_bonds=shake_b, shake_angles=shake_a,
                nbonds=len(topo["bonds"]), nangles=len(topo["angles"]), ndihedrals=len(dihedrals),
                nimpropers=len(topo["impropers"]), kspace=ks_lines, kspace_note=note, tilted=tilted, nmol=int(mol.max()),
                box=[lx, ly, lz, xy, xz, yz], gaps=gaps, warnings=warnings,
                reserved_water=dict(types=[tid[reserve["o"]], tid[reserve["h"]]], bond=bid[reserve["bond"]],
                                    angle=aid[reserve["angle"]]) if reserve else None)
