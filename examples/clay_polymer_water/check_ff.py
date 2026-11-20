"""Is the force field of clayff_oplsaa ClayFF and OPLS-AA?  The force field soapff writes with clayff_oplsaa for a start
structure is compared, atom by atom and term by term, with the force fields clayff and oplsaa write for their own atoms
(the mineral, water and ions with clayff, the molecules with oplsaa, each in the cell of the structure):

- the type and charge of every atom, and the LJ (epsilon, sigma) of every type, are the same;
- the bonds, angles, dihedrals and impropers (their atoms and coefficients) are those of the two put together;
- the LJ of every pair of types is mixed by the Lorentz-Berthelot rules (pair_modify mix arithmetic; clayff mixes so,
  oplsaa by geometric means: within the molecules the LJ is Szczerba and Kalinichev's, not OPLS-AA's);
- the 1-4 pairs are scaled by 0.5 (special_bonds of oplsaa; ClayFF has no 1-4 pair: no path of three of its bonds).

    python check_ff.py         (soapff)  ->  results/check_ff.json"""
import json
import os
import tempfile

import numpy as np
from ase.io import read

import soapff

HERE = os.path.dirname(os.path.abspath(__file__))
SECTIONS = {"Bonds": "bond", "Angles": "angle", "Dihedrals": "dihedral", "Impropers": "improper"}


def read_ff(d):
    """system.data and system.settings of an output folder -> type name and charge of every atom, LJ of every type name,
    {kind: sorted [(atoms, coefficients)]}, the pair_modify and special_bonds lines"""
    names, atoms, rows = {}, {}, {k: [] for k in SECTIONS.values()}
    sec = None
    for line in open(os.path.join(d, "system.data")):
        c = line.split("#")[0].split()
        if not c:
            continue
        if c[0] in ("Masses", "Atoms", *SECTIONS):
            sec = c[0]
            continue
        if sec == "Masses":
            names[int(c[0])] = line.split("#")[1].strip()
        elif sec == "Atoms":
            atoms[int(c[0])] = (names[int(c[2])], float(c[3]))
        elif sec in SECTIONS:
            rows[SECTIONS[sec]].append((int(c[1]), tuple(int(x) - 1 for x in c[2:])))
    coef = {k: {} for k in SECTIONS.values()}
    lj, lines = {}, {}
    for line in open(os.path.join(d, "system.settings")):
        c = line.split("#")[0].split()
        if not c:
            continue
        if c[0] == "pair_coeff":
            lj[names[int(c[1])]] = (float(c[3]), float(c[4]))
        elif c[0].endswith("_coeff") and c[0][:-6] in coef:
            coef[c[0][:-6]][int(c[1])] = tuple(float(x) for x in c[2:])
        elif c[0] in ("pair_modify", "special_bonds", "pair_style"):
            lines[c[0]] = " ".join(c[1:])
    ids = sorted(atoms)
    types = [atoms[i][0] for i in ids]
    q = np.array([atoms[i][1] for i in ids])
    terms = {}
    for k, r in rows.items():
        out = []
        for t, a in r:
            if k in ("bond", "angle", "dihedral") and a[0] > a[-1]:
                a = a[::-1]                     # the same term read from the other end
            out.append((a, coef[k][t]))
        terms[k] = sorted(out)
    return types, q, lj, terms, lines


def same(a, b, rel=1e-8):
    return len(a) == len(b) and all(abs(x - y) <= rel * max(1.0, abs(x), abs(y)) for x, y in zip(a, b))


def main(name):
    at = read(os.path.join(HERE, "structures", f"{name}_start.extxyz"))
    opls = set(soapff.library.load("oplsaa").types)
    part = {"oplsaa": np.nonzero(np.isin(at.arrays["fftype"], list(opls)))[0]}
    part["clayff"] = np.setdiff1d(np.arange(len(at)), part["oplsaa"])
    with tempfile.TemporaryDirectory() as d:
        soapff.build("clayff_oplsaa", at, os.path.join(d, "one"), log=lambda *a: None)
        T, Q, LJ, TERMS, LINES = read_ff(os.path.join(d, "one"))
        want_terms = {k: [] for k in TERMS}
        res = dict(structure=name, atoms=len(at), settings=LINES, parts={})
        for lib, idx in part.items():
            soapff.build(lib, at[idx], os.path.join(d, lib), log=lambda *a: None)
            t, q, lj, terms, lines = read_ff(os.path.join(d, lib))
            for k, r in terms.items():
                want_terms[k] += [(tuple(int(idx[i]) for i in a), c) for a, c in r]
            res["parts"][lib] = dict(
                atoms=len(idx), settings=lines,
                types_differ=int(sum(T[g] != x for g, x in zip(idx, t))),
                charges_max_difference=float(np.abs(Q[idx] - q).max()),
                lj_differ=sorted(x for x in lj if not same(LJ[x], lj[x])),
                terms={k: len(r) for k, r in terms.items()})
        res["terms"] = {}
        for k in TERMS:
            got = sorted((tuple(a[::-1]) if k != "improper" and a[0] > a[-1] else a, c) for a, c in TERMS[k])
            want = sorted((tuple(a[::-1]) if k != "improper" and a[0] > a[-1] else a, c) for a, c in want_terms[k])
            res["terms"][k] = dict(clayff_oplsaa=len(got), parts=len(want),
                                   same=len(got) == len(want) and all(x[0] == y[0] and same(x[1], y[1])
                                                                      for x, y in zip(got, want)))
    return res


if __name__ == "__main__":
    out = {n: main(n) for n in ("PEO", "PVA")}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "results", "check_ff.json"), "w"), indent=1)
    for n, r in out.items():
        print(n, r["settings"])
        for lib, p in r["parts"].items():
            print(" ", lib, {k: v for k, v in p.items() if k != "settings"}, p["settings"])
        print("  terms", r["terms"])
