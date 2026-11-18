"""LAMMPS forces (kcal/mol/A) on every atom with soapff's force field (library oplsaa), compared with the forces of
Foyer's OPLS-AA in OpenMM (forces_openmm.py).  LAMMPS is run from Python (the LAMMPS Python module).

    python forces_lammps.py TESTDIR openmm.npz WORK results/forces.json          (soapff, LAMMPS Python module)

- chains in vacuum: the system.data and the coefficients of system.settings with lj/cut/coul/cut 999 A and boundary f f f
  (no cutoff, as OpenMM's NoCutoff);
- amorphous cells: the force field as written (lj/cut/coul/long 12 A, geometric mixing, special_bonds lj/coul 0 0 0.5) with
  Ewald 1e-8 (OpenMM: PME 1e-8).
The forces are read from LAMMPS after run 0.  results/forces.json: per structure the number of atoms, the RMS of the
OpenMM force components, the RMS and the largest difference, and the correlation; WORK/forces.npz: both forces of every atom."""
import json
import os
import sys

import numpy as np
from lammps import lammps

import soapff

T, OMM, W, OUT = sys.argv[1:5]
omm = np.load(OMM)
quiet = dict(log=lambda *a: None)


def forces(s, commands):
    """the force on every atom (in the order of the atom ids) after run 0"""
    lmp = lammps(cmdargs=["-log", os.path.join(s.out_dir, "log.forces"), "-screen", "none"])
    lmp.commands_list(commands + ["thermo_style custom step pe ebond eangle edihed evdwl ecoul elong",
                                  "thermo_modify format float %.6f", "run 0"])
    n = lmp.extract_global("nlocal")                         # the arrays also hold the ghost atoms after that
    ids = np.array(lmp.numpy.extract_atom("id"))[:n].copy()
    f = np.array(lmp.numpy.extract_atom("f"))[:n].copy()
    lmp.close()
    return f[np.argsort(ids)]


res, pairs = {}, {}
for name in sorted(omm.files):
    f = os.path.join(T, name + ".extxyz")
    if name.startswith("cell_"):
        s = soapff.build("oplsaa", f, os.path.join(W, name), kspace="ewald 1.0e-8", **quiet)
        fl = forces(s, s.commands())
    else:
        s = soapff.build("oplsaa", f, os.path.join(W, name), **quiet)
        S = s.commands()
        first = lambda l: l.split()[0]                                          # noqa: E731
        styles = [l for l in S if first(l) in ("bond_style", "angle_style", "dihedral_style", "improper_style",
                                               "special_bonds")]
        coeffs = [l for l in S if first(l).endswith("_coeff")]
        mix = [l for l in S if first(l) == "pair_modify"][0].replace(" tail yes", "")
        data = [l for l in S if first(l) == "read_data"]
        fl = forces(s, ["units real", "atom_style full", "boundary f f f", "pair_style lj/cut/coul/cut 999.0", mix] + styles
                    + data + coeffs)
    fo = omm[name]
    assert fl.shape == fo.shape, name
    d = fl - fo
    k = int(np.argmax(np.linalg.norm(d, axis=1)))
    res[name] = dict(atoms=len(fo), rms_force=float(np.sqrt((fo ** 2).mean())), rms_diff=float(np.sqrt((d ** 2).mean())),
                     max_diff=float(np.abs(d).max()), max_diff_atom=k + 1, max_diff_type=s.types[k],
                     correlation=float(np.corrcoef(fl.ravel(), fo.ravel())[0, 1]),
                     dihedrals_without_parameters=s.report.get("dihedrals_without_parameters"))
    pairs[name + "_lammps"], pairs[name + "_openmm"] = fl.astype(np.float32), fo.astype(np.float32)
    r = res[name]
    print(f"{name:16s} {r['atoms']:6d} atoms  RMS F {r['rms_force']:7.3f}  RMS dF {r['rms_diff']:.2e}  max dF "
          f"{r['max_diff']:.2e} (atom {k + 1}, {s.types[k]})  r {r['correlation']:.12f}", flush=True)

os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=1)
np.savez_compressed(os.path.join(W, "forces.npz"), **pairs)
