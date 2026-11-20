"""Water on the basal (0 0 1) slabs (paper: test (4) and Fig. "interfaces"): the library slab of each of the 10 layered
minerals of the figure (the references of the clayff library, src/soapff/libraries/clayff/basal/; the second pyrophyllite
source is left out, as in the paper) with its
widest empty slice filled with water by packmol at 1 g/cm3 (soapff_cui.water: tolerance 2 A, seed 1; the structure.extxyz
kept is used again), the force field of soapff (soapff.build("clayff", ...)) compared with the right answer (the ClayFF types
of the bond-graph rule and Table 1 charges, soapff_cui.core.clayff_answer) and a LAMMPS minimisation with the LAMMPS Python
module (../clay_md.py).  The minimisation uses Ewald (1e-4) instead of the PPPM of system.settings, as in the paper (PPPM stopped at
the start of the minimisation on the LAMMPS used there).

    python run.py            (soapff and the LAMMPS Python module, pip install "lammps[mpi]", or LMP=/path/to/lmp; packmol)
                              ->  <mineral>/ (structure.extxyz; soapff's system.data, system.settings, charges.csv,
                                  report.json; the minimisation: its commands in.min, its log log.min), results.json"""
import json
import os
import shutil
import sys

import ase.io
import numpy as np

import soapff
from soapff_cui import core, water

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import clay_md  # noqa: E402  (examples/clay_md.py: the LAMMPS commands of the clay examples)

REFS = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "basal")
EWALD = "kspace_style ewald 1.0e-4"          # this test: Ewald 1e-4 instead of the PPPM of system.settings
MINERALS = ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
            "kaolinite", "lizerdite", "talc", "brucite"]


def minimize(d):
    """the two-stage minimisation (clay_md.minimize) of the force field in `d` with Ewald 1e-4 after system.settings, with the
    LAMMPS Python module (log: log.min); the same commands are written to in.min (lmp -in in.min in the folder).  Returns the
    stopping criteria and the energies"""
    s = soapff.System(d)
    L = open(clay_md.run(s, [EWALD] + clay_md.minimize(), "in.min", "log.min")).read().splitlines()
    k = [i for i, l in enumerate(L) if "Energy initial, next-to-last, final" in l]
    return dict(stop=[l.split("=")[1].strip() for l in L if "Stopping criterion" in l],
                energy_initial_final=[float(L[k[0] + 1].split()[0]), float(L[k[-1] + 1].split()[2])] if k else None)


if __name__ == "__main__":
    res = {}
    old = json.load(open(os.path.join(HERE, "results.json"))) if os.path.isfile(os.path.join(HERE, "results.json")) else {}
    for m in MINERALS:
        d = os.path.join(HERE, m)
        f = os.path.join(d, "structure.extxyz")
        at, q, types = core.read_lammps_data(os.path.join(REFS, m + "_001.data"))
        packed = {k: old.get(m, {}).get(k) for k in ("packmol_converged", "nearest_atom")}   # of the packmol run kept
        if not os.path.isfile(f):                                               # the structure: made once by packmol
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d)
            gap = water.gaps(at, min_gap=8.0, density=1.0)[0]                   # the vacuum above the slab
            st, rep = water.add_water(at, [(gap, gap["n_water"])], log=lambda *a: None)
            ase.io.write(f, st, format="extxyz")
            packed = dict(packmol_converged=rep[0]["converged"], nearest_atom=rep[0]["min_distance"])
        st = core.read_structure(f)
        s = soapff.build("clayff", f, d, log=lambda *a: None)
        t0, q0 = core.clayff_answer(st)
        res[m] = dict(slab_atoms=len(at), water=int((np.array(st.get_chemical_symbols()) == "O").sum()
                                                     - (np.array(at.get_chemical_symbols()) == "O").sum()),
                      atoms=len(st), **packed, **core.score(s.types, s.charges, t0, q0),
                      min_similarity=round(s.report["min_similarity_soap"], 4),
                      charge_sum=round(s.report["charge_sum"], 4), minimization=minimize(d))
        print(m, res[m], flush=True)
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
