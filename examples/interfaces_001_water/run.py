"""Water on the basal (0 0 1) slabs (paper: test (4) and Fig. "interfaces"): the library slab of each of the 10 layered
minerals of the figure (the references of the clayff library, src/soapff/libraries/clayff/basal/; the second pyrophyllite
source is left out, as in the paper) with its
widest empty slice filled with water by packmol at 1 g/cm3 (soapff_cui.water: tolerance 2 A, seed 1),
the force field of soapff (bundled clayff library, flexible SPC water) and a LAMMPS minimisation with the LAMMPS Python module
(../clay_md.py).  The minimisation uses Ewald (1e-4) instead of the PPPM of system.settings, as in the paper (PPPM stopped at
the start of the minimisation on the LAMMPS used there).

    python run.py            (with a Python where soapff and LAMMPS are installed: pip install "lammps[mpi]"; packmol)
                              ->  <mineral>/ (structure.extxyz; soapff's system.data, system.settings, charges.csv,
                                  report.json; the minimisation: its commands in.min, its log log.min), results.json"""
import json
import os
import shutil
import sys

import ase.io
from lammps import lammps

import soapff
from soapff_cui import core, water

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import clay_md  # noqa: E402  (examples/clay_md.py: the LAMMPS commands of the clay examples)

REFS = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "basal")
LIB = os.path.join(HERE, "..", "..", "src", "soapff", "libraries", "clayff", "prebuilt", "library_clayff_bulk.npz")
EWALD = "kspace_style ewald 1.0e-4"          # this test: Ewald 1e-4 instead of the PPPM of system.settings
MINERALS = ["montmorillonite_dry_Na", "muscovite_2M1", "muscovite_1M", "phlogopite_1M", "nanpingite", "pyrophyllite",
            "kaolinite", "lizerdite", "talc", "brucite"]


def minimize(d):
    """the two-stage minimisation (clay_md.minimize) of the force field in `d` with Ewald 1e-4 after system.settings, with the
    LAMMPS Python module (log: log.min); the same commands are written to in.min (lmp -in in.min in the folder).  Returns the
    stopping criteria and the energies"""
    s = soapff.System(d)
    cmds = [EWALD] + clay_md.minimize()
    clay_md.input_file(s, cmds, "in.min")
    lmp = lammps(cmdargs=["-log", os.path.join(d, "log.min"), "-screen", "none"])
    lmp.commands_list(s.commands() + cmds)
    lmp.close()
    L = open(os.path.join(d, "log.min")).read().splitlines()
    k = [i for i, l in enumerate(L) if "Energy initial, next-to-last, final" in l]
    return dict(stop=[l.split("=")[1].strip() for l in L if "Stopping criterion" in l],
                energy_initial_final=[float(L[k[0] + 1].split()[0]), float(L[k[-1] + 1].split()[2])] if k else None)


if __name__ == "__main__":
    res = {}
    for m in MINERALS:
        d = os.path.join(HERE, m)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        at, q, types = core.read_lammps_data(os.path.join(REFS, m + "_001.data"))
        gap = water.gaps(at, min_gap=8.0, density=1.0)[0]                       # the vacuum above the slab
        s, rep = water.add_water(at, [(gap, gap["n_water"])], log=lambda *a: None)
        ase.io.write(os.path.join(d, "structure.extxyz"), s, format="extxyz")
        r = core.make(os.path.join(d, "structure.extxyz"), d, LIB, strict=True, log=lambda *a: None)
        res[m] = dict(slab_atoms=len(at), water=rep[0]["n"], atoms=r["atoms"], packmol_converged=rep[0]["converged"],
                      nearest_atom=rep[0]["min_distance"], type_check_differ=r["type_check"]["differ"],
                      min_similarity=round(r["min_similarity"], 4),
                      charge_sum=round(r["charge_sum"], 4), minimization=minimize(d))
        print(m, res[m], flush=True)
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
