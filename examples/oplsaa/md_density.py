"""The melt densities at 413 K and 1 bar of the NPT MD with the files soapff writes (LAMMPS, md_lammps.py) and with Foyer's
OPLS-AA (OpenMM, md_openmm.py), and the measured densities of Fetters et al. (1994).

    python md_density.py WORK SKIP_PS results/md_density.json P1 [P2 ...]

The mean of the density after SKIP_PS (ps) and its standard error from block averages (5 blocks).  The measured values are
those of L. J. Fetters, D. J. Lohse, D. Richter, T. A. Witten, A. Zirkel, Macromolecules 27, 4639 (1994) at 413 K (as
listed in Table II of J. Luettmer-Strathmann, arXiv:cond-mat/9911250): high polymers, while the MD chains are
short (PE C100, PP 50 and PIB 40 repeat units), so the MD density is expected a little below them (more chain ends)."""
import json
import os
import sys

import numpy as np

W, SKIP, OUT = sys.argv[1], float(sys.argv[2]), sys.argv[3]
FETTERS = {"PE": 0.785, "PP": 0.791, "PIB": 0.849}          # g/cm3, 413 K
res = {}
for p in sys.argv[4:]:
    r = dict(fetters_1994=FETTERS.get(p))
    for code in ("lammps", "openmm"):
        f = os.path.join(W, f"{code}_{p}.csv")
        if not os.path.isfile(f):
            continue
        X = np.loadtxt(f, delimiter=",", skiprows=1, ndmin=2)
        x = X[X[:, 0] > SKIP]
        if len(x) < 10:
            continue
        b = np.array([v.mean() for v in np.array_split(x[:, 1], 5)])
        r[code] = dict(density=round(float(x[:, 1].mean()), 4), stderr=round(float(b.std(ddof=1) / np.sqrt(len(b))), 4),
                       temperature=round(float(x[:, 2].mean()), 1), from_ps=float(x[0, 0]), to_ps=float(x[-1, 0]))
    res[p] = r
    print(p, json.dumps(r))
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=1)
