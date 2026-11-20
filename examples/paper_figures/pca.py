"""figs/pca.png: the SOAP vectors of the O rows of the two tables on their first two principal components, coloured by
type (the type is not used to find the axes): (a) clayff, the step-3 rows of Na montmorillonite; (b) oplsaa, the O rows of
the 88 references of the 11 polymers.  Numbers in figdata/pca/numbers.json.

    PYTHONPATH=src python examples/paper_figures/pca.py      (Fig. 3 of the paper)"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import soapff  # noqa: E402
from soapff_cui import core, steps  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "figdata", "pca")
os.makedirs(OUTD, exist_ok=True)
os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)


def pca(X):
    Xc = X - X.mean(0)
    _, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return Xc @ Vt[:2].T, S ** 2 / (S ** 2).sum()


num = {}
lib = soapff.load_library("clayff")
T = steps.load_table(lib.table_path(log=lambda *a: None))
org = np.array([str(x) for x in T["s3_origin"]])
o = (T["s3_element"] == "O") & (np.char.find(org, "montmorillonite") >= 0)
P1, v1 = pca(T["s3_X"][o].astype(np.float64))
t1 = np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in T["s3_type"][o]])
num["clayff"] = dict(o_rows=int(o.sum()), explained=[round(float(v), 3) for v in v1[:2]],
                     rows_per_type={t: int((t1 == t).sum()) for t in sorted(set(t1))})

lib = soapff.load_library("oplsaa")
T = np.load(lib.table_path(log=lambda *a: None), allow_pickle=True)
o = T["element"].astype(str) == "O"
P2, v2 = pca(T["X"][o].astype(np.float64))
t2 = T["type"][o].astype(str)
num["oplsaa"] = dict(o_rows=int(o.sum()), explained=[round(float(v), 3) for v in v2[:2]],
                     rows_per_type={t: int((t2 == t).sum()) for t in sorted(set(t2))})
json.dump(num, open(os.path.join(OUTD, "numbers.json"), "w"), indent=1)
print(json.dumps(num))

plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})
fig, AX = plt.subplots(1, 2, figsize=(10.0, 4.0), gridspec_kw=dict(wspace=0.28))
PANELS = [(AX[0], P1, v1, t1, [("ob", "#0072B2", "ob"), ("obos", "#56B4E9", "obos"), ("obts", "#009E73", "obts"),
                               ("oh", "#E69F00", "oh"), ("ohs", "#D55E00", "ohs")],
           f"clayff: the {len(t1)} O rows of Na montmorillonite", "a"),
          (AX[1], P2, v2, t2, [("opls_154", "#E69F00", "opls_154 alcohol O"), ("opls_180", "#0072B2", "opls_180 ether O"),
                               ("opls_466", "#009E73", "opls_466 ester C=O"), ("opls_467", "#D55E00", "opls_467 ester O–R")],
           f"oplsaa: the {len(t2)} O rows (the 4 polymers with O)", "b")]
for ax, P, var, ty, styles, title, lab in PANELS:
    for t, c, name in styles:
        k = ty == t
        if k.any():
            ax.scatter(P[k, 0], P[k, 1], s=22, c=c, edgecolors="white", linewidths=0.4, label=f"{name} ({int(k.sum())})")
    ax.set_xlabel(f"PC 1 ({100 * var[0]:.0f}%)")
    ax.set_ylabel(f"PC 2 ({100 * var[1]:.0f}%)")
    ax.legend(fontsize=7.5, frameon=False, loc="best", markerscale=1.3)
    ax.text(-0.18, 1.03, lab, transform=ax.transAxes, weight="bold", fontsize=11)
    ax.set_title(title, fontsize=9)
fig.savefig(os.path.join(HERE, "figs", "pca.png"), dpi=250, bbox_inches="tight")
print("figs/pca.png")
