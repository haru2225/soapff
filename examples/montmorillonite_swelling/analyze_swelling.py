"""Swelling curve from run/swelling.dat (written by the swelling curve of swelling.py): for each water content,
d001 = lz averaged over the last 20 ps of its NPT (and the standard deviation), against the water content in g H2O / g clay
(clay = every atom that is not water, Na included); swelling and deswelling separately; an experimental curve can be drawn
with it from swelling_exp.csv in the folder (water_g_per_g, d001_A[, label]).
    python analyze_swelling.py [folder]     (default run)  ->  folder/swelling.csv, swelling.png"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def swelling(out_dir, dest):
    """swelling.dat (step added_waters phase lx ly lz xy xz yz pe temp, every 100 fs) -> one row per state (the
    consecutive lines with the same number of added waters; phase 1 = averaged): water molecules, g H2O / g clay, d001
    (lz) mean and standard deviation, cell a b c alpha beta gamma, energy; 'up' until the largest water content, then
    'down'.  Writes swelling.csv and swelling.png (with swelling_exp.csv if present)."""
    a = np.loadtxt(os.path.join(out_dir, "swelling.dat"), comments="#", ndmin=2)
    if not len(a):
        raise ValueError("swelling.dat is empty")
    rows = list(csv.DictReader(open(os.path.join(out_dir, "charges.csv"))))
    from soapff_cui._clayff.chem import MASS
    ty = [r.get("type", r.get("clayff_type")) for r in rows]          # clayff_type: charges.csv of soapff 1.0 (run/)
    clay_mass = sum(MASS[r["element"]] for r, t in zip(rows, ty) if t not in ("ow", "hw"))
    w0 = ty.count("ow")
    m_w = 2 * MASS["H"] + MASS["O"]
    states, start = [], 0
    nmax = a[:, 1].max()
    reached = False
    for k in range(1, len(a) + 1):
        if k == len(a) or a[k, 1] != a[start, 1]:
            blk = a[start:k]
            prod = blk[blk[:, 2] == 1]
            if len(prod):
                lx, ly, lz, xy, xz, yz = (prod[:, j] for j in range(3, 9))
                A = lx; B = np.sqrt(ly ** 2 + xy ** 2); C = np.sqrt(lz ** 2 + xz ** 2 + yz ** 2)
                al = np.degrees(np.arccos((xy * xz + ly * yz) / (B * C)))
                be = np.degrees(np.arccos(xz / C)); ga = np.degrees(np.arccos(xy / B))
                n = int(round(blk[0, 1])) + w0
                states.append(dict(waters=n, water_g_per_g=n * m_w / clay_mass, direction="down" if reached else "up",
                                   d001=float(lz.mean()), d001_sd=float(lz.std()), a=float(A.mean()), b=float(B.mean()),
                                   c=float(C.mean()), alpha=float(al.mean()), beta=float(be.mean()), gamma=float(ga.mean()),
                                   pe=float(prod[:, 9].mean()), T=float(prod[:, 10].mean()), frames=int(len(prod))))
            if a[start, 1] == nmax:
                reached = True
            start = k
    dt = np.median(np.diff(a[:, 0])) if len(a) > 1 else 1
    avg_steps = np.median([s_["frames"] for s_ in states]) * dt
    with open(os.path.join(dest, "swelling.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(states[0]))
        w.writeheader(); w.writerows(states)
    exp = []
    ef = os.path.join(out_dir, "swelling_exp.csv")
    if os.path.exists(ef):
        exp = list(csv.DictReader(open(ef)))
    plot_swelling(dest, states, exp, clay_mass)
    dt_fs = _time_step_fs(out_dir, 0.5)
    return dict(states=states, clay_mass_amu=clay_mass, water0=w0, avg_ps=float(avg_steps * dt_fs / 1000.0),
                experiment=os.path.basename(ef) if exp else None)


def plot_swelling(dest, states, exp, clay_mass):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    up = [r for r in states if r["direction"] == "up"]
    dn = [r for r in states if r["direction"] == "down"]
    ax.plot([r["water_g_per_g"] for r in up], [r["d001"] for r in up], "o-", color="#c0392b", ms=4, lw=1, label="MD (swelling)")
    if dn:
        ax.plot([r["water_g_per_g"] for r in dn], [r["d001"] for r in dn], "s--", color="#2c6fbb", mfc="none", ms=4, lw=1,
                label="MD (deswelling)")
    if exp:
        labs = sorted({e.get("label") or "experiment" for e in exp})
        for lab in labs:
            pts = [e for e in exp if (e.get("label") or "experiment") == lab]
            ax.plot([float(e["water_g_per_g"]) for e in pts], [float(e["d001_A"]) for e in pts], "^", color="#2a9d4b", ms=6,
                    label=lab)
    ax.set_xlabel("water content (g H$_2$O / g clay)")
    ax.set_ylabel("d$_{001}$ (Å)")
    sec = ax.secondary_xaxis("top", functions=(lambda x: x * clay_mass / 18.0153, lambda n: n * 18.0153 / clay_mass))
    sec.set_xlabel("H$_2$O per cell")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(dest, "swelling.png"), dpi=200)
    plt.close(fig)


def _time_step_fs(out_dir, default):
    """the time step of the run: the last `timestep` command of in.lammps (an input file), else of log.lammps (LAMMPS
    writes there every command, also those of the Python module), else the default."""
    for name in ("in.lammps", "log.lammps"):
        dt = None
        try:
            for line in open(os.path.join(out_dir, name)):
                w = line.split("#")[0].split()
                if len(w) == 2 and w[0] == "timestep":
                    try:
                        dt = float(w[1])
                    except ValueError:
                        pass
        except OSError:
            continue
        if dt is not None:
            return dt
    return default


if __name__ == "__main__":
    out = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(HERE, "run")
    r = swelling(out, out)
    for s in r["states"]:
        print(f"{s['waters']:3d} waters  {s['water_g_per_g']:.3f} g/g  d001 {s['d001']:.2f} +- {s['d001_sd']:.2f} A")
