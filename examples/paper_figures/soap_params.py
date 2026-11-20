"""How the typing depends on the SOAP settings: r_cut, n_max and l_max changed one at a time from those of the libraries
(r_cut 5 A, n_max 4, l_max 3; sigma 0.4 and the depth of 4 bonds kept).  For each setting, a library that extends the
bundled one with this [soap] (its SOAP table computed again), and the static test structures (no MD frame) typed with it
(strict=False, so that every atom gets the type of its most similar row):

  clayff       the 10 (001) interfaces with water (examples/interfaces_001_water) and the 6 water shapes
               (examples/interface_shapes); the answer: the types and charges of the default setting, which are those of the
               rule of ClayFF
  clayff_edge  the particle with edges, the particle twice and the particle with water and ions (examples/montmorillonite_edge);
               the answer: the default setting, which gives the types of the particle's files
  oplsaa       the 18 test structures of examples/oplsaa/test that are not MD frames (the chains, the copolymers and the
               amorphous cells); the answer: the types and charges of Foyer (fftype, ffcharge)

For each: the atoms with another type or charge than the answer, the atoms below the threshold 0.9 (soapff.build would stop),
the lowest similarity and the smallest margin (the similarity of the chosen row minus that of the most similar row with
another type or charge).

    SOAPFF_CACHE=WORK/cache PYTHONPATH=src python examples/paper_figures/soap_params.py run GROUP WORK   (GROUP: clayff, clayff_edge, oplsaa)
    python examples/paper_figures/soap_params.py summary WORK      -> figdata/soap_params/results.json
    python examples/paper_figures/soap_params.py figure          -> figs/soap_params.png (Fig. 6 of the paper)"""
import glob
import json
import os
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..")                       # the examples of the repository
BASE = dict(r_cut=5.0, n_max=4, l_max=3)
SCAN = dict(r_cut=[3.0, 4.0, 6.0, 8.0], n_max=[2, 3, 6, 8], l_max=[0, 1, 2, 4, 6])
STOPS = ("atoms are not in the library", "molecules have a charge that is not a whole number")


def settings():
    yield "default", dict(BASE)
    for k, vals in SCAN.items():
        for v in vals:
            yield f"{k}={v:g}", dict(BASE, **{k: v})


def structures(group):
    """(name, path or Atoms, answer types or None, answer charges or None)"""
    import ase.io
    if group == "clayff":
        for f in sorted(glob.glob(os.path.join(EX, "interfaces_001_water", "*", "structure.extxyz"))):
            yield "if_" + os.path.basename(os.path.dirname(f)), f, None, None
        for f in sorted(glob.glob(os.path.join(EX, "interface_shapes", "*", "structure.extxyz"))):
            yield "shape_" + os.path.basename(os.path.dirname(f)), f, None, None
    elif group == "clayff_edge":
        for n in ("clay_montmorillonite", "clay_montmorillonite_2", "particle_water"):
            yield n, os.path.join(EX, "montmorillonite_edge", n + ".extxyz"), None, None
    else:
        for f in sorted(glob.glob(os.path.join(EX, "oplsaa", "test", "*.extxyz"))):
            if os.path.basename(f).startswith("md"):            # the MD frames: left out
                continue
            at = ase.io.read(f)
            yield os.path.basename(f)[:-7], f, np.array(at.arrays["fftype"]), np.array(at.arrays["ffcharge"], float)


def library(group, name, soap, work):
    if name == "default":
        return group
    import soapff
    s = dict(soapff.load_library(group).soap, **soap)
    tag = name.replace("=", "")
    path = os.path.join(work, "libs", f"{group}_{tag}.txt")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(f"[library]\nname = {group}_{tag}\nextends = {group}\n\n[soap]\n"
                          + "".join(f"{k} = {v}\n" for k, v in s.items()))
    return path


def run(group, work):
    import soapff
    quiet = lambda *a: None                                                                  # noqa: E731
    out = os.path.join(work, "out", group)
    os.makedirs(out, exist_ok=True)
    tests = list(structures(group))
    answer = {}
    for name, soap in settings():
        f_json = os.path.join(out, name.replace("=", "") + ".json")
        if os.path.isfile(f_json):
            res = json.load(open(f_json))
            if name == "default":
                answer = {k: (np.array(v["types"]), np.array(v["charges"])) for k, v in res.pop("_atoms").items()}
            continue
        lib = library(group, name, soap, work)
        t0 = time.time()
        soapff.load_library(lib).table_path(log=quiet)
        res = dict(setting=name, soap=soap, table_seconds=round(time.time() - t0, 1), structures={})
        atoms_out = {}
        t1 = time.time()
        for sname, src, ft, fq in tests:
            d = os.path.join(work, "tmp", group)
            shutil.rmtree(d, ignore_errors=True)
            try:
                s = soapff.build(lib, src, d, strict=False, log=quiet)
            except soapff.TypingError as e:            # an atom without a type, or a force-field term without parameters
                n = len(ft) if ft is not None else len(answer[sname][0])
                res["structures"][sname] = dict(atoms=n, wrong_type=0, wrong_charge=0, below_threshold=0,
                                                min_similarity=np.nan, min_margin=None, other_stops=[str(e)[:300]],
                                                stopped=True)
                continue
            types, q = np.array(s.types), np.array(s.charges, float)
            rows = s.table()
            simk = "similarity_soap" if "similarity_soap" in rows[0] else "similarity"
            mk = "margin_soap" if "margin_soap" in rows[0] else "margin"
            sim = np.array([float(r[simk]) if r[simk] not in ("", None) else np.nan for r in rows])
            mar = np.array([float(r[mk]) if r.get(mk) not in ("", None) else np.nan for r in rows])
            if ft is None:
                if name == "default":
                    ft, fq = types, q
                else:
                    ft, fq = answer[sname]
            atoms_out[sname] = dict(types=[str(t) for t in types], charges=[float(x) for x in q])
            wt = types != ft
            wq = np.abs(q - fq) > 1e-6
            warn = [w for w in s.report.get("warnings", []) if any(x in w for x in STOPS)]
            r = dict(atoms=len(types), wrong_type=int(wt.sum()), wrong_charge=int(wq.sum()),
                     below_threshold=int(np.nansum(sim < 0.9)), min_similarity=float(np.nanmin(sim)),
                     min_margin=float(np.nanmin(mar)) if np.isfinite(mar).any() else None,
                     other_stops=warn[:1])
            if wt.any():
                pairs = {}
                for a, b in zip(ft[wt], types[wt]):
                    pairs[f"{a}->{b}"] = pairs.get(f"{a}->{b}", 0) + 1
                r["confused"] = pairs
            res["structures"][sname] = r
        res["typing_seconds"] = round(time.time() - t1, 1)
        if name == "default":
            answer = {k: (np.array(v["types"]), np.array(v["charges"])) for k, v in atoms_out.items()}
            res["_atoms"] = atoms_out
        tot = {k: sum(v[k] for v in res["structures"].values()) for k in ("atoms", "wrong_type", "wrong_charge",
                                                                          "below_threshold")}
        print(group, name, tot, f"table {res['table_seconds']} s, typing {res['typing_seconds']} s", flush=True)
        json.dump(res, open(f_json, "w"))


def summary(work):
    out = {}
    for group in ("clayff", "clayff_edge", "oplsaa"):
        for name, soap in settings():
            f = os.path.join(work, "out", group, name.replace("=", "") + ".json")
            if not os.path.isfile(f):
                continue
            r = json.load(open(f))
            S = {k: v for k, v in r["structures"].items() if not k.startswith("md")}     # static structures only
            row = dict(soap=soap, structures=len(S), atoms=sum(v["atoms"] for v in S.values()),
                       wrong_type=sum(v["wrong_type"] for v in S.values()),
                       wrong_charge=sum(v["wrong_charge"] for v in S.values()),
                       below_threshold=sum(v["below_threshold"] for v in S.values()),
                       structures_that_stop=sorted(k for k, v in S.items() if v["below_threshold"] or v["other_stops"]),
                       structures_stopped_by_bond_check=sorted(k for k, v in S.items() if v.get("stopped")),
                       min_similarity=round(float(np.nanmin([v["min_similarity"] for v in S.values()])), 4),
                       min_margin=round(min(v["min_margin"] for v in S.values() if v["min_margin"] is not None), 4),
                       table_seconds=r["table_seconds"], typing_seconds=r["typing_seconds"])
            conf = {}
            for v in S.values():
                for k, n in v.get("confused", {}).items():
                    conf[k] = conf.get(k, 0) + n
            if conf:
                row["confused"] = conf
            out.setdefault(group, {})[name] = row
    d = os.path.join(HERE, "figdata", "soap_params")
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "results.json"), "w"), indent=1)
    for g, rows in out.items():
        for n, r in rows.items():
            print(f"{g:12s} {n:10s} atoms {r['atoms']:6d} wrong type {r['wrong_type']:5d} charge {r['wrong_charge']:5d} "
                  f"below {r['below_threshold']:5d} stop {len(r['structures_that_stop']):2d} min sim {r['min_similarity']:.4f} "
                  f"min margin {r['min_margin']:+.4f}  {r.get('confused', '')}")


def figure():
    """figs/soap_params.png: against each SOAP setting, the atoms typed otherwise than the answer (top) and the structures on
    which soapff.build stops (bottom: an atom below the threshold, bonds that do not fit the type, or a molecule whose charge is
    not a whole number)"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "Hiragino Sans", "font.size": 9})
    R = json.load(open(os.path.join(HERE, "figdata", "soap_params", "results.json")))
    # categorical slots 1-3 of the dataviz palette (checked with its validate_palette.js), with a marker each
    style = {"clayff": ("#2a78d6", "o"), "clayff_edge": ("#eb6834", "s"), "oplsaa": ("#1baf7a", "^")}
    label = {"clayff": "clayff", "clayff_edge": "clayff_edge", "oplsaa": "oplsaa"}
    xlab = {"r_cut": r"$r_\mathrm{c}$ (Å)", "n_max": r"$n_{\max}$", "l_max": r"$l_{\max}$"}
    fig, AX = plt.subplots(2, 3, figsize=(7.2, 4.4), sharex="col", sharey="row")
    for col, key in enumerate(("r_cut", "n_max", "l_max")):
        xs = sorted(set([BASE[key]] + SCAN[key]))
        span = max(xs) - min(xs)
        for row, what in enumerate(("wrong_type", "stop")):
            ax = AX[row, col]
            ax.axvline(BASE[key], color="#bbbbbb", lw=1.0, ls="--", zorder=0)
            for k, g in enumerate(("clayff", "clayff_edge", "oplsaa")):
                rows = R[g]
                x, y = [], []
                for v in xs:
                    name = "default" if v == BASE[key] else f"{key}={v:g}"
                    if name in rows:
                        x.append(v + (k - 1) * 0.022 * span)
                        y.append(rows[name]["wrong_type"] if what == "wrong_type" else len(rows[name]["structures_that_stop"]))
                c, m = style[g]
                ax.plot(x, y, color=c, lw=1.4, marker=m, ms=5.5, mec="white", mew=0.6, zorder=3,
                        label=f"{label[g]}\n{rows['default']['structures']} structures, {rows['default']['atoms']:,} atoms")
            ax.grid(axis="y", alpha=0.25, lw=0.6)
            for s_ in ("top", "right"):
                ax.spines[s_].set_visible(False)
            if row == 0:
                ax.set_yscale("symlog", linthresh=1)
                ax.set_ylim(-0.3, 40)
                ax.set_yticks([0, 1, 10])
                ax.set_yticklabels(["0", "1", "10"])
            else:
                ax.set_ylim(-0.2, 2.2)
                ax.set_yticks([0, 1, 2])
                ax.set_xticks(xs)
                ax.set_xlabel(xlab[key])
    AX[0, 0].set_ylabel("atoms with\nanother type")
    AX[1, 0].set_ylabel("structures on which\nsoapff stops")
    AX[0, 1].text(BASE["n_max"], 28, "the libraries", ha="center", va="center", fontsize=7.5, color="#666",
                  bbox=dict(fc="white", ec="none", pad=0.5))
    h, l_ = AX[0, 0].get_legend_handles_labels()
    fig.legend(h, l_, loc="upper center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(os.path.join(HERE, "figs", "soap_params.png"), dpi=300)
    print("figs/soap_params.png")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], os.path.abspath(sys.argv[3]))
    elif sys.argv[1] == "figure":
        figure()
    else:
        summary(os.path.abspath(sys.argv[2]))
