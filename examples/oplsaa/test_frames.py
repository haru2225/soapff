"""Types and charges of the frames of the 413 K NPT melts, against Foyer's.  The test cells of test_typing.py were packed by
packmol and minimised; these frames are melts at 413 K and 1 bar, after md_prepare.py (OpenMM, 350 ps) and the NPT MD of
md_lammps.py (LAMMPS, one frame every 10 ps in WORK/lammps_<p>/traj.lammpstrj).

    python test_frames.py WORK results/frames.json P1 [P2 ...]     (soapff)  ->  also WORK/frames_<p>.npz (every atom)

Each frame (elements, coordinates and the box only) is given to soapff.build(oplsaa); its types and charges are compared with
those Foyer gave the start, and its bonds (made from the types) with Foyer's bonds of the start (WORK/start_<p>.json) (WORK/start_<p>.extxyz: the atoms keep their order in the data file and in the dump).  Each frame is
typed with strict=False, and whether soapff would have stopped on it (an atom below the threshold, a molecule whose charge is
not a whole number, an atom with no bond) is recorded: a frame with wrong atoms that does not stop is a force field written
with errors."""
import json
import os
import sys
import tempfile
import time

import ase.io
import numpy as np
from ase import Atoms

import soapff

W, OUT = os.path.abspath(sys.argv[1]), sys.argv[2]
res = {}
for p in sys.argv[3:]:
    start = ase.io.read(os.path.join(W, f"start_{p}.extxyz"))
    T, Q = np.array(start.arrays["fftype"]), np.array(start.arrays["ffcharge"], float)
    js = json.load(open(os.path.join(W, f"start_{p}.json")))      # Foyer's bonds of the start (the same in every frame)
    B, off = set(), 0
    for m in js["molecules"]:
        b = m["bonds"]
        B |= {(min(off + b[i], off + b[i + 1]), max(off + b[i], off + b[i + 1])) for i in range(0, len(b), 2)}
        off += len(m["elements"].split())
    frames = ase.io.read(os.path.join(W, f"lammps_{p}", "traj.lammpstrj"), index=":", format="lammps-dump-text")
    rows, per = [], []
    t0 = time.time()
    for k, fr in enumerate(frames):
        if fr.arrays.get("id") is not None:
            assert (np.asarray(fr.arrays["id"]) == np.arange(1, len(fr) + 1)).all()
        assert fr.get_chemical_symbols() == start.get_chemical_symbols(), f"{p} frame {k}: the atoms are not in the order of the start"
        at = Atoms(fr.get_chemical_symbols(), positions=fr.positions, cell=fr.cell, pbc=True)
        with tempfile.TemporaryDirectory(dir=W) as d:
            try:
                s = soapff.build("oplsaa", at, d, strict=False, log=lambda *a: None)
            except soapff.TypingError as e:                  # stops even with strict=False (an atom with no bond)
                per.append(dict(frame=k, wrong=None, stops=True, error=str(e)[:300]))
                continue
            tab = s.table()
            types, charges = np.array(s.types), np.array(s.charges, float)
            stops = any("not a whole number" in w or "not in the library" in w for w in s.report["warnings"])
            bonds = {(i, int(j) - 1) for i, r in enumerate(tab) for j in r["bonded_to"].split() if i < int(j) - 1}
        sim = np.array([float(r["similarity_soap"]) for r in tab])
        mar = np.array([float(r["margin_soap"]) for r in tab])
        ok = (types == T) & (np.abs(charges - Q) < 1e-6)
        rows.append((np.arange(1, len(at) + 1), np.full(len(at), k), sim, mar, ok))
        per.append(dict(frame=k, wrong=int((~ok).sum()), bonds_differing=len(bonds ^ B), stops=stops,
                        min_similarity=round(float(sim.min()), 4),
                        min_margin=round(float(mar.min()), 4)))
    a = [np.concatenate(x) for x in zip(*rows)]
    np.savez_compressed(os.path.join(W, f"frames_{p}.npz"), atom=a[0], frame=a[1], similarity=a[2], margin=a[3], right=a[4])
    res[p] = dict(frames=len(frames), atoms=len(start), atoms_typed=int(len(a[0])), wrong=int((~a[4]).sum()),
                  frames_that_stop=sum(r["stops"] for r in per),
                  frames_with_bonds_differing=sum(bool(r.get("bonds_differing")) for r in per),
                  frames_written_with_errors=sum(bool(r["wrong"]) and not r["stops"] for r in per),
                  frames_not_written=[r["frame"] for r in per if r["wrong"] is None],
                  min_similarity=round(float(a[2].min()), 4), min_margin=round(float(a[3].min()), 4),
                  margin_percentiles={str(q): round(float(np.percentile(a[3], q)), 4) for q in (0, 0.1, 1, 50)},
                  seconds_per_frame=round((time.time() - t0) / len(frames), 2), per_frame=per)
    print(p, json.dumps({k: v for k, v in res[p].items() if k != "per_frame"}), flush=True)
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=1)
