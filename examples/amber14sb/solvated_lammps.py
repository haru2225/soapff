"""soapff on the solvated ubiquitin: (1) the PDB file -> types and charges against OpenMM's; (2) the energy with PPPM against
OpenMM's PME; (3) the in.lammps soapff writes: minimisation and 1 ps of NVT with SHAKE water.
Ewald is used here, not PPPM: the Homebrew LAMMPS (22 Jul 2025, FFTW3 with threads) crashes when PPPM is set up a second time
(FFT3d destructor -> fftw_plan_awake), i.e. after the first run or minimize of an input
    python solvated_lammps.py . WORKDIR results/solvated.json      (after solvated.py, in this folder)"""
import json, os, subprocess, sys, time
import ase.io, numpy as np
import soapff
P, W, OUT = sys.argv[1:4]
res = {}
truth = ase.io.read(os.path.join(P, "ubq_water.extxyz"))
t0 = time.time()
s = soapff.build("amber14sb", os.path.join(P, "ubq_water.pdb"), os.path.join(W, "ubq_water"), md=True, relax_ps=0.0,
                 nvt_ps=1.0, dump_ps=0.5, kspace="ewald 1.0e-5", log=lambda *a: None)
tt, qq = np.array(s.types), s.charges
res["typing"] = dict(atoms=len(tt), seconds=round(time.time() - t0, 1), wrong_type=int((tt != truth.arrays["fftype"]).sum()),
                     wrong_charge=int((np.abs(qq - truth.arrays["ffcharge"]) > 1e-5).sum()), molecules=s.report["molecules"],
                     molecule_kinds=s.report["molecule_kinds"], net_charge=s.report["charge_sum"], sources=s.report["sources"])
print("typing", res["typing"], flush=True)
# energy: the extxyz coordinates (6 decimals, as OpenMM's) in the cell of the PDB, Ewald 1e-6, flexible water
at = truth.copy(); at.cell = ase.io.read(os.path.join(P, "ubq_water.pdb")).cell; at.pbc = True
e = soapff.build("amber14sb", at, os.path.join(W, "ubq_water_energy"), md=False, kspace="ewald 1.0e-6", log=lambda *a: None)
S = [l for l in open(e.files["settings"]).read().splitlines() if l.strip()]
inp = S + ["thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong", "thermo_modify format float %.6f", "run 0"]
open(os.path.join(e.out_dir, "in.energy"), "w").write("\n".join(inp) + "\n")
subprocess.run(["lmp_serial", "-in", "in.energy", "-log", "log.energy", "-screen", "none"], cwd=e.out_dir, check=True)
l = {k: float(v[-1]) for k, v in soapff.read_thermo(os.path.join(e.out_dir, "log.energy")).items()}
o = json.load(open(os.path.join(P, "results", "energies_solvated.json")))
lam = dict(bond=l["E_bond"], angle=l["E_angle"], torsion=l["E_dihed"] + l["E_impro"],
           nonbonded=l["E_vdwl"] + l["E_coul"] + l["E_long"], total=l["PotEng"])
omm = dict(bond=o["HarmonicBondForce"], angle=o["HarmonicAngleForce"], torsion=o["PeriodicTorsionForce"], nonbonded=o["NonbondedForce"],
           total=o["total"])
res["energy_ewald_vs_pme"] = dict(lammps=lam, openmm=omm, diff={k: lam[k] - omm[k] for k in omm})
print("energy", {k: f"{lam[k]:.3f} ({lam[k] - omm[k]:+.3f})" for k in omm}, flush=True)
# the in.lammps of soapff: minimisation, then NVT 300 K with SHAKE water
t0 = time.time()
log = s.run(np=1, lmp="lmp_serial", echo=False)
th = s.thermo()
res["run"] = dict(seconds=round(time.time() - t0, 1), steps=int(th["Step"][-1]), temp_last=float(th["Temp"][-1]),
                  pe_last=float(th["PotEng"][-1]), final_data=os.path.isfile(os.path.join(s.out_dir, "final.data")),
                  shake=[x for x in open(s.files["input"]).read().splitlines() if "shake" in x])
print("run", res["run"], flush=True)
json.dump(res, open(OUT, "w"), indent=1)
