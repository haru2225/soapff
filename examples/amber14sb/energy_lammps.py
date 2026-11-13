"""LAMMPS energies (kcal/mol) of the test structures with the files soapff writes (library amber14sb), no cutoff (lj/cut/coul/cut
999 A, boundary f f f), compared with OpenMM (energies_openmm.json).
    python energy_lammps.py TESTDIR energies_openmm.json WORKDIR OUT.json"""
import json, os, subprocess, sys
import soapff
T, ref, W, OUT = sys.argv[1:5]
omm = json.load(open(ref))
res = {}
for name, e in omm.items():
    s = soapff.build("amber14sb", os.path.join(T, name + ".extxyz"), os.path.join(W, name), md=False, log=lambda *a: None)
    S = open(s.files["settings"]).read().splitlines()
    styles = [l for l in S if l.split()[:1] and l.split()[0] in ("bond_style", "angle_style", "dihedral_style", "improper_style",
                                                               "special_bonds")]
    coeffs = [l for l in S if l.split()[:1] and l.split()[0].endswith("_coeff")]
    inp = ["units real", "atom_style full", "boundary f f f", "pair_style lj/cut/coul/cut 999.0", "pair_modify mix arithmetic"] \
        + styles + ["read_data system.data"] + coeffs + \
        ["thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong", "thermo_modify format float %.6f", "run 0"]
    open(os.path.join(s.out_dir, "in.energy"), "w").write("\n".join(inp) + "\n")
    out = subprocess.run(["lmp_serial", "-in", "in.energy", "-log", "log.energy", "-screen", "none"], cwd=s.out_dir,
                         capture_output=True, text=True)
    th = soapff.read_thermo(os.path.join(s.out_dir, "log.energy"))
    l = {k: float(v[-1]) for k, v in th.items()}
    lam = dict(bond=l["E_bond"], angle=l["E_angle"], torsion=l["E_dihed"] + l["E_impro"], nonbonded=l["E_vdwl"] + l["E_coul"],
               total=l["PotEng"])
    o = dict(bond=e["HarmonicBondForce"], angle=e["HarmonicAngleForce"], torsion=e["PeriodicTorsionForce"],
             nonbonded=e["NonbondedForce"], total=e["total"])
    res[name] = dict(atoms=s.report["atoms"], lammps=lam, openmm=o, dihedral=l["E_dihed"], improper=l["E_impro"],
                     diff={k: lam[k] - o[k] for k in o})
    print(f"{name:22s} " + "  ".join(f"{k} {lam[k]:10.4f} ({lam[k] - o[k]:+.4f})" for k in o), flush=True)
json.dump(res, open(OUT, "w"), indent=1)
