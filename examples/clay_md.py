"""The MD protocol of the clay examples, as LAMMPS commands for the LAMMPS Python module.

soapff writes the force field only (system.data, system.settings).  The runs of the clay examples are written here, in
LAMMPS: each function returns a list of commands, given to LAMMPS with lmp.commands_list(...), or written to an input
file for a run elsewhere (input_file).  They are the commands of the earlier in.lammps of these examples, so the results
kept in the folders are those of these runs.

    import soapff
    from lammps import lammps
    import clay_md

    s = soapff.build("clayff", "kaolinite_water.extxyz", "out")
    lmp = lammps(cmdargs=["-log", "out/log.lammps", "-screen", "none"])
    lmp.commands_list(s.commands())                     # the force field
    lmp.commands_list(clay_md.minimize())               # two-stage minimisation
    lmp.commands_list(clay_md.md(s, temp=300, nvt_ps=50, npt_ps=100, baro="aniso"))   # NVT, then NPT
"""
import os

TIMESTEP_FS = 0.5                     # the flexible SPC water of ClayFF
SEED = 4928459
ION_NAMES = ("Na", "K", "Ca", "Cs", "Ba", "Cl")
# NPT coupling: the three lengths independently (aniso), the full cell (tri), the length normal to the layers only (z),
# or the three together (iso).  z is the usual choice for slabs and interfaces, aniso or tri for bulk crystals.
BARO = {"iso": "iso", "aniso": "aniso", "tri": "tri", "z": "z"}


def thermo(every=1000):
    return ["thermo_style custom step temp press pe ke etotal vol density lx ly lz", f"thermo {every}"]


def minimize(thermo_every=1000):
    """the energy minimisation of a built structure (it may hold overlapping water): a small step limit first"""
    return thermo(thermo_every) + ["min_style cg", "min_modify dmax 0.01", "minimize 1.0e-4 1.0e-6 1000 10000",
                                   "min_modify dmax 0.1", "minimize 1.0e-6 1.0e-6 5000 50000"]


def msd_commands(s, every=200, out=None):
    """mean squared displacement of the water O and of each ion type, computed by LAMMPS (compute msd), written every
    `every` steps (0.1 ps at 0.5 fs) to msd_<group>.dat"""
    out = out or s.out_dir
    tid = s.type_ids
    groups = [("Ow", tid["ow"])] if "ow" in tid else []
    groups += [(n, tid[n]) for n in tid if n in ION_NAMES]
    L = []
    for g, t in groups:
        L += [f"group a_{g} type {t}", f"compute msd_{g} a_{g} msd",
              f"fix msd_{g} all ave/time 1 1 {every} c_msd_{g}[*] file {os.path.join(out, f'msd_{g}.dat')}"]
    return L


def md(s, temp=300.0, nvt_ps=50.0, npt_ps=0.0, press=1.0, baro="aniso", relax_ps=10.0, dump_ps=1.0,
       timestep_fs=TIMESTEP_FS, msd=True, out=None, seed=SEED):
    """after the minimisation: velocities at `temp` (K); NVT for relax_ps, not recorded (the step is reset to 0); then the
    production: NVT for nvt_ps with a trajectory (traj.lammpstrj: the ClayFF type name in the element column, wrapped
    positions with image flags, charges) and the mean squared displacements, and (npt_ps > 0) NPT at `press` (atm) with the
    coupling `baro`; the last structure is written to final.data.  The thermostat and barostat: Nose-Hoover, 100 fs and
    1000 fs.  Output files go to `out` (default: the folder of `s`)."""
    out = out or s.out_dir
    if baro not in BARO:
        raise ValueError(f"baro must be one of {list(BARO)}")
    steps = lambda ps: int(round(ps * 1000.0 / timestep_fs))
    T, P = f"{temp:g}", f"{press:g}"
    L = ["reset_timestep 0", f"velocity all create {T} {seed} mom yes rot no dist gaussian", f"timestep {timestep_fs:g}"]
    if npt_ps > 0 and baro == "tri":
        L.append("change_box all triclinic")           # before the dump (a LAMMPS rule)
    if relax_ps > 0:
        L += [f"fix relax all nvt temp {T} {T} 100.0", f"run {steps(relax_ps)}", "unfix relax", "reset_timestep 0"]
    L += [f"fix md all nvt temp {T} {T} 100.0",
          f"dump d all custom {max(1, steps(dump_ps))} {os.path.join(out, 'traj.lammpstrj')} id type element x y z ix iy iz q",
          "dump_modify d sort id element " + " ".join(s.type_names)]
    if msd:
        L += msd_commands(s, max(1, steps(0.1)), out)
    if nvt_ps > 0:
        L.append(f"run {steps(nvt_ps)}")
    if npt_ps > 0:
        L += ["unfix md", f"fix md all npt temp {T} {T} 100.0 {BARO[baro]} {P} {P} 1000.0", f"run {steps(npt_ps)}"]
    L.append(f"write_data {os.path.join(out, 'final.data')}")
    return L


def start(s, log="log.lammps", screen=False, cmdargs=()):
    """a LAMMPS instance (Python module) with the force field of the soapff System `s`; log in its folder"""
    from lammps import lammps
    lmp = lammps(cmdargs=["-log", os.path.join(s.out_dir, log)] + ([] if screen else ["-screen", "none"]) + list(cmdargs))
    lmp.commands_list(s.commands())
    return lmp


def input_file(s, commands, name="in.lammps"):
    """the same commands as a LAMMPS input in the folder of `s` (for a run on another computer: lmp -in in.lammps there);
    paths are made relative to the folder"""
    L = ["# the force field (soapff)", "include system.settings"]
    L += [c.replace(s.out_dir + os.sep, "") for c in commands]
    path = os.path.join(s.out_dir, name)
    open(path, "w").write("\n".join(L) + "\n")
    return path
