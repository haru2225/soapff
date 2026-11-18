"""The swelling curve of a clay in LAMMPS (Cygan et al. 2004): water molecules are added to the interlayer state by state,
each state is relaxed and run in NPT, and the cell is written to swelling.dat (d001 = lz; read by analyze_swelling.py,
function swelling).  This code was soapff_cui.core.swelling_schedule / swelling_lines; soapff now writes the force field
only, so it lives here, with the example that uses it (make_run.py).

  swelling_schedule(...)   the water counts of the states
  swelling_lines(...)      the LAMMPS commands as a block of an input file (loops with label / jump SELF), for in.lammps
                           and the lmp executable (as the run stored in run/)
  run(lmp, ...)            the same commands, given to the LAMMPS Python module; the loops are Python loops

The structure needs the water types (soapff_cui.core.reserve_water_types) and a molecule file water.mol of one water
(soapff_cui.core.water_molecule) in the folder `out`."""
import os

ION_NAMES = ("Na", "K", "Ca", "Cs", "Ba", "Cl")


def swelling_schedule(one_to=23, two_to=99, down=True, step_up2=2, step_down=2):
    """water molecules added to the starting structure, state by state: one at a time up to `one_to`, then `step_up2`
    at a time up to `two_to` (Cygan et al. 2004: 1 ... 23, 25 ... 99); then, if down, removed `step_down` at a time
    back to 1 (deswelling).  -> (up, down) lists of water counts (relative to the starting structure)."""
    up = list(range(1, one_to + 1)) + list(range(one_to + step_up2, two_to + 1, step_up2))
    dn = list(range(up[-1] - step_down, 0, -step_down)) if down else []
    return up, dn


def _blocks(type_ids, n_nonwater, n_water0, mol_max, temp, press, timestep_fs, tilted, up, down, eq_ps=30.0, avg_ps=20.0,
            seed=4928459, overlap=2.0, out=""):
    """the commands of the swelling curve in four blocks: head (once), insert (one try to add a water to the interlayer),
    state (the relaxation and the NPT of one water content), remove (the last added waters taken out, to v_sw_rem).
    Files (water.mol, swelling.dat, swelling_states.lammpstrj) are in the folder `out` ('' = the working folder of LAMMPS)."""
    f = lambda name: os.path.join(out, name) if out else name
    nonwater = {"ow", "hw"} | set(ION_NAMES)
    clay = " ".join(str(v) for k, v in sorted(type_ids.items(), key=lambda kv: kv[1]) if k not in nonwater)
    s_eq = int(round(eq_ps * 1000.0 / timestep_fs))
    s_avg = int(round(avg_ps * 1000.0 / timestep_fs))
    every = max(1, int(round(100.0 / timestep_fs)))
    # the cell: x, y, z and (tilted box) the in-plane tilt xy; the interlayer shear tilts xz, yz stay fixed: with water in
    # the interlayer the layers slide (Cygan et al. 2004: alpha and gamma drift), the tilts run away and fix npt stops
    # ("tilted box too far in one step"); d001 = lz does not depend on them
    baro = (f"x {{P}} {{P}} 1000.0 y {{P}} {{P}} 1000.0 z {{P}} {{P}} 1000.0 xy 0.0 0.0 1000.0 couple none" if tilted
            else "aniso {P} {P} 1000.0")
    T, P = f"{float(temp):g}", f"{float(press):g}"
    n_pre = int(round(2000.0 / timestep_fs))        # 2 ps
    state = [# the clay back to the middle of the cell in z (once per state; a lengthened box leaves it off centre)
             "displace_atoms all move 0 0 $(0.5*(zlo+zhi)-xcm(sw_clay,z)) units box",
             f"minimize 1.0e-4 1.0e-6 2000 20000",
             # the inserted water starts with large forces: a short damped run (limited displacement per step, Langevin)
             # before the barostat, which otherwise moves the cell too far in one step ("tilted box too far")
             f"fix sw_lim all nve/limit 0.1",
             f"fix sw_lan all langevin {T} {T} 50.0 {seed + 3} zero yes",
             f"run {n_pre}",
             "unfix sw_lan", "unfix sw_lim",
             f"fix sw_npt all npt temp {T} {T} 100.0 {baro.format(P=P)}",
             "variable sw_phase equal 0",
             f'fix sw_out all print {every} "${{sw_step}} ${{sw_nw}} ${{sw_phase}} ${{sw_lx}} ${{sw_ly}} ${{sw_lz}} ${{sw_xy}} '
             f'${{sw_xz}} ${{sw_yz}} ${{sw_pe}} ${{sw_temp}}" append {f("swelling.dat")} screen no',
             f"run {s_eq}",
             "variable sw_phase equal 1",
             f"run {s_avg}",
             "unfix sw_out", "unfix sw_npt",
             f"write_dump all custom {f('swelling_states.lammpstrj')} id type x y z modify append yes"]
    head = [f"# swelling curve (soapff; Cygan et al. 2004): water added to the interlayer ({len(up)} states), "
            + (f"then removed ({len(down)} states); " if down else "")
            + f"each state NPT {eq_ps:g} ps + {avg_ps:g} ps (recorded); d001 = lz in swelling.dat",
            f"velocity all create {T} {seed} mom yes rot no dist gaussian",
            f"timestep {timestep_fs:g}",
            f"shell rm -f {f('swelling_states.lammpstrj')}",
            f"molecule h2o {f('water.mol')}",
            f"group sw_clay type {clay}",
            f"variable sw_nw equal (atoms-{n_nonwater})/3-{n_water0}",
            "variable sw_step equal step", "variable sw_lx equal lx", "variable sw_ly equal ly", "variable sw_lz equal lz",
            "variable sw_xy equal xy", "variable sw_xz equal xz", "variable sw_yz equal yz", "variable sw_pe equal pe",
            "variable sw_temp equal temp",
            f'print "# step added_waters phase lx ly lz xy xz yz pe temp" file {f("swelling.dat")} screen no',
            "variable sw_k equal 0",
            "region sw_lo block INF INF INF INF EDGE 0.1 units box", "region sw_hi block INF INF INF INF 0.2 EDGE units box",
            "region sw_il union 2 sw_lo sw_hi"]
    # the interlayer is found from the clay: water is added between its top and the box top or its bottom and the box bottom
    # (the interlayer across the periodic boundary); when no water fits (create_atoms overlap), the box is lengthened by
    # 0.5 A along z
    insert = ["region sw_il delete", "region sw_lo delete", "region sw_hi delete",
              "variable sw_z1 equal $(bound(sw_clay,zmin)-0.3)", "variable sw_z2 equal $(bound(sw_clay,zmax)+0.3)",
              "region sw_lo block INF INF INF INF EDGE $((v_sw_z1>zlo+0.05)*v_sw_z1+(v_sw_z1<=zlo+0.05)*(zlo+0.05)) units box",
              "region sw_hi block INF INF INF INF $((v_sw_z2<zhi-0.05)*v_sw_z2+(v_sw_z2>=zhi-0.05)*(zhi-0.05)) EDGE units box",
              "region sw_il union 2 sw_lo sw_hi",
              "variable sw_k equal $(v_sw_k+1)",
              "variable sw_before equal $(v_sw_nw)",
              f"create_atoms 0 random 1 $({seed}+7*v_sw_k:%.0f) sw_il mol h2o $({seed}+11*v_sw_k:%.0f) overlap {overlap:g} maxtry 500 units box",
              'if "$(v_sw_nw) == ${sw_before}" then "change_box all z delta 0 0.5 units box"']
    # the last added waters (molecule IDs above mol_max) are removed down to v_sw_rem
    remove = [f"group sw_del molecule <> $({mol_max}+v_sw_rem+1:%.0f) $({mol_max}+v_sw_nw:%.0f)",
              "delete_atoms group sw_del mol yes compress no",
              "group sw_del delete"]
    return head, insert, state, remove


GIVE_UP = 5000          # insertion tries in all before the run stops


def swelling_lines(type_ids, n_nonwater, n_water0, mol_max, temp, press, timestep_fs, tilted, up, down,
                   eq_ps=30.0, avg_ps=20.0, seed=4928459, overlap=2.0, out=""):
    """LAMMPS commands of the swelling curve (Cygan et al. 2004), as a block of an input file: for each water content of `up`
    (then `down`), water molecules are added to (removed from) the interlayer, the structure is relaxed (minimize), then NPT
    (x, y, z and the in-plane tilt xy for a tilted box, else aniso; the interlayer shear tilts xz, yz fixed) for eq_ps and
    avg_ps more; every 100 fs (step, waters, phase, cell, energy) is appended to swelling.dat (phase 1 = the averaging
    part; d001 = lz against the water content).
    The interlayer is found from the clay: the clay is moved to the middle of the cell in z at each state (displace_atoms)
    and water is added between its top and the box top or its bottom and the box bottom (the interlayer across the periodic
    boundary).  When no water fits (create_atoms overlap), the box is lengthened by 0.5 A along z and the insertion tried
    again.  n_nonwater: atoms that are not water at the start; n_water0: water molecules at the start; mol_max: the largest
    molecule ID of system.data (added molecules get the IDs above it, removed are the last added).
    The loops use label / jump SELF: the block runs from an input file (lmp -in in.lammps, or lmp.file(...) of the Python
    module), not from lmp.commands_list (use run() there)."""
    head, insert, state, remove = _blocks(type_ids, n_nonwater, n_water0, mol_max, temp, press, timestep_fs, tilted, up, down,
                                          eq_ps, avg_ps, seed, overlap, out)
    L = head + [f"variable sw_add index {' '.join(map(str, up))}", "label sw_up", "label sw_ins"] + insert + [
        f'if "$(v_sw_k) > {GIVE_UP}" then "print \'swelling: water could not be added\'" "quit 1"',
        'if "$(v_sw_nw) < ${sw_add}" then "jump SELF sw_ins"'] + state + ["next sw_add", "jump SELF sw_up"]
    if down:
        L += [f"variable sw_rem index {' '.join(map(str, down))}", "label sw_dn"] + remove + state + [
            "next sw_rem", "jump SELF sw_dn"]
    return L


def run(lmp, type_ids, n_nonwater, n_water0, mol_max, temp, press, timestep_fs, tilted, up, down,
        eq_ps=30.0, avg_ps=20.0, seed=4928459, overlap=2.0, out="", log=print):
    """the swelling curve with the LAMMPS Python module (lmp: a lammps instance with the force field, minimised): the
    commands of swelling_lines, in the same order and with the same seeds, the loops written in Python.  Files go to the
    folder `out` (give an absolute path)."""
    head, insert, state, remove = _blocks(type_ids, n_nonwater, n_water0, mol_max, temp, press, timestep_fs, tilted, up, down,
                                          eq_ps, avg_ps, seed, overlap, out)
    lmp.commands_list(head)
    for n in up:
        while True:
            lmp.commands_list(insert)
            if lmp.extract_variable("sw_k") > GIVE_UP:
                raise RuntimeError("swelling: water could not be added")
            if lmp.extract_variable("sw_nw") >= n:
                break
        log(f"swelling: {n} waters added, box lz {lmp.get_thermo('lz'):.2f} A")
        lmp.commands_list(state)
    for n in down:
        lmp.command(f"variable sw_rem equal {n}")
        lmp.commands_list(remove)
        log(f"deswelling: {n} waters left")
        lmp.commands_list(state)
