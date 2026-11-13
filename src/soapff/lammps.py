"""LAMMPS files of the molecular engine (system.data, system.settings, in.lammps) and running LAMMPS (an executable, or the
LAMMPS Python module) in an output folder of soapff (either engine)."""
import os
import shutil
import subprocess
import sys

import numpy as np

from soapff_cui._clayff.ff import lammps_frame
from soapff_cui.core import kspace_setup, parse_boundary

TIMESTEP_FS = 1.0


def _num(x):
    return f"{x:.10g}"


def write(out_dir, at, nb, types, q, topo, lib, cutoff=None, kspace=None, boundary="p p p", pair_style=None,
          dihedral_style=None):
    """system.data and system.settings of a typed structure; returns the type ids and what in.lammps needs"""
    from .molecular import unwrapped
    lm = lib.lammps
    bnd, periodic = parse_boundary(boundary)
    el = at.get_chemical_symbols()
    tnames = sorted(set(types), key=lambda t: (lib.types[t]["element"], t))
    tid = {t: k + 1 for k, t in enumerate(tnames)}
    bkeys = sorted({k for _, k in topo["bonds"]})
    bid = {k: i + 1 for i, k in enumerate(bkeys)}
    akeys = sorted({k for _, k in topo["angles"]})
    aid = {k: i + 1 for i, k in enumerate(akeys)}
    bpar = dict(lib.bonds)
    apar = dict(lib.angles)
    dpar = dict(lib.dihedrals)
    ipar = dict(lib.impropers)
    dstyle = dihedral_style or lm.get("dihedral_style", "fourier")
    # dihedral types: one per row (fourier: all its terms in one type) or one per term (charmm, weight 0: the 1-4 pairs come
    # from special_bonds)
    dihedrals = []
    dtypes = {}
    for idx, key in topo["dihedrals"]:
        terms = [t for t in dpar[key] if t[0] != 0.0]     # terms with K = 0 add nothing (OpenMM leaves them out too)
        if not terms:
            continue
        if dstyle.split()[0] == "charmm":
            for K, n, d in terms:
                if abs(d - round(d)) > 1e-6:
                    raise ValueError(f"dihedral {'-'.join(key)}: phase {d} is not a whole degree; use dihedral_style fourier")
                tk = (key, K, n, int(round(d)))
                dtypes.setdefault(tk, len(dtypes) + 1)
                dihedrals.append((idx, dtypes[tk]))
        else:
            dtypes.setdefault(key, len(dtypes) + 1)
            dihedrals.append((idx, dtypes[key]))
    ikeys = sorted({k for _, k in topo["impropers"]})
    iid = {k: i + 1 for i, k in enumerate(ikeys)}
    # molecules and image flags (each molecule whole)
    from .molecular import molecules
    mol = molecules(nb) + 1
    cell = np.array(at.cell)
    pos_u = unwrapped(at, nb)
    frac = pos_u @ np.linalg.inv(cell)
    img = np.where(periodic, np.floor(frac), 0).astype(int)
    wrapped = at.copy()
    wrapped.positions = (frac - img) @ cell
    p, xl = lammps_frame(wrapped, wrap=False)
    lx, ly, lz, xy, xz, yz = p.get_lammps_prism()
    tilted = any(abs(v) > 1e-8 for v in (xy, xz, yz))
    L = [f"soapff: library {lib.name}, types and charges from the references (SOAP)", "",
         f"{len(at)} atoms", f"{len(tnames)} atom types", f"{len(topo['bonds'])} bonds", f"{len(bkeys)} bond types",
         f"{len(topo['angles'])} angles", f"{len(akeys)} angle types", f"{len(dihedrals)} dihedrals",
         f"{len(dtypes)} dihedral types", f"{len(topo['impropers'])} impropers", f"{len(ikeys)} improper types", "",
         f"0.0 {lx:.8f} xlo xhi", f"0.0 {ly:.8f} ylo yhi", f"0.0 {lz:.8f} zlo zhi"]
    if tilted:
        L.append(f"{xy:.8f} {xz:.8f} {yz:.8f} xy xz yz")
    L += ["", "Masses", ""] + [f"{tid[t]} {lib.types[t]['mass']:.5f}  # {t}" for t in tnames]
    L += ["", "Atoms  # full", ""] + [f"{i + 1} {mol[i]} {tid[types[i]]} {q[i]:.6f} {xl[i, 0]:.8f} {xl[i, 1]:.8f} {xl[i, 2]:.8f}"
                                      f" {img[i, 0]} {img[i, 1]} {img[i, 2]}" for i in range(len(at))]
    if topo["bonds"]:
        L += ["", "Bonds", ""] + [f"{k + 1} {bid[key]} {a + 1} {b + 1}" for k, ((a, b), key) in enumerate(topo["bonds"])]
    if topo["angles"]:
        L += ["", "Angles", ""] + [f"{k + 1} {aid[key]} {a + 1} {b + 1} {c + 1}"
                                   for k, ((a, b, c), key) in enumerate(topo["angles"])]
    if dihedrals:
        L += ["", "Dihedrals", ""] + [f"{k + 1} {t} {a + 1} {b + 1} {c + 1} {d + 1}" for k, ((a, b, c, d), t) in enumerate(dihedrals)]
    if topo["impropers"]:
        L += ["", "Impropers", ""] + [f"{k + 1} {iid[key]} {a + 1} {b + 1} {c + 1} {d + 1}"
                                      for k, ((a, b, c, d), key) in enumerate(topo["impropers"])]
    open(os.path.join(out_dir, "system.data"), "w").write("\n".join(L) + "\n")

    # settings
    ps = (pair_style or lm.get("pair_style", "lj/cut/coul/long 10.0")).split()
    if cutoff is not None:
        ps = ps[:1] + [str(cutoff)]
    long_range = ps[0].endswith("coul/long")
    ks_lines, note = [], None
    if long_range:
        coul, ks_lines, note = kspace_setup(kspace or lm.get("kspace_style", "pppm 1.0e-5"), bnd, periodic, tilted)
        if coul == "msm":
            ps[0] = ps[0].replace("coul/long", "coul/msm")
    S = [f"# soapff settings: library {lib.name} ({lib.library.get('reference', lib.description)}); types and charges in "
         f"system.data from the references of the library", "units real", "atom_style full", f"boundary {' '.join(bnd)}",
         f"pair_style {' '.join(ps)}"]
    if lm.get("pair_modify"):
        S.append(f"pair_modify {lm['pair_modify']}")
    S += [f"bond_style {lm.get('bond_style', 'harmonic')}", f"angle_style {lm.get('angle_style', 'harmonic')}"]
    if dtypes:
        S.append(f"dihedral_style {dstyle}")
    if ikeys:
        S.append(f"improper_style {lm.get('improper_style', 'cvff')}")
    S += [f"special_bonds {lm.get('special_bonds', 'lj/coul 0.0 0.0 0.0')}", "read_data system.data"]
    for t in tnames:
        x = lib.types[t]
        S.append(f"pair_coeff {tid[t]} {tid[t]} {_num(x['epsilon'])} {_num(x['sigma'])}  # {t}")
    for k in bkeys:
        K, r0 = bpar[k]
        S.append(f"bond_coeff {bid[k]} {_num(K)} {_num(r0)}  # {'-'.join(k)}")
    for k in akeys:
        K, th = apar[k]
        S.append(f"angle_coeff {aid[k]} {_num(K)} {_num(th)}  # {'-'.join(k)}")
    for k, i in sorted(dtypes.items(), key=lambda kv: kv[1]):
        if dstyle.split()[0] == "charmm":
            key, K, n, d = k
            S.append(f"dihedral_coeff {i} {_num(K)} {n} {d} 0.0  # {'-'.join(key)}")
        else:
            terms = [t for t in dpar[k] if t[0] != 0.0]
            S.append(f"dihedral_coeff {i} {len(terms)} " + " ".join(f"{_num(K)} {n} {_num(d)}" for K, n, d in terms)
                     + f"  # {'-'.join(k)}")
    for k in ikeys:
        K, d, n = ipar[k]
        S.append(f"improper_coeff {iid[k]} {_num(K)} {int(round(d))} {int(round(n))}  # {'-'.join(k)}")
    S += ks_lines
    open(os.path.join(out_dir, "system.settings"), "w").write("\n".join(S) + "\n")
    # SHAKE: the bonds and angles whose atoms are all of the rigid types
    rigid = set(lm.get("rigid", "").split())
    shake_b = [bid[k] for k in bkeys if rigid and set(k) <= rigid]
    shake_a = [aid[k] for k in akeys if rigid and set(k) <= rigid]
    return dict(type_ids=tid, elements=[lib.types[t]["element"] for t in tnames], shake_bonds=shake_b, shake_angles=shake_a,
                nbonds=len(topo["bonds"]), nangles=len(topo["angles"]), ndihedrals=len(dihedrals),
                nimpropers=len(topo["impropers"]), kspace=ks_lines, kspace_note=note, tilted=tilted, nmol=int(mol.max()),
                box=[lx, ly, lz, xy, xz, yz])


def input_file(lib, info, md=True, temp=300.0, relax_ps=10.0, nvt_ps=50.0, npt_ps=0.0, press=1.0, baro="aniso", dump_ps=1.0,
               timestep_fs=None, seed=4928459):
    """in.lammps: minimisation, then (md) a relaxation NVT (not recorded), NVT and optionally NPT, with a trajectory"""
    ts = float(timestep_fs or lib.lammps.get("timestep", TIMESTEP_FS))
    steps = lambda ps: int(round(ps * 1000.0 / ts))
    I = [f"# soapff: library {lib.name}.  lmp -in in.lammps", "include system.settings", "neighbor 2.0 bin",
         "neigh_modify every 1 delay 0 check yes",
         "thermo_style custom step temp press pe ebond eangle edihed eimp evdwl ecoul elong vol density", "thermo 100",
         "# minimisation (built structures may have close contacts)", "min_style cg", "minimize 1.0e-4 1.0e-6 1000 10000",
         "reset_timestep 0"]
    if md:
        # SHAKE before the velocities: velocity create then counts the degrees of freedom SHAKE removes (3 per rigid
        # water), so the start is at `temp` (with the velocities first it reads about 1.4 x temp once SHAKE is on)
        if info["shake_bonds"] or info["shake_angles"]:
            I.append("fix rigid all shake 1.0e-4 20 0" + (" b " + " ".join(map(str, info["shake_bonds"])) if info["shake_bonds"] else "")
                     + (" a " + " ".join(map(str, info["shake_angles"])) if info["shake_angles"] else "")
                     + "   # rigid water (library: rigid = " + lib.lammps.get("rigid", "") + ")")
        I += [f"velocity all create {temp:g} {seed} mom yes rot no dist gaussian", f"timestep {ts:g}"]
        tdamp = 100.0 * ts
        if relax_ps > 0:
            I += ["# relaxation, not recorded", f"fix relax all nvt temp {temp:g} {temp:g} {tdamp:g}", f"run {steps(relax_ps)}",
                  "unfix relax", "reset_timestep 0", "# production"]
        every = max(1, steps(dump_ps))
        I += [f"fix md all nvt temp {temp:g} {temp:g} {tdamp:g}",
              f"dump d all custom {every} traj.lammpstrj id mol type element x y z ix iy iz q",
              "dump_modify d sort id element " + " ".join(info["elements"])]
        if nvt_ps > 0:
            I.append(f"run {steps(nvt_ps)}")
        if npt_ps > 0:
            if baro == "tri":
                I.insert(I.index(f"fix md all nvt temp {temp:g} {temp:g} {tdamp:g}"), "change_box all triclinic")
            I += ["unfix md", f"fix md all npt temp {temp:g} {temp:g} {tdamp:g} {baro} {press:g} {press:g} {1000.0 * ts:g}",
                  f"run {steps(npt_ps)}"]
        I.append("write_data final.data")
    return "\n".join(I) + "\n"


# ------------------------------------------------------------------------------------------------- running LAMMPS
def find_lammps():
    """(lmp, mpirun): $SOAPFF_LMP / $SOAPFF_MPIRUN, else lmp / lmp_mpi / lmp_serial next to this Python, then on the PATH; the
    mpirun (or mpiexec) next to lmp is preferred (same MPI)"""
    lmp = os.environ.get("SOAPFF_LMP")
    mpirun = os.environ.get("SOAPFF_MPIRUN")
    if lmp:
        return shutil.which(lmp) or lmp, mpirun
    dirs = [os.path.dirname(sys.executable)] + os.environ.get("PATH", "").split(os.pathsep)

    def first(names):
        for d in dirs:
            for n in names:
                p = os.path.join(d, n)
                if d and os.path.isfile(p) and os.access(p, os.X_OK):
                    return p
        return None
    lmp = first(["lmp", "lmp_mpi", "lmp_serial"])
    if not mpirun and lmp and not lmp.endswith("lmp_serial"):
        mpirun = next((os.path.join(os.path.dirname(lmp), n) for n in ("mpirun", "mpiexec")
                       if os.path.isfile(os.path.join(os.path.dirname(lmp), n))), None) or first(["mpirun", "mpiexec"])
    return lmp, mpirun


def run(out_dir, np_=1, lmp=None, mpirun=None, input="in.lammps", log="log.lammps", extra=(), echo=True):
    """run LAMMPS on `input` in out_dir (np_ > 1: with mpirun).  Returns the path of the log; raises when LAMMPS stops with an
    error (the last lines of the log are in the message)."""
    lmp0, mpi0 = find_lammps()
    lmp = lmp or lmp0
    if not lmp:
        raise FileNotFoundError("LAMMPS not found: give lmp='/path/to/lmp' or set $SOAPFF_LMP")
    cmd = [lmp, "-in", input, "-log", log] + list(extra)
    if np_ and np_ > 1:
        mpirun = mpirun or mpi0 or shutil.which("mpirun")
        if not mpirun:
            raise FileNotFoundError("mpirun not found: give mpirun= or run with np=1")
        cmd = [mpirun, "-np", str(int(np_))] + cmd
    with subprocess.Popen(cmd, cwd=out_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1) as p:
        tail = []
        for line in p.stdout:
            tail = (tail + [line])[-30:]
            if echo:
                sys.stdout.write(line)
        p.wait()
    if p.returncode != 0:
        raise RuntimeError(f"LAMMPS stopped with code {p.returncode} ({' '.join(cmd)} in {out_dir}):\n" + "".join(tail[-15:]))
    return os.path.join(out_dir, log)


def session(out_dir, input="in.lammps", cmdargs=("-screen", "none"), run_input=True):
    """the LAMMPS Python module (from lammps import lammps) started in out_dir, with `input` run when run_input; returns the
    lammps object (lmp.command(...), lmp.get_thermo('pe'), lmp.extract_atom(...))"""
    try:
        from lammps import lammps
    except ImportError as e:
        raise ImportError("the LAMMPS Python module is not installed (pip install lammps, or conda install -c conda-forge "
                          "lammps); soapff.run() uses the lmp executable instead") from e
    here = os.getcwd()
    out_dir = os.path.abspath(out_dir)
    os.chdir(out_dir)                      # the input reads system.settings / system.data from its folder
    try:
        lmp = lammps(cmdargs=list(cmdargs) + ["-log", os.path.join(out_dir, "log.lammps")])
        if run_input:
            lmp.file(input)
    finally:
        os.chdir(here)
    return lmp
