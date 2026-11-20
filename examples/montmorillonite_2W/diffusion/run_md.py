"""The MD of the 2-layer Na-montmorillonite hydrate (../initial.xyz, 1,554 atoms; the force field built by soapff with the
bundled clayff library: ClayFF + flexible SPC water), written in LAMMPS and run with its Python module
(pip install "lammps[mpi]"), or with an lmp executable (--lmp).  ../initial.xyz -> start/: minimisation, 10 ps NVT
relaxation, 10 ps NVT, 100 ps NPT along z (1 atm), 300 K, 0.5 fs (../../clay_md.py), trajectory every 10 ps, compute msd of
the water O and Na+ every 0.1 ps.  Its last frame (start/final.data) is one of the MD frames of ../../clayff/water_md_frames.py.
start/ gets in.lammps (include system.settings) with the same commands: lmp -in in.lammps runs it again.

    python run_md.py [--lmp "mpirun -np 3 lmp"]
        --lmp CMD: run in.lammps with this command (CMD -in in.lammps -log log.lammps -screen none) instead of the LAMMPS
                   Python module in this process (serial)"""
import argparse
import os
import subprocess
import sys

import soapff

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))           # clay_md.py
import clay_md  # noqa: E402

START = os.path.join(HERE, "start")
TIMESTEP_FS = 0.5
# start/: the MD of clay_md.md
START_MD = dict(temp=300.0, relax_ps=10.0, nvt_ps=10.0, npt_ps=100.0, press=1.0, baro="z", dump_ps=10.0)


def execute(s, force_field, commands, lmp_cmd=None):
    """the commands in the folder of `s`: written to in.lammps, then run with the LAMMPS Python module (force field first)
    or with the executable lmp_cmd"""
    clay_md.input_file(s, commands)
    if lmp_cmd:
        env = dict(os.environ, OMP_NUM_THREADS="1", HWLOC_COMPONENTS="-gl")
        subprocess.run(f"{lmp_cmd} -in in.lammps -log log.lammps -screen none", shell=True, cwd=s.out_dir, env=env, check=True)
        return
    from lammps import lammps
    lmp = lammps(cmdargs=["-log", os.path.join(s.out_dir, "log.lammps"), "-screen", "none"])
    try:
        lmp.commands_list(force_field + commands)
    finally:
        lmp.close()


def start(lmp_cmd=None):
    s = soapff.build("clayff", os.path.join(HERE, "..", "initial.xyz"), START)
    execute(s, s.commands(), clay_md.minimize() + clay_md.md(s, timestep_fs=TIMESTEP_FS, **START_MD), lmp_cmd)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="the MD of the 2-layer hydrate (start/)")
    ap.add_argument("--lmp", default=None, help='run in.lammps with this command, e.g. "mpirun -np 3 lmp"')
    start(ap.parse_args().lmp)
