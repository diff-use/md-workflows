"""Re-solvate after equilibration and run a second round of min+equil.

Corresponds to run_all.sh line:
  bash scripts/resolvate.sh
"""

import re
import subprocess
from pathlib import Path


def run(ntmpi: int = 8, ntomp: int = 1):
    maxsol = _compute_maxsol()
    _resolvate(maxsol)
    _update_topology(maxsol)
    _minimize(ntmpi, ntomp)
    _equilibrate(ntmpi, ntomp)


def _compute_maxsol() -> int:
    """Do a trial solvation to figure out 25 % fill."""
    result = subprocess.run(
        ["gmx", "solvate",
         "-cp", "md_equil.gro",
         "-cs", "waterbox/water_equil.gro",
         "-o", "tmp.pdb"],
        capture_output=True, text=True,
    )
    log_text = result.stdout + result.stderr
    with open("gmx_solvate.log", "w") as fh:
        fh.write(log_text)

    for line in log_text.splitlines():
        m = re.search(r"Number of solvent molecules:\s+(\d+)", line)
        if m:
            return int(int(m.group(1)) * 0.25)
    raise RuntimeError("Could not parse solvent molecule count from gmx solvate output")


def _resolvate(maxsol: int):
    result = subprocess.run(
        ["gmx", "solvate",
         "-cp", "md_equil.gro",
         "-cs", "waterbox/water_equil.gro",
         "-o", "md_resolv.pdb",
         "-maxsol", str(maxsol)],
        capture_output=True, text=True, check=True,
    )
    with open("gmx_resolvate.log", "w") as fh:
        fh.write(result.stdout + result.stderr)


def _update_topology(maxsol: int):
    with open("md_model_posre.top", "a") as fh:
        fh.write(f"WAT {maxsol}\n")


def _minimize(ntmpi: int, ntomp: int):
    artifacts_dir = Path("artifacts")
    subprocess.run([
        "gmx", "grompp",
        "-f", str(artifacts_dir / "min.mdp"),
        "-c", "md_resolv.pdb",
        "-o", "md_resolv_min.tpr",
        "-p", "md_model_posre.top",
    ], capture_output=True, text=True, check=True)

    subprocess.run([
        "gmx", "mdrun",
        "-ntmpi", str(ntmpi),
        "-ntomp", str(ntomp),
        "-deffnm", "md_resolv_min",
        "-v",
    ], check=True)


def _equilibrate(ntmpi: int, ntomp: int):
    artifacts_dir = Path("artifacts")
    subprocess.run([
        "gmx", "grompp",
        "-f", str(artifacts_dir / "equil.mdp"),
        "-c", "md_resolv_min.gro",
        "-o", "md_resolv_equil.tpr",
        "-p", "md_model_posre.top",
        "-r", "md_model.pdb",
        "-maxwarn", "2",
    ], capture_output=True, text=True, check=True)

    subprocess.run([
        "gmx", "mdrun",
        "-ntmpi", str(ntmpi),
        "-ntomp", str(ntomp),
        "-deffnm", "md_resolv_equil",
        "-v",
    ], check=True)


if __name__ == "__main__":
    run()
