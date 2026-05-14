"""Re-solvate after equilibration and run a second round of min+equil.

Corresponds to run_all.sh line:
  bash scripts/resolvate.sh
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .pdb_file_processing import resolve_artifacts_dir


def run(
    ntmpi: int = 8,
    ntomp: int = 1,
    *,
    work_dir: Path | None = None,
    artifacts_dir: Path | None = None,
    project_root: Path | None = None,
):
    """Resolvate under ``work_dir`` using project ``artifacts/`` MDP templates."""
    wd = (work_dir or Path.cwd()).resolve()
    ad = resolve_artifacts_dir(wd, project_root=project_root, artifacts_dir=artifacts_dir)
    maxsol = _compute_maxsol(wd)
    _resolvate(wd, maxsol)
    _update_topology(wd, maxsol)
    _minimize(wd, ad, ntmpi, ntomp)
    _equilibrate(wd, ad, ntmpi, ntomp)


def _compute_maxsol(wd: Path) -> int:
    """Do a trial solvation to figure out 25 % fill."""
    result = subprocess.run(
        ["gmx", "solvate",
         "-cp", str(wd / "md_equil.gro"),
         "-cs", str(wd / "waterbox" / "water_equil.gro"),
         "-o", str(wd / "tmp.pdb")],
        capture_output=True,
        text=True,
        cwd=str(wd),
    )
    log_text = result.stdout + result.stderr
    with open(wd / "gmx_solvate.log", "w", encoding="utf-8", errors="replace") as fh:
        fh.write(log_text)

    for line in log_text.splitlines():
        m = re.search(r"Number of solvent molecules:\s+(\d+)", line)
        if m:
            return int(int(m.group(1)) * 0.25)
    raise RuntimeError("Could not parse solvent molecule count from gmx solvate output")


def _resolvate(wd: Path, maxsol: int):
    result = subprocess.run(
        ["gmx", "solvate",
         "-cp", str(wd / "md_equil.gro"),
         "-cs", str(wd / "waterbox" / "water_equil.gro"),
         "-o", str(wd / "md_resolv.pdb"),
         "-maxsol", str(maxsol)],
        capture_output=True,
        text=True,
        cwd=str(wd),
        check=True,
    )
    with open(wd / "gmx_resolvate.log", "w", encoding="utf-8", errors="replace") as fh:
        fh.write(result.stdout + result.stderr)


def _update_topology(wd: Path, maxsol: int):
    with open(wd / "md_model_posre.top", "a", encoding="utf-8", errors="replace") as fh:
        fh.write(f"WAT {maxsol}\n")


def _minimize(wd: Path, artifacts_dir: Path, ntmpi: int, ntomp: int):
    subprocess.run(
        [
            "gmx",
            "grompp",
            "-f",
            str(artifacts_dir / "min.mdp"),
            "-c",
            str(wd / "md_resolv.pdb"),
            "-o",
            str(wd / "md_resolv_min.tpr"),
            "-p",
            str(wd / "md_model_posre.top"),
        ],
        capture_output=True,
        text=True,
        cwd=str(wd),
        check=True,
    )

    subprocess.run(
        [
            "gmx",
            "mdrun",
            "-ntmpi",
            str(ntmpi),
            "-ntomp",
            str(ntomp),
            "-deffnm",
            "md_resolv_min",
            "-v",
        ],
        cwd=str(wd),
        check=True,
    )


def _equilibrate(wd: Path, artifacts_dir: Path, ntmpi: int, ntomp: int):
    subprocess.run(
        [
            "gmx",
            "grompp",
            "-f",
            str(artifacts_dir / "equil.mdp"),
            "-c",
            str(wd / "md_resolv_min.gro"),
            "-o",
            str(wd / "md_resolv_equil.tpr"),
            "-p",
            str(wd / "md_model_posre.top"),
            "-r",
            str(wd / "md_model.pdb"),
            "-maxwarn",
            "2",
        ],
        capture_output=True,
        text=True,
        cwd=str(wd),
        check=True,
    )

    subprocess.run(
        [
            "gmx",
            "mdrun",
            "-ntmpi",
            str(ntmpi),
            "-ntomp",
            str(ntomp),
            "-deffnm",
            "md_resolv_equil",
            "-v",
        ],
        cwd=str(wd),
        check=True,
    )


if __name__ == "__main__":
    run()
