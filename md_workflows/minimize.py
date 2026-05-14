"""Energy-minimize the solvated MD model.

Corresponds to run_all.sh line:
  bash scripts/minimize.sh
"""

import subprocess
from pathlib import Path

from .pdb_file_processing import resolve_artifacts_dir


def run(
    ntomp: int = 26,
    *,
    work_dir: Path | None = None,
    artifacts_dir: Path | None = None,
    project_root: Path | None = None,
):
    """Minimize under ``work_dir`` using ``<project_root>/artifacts/min.mdp``."""
    wd = (work_dir or Path.cwd()).resolve()
    ad = resolve_artifacts_dir(wd, project_root=project_root, artifacts_dir=artifacts_dir)
    subprocess.run(
        [
            "gmx",
            "grompp",
            "-f",
            str(ad / "min.mdp"),
            "-c",
            str(wd / "md_model.pdb"),
            "-o",
            str(wd / "md_min.tpr"),
            "-p",
            str(wd / "md_model.top"),
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
            "1",
            "-ntomp",
            str(ntomp),
            "-deffnm",
            "md_min",
            "-v",
        ],
        cwd=str(wd),
        check=True,
    )


if __name__ == "__main__":
    run()
