"""Run the MD pipeline in the order defined by ``scripts/run_all.sh``.

All protein/crystal/solvation inputs and outputs use the per-run directory
``{project_root}/{PDB}/{PDB}_{timestamp}/`` created by ``param_prot``. MDP
templates are read only from ``{project_root}/artifacts/`` at the project root
(not inside the run directory). The entry coordinate file is downloaded from
RCSB into the run directory (or copied from ``{project_root}/{PDB}.pdb`` when
present).

Pipeline stages:

1. ``param_prot``
2. ``make_crystal``
3. ``make_waterbox``
4. ``solvate``
5. ``minimize``
6. ``equilibrate``
7. ``resolvate``

Adjust defaults via CLI flags below.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..equilibrate import run as run_equilibrate
from ..make_crystal import run as run_make_crystal
from ..make_waterbox import run as run_make_waterbox
from ..minimize import run as run_minimize
from ..param_prot import run as run_param_prot
from ..resolvate import run as run_resolvate
from ..solvate import run as run_solvate


def main(
    *,
    ntomp: int = 26,
    param_pdb_id: str = "6B8X",
    crystal_ix: int = 1,
    crystal_iy: int | None = None,
    crystal_iz: int | None = None,
    resolv_ntmpi: int = 8,
    resolv_ntomp: int = 1,
    project_root: Path | None = None,
) -> Path:
    """Execute workflow stages in ``run_all.sh`` order.

    Returns the per-run directory ``{project_root}/{PDB}/{PDB}_{timestamp}/``
    where all pipeline outputs are written.
    """
    root = (project_root or Path.cwd()).resolve()
    artifacts_dir = root / "artifacts"
    work_dir = run_param_prot(pdb_id=param_pdb_id, project_root=root)
    run_make_crystal(ix=crystal_ix, iy=crystal_iy, iz=crystal_iz, work_dir=work_dir)
    run_make_waterbox(
        ntomp=ntomp, work_dir=work_dir, artifacts_dir=artifacts_dir, project_root=root
    )
    run_solvate(work_dir=work_dir)
    run_minimize(ntomp=ntomp, work_dir=work_dir, artifacts_dir=artifacts_dir, project_root=root)
    run_equilibrate(ntomp=ntomp, work_dir=work_dir, artifacts_dir=artifacts_dir, project_root=root)
    run_resolvate(
        ntmpi=resolv_ntmpi,
        ntomp=resolv_ntomp,
        work_dir=work_dir,
        artifacts_dir=artifacts_dir,
        project_root=root,
    )
    return work_dir


def _cli() -> None:
    parser = argparse.ArgumentParser(description="Run full pipeline per scripts/run_all.sh")
    parser.add_argument(
        "--ntomp",
        type=int,
        default=26,
        help="OpenMP threads for GROMACS steps (waterbox, minimize, equilibrate)",
    )
    parser.add_argument(
        "--param-pdb-id",
        default="6B8X",
        help="PDB ID passed to param_prot (Coordinates file should be <ID>.pdb in project root)",
    )
    parser.add_argument("--ix", type=int, default=1, help="make_crystal supercell replication (x; also y/z if omitted)")
    parser.add_argument("--iy", type=int, default=None, help="make_crystal y replication (optional)")
    parser.add_argument("--iz", type=int, default=None, help="make_crystal z replication (optional)")
    parser.add_argument("--resolv-ntmpi", type=int, default=8, help="resolvate gmx mdrun -ntmpi")
    parser.add_argument("--resolv-ntomp", type=int, default=1, help="resolvate gmx mdrun -ntomp")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=None,
        help="Project root (directory containing artifacts/ and optional <PDB>.pdb); default cwd",
    )

    args = parser.parse_args()
    wd = main(
        ntomp=args.ntomp,
        param_pdb_id=args.param_pdb_id,
        crystal_ix=args.ix,
        crystal_iy=args.iy,
        crystal_iz=args.iz,
        resolv_ntmpi=args.resolv_ntmpi,
        resolv_ntomp=args.resolv_ntomp,
        project_root=args.project_root,
    )
    print(f"Pipeline outputs: {wd}")


if __name__ == "__main__":
    _cli()
