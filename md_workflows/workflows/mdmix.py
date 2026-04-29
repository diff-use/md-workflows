"""Run the MD pipeline in the order defined by ``scripts/run_all.sh``.

Matches the current shell script:

1. ``run_params_gaussian`` (under ``ligand/``, same as ``cd ligand && bash ../run_params_gaussian.sh``)
2. ``param_prot``
3. ``make_crystal``
4. ``make_waterbox``
5. ``solvate``
6. ``minimize``
7. ``equilibrate``
8. ``resolvate``

Adjust defaults via CLI flags below.
"""

from __future__ import annotations

import argparse

from ..equilibrate import run as run_equilibrate
from ..make_crystal import run as run_make_crystal
from ..make_waterbox import run as run_make_waterbox
from ..minimize import run as run_minimize
from ..param_prot import run as run_param_prot
from ..resolvate import run as run_resolvate
from ..run_params_gaussian import run as run_run_params_gaussian
from ..solvate import run as run_solvate


def main(
    *,
    ntomp: int = 26,
    param_pdb_id: str = "6B8X",
    crystal_ix: int = 1,
    crystal_iy: int | None = None,
    crystal_iz: int | None = None,
    chimerax_exec: str = "chimerax",
    ligand_resn: str = "AR6",
    g16root: str = "/Users/mewall/packages",
    gaussian_nproc: int = 8,
    gaussian_pdb_id: str | None = None,
    resolv_ntmpi: int = 8,
    resolv_ntomp: int = 1,
) -> None:
    """Execute workflow stages in ``run_all.sh`` order."""

    run_run_params_gaussian(
        resn=ligand_resn,
        g16root=g16root,
        nproc=gaussian_nproc,
        pdb_id=gaussian_pdb_id,
    )
    run_param_prot(pdb_id=param_pdb_id)
    run_make_crystal(ix=crystal_ix, iy=crystal_iy, iz=crystal_iz, chimerax_exec=chimerax_exec)
    run_make_waterbox(ntomp=ntomp)
    run_solvate()
    run_minimize(ntomp=ntomp)
    run_equilibrate(ntomp=ntomp)
    run_resolvate(ntmpi=resolv_ntmpi, ntomp=resolv_ntomp)


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
        help="PDB ID passed to param_prot (Coordinates file should be <ID>.pdb in cwd)",
    )
    parser.add_argument("--ix", type=int, default=1, help="make_crystal supercell replication (x; also y/z if omitted)")
    parser.add_argument("--iy", type=int, default=None, help="make_crystal y replication (optional)")
    parser.add_argument("--iz", type=int, default=None, help="make_crystal z replication (optional)")
    parser.add_argument("--chimerax-exec", default="chimerax", help="ChimeraX executable for make_crystal")

    parser.add_argument("--ligand-resn", default="AR6", help="Ligand residue name for run_params_gaussian")
    parser.add_argument("--g16root", default="/Users/mewall/packages", help="Gaussian install root")
    parser.add_argument("--gaussian-nproc", type=int, default=8, help="Threads for Gaussian steps")
    parser.add_argument(
        "--gaussian-pdb-id",
        default=None,
        help="Optional RCSB PDB ID: fetch legacy PDB and list hetero residues before Gaussian",
    )

    parser.add_argument("--resolv-ntmpi", type=int, default=8, help="resolvate gmx mdrun -ntmpi")
    parser.add_argument("--resolv-ntomp", type=int, default=1, help="resolvate gmx mdrun -ntomp")

    args = parser.parse_args()
    main(
        ntomp=args.ntomp,
        param_pdb_id=args.param_pdb_id,
        crystal_ix=args.ix,
        crystal_iy=args.iy,
        crystal_iz=args.iz,
        chimerax_exec=args.chimerax_exec,
        ligand_resn=args.ligand_resn,
        g16root=args.g16root,
        gaussian_nproc=args.gaussian_nproc,
        gaussian_pdb_id=args.gaussian_pdb_id,
        resolv_ntmpi=args.resolv_ntmpi,
        resolv_ntomp=args.resolv_ntomp,
    )


if __name__ == "__main__":
    _cli()
