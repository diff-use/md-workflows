"""CLI entrypoints for workflow commands."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import equilibrate
from . import make_crystal
from . import make_waterbox
from . import minimize
from . import param_prot
from . import resolvate
from . import run_params_gaussian
from . import solvate


def _single_command_cli(command: str) -> None:
    parser = argparse.ArgumentParser(description=f"Run {command} workflow")

    if command == "param_prot":
        parser.add_argument("--pdb-id", default="6B8X")
        parser.add_argument(
            "--project-root",
            type=Path,
            default=None,
            help="Directory containing artifacts/ and optional <PDB>.pdb (default: cwd)",
        )
    elif command == "make_crystal":
        parser.add_argument("--ix", type=int, default=1)
        parser.add_argument("--iy", type=int, default=None)
        parser.add_argument("--iz", type=int, default=None)
        parser.add_argument(
            "--work-dir",
            type=Path,
            default=None,
            help="Run directory with prot.pdb, pdb_clean.pdb, etc. (default: cwd)",
        )
    elif command == "make_waterbox":
        parser.add_argument("--ntomp", type=int, default=26)
        parser.add_argument("--work-dir", type=Path, default=None)
        parser.add_argument(
            "--artifacts-dir",
            type=Path,
            default=None,
            help="Override MDP directory (default: <project-root>/artifacts)",
        )
        parser.add_argument(
            "--project-root",
            type=Path,
            default=None,
            help="Project root containing artifacts/ (default: cwd; inferred from --work-dir if possible)",
        )
    elif command == "solvate":
        parser.add_argument("--work-dir", type=Path, default=None)
    elif command == "minimize":
        parser.add_argument("--ntomp", type=int, default=26)
        parser.add_argument("--work-dir", type=Path, default=None)
        parser.add_argument("--artifacts-dir", type=Path, default=None)
        parser.add_argument("--project-root", type=Path, default=None)
    elif command == "equilibrate":
        parser.add_argument("--ntomp", type=int, default=26)
        parser.add_argument("--work-dir", type=Path, default=None)
        parser.add_argument("--artifacts-dir", type=Path, default=None)
        parser.add_argument("--project-root", type=Path, default=None)
    elif command == "resolvate":
        parser.add_argument("--ntmpi", type=int, default=8)
        parser.add_argument("--ntomp", type=int, default=1)
        parser.add_argument("--work-dir", type=Path, default=None)
        parser.add_argument("--artifacts-dir", type=Path, default=None)
        parser.add_argument("--project-root", type=Path, default=None)
    elif command == "run_params_gaussian":
        parser.add_argument("--g16root", default="/Users/mewall/packages")
        parser.add_argument("--nproc", type=int, default=8)
    else:
        raise ValueError(f"Unknown command: {command}")

    args = parser.parse_args(sys.argv[1:])

    if command == "param_prot":
        root = (args.project_root or Path.cwd()).resolve()
        wd = param_prot.run(pdb_id=args.pdb_id, project_root=root)
        print(wd)
    elif command == "make_crystal":
        wd = (args.work_dir or Path.cwd()).resolve()
        make_crystal.run(ix=args.ix, iy=args.iy, iz=args.iz, work_dir=wd)
    elif command == "make_waterbox":
        wd = (args.work_dir or Path.cwd()).resolve()
        pr = args.project_root.resolve() if args.project_root else None
        make_waterbox.run(
            ntomp=args.ntomp,
            work_dir=wd,
            artifacts_dir=args.artifacts_dir,
            project_root=pr,
        )
    elif command == "solvate":
        wd = (args.work_dir or Path.cwd()).resolve()
        solvate.run(work_dir=wd)
    elif command == "minimize":
        wd = (args.work_dir or Path.cwd()).resolve()
        pr = args.project_root.resolve() if args.project_root else None
        minimize.run(
            ntomp=args.ntomp,
            work_dir=wd,
            artifacts_dir=args.artifacts_dir,
            project_root=pr,
        )
    elif command == "equilibrate":
        wd = (args.work_dir or Path.cwd()).resolve()
        pr = args.project_root.resolve() if args.project_root else None
        equilibrate.run(
            ntomp=args.ntomp,
            work_dir=wd,
            artifacts_dir=args.artifacts_dir,
            project_root=pr,
        )
    elif command == "resolvate":
        wd = (args.work_dir or Path.cwd()).resolve()
        pr = args.project_root.resolve() if args.project_root else None
        resolvate.run(
            ntmpi=args.ntmpi,
            ntomp=args.ntomp,
            work_dir=wd,
            artifacts_dir=args.artifacts_dir,
            project_root=pr,
        )
    elif command == "run_params_gaussian":
        run_params_gaussian.run(
            g16root=args.g16root,
            nproc=args.nproc,
        )


def param_prot_cli() -> None:
    _single_command_cli("param_prot")


def make_crystal_cli() -> None:
    _single_command_cli("make_crystal")


def make_waterbox_cli() -> None:
    _single_command_cli("make_waterbox")


def solvate_cli() -> None:
    _single_command_cli("solvate")


def minimize_cli() -> None:
    _single_command_cli("minimize")


def equilibrate_cli() -> None:
    _single_command_cli("equilibrate")


def resolvate_cli() -> None:
    _single_command_cli("resolvate")


def run_params_gaussian_cli() -> None:
    _single_command_cli("run_params_gaussian")
