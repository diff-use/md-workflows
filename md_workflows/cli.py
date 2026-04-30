"""CLI entrypoints for workflow commands."""

from __future__ import annotations

import argparse
import sys

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
    elif command == "make_crystal":
        parser.add_argument("--ix", type=int, default=1)
        parser.add_argument("--iy", type=int, default=None)
        parser.add_argument("--iz", type=int, default=None)
        parser.add_argument("--chimerax-exec", default="chimerax")
    elif command == "make_waterbox":
        parser.add_argument("--ntomp", type=int, default=26)
    elif command == "solvate":
        pass
    elif command == "minimize":
        parser.add_argument("--ntomp", type=int, default=26)
    elif command == "equilibrate":
        parser.add_argument("--ntomp", type=int, default=26)
    elif command == "resolvate":
        parser.add_argument("--ntmpi", type=int, default=8)
        parser.add_argument("--ntomp", type=int, default=1)
    elif command == "run_params_gaussian":
        parser.add_argument("--g16root", default="/Users/mewall/packages")
        parser.add_argument("--nproc", type=int, default=8)
        parser.add_argument(
            "--pdb-id",
            required=True,
            help=(
                "RCSB PDB ID: download legacy .pdb if missing, detect ligand resn, "
                "and create <resn>.pdb if missing."
            ),
        )
    else:
        raise ValueError(f"Unknown command: {command}")

    args = parser.parse_args(sys.argv[1:])

    if command == "param_prot":
        param_prot.run(pdb_id=args.pdb_id)
    elif command == "make_crystal":
        make_crystal.run(ix=args.ix, iy=args.iy, iz=args.iz, chimerax_exec=args.chimerax_exec)
    elif command == "make_waterbox":
        make_waterbox.run(ntomp=args.ntomp)
    elif command == "solvate":
        solvate.run()
    elif command == "minimize":
        minimize.run(ntomp=args.ntomp)
    elif command == "equilibrate":
        equilibrate.run(ntomp=args.ntomp)
    elif command == "resolvate":
        resolvate.run(ntmpi=args.ntmpi, ntomp=args.ntomp)
    elif command == "run_params_gaussian":
        run_params_gaussian.run(
            g16root=args.g16root,
            nproc=args.nproc,
            pdb_id=args.pdb_id,
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
