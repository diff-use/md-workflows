"""Build a crystal supercell from the protein structure.

Corresponds to run_all.sh line:
  bash scripts/make_crystal.sh 1
"""

import os
import re
import subprocess
from pathlib import Path


def run(
    ix: int = 1,
    iy: int | None = None,
    iz: int | None = None,
    chimerax_exec: str = "chimerax",
):
    """
    Args:
        ix: Number of unit-cell replications in x (also used for y, z if iy/iz
            are not given).
        iy: Optional separate y replication count.
        iz: Optional separate z replication count.
        chimerax_exec: Path to the ChimeraX executable.
    """
    if iy is None:
        iy = ix
    if iz is None:
        iz = ix

    workdir = Path.cwd()

    subprocess.run(["pdb4amber", "-i", "prot.pdb", "-o", "prot_dry.pdb", "--dry"], check=True)

    _prepend_cryst1("pdb_clean.pdb", "prot_dry.pdb")
    _expand_unit_cell(workdir, chimerax_exec)
    _set_p1_spacegroup("prot_dry.pdb", "prot_dry_cell.pdb")
    _propagate_crystal(ix, iy, iz)


def _prepend_cryst1(source_pdb: str, target_pdb: str):
    """Copy the CRYST1 line from source_pdb to the top of target_pdb,
    removing any Na+/Cl-/existing CRYST1 lines from target_pdb."""
    cryst1 = ""
    with open(source_pdb) as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                cryst1 = line
                break

    with open(target_pdb) as fh:
        lines = fh.readlines()

    filtered = [
        l for l in lines
        if "Na+" not in l and "Cl-" not in l and not l.startswith("CRYST1")
    ]

    with open(target_pdb, "w") as fh:
        fh.write(cryst1)
        fh.writelines(filtered)


def _expand_unit_cell(workdir: Path, chimerax_exec: str):
    """Use ChimeraX to expand the unit cell."""
    cxc_script = f"""\
open {workdir}/prot_dry.pdb
changechains #1 A
unitcell #1
combine #2
save {workdir}/prot_dry_cell.pdb #3
quit
"""
    with open("expand.cxc", "w") as fh:
        fh.write(cxc_script)
    subprocess.run([chimerax_exec, "--nogui", "expand.cxc"], check=True)


def _set_p1_spacegroup(dry_pdb: str, cell_pdb: str):
    """Rewrite the CRYST1 line with P 1 spacegroup and prepend to cell PDB."""
    with open(dry_pdb) as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                cryst1_p1 = line[:55] + "P 1\n"
                break

    with open("cryst1_p1.pdb", "w") as fh:
        fh.write(cryst1_p1)

    with open(cell_pdb) as fh:
        cell_content = fh.read()

    with open(cell_pdb, "w") as fh:
        fh.write(cryst1_p1)
        fh.write(cell_content)


def _propagate_crystal(ix: int, iy: int, iz: int):
    """Use PropPDB to replicate the unit cell, or just copy if 0."""
    if ix > 0 or iy > 0 or iz > 0:
        subprocess.run([
            "PropPDB", "-p", "prot_dry_cell.pdb", "-o", "xtal.pdb",
            "-ix", str(ix), "-iy", str(iy), "-iz", str(iz),
        ], check=True)
    else:
        import shutil
        shutil.copy("prot_dry_cell.pdb", "xtal.pdb")


if __name__ == "__main__":
    run()
