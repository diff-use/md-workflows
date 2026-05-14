"""Solvate the crystal with water and add ions.

Corresponds to run_all.sh line:
  bash scripts/solvate.sh
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path


def run(*, work_dir: Path | None = None):
    """Read/write all files under ``work_dir`` (defaults to cwd)."""
    wd = (work_dir or Path.cwd()).resolve()

    ncopies = _count_copies(wd)
    _solvate_crystal(wd)
    nwat_initial = _read_solvate_nwat(wd)
    _write_topology_header(wd, ncopies)
    ions_pos, ions_neg = _count_ions(wd)
    net_ion_charge = (ions_pos - ions_neg) * ncopies
    print(f"Positive ions, negative ions, net charge: {ions_pos} {ions_neg} {net_ion_charge}")

    nna, ncl = _compute_ion_counts(nwat_initial, net_ion_charge)
    _insert_ions(wd, ncl, nna)
    nwat = _count_final_water(wd)
    _finalize_topology(wd, ncopies, nwat, ncl, nna)

    shutil.copy(wd / "xtal_solv_cl_na.pdb", wd / "md_model.pdb")


def _count_copies(wd: Path) -> int:
    nats_one = 0
    with open(wd / "prot_dry.pdb", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                nats_one += 1

    nats_cell = 0
    with open(wd / "xtal.pdb", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                nats_cell += 1

    ncopies = nats_cell // nats_one
    print(f"There are {ncopies} copies of the asymmetric unit")
    return ncopies


def _solvate_crystal(wd: Path):
    # Match scripts/solvate.sh: >& gmx_solvate.log so _read_solvate_nwat() can parse it.
    with open(wd / "gmx_solvate.log", "w", encoding="utf-8", errors="replace") as log_fh:
        subprocess.run(
            [
                "gmx",
                "solvate",
                "-cp",
                str(wd / "xtal.pdb"),
                "-cs",
                str(wd / "waterbox" / "water_equil.gro"),
                "-o",
                str(wd / "xtal_solv.pdb"),
            ],
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            cwd=str(wd),
            check=True,
        )


def _read_solvate_nwat(wd: Path) -> int:
    with open(wd / "gmx_solvate.log", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.search(r"Output configuration contains\s+(\d+)", line)
            if m:
                return int(m.group(1))
    raise RuntimeError("Could not parse water count from gmx_solvate.log")


def _write_topology_header(wd: Path, ncopies: int):
    with open(wd / "prot.top", encoding="utf-8", errors="replace") as fh:
        header_lines = []
        for line in fh:
            if "molecules" in line.lower():
                break
            header_lines.append(line)

    with open(wd / "md_model.top", "w", encoding="utf-8", errors="replace") as fh:
        fh.writelines(header_lines)


def _count_ions(wd: Path) -> tuple[int, int]:
    ions_pos = 0
    ions_neg = 0
    with open(wd / "prot.pdb", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("HETATM"):
                if "Na+" in line:
                    ions_pos += 1
                elif "Cl-" in line:
                    ions_neg += 1
    return ions_pos, ions_neg


def _compute_ion_counts(nwat: int, net_ion_charge: int) -> tuple[int, int]:
    """Compute Na+ and Cl- counts for ~0.1 M ionic strength plus neutralization."""
    if net_ion_charge >= 0:
        ncl = nwat * 0.1 // 55
        nna = ncl + net_ion_charge
    else:
        nna = nwat * 0.1 // 55
        ncl = nna - net_ion_charge
    return int(nna), int(ncl)


def _insert_ions(wd: Path, ncl: int, nna: int):
    subprocess.run(
        [
            "gmx",
            "insert-molecules",
            "-f",
            str(wd / "xtal_solv.pdb"),
            "-ci",
            str(wd / "Cl-.pdb"),
            "-o",
            str(wd / "xtal_solv_cl.pdb"),
            "-replace",
            "SOL",
            "-nmol",
            str(ncl),
        ],
        cwd=str(wd),
        check=True,
    )

    subprocess.run(
        [
            "gmx",
            "insert-molecules",
            "-f",
            str(wd / "xtal_solv_cl.pdb"),
            "-ci",
            str(wd / "Na+.pdb"),
            "-o",
            str(wd / "xtal_solv_cl_na.pdb"),
            "-replace",
            "SOL",
            "-nmol",
            str(nna),
        ],
        cwd=str(wd),
        check=True,
    )


def _count_final_water(wd: Path) -> int:
    wat_atoms = 0
    with open(wd / "xtal_solv_cl_na.pdb", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "WAT" in line:
                wat_atoms += 1
    return wat_atoms // 3


def _finalize_topology(wd: Path, ncopies: int, nwat: int, ncl: int, nna: int):
    topcopy = "system1              1\n"
    topall = topcopy * ncopies

    with open(wd / "md_model.top", "a", encoding="utf-8", errors="replace") as fh:
        fh.write("[ molecules ]\n")
        fh.write("; Compound       #mols\n")
        fh.write(topall)
        fh.write(f"WAT              {nwat}\n")
        fh.write(f"Cl-              {ncl}\n")
        fh.write(f"Na+              {nna}\n")


if __name__ == "__main__":
    run()
