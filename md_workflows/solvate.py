"""Solvate the crystal with water and add ions.

Corresponds to run_all.sh line:
  bash scripts/solvate.sh
"""

import re
import subprocess
from pathlib import Path


def run():
    ncopies = _count_copies()
    _solvate_crystal()
    nwat_initial = _read_solvate_nwat()
    _write_topology_header(ncopies)
    ions_pos, ions_neg = _count_ions()
    net_ion_charge = (ions_pos - ions_neg) * ncopies
    print(f"Positive ions, negative ions, net charge: {ions_pos} {ions_neg} {net_ion_charge}")

    nna, ncl = _compute_ion_counts(nwat_initial, net_ion_charge)
    _insert_ions(ncl, nna)
    nwat = _count_final_water()
    _finalize_topology(ncopies, nwat, ncl, nna)

    import shutil
    shutil.copy("xtal_solv_cl_na.pdb", "md_model.pdb")


def _count_copies() -> int:
    nats_one = 0
    with open("prot_dry.pdb") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                nats_one += 1

    nats_cell = 0
    with open("xtal.pdb") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                nats_cell += 1

    ncopies = nats_cell // nats_one
    print(f"There are {ncopies} copies of the asymmetric unit")
    return ncopies


def _solvate_crystal():
    subprocess.run([
        "gmx", "solvate",
        "-cp", "xtal.pdb",
        "-cs", "waterbox/water_equil.gro",
        "-o", "xtal_solv.pdb",
    ], capture_output=True, text=True, check=True)


def _read_solvate_nwat() -> int:
    with open("gmx_solvate.log") as fh:
        for line in fh:
            m = re.search(r"Output configuration contains\s+(\d+)", line)
            if m:
                return int(m.group(1))
    raise RuntimeError("Could not parse water count from gmx_solvate.log")


def _write_topology_header(ncopies: int):
    with open("prot.top") as fh:
        header_lines = []
        for line in fh:
            if "molecules" in line.lower():
                break
            header_lines.append(line)

    with open("md_model.top", "w") as fh:
        fh.writelines(header_lines)


def _count_ions() -> tuple[int, int]:
    ions_pos = 0
    ions_neg = 0
    with open("prot.pdb") as fh:
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


def _insert_ions(ncl: int, nna: int):
    subprocess.run([
        "gmx", "insert-molecules",
        "-f", "xtal_solv.pdb",
        "-ci", "Cl-.pdb",
        "-o", "xtal_solv_cl.pdb",
        "-replace", "SOL",
        "-nmol", str(ncl),
    ], check=True)

    subprocess.run([
        "gmx", "insert-molecules",
        "-f", "xtal_solv_cl.pdb",
        "-ci", "Na+.pdb",
        "-o", "xtal_solv_cl_na.pdb",
        "-replace", "SOL",
        "-nmol", str(nna),
    ], check=True)


def _count_final_water() -> int:
    wat_atoms = 0
    with open("xtal_solv_cl_na.pdb") as fh:
        for line in fh:
            if "WAT" in line:
                wat_atoms += 1
    return wat_atoms // 3


def _finalize_topology(ncopies: int, nwat: int, ncl: int, nna: int):
    topcopy = "system1              1\n"
    topall = topcopy * ncopies

    with open("md_model.top", "a") as fh:
        fh.write("[ molecules ]\n")
        fh.write("; Compound       #mols\n")
        fh.write(topall)
        fh.write(f"WAT              {nwat}\n")
        fh.write(f"Cl-              {ncl}\n")
        fh.write(f"Na+              {nna}\n")


if __name__ == "__main__":
    run()
