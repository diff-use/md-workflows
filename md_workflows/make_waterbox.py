"""Build an equilibrated bulk-water box matching the crystal unit cell.

Corresponds to run_all.sh line:
  bash scripts/make_waterbox.sh
"""

import re
import subprocess
from pathlib import Path


def run(ntomp: int = 26):
    workdir = Path.cwd()
    wb_dir = workdir / "waterbox"
    wb_dir.mkdir(parents=True, exist_ok=True)

    _extract_cryst1(workdir, wb_dir)
    _create_box_pdb(workdir, wb_dir)
    _insert_water(workdir, wb_dir)
    _expand_waterbox(workdir, wb_dir)
    nwat = _count_water(wb_dir)
    _write_topology(workdir, wb_dir, nwat)
    _minimize_waterbox(workdir, wb_dir, ntomp)
    _equilibrate_waterbox(workdir, wb_dir, ntomp)


def _extract_cryst1(workdir: Path, wb_dir: Path):
    with open(workdir / "xtal.pdb") as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                with open(wb_dir / "cryst1_xtal.pdb", "w") as out:
                    out.write(line)
                return


def _create_box_pdb(workdir: Path, wb_dir: Path):
    """Create a box.pdb with CRYST1 dimensions scaled from Angstroms to nm."""
    with open(workdir / "xtal.pdb") as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                a = float(line[6:15]) / 10.0
                b = float(line[15:24]) / 10.0
                c = float(line[24:33]) / 10.0
                alpha = float(line[33:40])
                beta = float(line[40:47])
                gamma = float(line[47:54])
                cryst1 = f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f}\n"
                with open(wb_dir / "box.pdb", "w") as out:
                    out.write(cryst1)
                return


def _insert_water(workdir: Path, wb_dir: Path):
    subprocess.run([
        "gmx", "insert-molecules",
        "-f", str(wb_dir / "box.pdb"),
        "-ci", str(workdir / "WAT.pdb"),
        "-conc", "58.0",
        "-o", str(wb_dir / "box_solv.pdb"),
    ], capture_output=True, text=True, cwd=str(wb_dir), check=True)


def _expand_waterbox(workdir: Path, wb_dir: Path):
    subprocess.run([
        "PropPDB",
        "-p", str(wb_dir / "box_solv.pdb"),
        "-o", str(wb_dir / "box_solv_expand.pdb"),
        "-ix", "10", "-iy", "10", "-iz", "10",
    ], check=True)

    with open(wb_dir / "cryst1_xtal.pdb") as fh:
        cryst1 = fh.read()

    with open(wb_dir / "box_solv_expand.pdb") as fh:
        lines = [
            l for l in fh
            if not l.startswith(("CRYST1", "HEADER"))
        ]

    with open(wb_dir / "box_solv_expand.pdb", "w") as fh:
        fh.write(cryst1)
        fh.writelines(lines)


def _count_water(wb_dir: Path) -> int:
    log = wb_dir / "insert-molecules.log"
    if log.exists():
        with open(log) as fh:
            for line in fh:
                m = re.search(r"Output configuration contains\s+(\d+)", line)
                if m:
                    return int(m.group(1)) * 1000

    with open(wb_dir / "box_solv.pdb") as fh:
        wat_atoms = sum(1 for l in fh if "WAT" in l)
    return wat_atoms // 3


def _write_topology(workdir: Path, wb_dir: Path, nwat: int):
    with open(workdir / "prot.top") as fh:
        header_lines = []
        for line in fh:
            if "molecules" in line.lower():
                break
            header_lines.append(line)

    with open(wb_dir / "waterbox.top", "w") as fh:
        fh.writelines(header_lines)
        fh.write("[ molecules ]\n")
        fh.write("; Compound       #mols\n")
        fh.write(f"WAT              {nwat}\n")


def _minimize_waterbox(workdir: Path, wb_dir: Path, ntomp: int):
    subprocess.run([
        "gmx", "grompp",
        "-f", str(workdir / "min_water.mdp"),
        "-c", str(wb_dir / "box_solv_expand.pdb"),
        "-o", str(wb_dir / "water_min.tpr"),
        "-p", str(wb_dir / "waterbox.top"),
    ], cwd=str(wb_dir), check=True)

    subprocess.run([
        "gmx", "mdrun",
        "-ntmpi", "1", "-ntomp", str(ntomp),
        "-deffnm", "water_min", "-v",
    ], cwd=str(wb_dir), check=True)


def _equilibrate_waterbox(workdir: Path, wb_dir: Path, ntomp: int):
    subprocess.run([
        "gmx", "grompp",
        "-f", str(workdir / "equil_water.mdp"),
        "-c", str(wb_dir / "water_min.gro"),
        "-o", str(wb_dir / "water_equil.tpr"),
        "-p", str(wb_dir / "waterbox.top"),
    ], cwd=str(wb_dir), check=True)

    subprocess.run([
        "gmx", "mdrun",
        "-ntmpi", "1", "-ntomp", str(ntomp),
        "-deffnm", "water_equil", "-v",
    ], cwd=str(wb_dir), check=True)


if __name__ == "__main__":
    run()
