"""Set up positional restraints and run NVT equilibration.

Corresponds to run_all.sh line:
  bash scripts/equilibrate.sh
"""

from __future__ import annotations

import glob
import os
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
    """Equilibrate under ``work_dir`` using ``<project_root>/artifacts/equil.mdp``."""
    wd = (work_dir or Path.cwd()).resolve()
    ad = resolve_artifacts_dir(wd, project_root=project_root, artifacts_dir=artifacts_dir)
    _extract_first_copy(wd)
    chain_files = _split_chains(wd)
    _generate_restraints(wd, chain_files)
    _build_restrained_topology(wd, chain_files)

    subprocess.run(
        [
            "gmx",
            "grompp",
            "-f",
            str(ad / "equil.mdp"),
            "-c",
            str(wd / "md_min.gro"),
            "-o",
            str(wd / "md_equil.tpr"),
            "-p",
            str(wd / "md_model_posre.top"),
            "-r",
            str(wd / "md_model.pdb"),
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
            "-nb",
            "gpu",
            "-pme", 
            "gpu", 
            "-bonded",
            "gpu",
            "-dlb", 
            "no",
            "-notunepme",
            "-deffnm",
            "md_equil",
            "-v",
        ],
        cwd=str(wd),
        check=True,
    )

def _extract_first_copy(wd: Path):
    """Extract the first copy of the asymmetric unit (up to first GOL residue)."""
    with open(wd / "pdb_clean.pdb", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    lines = [l for l in lines if not l.startswith("JRNL")]

    kept = []
    found_gol = False
    gol_next_idx = None
    for i, line in enumerate(lines):
        if "GOL" in line:
            found_gol = True
            kept.append(line)
        elif found_gol and gol_next_idx is None:
            gol_next_idx = i
            break
        else:
            kept.append(line)

    with open(wd / "first_copy.pdb", "w", encoding="utf-8", errors="replace") as fh:
        fh.writelines(kept)

    prot_lines = [
        l for l in kept
        if (l.startswith(("ATOM", "HETATM", "TER")) and "GOL" not in l)
    ]
    with open(wd / "first_copy_prot.pdb", "w", encoding="utf-8", errors="replace") as fh:
        fh.writelines(prot_lines)


def _split_chains(wd: Path) -> list[str]:
    """Split first_copy_prot.pdb at TER records into part00, part01, etc."""
    for f in glob.glob(str(wd / "part??")):
        os.remove(f)

    with open(wd / "first_copy_prot.pdb", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    # Split after each TER line (PDB chain end). Python re forbids variable-width lookbehinds
    # like (?<=TER.*\n).
    chunks: list[str] = []
    buf: list[str] = []
    for line in lines:
        buf.append(line)
        if line.startswith("TER"):
            chunks.append("".join(buf))
            buf = []
    if buf:
        chunks.append("".join(buf))
    chunks = [c for c in chunks if c.strip()]

    files = []
    for i, chunk in enumerate(chunks):
        fname = f"part{i:02d}"
        with open(wd / fname, "w", encoding="utf-8", errors="replace") as fh:
            fh.write(chunk)
        files.append(fname)
    return files


def _generate_restraints(wd: Path, chain_files: list[str]):
    """Run pdb4amber + gmx genrestr for each chain fragment."""
    for f in chain_files:
        subprocess.run(
            ["pdb4amber", "-i", f, "-o", f"{f}_amber.pdb"],
            cwd=str(wd),
            check=True,
        )
        subprocess.run(
            ["gmx", "genrestr", "-fc", "209.2", "209.2", "209.2",
             "-f", f"{f}_amber.pdb", "-o", f"posre_{f}.itp"],
            input="Protein-H\nq\n",
            text=True,
            cwd=str(wd),
            check=True,
        )


def _build_restrained_topology(wd: Path, chain_files: list[str]):
    """Insert #ifdef POSRES_partXX blocks into md_model_posre.top."""
    import shutil

    shutil.copy(wd / "md_model.top", wd / "md_model_posre.top")

    for molnum_offset, f in enumerate(chain_files):
        target_moltype_count = molnum_offset + 2  # 1-indexed, skip first

        with open(wd / "md_model_posre.top", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()

        new_lines = []
        cnt = 0
        for line in lines:
            if "moleculetype" in line.lower():
                cnt += 1
                if cnt == target_moltype_count:
                    posre_block = (
                        f'#ifdef POSRES_{f}\n'
                        f'#include "posre_{f}.itp"\n'
                        f'#endif\n\n'
                    )
                    new_lines.append(posre_block)
            new_lines.append(line)

        with open(wd / "md_model_posre.top", "w", encoding="utf-8", errors="replace") as fh:
            fh.writelines(new_lines)

#gmx mdrun -ntmpi 1 -ntomp 26 -nb gpu -pme gpu -bonded gpu -dlb no -notunepme -deffnm md_equil -v
if __name__ == "__main__":
    run()
