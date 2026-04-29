"""Set up positional restraints and run NPT equilibration.

Corresponds to run_all.sh line:
  bash scripts/equilibrate.sh
"""

import glob
import os
import re
import subprocess
import tempfile
from pathlib import Path


def run(ntomp: int = 26):
    artifacts_dir = Path("artifacts")
    _extract_first_copy()
    chain_files = _split_chains()
    _generate_restraints(chain_files)
    _build_restrained_topology(chain_files)

    subprocess.run([
        "gmx", "grompp",
        "-f", str(artifacts_dir / "equil.mdp"),
        "-c", "md_min.gro",
        "-o", "md_equil.tpr",
        "-p", "md_model_posre.top",
        "-r", "md_model.pdb",
    ], capture_output=True, text=True, check=True)

    subprocess.run([
        "gmx", "mdrun",
        "-ntmpi", "1",
        "-ntomp", str(ntomp),
        "-deffnm", "md_equil",
        "-v",
    ], check=True)


def _extract_first_copy():
    """Extract the first copy of the asymmetric unit (up to first GOL residue)."""
    with open("pdb_clean.pdb") as fh:
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

    with open("first_copy.pdb", "w") as fh:
        fh.writelines(kept)

    prot_lines = [
        l for l in kept
        if (l.startswith(("ATOM", "HETATM", "TER")) and "GOL" not in l)
    ]
    with open("first_copy_prot.pdb", "w") as fh:
        fh.writelines(prot_lines)


def _split_chains() -> list[str]:
    """Split first_copy_prot.pdb at TER records into part00, part01, etc."""
    for f in glob.glob("part??"):
        os.remove(f)

    with open("first_copy_prot.pdb") as fh:
        content = fh.read()

    chunks = re.split(r"(?<=TER.*\n)", content)
    chunks = [c for c in chunks if c.strip()]

    files = []
    for i, chunk in enumerate(chunks):
        fname = f"part{i:02d}"
        with open(fname, "w") as fh:
            fh.write(chunk)
        files.append(fname)
    return files


def _generate_restraints(chain_files: list[str]):
    """Run pdb4amber + gmx genrestr for each chain fragment."""
    for f in chain_files:
        subprocess.run(["pdb4amber", "-i", f, "-o", f"{f}_amber.pdb"], check=True)
        subprocess.run(
            ["gmx", "genrestr", "-fc", "209.2", "209.2", "209.2",
             "-f", f"{f}_amber.pdb", "-o", f"posre_{f}.itp"],
            input="Protein-H\nq\n",
            text=True,
            check=True,
        )


def _build_restrained_topology(chain_files: list[str]):
    """Insert #ifdef POSRES_partXX blocks into md_model_posre.top."""
    import shutil
    shutil.copy("md_model.top", "md_model_posre.top")

    for molnum_offset, f in enumerate(chain_files):
        target_moltype_count = molnum_offset + 2  # 1-indexed, skip first

        with open("md_model_posre.top") as fh:
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

        with open("md_model_posre.top", "w") as fh:
            fh.writelines(new_lines)


if __name__ == "__main__":
    run()
