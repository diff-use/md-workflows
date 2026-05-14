"""Parameterize the protein: clean PDB, build Amber topology, convert to GROMACS.

Creates a per-run directory ``{project_root}/{PDB}/{PDB}_{timestamp}/``. All
intermediate and output files from this step go there. The RCSB entry PDB (or a
copy of ``{PDB}.pdb`` from ``project_root``) is placed in that directory.

Corresponds to run_all.sh line:
  bash scripts/param_prot.sh 6B8X
"""

from __future__ import annotations

import shutil
import subprocess
import textwrap
from datetime import datetime
from pathlib import Path

from .pdb_file_processing import ensure_entry_pdb_file


def _make_run_directory(pdb_id: str, project_root: Path) -> Path:
    """Create ``{project_root}/{PDB}/{PDB}_{timestamp}/`` and return the child path."""
    nid = pdb_id.strip()
    if len(nid) != 4 or not nid.isalnum():
        raise ValueError(f"Expected a valid 4-letter PDB ID, got {pdb_id!r}")
    nid_u = nid.upper()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    parent = project_root / nid_u
    parent.mkdir(parents=True, exist_ok=True)
    run_dir = parent / f"{nid_u}_{stamp}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir.resolve()


def run(pdb_id: str = "6B8X", *, project_root: Path | None = None) -> Path:
    """Create run directory, parameterize protein; return path to the run directory.

    All outputs are written under ``{project_root}/{PDB}/{PDB}_{timestamp}/``.
    ``project_root`` defaults to the current working directory.
    """
    root = (project_root or Path.cwd()).resolve()
    work_dir = _make_run_directory(pdb_id, root)
    _clean_pdb(pdb_id, work_dir, root)
    _initial_solvation(work_dir)
    _amber_to_gromacs(work_dir)
    _extract_solvent_pdbs(work_dir)
    _final_protein_topology(work_dir)
    _amber_to_gromacs(work_dir)
    return work_dir


def _clean_pdb(pdb_id: str, work_dir: Path, project_root: Path):
    """Strip REMARK/KEYWDS, keep CRYST1/ATOM/HETATM/TER/END lines."""
    nid = pdb_id.strip().upper()
    src = project_root / f"{nid}.pdb"
    dst = work_dir / f"{nid}.pdb"
    if src.is_file() and not dst.exists():
        shutil.copy2(src, dst)

    pdb_file = str(ensure_entry_pdb_file(pdb_id, work_dir))
    with open(pdb_file) as fh:
        lines = fh.readlines()

    kept = []
    for line in lines:
        if line.startswith(("REMARK", "KEYWDS")):
            continue
        if line.startswith(("CRYST1", "ATOM", "HETATM", "TER", "END")):
            kept.append(line)

    with open(work_dir / "pdb_clean.pdb", "w") as fh:
        fh.writelines(kept)

    subprocess.run(
        [
            "pdb4amber",
            "-i",
            str(work_dir / "pdb_clean.pdb"),
            "--prot",
            "-o",
            str(work_dir / "pdb_clean_amber.pdb"),
        ],
        cwd=str(work_dir),
        check=True,
    )


def _initial_solvation(work_dir: Path):
    """Run tleap with minimal solvation to obtain solvent PDB templates."""
    tleap_input = textwrap.dedent("""\
        source leaprc.protein.ff19SB
        source leaprc.DNA.OL15
        source leaprc.RNA.OL3
        source leaprc.water.spceb
        source leaprc.gaff2
        p = loadpdb pdb_clean_amber.pdb
        x = combine{p}
        addions2 x Cl- 0
        addions2 x Na+ 0
        addions2 x Na+ 1
        addions2 x Cl- 1
        solvateBox x SPCBOX 1.
        set default PBradii mbondi3
        set default nocenter on
        saveAmberParm x prot.parm7 prot.rst7
        quit
    """)
    with open(work_dir / "tleap_temp.in", "w") as fh:
        fh.write(tleap_input)
    subprocess.run(["tleap", "-f", "tleap_temp.in"], cwd=str(work_dir), check=True)

    for f in ["prot.top", "prot.pdb"]:
        (work_dir / f).unlink(missing_ok=True)


def _amber_to_gromacs(work_dir: Path):
    """Convert Amber parm7/rst7 to GROMACS topology and PDB."""
    script = textwrap.dedent("""\
        import parmed as pmd
        parm = pmd.load_file("prot.parm7", "prot.rst7")
        parm.save("prot.top")
        parm.save("prot.pdb")
    """)
    with open(work_dir / "amber_to_gromacs.py", "w") as fh:
        fh.write(script)
    subprocess.run(["python", "amber_to_gromacs.py"], cwd=str(work_dir), check=True)


def _extract_solvent_pdbs(work_dir: Path):
    """Pull single Na+, Cl-, and WAT coordinate templates from the solvated PDB."""
    with open(work_dir / "prot.pdb") as fh:
        lines = fh.readlines()

    hetatm_lines = [l for l in lines if l.startswith("HETATM")]

    na_lines = [l for l in hetatm_lines if "Na+" in l]
    cl_lines = [l for l in hetatm_lines if "Cl-" in l]
    wat_lines = [l for l in hetatm_lines if "WAT" in l]

    if na_lines:
        with open(work_dir / "Na+.pdb", "w") as fh:
            fh.write(na_lines[0])
    if cl_lines:
        with open(work_dir / "Cl-.pdb", "w") as fh:
            fh.write(cl_lines[0])
    if wat_lines:
        with open(work_dir / "WAT.pdb", "w") as fh:
            fh.writelines(wat_lines[:3])


def _final_protein_topology(work_dir: Path):
    """Build the final protein topology with one water + ions."""
    tleap_input = textwrap.dedent("""\
        source leaprc.protein.ff19SB
        source leaprc.DNA.OL15
        source leaprc.RNA.OL3
        source leaprc.water.spceb
        source leaprc.gaff2
        p = loadpdb pdb_clean_amber.pdb
        w = loadpdb WAT.pdb
        x = combine{p w}
        addions2 x Na+ 0
        addions2 x Cl- 0
        addions2 x Na+ 1
        addions2 x Cl- 1
        set default PBradii mbondi3
        set default nocenter on
        saveAmberParm x prot.parm7 prot.rst7
        quit
    """)
    with open(work_dir / "tleap_prot.in", "w") as fh:
        fh.write(tleap_input)
    subprocess.run(["tleap", "-f", "tleap_prot.in"], cwd=str(work_dir), check=True)

    for f in ["prot.top", "prot.pdb"]:
        (work_dir / f).unlink(missing_ok=True)


if __name__ == "__main__":
    wd = run()
    print(f"Run directory: {wd}")
