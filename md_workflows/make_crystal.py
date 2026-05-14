"""Build a crystal supercell from the protein structure.

Corresponds to run_all.sh line:
  bash scripts/make_crystal.sh 1
"""

import shutil
import string
import subprocess
from pathlib import Path

import gemmi


def _chain_name_iter():
    """Yield A..Z, AA..ZZ for unique PDB chain identifiers."""
    letters = string.ascii_uppercase
    for c in letters:
        yield c
    for c1 in letters:
        for c2 in letters:
            yield c1 + c2


def run(
    ix: int = 1,
    iy: int | None = None,
    iz: int | None = None,
    *,
    work_dir: Path | None = None,
):
    """
    Args:
        ix: Number of unit-cell replications in x (also used for y, z if iy/iz
            are not given).
        iy: Optional separate y replication count.
        iz: Optional separate z replication count.
        work_dir: Directory containing inputs/outputs (defaults to cwd).
    """
    wd = work_dir or Path.cwd()
    wd = wd.resolve()

    if iy is None:
        iy = ix
    if iz is None:
        iz = ix

    subprocess.run(
        [
            "pdb4amber",
            "-i",
            str(wd / "prot.pdb"),
            "-o",
            str(wd / "prot_dry.pdb"),
            "--dry",
        ],
        cwd=str(wd),
        check=True,
    )

    _prepend_cryst1(wd / "pdb_clean.pdb", wd / "prot_dry.pdb")
    _expand_unit_cell(wd / "prot_dry.pdb", wd / "prot_dry_cell.pdb")
    _set_p1_spacegroup(wd / "prot_dry.pdb", wd / "prot_dry_cell.pdb")
    _propagate_crystal(ix, iy, iz, wd)


def _prepend_cryst1(source_pdb: Path, target_pdb: Path):
    """Copy the CRYST1 line from source_pdb to the top of target_pdb,
    removing any Na+/Cl-/existing CRYST1 lines from target_pdb."""
    cryst1 = ""
    with open(source_pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                cryst1 = line
                break

    with open(target_pdb, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    filtered = [
        l for l in lines
        if "Na+" not in l and "Cl-" not in l and not l.startswith("CRYST1")
    ]

    with open(target_pdb, "w", encoding="utf-8", errors="replace") as fh:
        fh.write(cryst1)
        fh.writelines(filtered)


def _expand_unit_cell(input_pdb: Path, output_pdb: Path):
    """Expand symmetry mates into a P1 unit cell.

    Uses gemmi only for the stable parts of its API (UnitCell, spacegroup
    lookup, Op application, fractional/orthogonal conversion) and rewrites
    PDB ATOM/HETATM lines directly. This avoids gemmi.Model / gemmi.Residue
    / gemmi.Atom construction differences across gemmi versions.

    The output file does not include a CRYST1 record; ``_set_p1_spacegroup``
    runs after this and prepends the correct P1 CRYST1 line.
    """
    with open(input_pdb, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    cryst1_line = next((line for line in lines if line.startswith("CRYST1")), None)
    if cryst1_line is None:
        raise RuntimeError(f"No CRYST1 record found in {input_pdb!s}")

    a = float(cryst1_line[6:15])
    b = float(cryst1_line[15:24])
    c = float(cryst1_line[24:33])
    alpha = float(cryst1_line[33:40])
    beta = float(cryst1_line[40:47])
    gamma = float(cryst1_line[47:54])
    sg_field = cryst1_line[55:66].strip() if len(cryst1_line) >= 66 else cryst1_line[55:].strip()

    cell = gemmi.UnitCell(a, b, c, alpha, beta, gamma)
    spacegroup = gemmi.find_spacegroup_by_name(sg_field) if sg_field else None
    if spacegroup is None:
        raise RuntimeError(
            f"Could not look up spacegroup '{sg_field}' from CRYST1 in {input_pdb!s}"
        )

    operations = list(spacegroup.operations())
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]

    chain_ids = _chain_name_iter()
    out_lines: list[str] = []

    for op in operations:
        chain_id = next(chain_ids)
        if len(chain_id) != 1:
            raise RuntimeError(
                "Too many symmetry copies for single-character PDB chain ids "
                f"(next would be {chain_id!r})"
            )
        for line in atom_lines:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            frac = cell.fractionalize(gemmi.Position(x, y, z))
            triple = op.apply_to_xyz([frac.x, frac.y, frac.z])
            new_pos = cell.orthogonalize(
                gemmi.Fractional(triple[0], triple[1], triple[2])
            )
            tail = line[54:] if line.endswith("\n") else line[54:] + "\n"
            new_line = (
                line[:21]
                + chain_id
                + line[22:30]
                + f"{new_pos.x:8.3f}{new_pos.y:8.3f}{new_pos.z:8.3f}"
                + tail
            )
            out_lines.append(new_line)

    with open(output_pdb, "w", encoding="utf-8", errors="replace") as fh:
        fh.writelines(out_lines)
        fh.write("END\n")


def _set_p1_spacegroup(dry_pdb: Path, cell_pdb: Path):
    """Rewrite the CRYST1 line with P 1 spacegroup and prepend to cell PDB."""
    with open(dry_pdb, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                cryst1_p1 = line[:55] + "P 1\n"
                break

    wd = dry_pdb.parent
    with open(wd / "cryst1_p1.pdb", "w", encoding="utf-8", errors="replace") as fh:
        fh.write(cryst1_p1)

    with open(cell_pdb, encoding="utf-8", errors="replace") as fh:
        cell_content = fh.read()

    with open(cell_pdb, "w", encoding="utf-8", errors="replace") as fh:
        fh.write(cryst1_p1)
        fh.write(cell_content)


def _propagate_crystal(ix: int, iy: int, iz: int, work_dir: Path):
    """Use PropPDB to replicate the unit cell, or just copy if 0."""
    if ix > 0 or iy > 0 or iz > 0:
        subprocess.run(
            [
                "PropPDB",
                "-p",
                str(work_dir / "prot_dry_cell.pdb"),
                "-o",
                str(work_dir / "xtal.pdb"),
                "-ix",
                str(ix),
                "-iy",
                str(iy),
                "-iz",
                str(iz),
            ],
            cwd=str(work_dir),
            check=True,
        )
    else:
        shutil.copy(work_dir / "prot_dry_cell.pdb", work_dir / "xtal.pdb")


if __name__ == "__main__":
    run()
