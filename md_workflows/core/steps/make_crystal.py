"""Build a crystal supercell from the parameterized protein.

Corresponds to ``make_crystal.sh``: dry the protein, restore the CRYST1 cell, expand
the asymmetric unit to the full unit cell with ChimeraX, set the spacegroup to P1, then
replicate into a supercell with PropPDB.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from ..config import CrystalParams, SystemConfig
from ..gmx import checksums
from ..results import StepInputs, StepResult, StepStatus

STEP = "make_crystal"


class MakeCrystalInputs(StepInputs):
    prot_pdb: Path
    pdb_clean: Path

    def consumed_paths(self) -> list[Path]:
        return [self.prot_pdb, self.pdb_clean]


class MakeCrystalResult(StepResult):
    @property
    def xtal(self) -> Path:
        return self.output("xtal")


def resolve_inputs(workdir: Path, cfg: SystemConfig) -> MakeCrystalInputs:
    workdir = Path(workdir)
    return MakeCrystalInputs(
        workdir=workdir,
        prot_pdb=workdir / "prot.pdb",
        pdb_clean=workdir / "pdb_clean.pdb",
    )


def check_inputs(inputs: MakeCrystalInputs) -> None:
    inputs.check_exists(STEP)


def make_crystal(
    inputs: MakeCrystalInputs,
    params: CrystalParams,
    *,
    resume: bool = False,
) -> MakeCrystalResult:
    from ..gmx import run_tool  # local import keeps module import light for tests

    check_inputs(inputs)
    wd = inputs.workdir
    ix = params.ix
    iy = params.iy if params.iy is not None else ix
    iz = params.iz if params.iz is not None else ix

    prot_dry = wd / "prot_dry.pdb"
    prot_dry_cell = wd / "prot_dry_cell.pdb"
    xtal = wd / "xtal.pdb"

    outputs = {"prot_dry": prot_dry, "prot_dry_cell": prot_dry_cell, "xtal": xtal}
    consumed = {"prot_pdb": inputs.prot_pdb, "pdb_clean": inputs.pdb_clean}

    if resume and xtal.exists():
        return MakeCrystalResult(
            step=STEP,
            status=StepStatus.SKIPPED,
            workdir=wd,
            outputs=outputs,
            params=params.model_dump(),
            input_checksums=checksums(consumed),
        )

    # 1. Dry the protein (strip waters/ions) -> prot_dry.pdb
    run_tool(
        ["pdb4amber", "-i", str(inputs.prot_pdb), "-o", str(prot_dry), "--dry"],
        tool="pdb4amber",
        cwd=wd,
        log_path=wd / "pdb4amber_dry.log",
    )

    # 2. Restore the CRYST1 cell from pdb_clean; drop stray ions / any wrong CRYST1.
    _prepend_cryst1(inputs.pdb_clean, prot_dry)

    # 3. Expand the asymmetric unit to the full unit cell with ChimeraX.
    _expand_unit_cell(wd, prot_dry, prot_dry_cell, params.chimerax_exec)

    # 4. Rewrite the CRYST1 spacegroup to P1 and prepend to the cell PDB.
    _set_p1_spacegroup(prot_dry, prot_dry_cell)

    # 5. Replicate the P1 cell into the requested supercell.
    if ix > 0 or iy > 0 or iz > 0:
        run_tool(
            [
                "PropPDB",
                "-p",
                str(prot_dry_cell),
                "-o",
                str(xtal),
                "-ix",
                str(ix),
                "-iy",
                str(iy),
                "-iz",
                str(iz),
            ],
            tool="PropPDB",
            cwd=wd,
            log_path=wd / "proppdb.log",
        )
    else:
        shutil.copy(prot_dry_cell, xtal)

    return MakeCrystalResult(
        step=STEP,
        status=StepStatus.COMPLETED,
        workdir=wd,
        outputs=outputs,
        params=params.model_dump(),
        input_checksums=checksums(consumed),
        metrics={"ix": ix, "iy": iy, "iz": iz},
    )


def _prepend_cryst1(source_pdb: Path, target_pdb: Path) -> None:
    """Copy CRYST1 from source to the top of target, dropping Na+/Cl-/old CRYST1."""
    cryst1 = ""
    with open(source_pdb) as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                cryst1 = line
                break

    with open(target_pdb) as fh:
        lines = fh.readlines()
    filtered = [
        line
        for line in lines
        if "Na+" not in line and "Cl-" not in line and not line.startswith("CRYST1")
    ]
    with open(target_pdb, "w") as fh:
        fh.write(cryst1)
        fh.writelines(filtered)


def _expand_unit_cell(workdir: Path, dry_pdb: Path, cell_pdb: Path, chimerax_exec: str) -> None:
    from ..gmx import run_tool

    cxc = workdir / "expand.cxc"
    cxc.write_text(
        f"open {dry_pdb}\nchangechains #1 A\nunitcell #1\ncombine #2\nsave {cell_pdb} #3\nquit\n"
    )
    run_tool(
        [chimerax_exec, "--offscreen", "--nogui", str(cxc)],
        tool="ChimeraX",
        cwd=workdir,
        log_path=workdir / "chimerax_expand.log",
    )


def _set_p1_spacegroup(dry_pdb: Path, cell_pdb: Path) -> None:
    """Rewrite the CRYST1 line with a P 1 spacegroup and prepend to the cell PDB."""
    cryst1_p1 = ""
    with open(dry_pdb) as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                cryst1_p1 = line[:55] + "P 1\n"
                break

    with open(cell_pdb) as fh:
        cell_content = fh.read()
    with open(cell_pdb, "w") as fh:
        fh.write(cryst1_p1)
        fh.write(cell_content)
