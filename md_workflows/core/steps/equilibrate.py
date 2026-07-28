"""Set up position restraints and run the NPT equilibration.

Corresponds to ``equilibrate.sh``: isolate one asymmetric-unit copy, split it into
chains, generate per-chain position restraints, weave ``#ifdef POSRES_*`` blocks into
the topology, then equilibrate.

Ligand handling is deliberately minimal here (full support is a later feature):

* ``chain_split_mode="ignore_ligand"`` (default) restrains protein chains only (taylor
  behavior).
* ``chain_split_mode="split_intelligently"`` isolates copy 1 up to the ligand before
  splitting (original behavior, generalized to any ``ligand_resname``).
* ``handle_ligand_restraint`` (default ``False``) additionally restrains the ligand via
  a ``make_ndx`` group — provisional and off by default.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from ...pdb_file_processing import find_ligands_in_legacy_pdb_file
from ..config import EquilibrateParams, GromacsRunProfile, SystemConfig
from ..exceptions import MDWorkflowError
from ..gmx import checksums, gmx_version, run_equil, run_tool
from ..paths import under
from ..results import StepInputs, StepResult, StepStatus

STEP = "equilibrate"
DEFFNM = "md_equil"


class EquilibrateInputs(StepInputs):
    pdb_clean: Path
    model_top: Path
    model_pdb: Path
    min_gro: Path
    equil_mdp: Path

    def consumed_paths(self) -> list[Path]:
        return [self.pdb_clean, self.model_top, self.model_pdb, self.min_gro, self.equil_mdp]


class EquilibrateResult(StepResult):
    @property
    def gro(self) -> Path:
        return self.output("gro")

    @property
    def edr(self) -> Path:
        return self.output("edr")

    @property
    def posre_top(self) -> Path:
        return self.output("posre_top")


def resolve_inputs(workdir: Path, cfg: SystemConfig) -> EquilibrateInputs:
    workdir = Path(workdir)
    mdp_dir = under(workdir, cfg.mdp_dir)
    return EquilibrateInputs(
        workdir=workdir,
        pdb_clean=workdir / "pdb_clean.pdb",
        model_top=workdir / "md_model.top",
        model_pdb=workdir / "md_model.pdb",
        min_gro=workdir / "md_min.gro",
        equil_mdp=mdp_dir / cfg.equilibrate.equil_mdp,
    )


def check_inputs(inputs: EquilibrateInputs) -> None:
    inputs.check_exists(STEP)


def equilibrate(
    inputs: EquilibrateInputs,
    params: EquilibrateParams,
    profile: GromacsRunProfile,
    *,
    resume: bool = False,
) -> EquilibrateResult:
    check_inputs(inputs)
    wd = inputs.workdir

    posre_top = wd / "md_model_posre.top"
    tpr = wd / f"{DEFFNM}.tpr"
    gro = wd / f"{DEFFNM}.gro"
    edr = wd / f"{DEFFNM}.edr"
    outputs = {"posre_top": posre_top, "tpr": tpr, "gro": gro, "edr": edr}
    consumed = {
        "pdb_clean": inputs.pdb_clean,
        "model_top": inputs.model_top,
        "model_pdb": inputs.model_pdb,
        "min_gro": inputs.min_gro,
        "equil_mdp": inputs.equil_mdp,
    }

    if resume and gro.exists() and posre_top.exists():
        return EquilibrateResult(
            step=STEP,
            status=StepStatus.SKIPPED,
            workdir=wd,
            outputs=outputs,
            params=params.model_dump(),
            input_checksums=checksums(consumed),
        )

    _clean_stale(wd)
    _build_first_copy(wd, inputs.pdb_clean, params)
    chain_files = _split_chains(wd)
    _generate_restraints(wd, chain_files, params)

    shutil.copy(inputs.model_top, posre_top)
    _weave_posre_includes(posre_top, chain_files)

    if params.handle_ligand_restraint:
        _add_ligand_restraint(wd, inputs.pdb_clean, posre_top, params, len(chain_files))

    run = run_equil(
        inputs.min_gro, posre_top, inputs.equil_mdp, inputs.model_pdb, DEFFNM, profile, wd
    )

    return EquilibrateResult(
        step=STEP,
        status=StepStatus.COMPLETED,
        workdir=wd,
        outputs={"posre_top": posre_top, "tpr": run.tpr, "gro": run.gro, "edr": run.edr},
        params=params.model_dump(),
        input_checksums=checksums(consumed),
        run_profiles_used={"equil": profile},
        tool_versions={"gmx": gmx_version(profile.gmx_bin)},
        metrics={"n_chains": len(chain_files)},
        log_paths=[run.log],
    )


def _clean_stale(workdir: Path) -> None:
    """Remove leftover chain fragments and restraint files from a prior run."""
    for pattern in ("part??", "part??_amber.pdb", "posre_part??.itp"):
        for f in workdir.glob(pattern):
            f.unlink()


def _build_first_copy(workdir: Path, pdb_clean: Path, params: EquilibrateParams) -> None:
    """Write ``first_copy_prot.pdb`` — the protein atoms to restrain."""
    lines = [
        ln for ln in pdb_clean.read_text().splitlines(keepends=True) if not ln.startswith("JRNL")
    ]

    if params.chain_split_mode == "ignore_ligand":
        prot = [ln for ln in lines if ln.startswith(("ATOM", "TER"))]
        (workdir / "first_copy_prot.pdb").write_text("".join(prot))
        return

    # split_intelligently: isolate copy 1 up to the first ligand residue, then drop it.
    resn = _resolve_ligand_resname(params, pdb_clean)
    kept: list[str] = []
    found = False
    for line in lines:
        if resn in line:
            found = True
            kept.append(line)
        elif found:
            break
        else:
            kept.append(line)
    prot = [ln for ln in kept if ln.startswith(("ATOM", "HETATM", "TER")) and resn not in ln]
    (workdir / "first_copy_prot.pdb").write_text("".join(prot))


def _split_chains(workdir: Path) -> list[str]:
    """Split ``first_copy_prot.pdb`` after each TER into part00, part01, ... files."""
    lines = (workdir / "first_copy_prot.pdb").read_text().splitlines(keepends=True)
    chunks: list[str] = []
    buf: list[str] = []
    for line in lines:
        buf.append(line)
        if line.startswith("TER"):
            chunks.append("".join(buf))
            buf = []
    if buf:
        chunks.append("".join(buf))

    files: list[str] = []
    for i, chunk in enumerate(c for c in chunks if c.strip()):
        # Skip fragments with no atoms (empty-part guard, taylor).
        if not any(ln.startswith("ATOM") for ln in chunk.splitlines()):
            continue
        name = f"part{i:02d}"
        (workdir / name).write_text(chunk)
        files.append(name)
    return files


def _generate_restraints(workdir: Path, chain_files: list[str], params: EquilibrateParams) -> None:
    fc = str(params.restraint_fc)
    for f in chain_files:
        run_tool(
            ["pdb4amber", "-i", f, "-o", f"{f}_amber.pdb"],
            tool="pdb4amber",
            cwd=workdir,
            log_path=workdir / f"pdb4amber_{f}.log",
        )
        run_tool(
            ["gmx", "genrestr", "-fc", fc, fc, fc, "-f", f"{f}_amber.pdb", "-o", f"posre_{f}.itp"],
            tool="gmx genrestr",
            cwd=workdir,
            input=f"{params.restraint_group}\nq\n",
            log_path=workdir / f"genrestr_{f}.log",
        )


def _weave_posre_includes(posre_top: Path, chain_files: list[str]) -> None:
    """Insert ``#ifdef POSRES_partXX`` include blocks before each chain's moleculetype."""
    for offset, f in enumerate(chain_files):
        target = offset + 2  # 1-indexed, skipping the first moleculetype
        lines = posre_top.read_text().splitlines(keepends=True)
        out: list[str] = []
        count = 0
        for line in lines:
            if "moleculetype" in line.lower():
                count += 1
                if count == target:
                    out.append(f'#ifdef POSRES_{f}\n#include "posre_{f}.itp"\n#endif\n\n')
            out.append(line)
        posre_top.write_text("".join(out))


def _add_ligand_restraint(
    workdir: Path, pdb_clean: Path, posre_top: Path, params: EquilibrateParams, n_parts: int
) -> None:
    """Provisional ligand restraint via a make_ndx group (off by default)."""
    resn = _resolve_ligand_resname(params, pdb_clean)
    f = f"part{n_parts:02d}"
    lig_lines = [
        ln
        for ln in pdb_clean.read_text().splitlines(keepends=True)
        if ln.startswith(("ATOM", "HETATM")) and resn in ln
    ]
    if not lig_lines:
        return
    (workdir / f).write_text("".join(lig_lines))
    run_tool(["pdb4amber", "-i", f, "-o", f"{f}_amber.pdb"], tool="pdb4amber", cwd=workdir)
    run_tool(
        ["gmx", "make_ndx", "-f", f"{f}_amber.pdb", "-o", f"{f}_amber.ndx"],
        tool="gmx make_ndx",
        cwd=workdir,
        input=f"! a H*\nname 3 {resn}-H\nq\n",
    )
    fc = str(params.restraint_fc)
    run_tool(
        [
            "gmx",
            "genrestr",
            "-fc",
            fc,
            fc,
            fc,
            "-f",
            f"{f}_amber.pdb",
            "-o",
            f"posre_{f}.itp",
            "-n",
            f"{f}_amber.ndx",
        ],
        tool="gmx genrestr",
        cwd=workdir,
        input=f"{resn}-H\nq\n",
    )


def _resolve_ligand_resname(params: EquilibrateParams, pdb_clean: Path) -> str:
    if params.ligand_resname:
        return params.ligand_resname
    hits = find_ligands_in_legacy_pdb_file(pdb_clean)
    unique = sorted({h["resname"] for h in hits})
    if len(unique) == 1:
        return unique[0]
    raise MDWorkflowError(
        "equilibrate: ligand_resname is required for this mode but could not be "
        f"auto-detected (found {unique or 'none'}); set equilibrate.ligand_resname."
    )
