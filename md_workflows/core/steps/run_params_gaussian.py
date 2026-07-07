"""Ligand parameterization with Gaussian + AmberTools (standalone utility).

Corresponds to ``run_params_gaussian.sh`` (run under ``ligand/``). This is kept off the
main prep chain — wiring the resulting ligand parameters into the protein topology is
part of the deferred full ligand-handling feature.

Ported to be cwd-independent: instead of ``os.chdir(ligand)`` + mutating ``os.environ``,
every external tool runs with ``cwd=lig_dir`` and an explicit ``env`` carrying
``g16root``/``OMP_NUM_THREADS``.

Requires ``<resn>.pdb`` (the ligand coordinates) to already exist in ``ligand/``; the
residue name is auto-detected from the RCSB legacy PDB.
"""

from __future__ import annotations

import os
from pathlib import Path

from ...pdb_file_processing import prepare_pdb_and_resn_files
from ..config import GaussianParams, SystemConfig
from ..exceptions import MDWorkflowError
from ..gmx import checksums, run_tool
from ..results import StepInputs, StepResult, StepStatus

STEP = "run_params_gaussian"


class RunParamsGaussianInputs(StepInputs):
    pdb_id: str

    def consumed_paths(self) -> list[Path]:
        return []  # legacy PDB is downloaded; <resn>.pdb is checked after detection


class RunParamsGaussianResult(StepResult):
    @property
    def frcmod(self) -> Path:
        return self.output("frcmod")

    @property
    def mol2(self) -> Path:
        return self.output("mol2")


def resolve_inputs(workdir: Path, cfg: SystemConfig) -> RunParamsGaussianInputs:
    return RunParamsGaussianInputs(workdir=Path(workdir), pdb_id=cfg.pdb_id)


def check_inputs(inputs: RunParamsGaussianInputs) -> None:
    inputs.check_exists(STEP)


def run_params_gaussian(
    inputs: RunParamsGaussianInputs,
    params: GaussianParams,
    *,
    resume: bool = False,
) -> RunParamsGaussianResult:
    if not params.g16root:
        raise MDWorkflowError("run_params_gaussian: gaussian.g16root must be set")
    check_inputs(inputs)

    lig_dir = inputs.workdir / "ligand"
    lig_dir.mkdir(parents=True, exist_ok=True)
    resn, _ = prepare_pdb_and_resn_files(lig_dir=lig_dir, pdb_id=inputs.pdb_id)

    ligand_pdb = lig_dir / f"{resn}.pdb"
    if not ligand_pdb.exists():
        raise MDWorkflowError(
            f"run_params_gaussian: expected ligand coordinates {ligand_pdb} (residue "
            f"{resn}); provide it before running."
        )

    outputs = {
        "mol2": lig_dir / f"{resn}_resp.mol2",
        "frcmod": lig_dir / f"{resn}_resp.frcmod",
        "parm7": lig_dir / f"{resn}_resp.parm7",
        "rst7": lig_dir / f"{resn}_resp.rst7",
        "pdb": lig_dir / f"{resn}_resp.pdb",
    }
    if resume and outputs["frcmod"].exists():
        return RunParamsGaussianResult(
            step=STEP,
            status=StepStatus.SKIPPED,
            workdir=inputs.workdir,
            outputs=outputs,
            params=params.model_dump(),
        )

    env = _g16_env(params.g16root, params.nproc)

    run_tool(
        [
            "antechamber",
            "-fi",
            "pdb",
            "-fo",
            "gcrt",
            "-i",
            f"{resn}.pdb",
            "-o",
            f"{resn}.gau",
            "-nc",
            str(params.net_charge),
            "-m",
            "1",
        ],
        tool="antechamber",
        cwd=lig_dir,
        env=env,
        log_path=lig_dir / "antechamber_gcrt.log",
    )
    _patch_gaussian_input(lig_dir, resn, params.nproc, params.method)
    _run_gaussian_opt(lig_dir, resn, env)
    _run_gaussian_esp(lig_dir, resn, params.method, env)
    _process_resp_charges(lig_dir, resn, env)
    _build_amber_lib(lig_dir, resn, env)

    return RunParamsGaussianResult(
        step=STEP,
        status=StepStatus.COMPLETED,
        workdir=inputs.workdir,
        outputs=outputs,
        params=params.model_dump(),
        input_checksums=checksums({"ligand_pdb": ligand_pdb}),
        metrics={"net_charge": params.net_charge},
    )


def _g16_env(g16root: str, nproc: int) -> dict[str, str]:
    env = dict(os.environ)
    env["g16root"] = g16root
    env["OMP_NUM_THREADS"] = str(nproc)
    return env


def _patch_gaussian_input(lig_dir: Path, resn: str, nproc: int, method: str) -> None:
    """Insert NprocShared and switch the HF header to the requested method."""
    lines = (lig_dir / f"{resn}.gau").read_text().splitlines(keepends=True)
    patched: list[str] = []
    for line in lines:
        patched.append(line)
        if "Link" in line:
            patched.append(f"%NprocShared={nproc}\n")
    text = "".join(patched).replace("#HF", f"#{method}")
    (lig_dir / "tmp").write_text(text)


def _run_gaussian_opt(lig_dir: Path, resn: str, env: dict[str, str]) -> None:
    """Run the Gaussian geometry optimization if the input changed or is new."""
    gau = lig_dir / f"{resn}.gau"
    log = lig_dir / f"{resn}.log"
    tmp = lig_dir / "tmp"
    if log.exists() and gau.read_text() == tmp.read_text():
        tmp.unlink(missing_ok=True)
        return
    os.replace(tmp, gau)
    run_tool(["g16", f"{resn}.gau"], tool="g16 (opt)", cwd=lig_dir, env=env)


def _run_gaussian_esp(lig_dir: Path, resn: str, method: str, env: dict[str, str]) -> None:
    """Run the Gaussian ESP / CHELPG calculation if the input changed or is new."""
    text = (lig_dir / f"{resn}.gau").read_text()
    text = text.replace("opt", "pop(chelpg,regular)").replace("molecule", "grid")
    resp_gau = lig_dir / f"{resn}_resp.gau"
    resp_log = lig_dir / f"{resn}_resp.log"
    if resp_log.exists() and resp_gau.exists() and resp_gau.read_text() == text:
        return
    resp_gau.write_text(text)
    run_tool(["g16", f"{resn}_resp.gau"], tool="g16 (esp)", cwd=lig_dir, env=env)


def _process_resp_charges(lig_dir: Path, resn: str, env: dict[str, str]) -> None:
    """Derive RESP charges, graft original coordinates, and build GAFF parameters."""
    run_tool(
        [
            "antechamber",
            "-fi",
            "gout",
            "-i",
            f"{resn}_resp.log",
            "-cf",
            f"{resn}_resp.crg",
            "-c",
            "resp",
            "-o",
            f"{resn}_gauss.ac",
            "-fo",
            "ac",
            "-rn",
            resn,
        ],
        tool="antechamber (resp)",
        cwd=lig_dir,
        env=env,
    )
    run_tool(
        [
            "antechamber",
            "-fi",
            "gout",
            "-i",
            f"{resn}_resp.log",
            "-o",
            f"{resn}_gauss.pdb",
            "-fo",
            "pdb",
            "-rn",
            resn,
        ],
        tool="antechamber (pdb)",
        cwd=lig_dir,
        env=env,
    )

    coords = _extract_coords_from_pdb(lig_dir / f"{resn}.pdb")
    _graft_coords_to_ac(lig_dir / f"{resn}_gauss.ac", coords, lig_dir / f"{resn}_resp.ac")
    _correct_charge(lig_dir / f"{resn}_resp.ac")

    run_tool(
        [
            "antechamber",
            "-fi",
            "ac",
            "-i",
            f"{resn}_resp.ac",
            "-fo",
            "mol2",
            "-o",
            f"{resn}_resp.mol2",
            "-rn",
            resn,
        ],
        tool="antechamber (mol2)",
        cwd=lig_dir,
        env=env,
    )
    run_tool(
        ["atomtype", "-i", f"{resn}_resp.ac", "-o", f"{resn}_resp_gaff.ac", "-p", "gaff"],
        tool="atomtype",
        cwd=lig_dir,
        env=env,
    )
    run_tool(
        ["prepgen", "-i", f"{resn}_resp_gaff.ac", "-o", f"{resn}_resp_gaff.prepc", "-f", "car"],
        tool="prepgen",
        cwd=lig_dir,
        env=env,
    )
    run_tool(
        ["parmchk2", "-i", f"{resn}_resp_gaff.prepc", "-o", f"{resn}_resp.frcmod", "-f", "prepc"],
        tool="parmchk2",
        cwd=lig_dir,
        env=env,
    )


def _build_amber_lib(lig_dir: Path, resn: str, env: dict[str, str]) -> None:
    import textwrap

    tleap_input = textwrap.dedent(f"""\
        source leaprc.protein.ff19SB
        source leaprc.gaff
        loadamberparams {resn}_resp.frcmod
        loadamberprep {resn}_resp_gaff.prepc
        lig = loadmol2 {resn}_resp.mol2
        savepdb lig {resn}_resp.pdb
        saveamberparm lig {resn}_resp.parm7 {resn}_resp.rst7
        quit
    """)
    (lig_dir / "tleap_lig.in").write_text(tleap_input)
    run_tool(["tleap", "-f", "tleap_lig.in"], tool="tleap (ligand)", cwd=lig_dir, env=env)


def _extract_coords_from_pdb(pdb_file: Path) -> list[str]:
    coords: list[str] = []
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                coords.append(line[30:53])
    return coords


def _graft_coords_to_ac(ac_file: Path, coords: list[str], out_file: Path) -> None:
    anum = 0
    with open(ac_file) as fh, open(out_file, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                out.write(line[:30] + coords[anum] + line[53:])
                anum += 1
            else:
                out.write(line)


def _correct_charge(ac_file: Path) -> None:
    """Adjust the largest-magnitude charge so the total is an exact integer."""
    lines = ac_file.read_text().splitlines(keepends=True)
    charges = [float(line[54:63]) for line in lines if line.startswith("ATOM")]
    if not charges:
        return
    correction = round(sum(charges)) - sum(charges)
    max_idx = max(range(len(charges)), key=lambda i: abs(charges[i]))
    new_charge = charges[max_idx] + correction

    out: list[str] = []
    anum = 0
    for line in lines:
        if line.startswith("ATOM"):
            if anum == max_idx:
                out.append(f"{line[:54]}{new_charge:9.6f}{line[63:]}")
            else:
                out.append(line)
            anum += 1
        else:
            out.append(line)
    ac_file.write_text("".join(out))
