"""Adaptive resolvation to bring the equilibrated model to the target pressure.

Corresponds to the taylor ``resolvate_1.sh`` -> ``resolvate_2.sh`` -> ``resolvate_final.sh``
sequence. The reusable :func:`resolvation_cycle` (solvate -> manage WAT topology ->
minimize -> equilibrate) is invoked up to three times; the branching between them is the
dynamic control flow that would become a Prefect ``@flow`` body unchanged:

1. run 1: add ``trial_fraction`` of a trial solvation, minimize + equilibrate.
2. read run-1 pressure. Case A (under target): add more water (run-2 probe) and
   interpolate between run 1 and run 2. Case B (overshot): interpolate between the
   pre-resolvation equilibration and run 1.
3. final: apply the interpolated ``dNw`` — positive => add to the chosen reference,
   negative => restart from ``md_equil`` with fewer waters, zero/converged => keep run 1.

WAT topology is managed as ordered blocks (``manage_wat_topology``): each stage appends
its block; the final stage drops the discarded run-2 probe block first.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from pydantic import BaseModel

from ..config import GromacsRunProfile, ResolvateParams, SystemConfig
from ..gmx import (
    add_water,
    checksums,
    gmx_version,
    manage_wat_topology,
    read_mean_pressure,
    run_equil,
    run_min,
)
from ..paths import under
from ..results import StepInputs, StepResult, StepStatus
from .pressure_interpolate import pressure_interpolate

STEP = "resolvate"
FINAL_PREFIX = "md_resolv_final"


class ResolvateInputs(StepInputs):
    md_equil_gro: Path
    md_equil_edr: Path
    water_equil: Path
    posre_top: Path
    model_pdb: Path
    min_mdp: Path
    equil_mdp: Path

    def consumed_paths(self) -> list[Path]:
        return [
            self.md_equil_gro,
            self.md_equil_edr,
            self.water_equil,
            self.posre_top,
            self.model_pdb,
            self.min_mdp,
            self.equil_mdp,
        ]


class ResolvationResult(BaseModel):
    """Result of one resolvation cycle (internal, not a StepResult)."""

    prefix: str
    pdb: Path
    min_gro: Path
    equil_gro: Path
    equil_edr: Path
    nwat_added: int


class ResolvateResult(StepResult):
    @property
    def final_equil_gro(self) -> Path:
        return self.output("final_equil_gro")

    @property
    def final_equil_edr(self) -> Path:
        return self.output("final_equil_edr")


def resolve_inputs(workdir: Path, cfg: SystemConfig) -> ResolvateInputs:
    workdir = Path(workdir)
    mdp_dir = under(workdir, cfg.mdp_dir)
    return ResolvateInputs(
        workdir=workdir,
        md_equil_gro=workdir / "md_equil.gro",
        md_equil_edr=workdir / "md_equil.edr",
        water_equil=workdir / "waterbox" / "water_equil.gro",
        posre_top=workdir / "md_model_posre.top",
        model_pdb=workdir / "md_model.pdb",
        min_mdp=mdp_dir / cfg.minimize.min_mdp,
        equil_mdp=mdp_dir / cfg.equilibrate.equil_mdp,
    )


def check_inputs(inputs: ResolvateInputs) -> None:
    inputs.check_exists(STEP)


def resolvation_cycle(
    *,
    input_gro: Path,
    maxsol: int,
    prefix: str,
    ref_pdb: Path,
    reservoir: Path,
    top: Path,
    min_mdp: Path,
    equil_mdp: Path,
    min_profile: GromacsRunProfile,
    equil_profile: GromacsRunProfile,
    workdir: Path,
    wat_mode: str = "append",
) -> ResolvationResult:
    """One resolvation cycle: solvate -> WAT topology -> minimize -> equilibrate."""
    workdir = Path(workdir)
    aw = add_water(
        input_gro,
        reservoir,
        maxsol,
        workdir / f"{prefix}.pdb",
        workdir,
        gmx_bin=min_profile.gmx_bin,
        log_path=workdir / f"gmx_{prefix}.log",
    )
    manage_wat_topology(top, aw.nwat_added, mode=wat_mode)
    mn = run_min(aw.out_pdb, top, min_mdp, f"{prefix}_min", min_profile, workdir)
    eq = run_equil(
        mn.gro, top, equil_mdp, ref_pdb, f"{prefix}_equil", equil_profile, workdir, maxwarn=2
    )
    return ResolvationResult(
        prefix=prefix,
        pdb=aw.out_pdb,
        min_gro=mn.gro,
        equil_gro=eq.gro,
        equil_edr=eq.edr,
        nwat_added=aw.nwat_added,
    )


def resolvate(
    inputs: ResolvateInputs,
    params: ResolvateParams,
    min_profile: GromacsRunProfile,
    equil_profile: GromacsRunProfile,
    *,
    resume: bool = False,
) -> ResolvateResult:
    check_inputs(inputs)
    wd = inputs.workdir
    gmx_bin = min_profile.gmx_bin

    final_equil_gro = wd / f"{FINAL_PREFIX}_equil.gro"
    final_equil_edr = wd / f"{FINAL_PREFIX}_equil.edr"
    outputs = {
        "posre_top": inputs.posre_top,
        "final_equil_gro": final_equil_gro,
        "final_equil_edr": final_equil_edr,
    }
    consumed = {
        "md_equil_gro": inputs.md_equil_gro,
        "md_equil_edr": inputs.md_equil_edr,
        "water_equil": inputs.water_equil,
        "posre_top": inputs.posre_top,
        "model_pdb": inputs.model_pdb,
        "min_mdp": inputs.min_mdp,
        "equil_mdp": inputs.equil_mdp,
    }

    if resume and final_equil_gro.exists():
        return ResolvateResult(
            step=STEP,
            status=StepStatus.SKIPPED,
            workdir=wd,
            outputs=outputs,
            params=params.model_dump(),
            input_checksums=checksums(consumed),
        )

    def cycle(input_gro: Path, maxsol: int, prefix: str, ref_pdb: Path, wat_mode: str):
        return resolvation_cycle(
            input_gro=input_gro,
            maxsol=maxsol,
            prefix=prefix,
            ref_pdb=ref_pdb,
            reservoir=inputs.water_equil,
            top=inputs.posre_top,
            min_mdp=inputs.min_mdp,
            equil_mdp=inputs.equil_mdp,
            min_profile=min_profile,
            equil_profile=equil_profile,
            workdir=wd,
            wat_mode=wat_mode,
        )

    # Trial solvation just sizes run 1 (no topology change).
    trial = add_water(
        inputs.md_equil_gro,
        inputs.water_equil,
        None,
        wd / "tmp_resolv_trial.pdb",
        wd,
        gmx_bin=gmx_bin,
        log_path=wd / "gmx_resolv_trial.log",
    )
    maxsol1 = int(trial.nwat_added * params.trial_fraction)

    r1 = cycle(inputs.md_equil_gro, maxsol1, "md_resolv1", inputs.model_pdb, "append")
    p1 = read_mean_pressure(r1.equil_edr, wd, gmx_bin=gmx_bin)

    if p1 < params.target_bar:  # Case A: still under-solvated -> probe + bracket r1<->r2
        maxsol2 = round(r1.nwat_added * params.scale_add)
        r2 = cycle(r1.equil_gro, maxsol2, "md_resolv2", inputs.model_pdb, "append")
        itp = pressure_interpolate(
            r1.pdb,
            r1.equil_edr,
            r2.pdb,
            r2.equil_edr,
            wd,
            target=params.target_bar,
            pressure_tol=params.pressure_tol,
            max_dNw_factor=params.max_dNw_factor,
            ref_gro=r1.equil_gro,
            gmx_bin=gmx_bin,
        )
        case = "A"
    else:  # Case B: overshot on run 1 -> bracket md_equil<->r1
        itp = pressure_interpolate(
            inputs.md_equil_gro,
            inputs.md_equil_edr,
            r1.pdb,
            r1.equil_edr,
            wd,
            target=params.target_bar,
            pressure_tol=params.pressure_tol,
            max_dNw_factor=params.max_dNw_factor,
            ref_gro=inputs.md_equil_gro,
            gmx_bin=gmx_bin,
        )
        case = "B"

    _write_interpolation_report(wd / "interpolation.txt", case, p1, itp)

    if itp.converged or itp.dNw == 0:
        # Already at target: run 1 is the final state; publish it under the final names.
        shutil.copy(r1.equil_gro, final_equil_gro)
        shutil.copy(r1.equil_edr, final_equil_edr)
        final_added = 0
    elif itp.dNw > 0:
        cycle(itp.ref_gro, itp.dNw, FINAL_PREFIX, inputs.model_pdb, "drop_last_then_append")
        final_added = itp.dNw
    else:  # dNw < 0: cannot remove water by solvating -> restart from md_equil with fewer
        cycle(
            inputs.md_equil_gro,
            max(0, itp.Nw_target),
            FINAL_PREFIX,
            inputs.model_pdb,
            "drop_last_then_append",
        )
        final_added = itp.Nw_target

    return ResolvateResult(
        step=STEP,
        status=StepStatus.COMPLETED,
        workdir=wd,
        outputs=outputs,
        params=params.model_dump(),
        input_checksums=checksums(consumed),
        run_profiles_used={"resolv_min": min_profile, "resolv_equil": equil_profile},
        tool_versions={"gmx": gmx_version(gmx_bin)},
        metrics={
            "trial_nwat": trial.nwat_added,
            "maxsol1": maxsol1,
            "p1_bar": p1,
            "dNw": itp.dNw,
            "Nw_target": itp.Nw_target,
            "final_added": final_added,
            "converged": int(itp.converged),
        },
        log_paths=[wd / "interpolation.txt"],
    )


def _write_interpolation_report(path: Path, case: str, p1: float, itp) -> None:
    lines = [
        f"case: {case}",
        f"run1 mean pressure: {p1:.1f} bar",
        f"P1={itp.P1:.1f} bar  Nw1={itp.Nw1}",
        f"P2={itp.P2:.1f} bar  Nw2={itp.Nw2}",
        f"dNw={itp.dNw}  Nw_target={itp.Nw_target}  converged={itp.converged}",
    ]
    if itp.warning:
        lines.append(f"WARNING: {itp.warning}")
    path.write_text("\n".join(lines) + "\n")
