"""The standard crystalline-MD prep pipeline, as plain function composition.

``param_prot -> make_crystal -> make_waterbox -> solvate -> minimize -> equilibrate ->
resolvate`` (the last stage runs the adaptive pressure loop). Every step reads/writes
the shared ``workdir`` under fixed conventional names, so each stage's inputs are
resolved from that workdir; with ``resume=True`` a stage whose durable outputs already
exist is skipped.

This function is the literal thing a Prefect ``@flow`` would wrap later.
"""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel

from ..config import SystemConfig
from ..results import StepResult
from ..steps import (
    equilibrate as equilibrate_step,
)
from ..steps import (
    make_crystal as make_crystal_step,
)
from ..steps import (
    make_waterbox as make_waterbox_step,
)
from ..steps import (
    minimize as minimize_step,
)
from ..steps import (
    param_prot as param_prot_step,
)
from ..steps import (
    resolvate as resolvate_step,
)
from ..steps import (
    solvate as solvate_step,
)


class PipelineResult(BaseModel):
    workdir: Path
    steps: list[StepResult]

    @property
    def final(self) -> StepResult:
        return self.steps[-1]


def standard_md_pipeline(
    workdir: Path,
    cfg: SystemConfig,
    *,
    resume: bool = False,
) -> PipelineResult:
    """Run the full prep pipeline in ``workdir`` under ``cfg`` and return every result.

    Inputs are resolved (and guarded, inside each step) immediately before that step
    runs — validation stays interleaved with execution because each step consumes the
    previous step's outputs, so it can't be hoisted ahead of the run.
    """
    workdir = Path(workdir)

    pp_in = param_prot_step.resolve_inputs(workdir, cfg)
    pp = param_prot_step.param_prot(pp_in, resume=resume)

    mc_in = make_crystal_step.resolve_inputs(workdir, cfg)
    mc = make_crystal_step.make_crystal(mc_in, cfg.crystal, resume=resume)

    wb_in = make_waterbox_step.resolve_inputs(workdir, cfg)
    wb = make_waterbox_step.make_waterbox(
        wb_in,
        cfg.waterbox,
        cfg.profile("waterbox_min"),
        cfg.profile("waterbox_equil"),
        resume=resume,
    )

    sv_in = solvate_step.resolve_inputs(workdir, cfg)
    sv = solvate_step.solvate(sv_in, cfg.solvate, resume=resume)

    mn_in = minimize_step.resolve_inputs(workdir, cfg)
    mn = minimize_step.minimize(mn_in, cfg.profile("min"), resume=resume)

    eq_in = equilibrate_step.resolve_inputs(workdir, cfg)
    eq = equilibrate_step.equilibrate(eq_in, cfg.equilibrate, cfg.profile("equil"), resume=resume)

    rv_in = resolvate_step.resolve_inputs(workdir, cfg)
    rv = resolvate_step.resolvate(
        rv_in,
        cfg.resolvate,
        cfg.profile("resolv_min"),
        cfg.profile("resolv_equil"),
        resume=resume,
    )

    return PipelineResult(workdir=workdir, steps=[pp, mc, wb, sv, mn, eq, rv])
