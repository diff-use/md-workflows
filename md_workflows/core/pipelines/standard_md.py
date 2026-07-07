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
    """Run the full prep pipeline in ``workdir`` under ``cfg`` and return every result."""
    workdir = Path(workdir)

    pp = param_prot_step.param_prot(param_prot_step.resolve_inputs(workdir, cfg), resume=resume)
    mc = make_crystal_step.make_crystal(
        make_crystal_step.resolve_inputs(workdir, cfg), cfg.crystal, resume=resume
    )
    wb = make_waterbox_step.make_waterbox(
        make_waterbox_step.resolve_inputs(workdir, cfg),
        cfg.waterbox,
        cfg.profile("waterbox_min"),
        cfg.profile("waterbox_equil"),
        resume=resume,
    )
    sv = solvate_step.solvate(solvate_step.resolve_inputs(workdir, cfg), cfg.solvate, resume=resume)
    mn = minimize_step.minimize(
        minimize_step.resolve_inputs(workdir, cfg), cfg.profile("min"), resume=resume
    )
    eq = equilibrate_step.equilibrate(
        equilibrate_step.resolve_inputs(workdir, cfg),
        cfg.equilibrate,
        cfg.profile("equil"),
        resume=resume,
    )
    rv = resolvate_step.resolvate(
        resolvate_step.resolve_inputs(workdir, cfg),
        cfg.resolvate,
        cfg.profile("resolv_min"),
        cfg.profile("resolv_equil"),
        resume=resume,
    )

    return PipelineResult(workdir=workdir, steps=[pp, mc, wb, sv, mn, eq, rv])
