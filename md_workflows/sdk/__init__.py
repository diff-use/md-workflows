"""Public Python API for md-workflows.

Thin, curated re-exports of the pure ``core`` functions plus a couple of ergonomic
wrappers. Import from here for a stable surface; ``core`` internals may move.

    from md_workflows.sdk import run_standard_md, standard_md_pipeline, SystemConfig
"""

from __future__ import annotations

from ..core.config import (
    CrystalParams,
    EquilibrateParams,
    GaussianParams,
    GromacsRunProfile,
    MinimizeParams,
    ResolvateParams,
    SolvateParams,
    SystemConfig,
    WaterboxParams,
)
from ..core.exceptions import (
    GmxOutputParseError,
    MDWorkflowError,
    MissingInputError,
    StepToolError,
)
from ..core.pipelines.standard_md import PipelineResult, standard_md_pipeline
from ..core.results import StepResult, StepStatus
from ..core.steps.equilibrate import equilibrate
from ..core.steps.make_crystal import make_crystal
from ..core.steps.make_waterbox import make_waterbox
from ..core.steps.minimize import minimize
from ..core.steps.param_prot import param_prot
from ..core.steps.pressure_interpolate import pressure_interpolate
from ..core.steps.resolvate import resolvate, resolvation_cycle
from ..core.steps.run_params_gaussian import run_params_gaussian
from ..core.steps.solvate import solvate
from .convenience import run_standard_md

__all__ = [
    # pipeline + convenience
    "standard_md_pipeline",
    "PipelineResult",
    "run_standard_md",
    # steps
    "param_prot",
    "make_crystal",
    "make_waterbox",
    "solvate",
    "minimize",
    "equilibrate",
    "resolvate",
    "resolvation_cycle",
    "pressure_interpolate",
    "run_params_gaussian",
    # config + results
    "SystemConfig",
    "GromacsRunProfile",
    "CrystalParams",
    "WaterboxParams",
    "SolvateParams",
    "MinimizeParams",
    "EquilibrateParams",
    "ResolvateParams",
    "GaussianParams",
    "StepResult",
    "StepStatus",
    # exceptions
    "MDWorkflowError",
    "MissingInputError",
    "StepToolError",
    "GmxOutputParseError",
]
