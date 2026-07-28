"""Typed inputs and results shared by every core step.

These are Pydantic models so results are serializable (JSON round-trip) and carry
enough metadata — parameters used, input checksums, tool versions — to later prove a
CLI run and an orchestrated run are scientifically equivalent. Per-step ``*Inputs`` and
``*Result`` subclasses live next to each step in ``core/steps``.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from .config import GromacsRunProfile
from .exceptions import MissingInputError


class StepStatus(str, Enum):
    COMPLETED = "completed"
    SKIPPED = "skipped"  # outputs already present and inputs unchanged (--resume)
    FAILED = "failed"


class StepInputs(BaseModel):
    """Explicit, typed inputs for one step. Subclasses add the concrete file fields.

    Core steps take one of these instead of assuming filenames in ``cwd``.
    """

    model_config = ConfigDict(extra="forbid")

    workdir: Path

    def consumed_paths(self) -> list[Path]:
        """Files that must exist before the step runs. Overridden per step."""
        return []

    def check_exists(self, step: str) -> None:
        """Raise :class:`MissingInputError` if any required input is absent."""
        missing = [p for p in self.consumed_paths() if not Path(p).exists()]
        if missing:
            raise MissingInputError(step, missing)


class StepResult(BaseModel):
    """Typed result of a step: what it produced plus reproducibility metadata."""

    model_config = ConfigDict(extra="forbid")

    step: str
    status: StepStatus = StepStatus.COMPLETED
    workdir: Path
    outputs: dict[str, Path] = Field(default_factory=dict)
    params: dict[str, Any] = Field(default_factory=dict)
    run_profiles_used: dict[str, GromacsRunProfile] = Field(default_factory=dict)
    input_checksums: dict[str, str] = Field(default_factory=dict)
    tool_versions: dict[str, str] = Field(default_factory=dict)
    metrics: dict[str, float | int] = Field(default_factory=dict)
    log_paths: list[Path] = Field(default_factory=list)

    def output(self, name: str) -> Path:
        """Return a named output path, or raise KeyError with a helpful message."""
        try:
            return self.outputs[name]
        except KeyError as exc:
            raise KeyError(
                f"{self.step}: no output named {name!r}; have {sorted(self.outputs)}"
            ) from exc
