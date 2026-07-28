"""Ergonomic wrappers over the core pipeline for SDK users."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..core.config import SystemConfig
from ..core.pipelines.standard_md import PipelineResult, standard_md_pipeline


def run_standard_md(
    pdb_id: str,
    workdir: str | Path,
    *,
    config: str | Path | None = None,
    resume: bool = False,
    **overrides: Any,
) -> PipelineResult:
    """Run the standard MD prep pipeline for ``pdb_id`` in ``workdir``.

    ``config`` is an optional YAML/TOML file; keyword ``overrides`` (deep-merged, so
    nested groups like ``crystal={"ix": 5}`` work) win over it. Example::

        run_standard_md("4LZT", "runs/4lzt", crystal={"ix": 5}, resume=True)
    """
    cfg = SystemConfig.load(config, overrides={"pdb_id": pdb_id, **overrides})
    return standard_md_pipeline(Path(workdir), cfg, resume=resume)
