"""Two-state linear interpolation for the resolvation water correction.

Port of the taylor ``pressure_interpolate.py`` (Wych & Wall MDPreparation) as a pure
function. Given two ``(water count, mean pressure)`` states it solves for the water
delta that brings the system to the target pressure:

    dNw = round((target - P1) * (Nw2 - Nw1) / (P2 - P1))

The script's ``argparse``/``print`` command-suggestion UX is dropped — that belongs to
the CLI/flow layer — and errors are raised rather than ``sys.exit``.

The numeric solve (:func:`interpolate`) is separated from the gmx I/O
(:func:`pressure_interpolate`) so it is unit-testable without GROMACS.

Two robustness guards beyond the original script (both use the same rule: if the system
is already within ``pressure_tol`` of target, report convergence with ``dNw=0``;
otherwise raise, because the linear model cannot be trusted to reach target):

* **Zero / near-zero slope** (``P2 ≈ P1``): adding water did not move the pressure, so
  the slope is undefined — the original ``P1 == P2`` guard, generalized to an epsilon.
* **Absurd ``dNw``**: a tiny (but nonzero) slope can produce an impossible water count;
  reject any ``|dNw|`` beyond ``max_dNw_factor * Nw1``.
"""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel

from ..exceptions import MDWorkflowError
from ..gmx import count_waters, read_mean_pressure

_SLOPE_EPS = 1e-9  # bar; below this |P2-P1| is treated as zero slope


class InterpolationResult(BaseModel):
    dNw: int
    Nw_target: int
    ref_gro: Path
    P1: float
    P2: float
    Nw1: int
    Nw2: int
    converged: bool = False  # True when already within pressure_tol of target
    warning: str | None = None


def interpolate(
    p1: float,
    p2: float,
    nw1: int,
    nw2: int,
    ref_gro: Path,
    *,
    target: float = 1.0,
    pressure_tol: float = 100.0,
    max_dNw_factor: float = 1.0,
) -> InterpolationResult:
    """Pure numeric solve for the water delta (no gmx I/O)."""
    ref = Path(ref_gro)
    near_target = abs(p1 - target) <= pressure_tol

    def _converged(reason: str) -> InterpolationResult:
        return _result(0, nw1, ref, p1, p2, nw1, nw2, converged=True, warning=reason)

    # Guard 1: zero / near-zero slope — linear model is undefined.
    if abs(p2 - p1) < _SLOPE_EPS:
        if near_target:
            return _converged(
                f"pressure insensitive to water (P1={p1:.1f}≈P2={p2:.1f} bar) but within "
                f"{pressure_tol:.0f} bar of target — keeping current solvation."
            )
        raise MDWorkflowError(
            f"pressure_interpolate: zero slope (P1={p1:.1f}≈P2={p2:.1f} bar) and pressure "
            f"is {abs(p1 - target):.0f} bar from target ({target:.1f}); cannot interpolate."
        )

    dNw = int(round((target - p1) * (nw2 - nw1) / (p2 - p1)))

    # Guard 2: sanity-bound the correction against tiny-slope blow-ups.
    if abs(dNw) > max_dNw_factor * nw1:
        if near_target:
            return _converged(
                f"computed dNw={dNw} exceeds sanity bound ({max_dNw_factor}*{nw1}) but "
                f"pressure is within {pressure_tol:.0f} bar of target — keeping current."
            )
        raise MDWorkflowError(
            f"pressure_interpolate: computed dNw={dNw} exceeds sanity bound "
            f"({max_dNw_factor}*{nw1} waters); the two-state model is unreliable here."
        )

    warning = None
    if (p1 - target) * (p2 - target) > 0:
        warning = (
            f"extrapolation: P1={p1:.1f} and P2={p2:.1f} bar are on the same side of "
            f"target={target:.1f} bar; result may be inaccurate."
        )
    return _result(dNw, nw1 + dNw, ref, p1, p2, nw1, nw2, converged=False, warning=warning)


def pressure_interpolate(
    gro1: Path,
    edr1: Path,
    gro2: Path,
    edr2: Path,
    workdir: Path,
    *,
    target: float = 1.0,
    pressure_tol: float = 100.0,
    max_dNw_factor: float = 1.0,
    ref_gro: Path | None = None,
    gmx_bin: str = "gmx",
) -> InterpolationResult:
    """Read the two states from gmx and interpolate the water delta to ``target``."""
    p1 = read_mean_pressure(edr1, workdir, gmx_bin=gmx_bin)
    p2 = read_mean_pressure(edr2, workdir, gmx_bin=gmx_bin)
    nw1 = count_waters(gro1)
    nw2 = count_waters(gro2)
    return interpolate(
        p1,
        p2,
        nw1,
        nw2,
        Path(ref_gro) if ref_gro is not None else Path(gro1),
        target=target,
        pressure_tol=pressure_tol,
        max_dNw_factor=max_dNw_factor,
    )


def _result(
    dNw: int,
    nw_target: int,
    ref: Path,
    p1: float,
    p2: float,
    nw1: int,
    nw2: int,
    *,
    converged: bool,
    warning: str | None,
) -> InterpolationResult:
    return InterpolationResult(
        dNw=dNw,
        Nw_target=nw_target,
        ref_gro=ref,
        P1=p1,
        P2=p2,
        Nw1=nw1,
        Nw2=nw2,
        converged=converged,
        warning=warning,
    )
