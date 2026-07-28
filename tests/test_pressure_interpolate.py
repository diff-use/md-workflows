from pathlib import Path

import pytest

from md_workflows.core.exceptions import MDWorkflowError
from md_workflows.core.steps.pressure_interpolate import interpolate

REF = Path("ref.gro")


def test_normal_interpolation():
    # dNw = round((1-(-500))*(1200-1000)/(300-(-500))) = round(125.25) = 125
    r = interpolate(-500.0, 300.0, 1000, 1200, REF, target=1.0)
    assert r.dNw == 125
    assert r.Nw_target == 1125
    assert not r.converged
    assert r.warning is None


def test_extrapolation_same_side_warns_and_can_be_negative():
    r = interpolate(400.0, 900.0, 1000, 1200, REF, target=1.0)
    assert r.dNw < 0
    assert r.warning is not None


def test_zero_slope_near_target_converges():
    r = interpolate(50.0, 50.0, 1000, 1000, REF, target=1.0, pressure_tol=100.0)
    assert r.dNw == 0
    assert r.converged


def test_zero_slope_far_from_target_raises():
    with pytest.raises(MDWorkflowError):
        interpolate(-800.0, -800.0, 1000, 1050, REF, target=1.0, pressure_tol=100.0)


def test_absurd_dNw_far_from_target_raises():
    with pytest.raises(MDWorkflowError):
        interpolate(-5000.0, -4999.999, 1000, 1001, REF, target=1.0, max_dNw_factor=1.0)


def test_absurd_dNw_near_target_converges():
    r = interpolate(
        50.0, 50.001, 1000, 1001, REF, target=1.0, pressure_tol=100.0, max_dNw_factor=1.0
    )
    assert r.dNw == 0
    assert r.converged
