from pathlib import Path

import pytest

from md_workflows.core.config import SystemConfig
from md_workflows.core.exceptions import MissingInputError
from md_workflows.core.results import StepStatus
from md_workflows.core.steps import (
    equilibrate,
    make_crystal,
    make_waterbox,
    minimize,
    param_prot,
    resolvate,
    solvate,
)

# Steps with required pre-existing inputs (param_prot / gaussian download or detect theirs).
STEPS_WITH_REQUIRED_INPUTS = [
    make_crystal,
    make_waterbox,
    solvate,
    minimize,
    equilibrate,
    resolvate,
]


@pytest.mark.parametrize("mod", STEPS_WITH_REQUIRED_INPUTS)
def test_resolve_inputs_are_under_workdir(mod, tmp_path: Path):
    cfg = SystemConfig()
    inputs = mod.resolve_inputs(tmp_path, cfg)
    assert inputs.workdir == tmp_path
    consumed = inputs.consumed_paths()
    assert consumed, f"{mod.STEP} should declare required inputs"
    for p in consumed:
        assert tmp_path in p.parents, f"{p} not under workdir"


@pytest.mark.parametrize("mod", STEPS_WITH_REQUIRED_INPUTS)
def test_check_inputs_raises_listing_all_missing(mod, tmp_path: Path):
    cfg = SystemConfig()
    inputs = mod.resolve_inputs(tmp_path, cfg)
    with pytest.raises(MissingInputError) as exc:
        mod.check_inputs(inputs)
    assert exc.value.step == mod.STEP
    assert len(exc.value.missing) == len(inputs.consumed_paths())


def test_param_prot_has_no_required_files(tmp_path: Path):
    cfg = SystemConfig()
    inputs = param_prot.resolve_inputs(tmp_path, cfg)
    assert inputs.consumed_paths() == []
    # guard does not raise when nothing is required
    assert param_prot.check_inputs(inputs) is None


def test_minimize_resume_skips_without_running_gmx(tmp_path: Path):
    (tmp_path / "artifacts").mkdir()
    for name in ("md_model.pdb", "md_model.top"):
        (tmp_path / name).write_text("x")
    (tmp_path / "artifacts" / "min.mdp").write_text("integrator=steep\n")
    (tmp_path / "md_min.gro").write_text("fake existing output")

    cfg = SystemConfig()
    inputs = minimize.resolve_inputs(tmp_path, cfg)
    result = minimize.minimize(inputs, cfg.profile("min"), resume=True)

    assert result.status == StepStatus.SKIPPED
    assert result.outputs["gro"] == tmp_path / "md_min.gro"
    assert set(result.input_checksums) == {"model_pdb", "model_top", "min_mdp"}


def test_mdp_dir_resolves_under_workdir(tmp_path: Path):
    cfg = SystemConfig()  # mdp_dir default "artifacts" (relative)
    inputs = minimize.resolve_inputs(tmp_path, cfg)
    assert inputs.min_mdp == tmp_path / "artifacts" / "min.mdp"


def test_absolute_mdp_dir_is_honored(tmp_path: Path):
    cfg = SystemConfig.load(overrides={"mdp_dir": "/opt/mdps"})
    inputs = minimize.resolve_inputs(tmp_path, cfg)
    assert inputs.min_mdp == Path("/opt/mdps/min.mdp")
