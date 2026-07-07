from pathlib import Path

from typer.testing import CliRunner

from md_workflows.cli import app

runner = CliRunner()


def test_help_lists_all_commands():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    for cmd in [
        "param-prot",
        "make-crystal",
        "make-waterbox",
        "solvate",
        "minimize",
        "equilibrate",
        "resolvate",
        "run-params-gaussian",
        "run-pipeline",
    ]:
        assert cmd in result.output


def test_subcommand_help_exit_zero():
    assert runner.invoke(app, ["run-pipeline", "--help"]).exit_code == 0


def test_missing_input_exits_2(tmp_path: Path):
    result = runner.invoke(app, ["--workdir", str(tmp_path), "minimize"])
    assert result.exit_code == 2
    assert "missing required input" in result.output


def test_bad_config_exits_1(tmp_path: Path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("a: b: c: [\n")
    result = runner.invoke(app, ["--workdir", str(tmp_path), "--config", str(bad), "minimize"])
    assert result.exit_code == 1
    assert "config error" in result.output


def test_valid_config_parses_then_reports_missing_inputs(tmp_path: Path):
    cfg = tmp_path / "ok.yaml"
    cfg.write_text("pdb_id: 4LZT\ncrystal:\n  ix: 5\n")
    result = runner.invoke(app, ["--workdir", str(tmp_path), "--config", str(cfg), "minimize"])
    # config parsed fine; fails only on missing inputs (exit 2)
    assert result.exit_code == 2
