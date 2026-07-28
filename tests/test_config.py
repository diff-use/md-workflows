from pathlib import Path

import pytest

from md_workflows.core.config import RUN_PROFILE_KEYS, GromacsRunProfile, SystemConfig


def test_default_profiles_render_taylor_flags():
    cfg = SystemConfig()
    assert cfg.profile("min").to_mdrun_flags() == [
        "-ntmpi",
        "1",
        "-ntomp",
        "16",
        "-nb",
        "gpu",
        "-pme",
        "cpu",
        "-bonded",
        "cpu",
        "-notunepme",
    ]
    assert cfg.profile("equil").to_mdrun_flags() == [
        "-ntmpi",
        "1",
        "-ntomp",
        "16",
        "-nb",
        "gpu",
        "-pme",
        "gpu",
        "-bonded",
        "gpu",
        "-notunepme",
    ]


def test_to_mdrun_flags_omits_none_and_handles_tunepme():
    p = GromacsRunProfile(ntmpi=None, ntomp=8, nb="gpu", pme=None, bonded=None, tunepme=True)
    assert p.to_mdrun_flags() == ["-ntomp", "8", "-nb", "gpu"]


def test_overrides_win_and_merge_into_defaults():
    cfg = SystemConfig.load(
        overrides={"pdb_id": "4LZT", "crystal": {"ix": 5}, "run_profiles": {"min": {"ntomp": 26}}}
    )
    assert cfg.pdb_id == "4LZT"
    assert cfg.crystal.ix == 5
    # overriding one field of one profile keeps that profile's other defaults ...
    assert "-pme" in cfg.profile("min").to_mdrun_flags()
    assert cfg.profile("min").pme == "cpu"
    assert cfg.profile("min").ntomp == 26
    # ... and leaves the other profiles present
    assert set(cfg.run_profiles) == set(RUN_PROFILE_KEYS)
    assert cfg.profile("equil").pme == "gpu"


def test_yaml_round_trip(tmp_path: Path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("pdb_id: 6LYZ\nwaterbox:\n  nc_scale: 3\n")
    cfg = SystemConfig.load(cfg_file)
    assert cfg.pdb_id == "6LYZ"
    assert cfg.waterbox.nc_scale == 3


def test_flag_overrides_beat_file(tmp_path: Path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("pdb_id: 6LYZ\n")
    cfg = SystemConfig.load(cfg_file, overrides={"pdb_id": "4LZT"})
    assert cfg.pdb_id == "4LZT"


def test_bad_yaml_raises_valueerror(tmp_path: Path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("a: b: c: [\n")
    with pytest.raises(ValueError):
        SystemConfig.load(bad)


def test_unknown_format_raises(tmp_path: Path):
    weird = tmp_path / "c.ini"
    weird.write_text("x=1\n")
    with pytest.raises(ValueError):
        SystemConfig.load(weird)
