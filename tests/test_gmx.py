from pathlib import Path

import pytest

from md_workflows.core.gmx import checksums, count_waters, manage_wat_topology, sha256_file


def test_count_waters_pdb(tmp_path: Path):
    pdb = tmp_path / "x.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1       0.0   0.0   0.0\n"
        "HETATM    2  O   WAT A   2       1.0   0.0   0.0\n"
        "HETATM    3  H1  WAT A   2       1.1   0.0   0.0\n"
        "HETATM    4  H2  WAT A   2       1.2   0.0   0.0\n"
        "HETATM    5  O   WAT A   3       2.0   0.0   0.0\n"
        "HETATM    6  H1  WAT A   3       2.1   0.0   0.0\n"
        "HETATM    7  H2  WAT A   3       2.2   0.0   0.0\n"
    )
    assert count_waters(pdb) == 2


def test_count_waters_gro(tmp_path: Path):
    gro = tmp_path / "y.gro"
    gro.write_text(
        "title\n 4\n"
        "    1ALA      N    1   0.000   0.000   0.000\n"
        "    2SOL     OW    2   1.000   0.000   0.000\n"
        "    2SOL    HW1    3   1.100   0.000   0.000\n"
        "    2SOL    HW2    4   1.200   0.000   0.000\n"
        "   2.0 2.0 2.0\n"
    )
    assert count_waters(gro) == 1


def _top(tmp_path: Path) -> Path:
    top = tmp_path / "m.top"
    top.write_text("[ molecules ]\n; Compound  #mols\nsystem1 1\nWAT 100\nCl- 5\nNa+ 7\n")
    return top


def test_manage_wat_topology_append_is_idempotent(tmp_path: Path):
    top = _top(tmp_path)
    manage_wat_topology(top, 863, mode="append")
    manage_wat_topology(top, 863, mode="append")  # re-run must not double-add
    assert top.read_text().count("WAT 863") == 1
    assert top.read_text().rstrip().endswith("WAT 863")


def test_manage_wat_topology_multi_block_sequence(tmp_path: Path):
    top = _top(tmp_path)
    manage_wat_topology(top, 863, mode="append")  # run 1
    manage_wat_topology(top, 1294, mode="append")  # run 2 probe
    manage_wat_topology(top, 412, mode="drop_last_then_append")  # final drops probe
    wat_lines = [ln for ln in top.read_text().splitlines() if ln.startswith("WAT")]
    # initial 100 kept, run-1 863 kept, probe 1294 dropped, final 412 appended
    assert wat_lines == ["WAT 100", "WAT 863", "WAT 412"]


def test_manage_wat_topology_requires_molecules_section(tmp_path: Path):
    top = tmp_path / "bad.top"
    top.write_text("[ atoms ]\n")
    with pytest.raises(ValueError):
        manage_wat_topology(top, 10)


def test_checksums_skips_missing(tmp_path: Path):
    f = tmp_path / "a.txt"
    f.write_text("hello")
    result = checksums({"a": f, "missing": tmp_path / "nope"})
    assert set(result) == {"a"}
    assert result["a"] == sha256_file(f)
    assert len(result["a"]) == 64
