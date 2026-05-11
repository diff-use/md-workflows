"""Unit tests for md_workflows.workflows.mdmx module."""

from unittest.mock import Mock, patch
from md_workflows.workflows import mdmx


class TestMain:
    """Test the main() workflow orchestration function."""

    @patch("md_workflows.workflows.mdmx.run_resolvate")
    @patch("md_workflows.workflows.mdmx.run_equilibrate")
    @patch("md_workflows.workflows.mdmx.run_minimize")
    @patch("md_workflows.workflows.mdmx.run_solvate")
    @patch("md_workflows.workflows.mdmx.run_make_waterbox")
    @patch("md_workflows.workflows.mdmx.run_make_crystal")
    @patch("md_workflows.workflows.mdmx.run_param_prot")
    def test_main_with_defaults(
        self,
        mock_param_prot: Mock,
        mock_make_crystal: Mock,
        mock_make_waterbox: Mock,
        mock_solvate: Mock,
        mock_minimize: Mock,
        mock_equilibrate: Mock,
        mock_resolvate: Mock,
    ):
        """Test main() executes all workflow steps in correct order with default parameters."""
        mdmx.main(
            min_mdp="tests/test_artifacts/min.mdp",
            min_water_mdp="tests/test_artifacts/min_water.mdp",
            equil_mdp="tests/test_artifacts/equil.mdp",
            equil_water_mdp="tests/test_artifacts/equil_water.mdp",
        )

        # Verify all steps were called exactly once
        mock_param_prot.assert_called_once_with(pdb_id="6B8X")
        mock_make_crystal.assert_called_once_with(
            ix=1, iy=None, iz=None, chimerax_exec="/usr/bin/chimerax-daily"
        )
        mock_make_waterbox.assert_called_once_with(
            ntomp=26,
            min_water_mdp="tests/test_artifacts/min_water.mdp",
            equil_water_mdp="tests/test_artifacts/equil_water.mdp",
        )
        mock_solvate.assert_called_once_with()
        mock_minimize.assert_called_once_with(ntomp=26, min_mdp="tests/test_artifacts/min.mdp")
        mock_equilibrate.assert_called_once_with(ntomp=26, equil_mdp="tests/test_artifacts/equil.mdp")
        mock_resolvate.assert_called_once_with(ntmpi=8, ntomp=1)

        # Verify call order matches the documented pipeline order
        manager = Mock()
        manager.attach_mock(mock_param_prot, "param_prot")
        manager.attach_mock(mock_make_crystal, "make_crystal")
        manager.attach_mock(mock_make_waterbox, "make_waterbox")
        manager.attach_mock(mock_solvate, "solvate")
        manager.attach_mock(mock_minimize, "minimize")
        manager.attach_mock(mock_equilibrate, "equilibrate")
        manager.attach_mock(mock_resolvate, "resolvate")

        expected_calls = [
            ("param_prot", (), {"pdb_id": "6B8X"}),
            ("make_crystal", (), {"ix": 1, "iy": None, "iz": None, "chimerax_exec": "/usr/bin/chimerax-daily"}),
            ("make_waterbox", (), {"ntomp": 26, "min_water_mdp": "tests/test_artifacts/min_water.mdp", "equil_water_mdp": "tests/test_artifacts/equil_water.mdp"}),
            ("solvate", (), {}),
            ("minimize", (), {"ntomp": 26, "min_mdp": "tests/test_artifacts/min.mdp"}),
            ("equilibrate", (), {"ntomp": 26, "equil_mdp": "tests/test_artifacts/equil.mdp"}),
            ("resolvate", (), {"ntmpi": 8, "ntomp": 1}),
        ]

        actual_calls = [(name, call.args, call.kwargs) for name, call in manager.method_calls]
        assert actual_calls == expected_calls

    @patch("md_workflows.workflows.mdmx.run_resolvate")
    @patch("md_workflows.workflows.mdmx.run_equilibrate")
    @patch("md_workflows.workflows.mdmx.run_minimize")
    @patch("md_workflows.workflows.mdmx.run_solvate")
    @patch("md_workflows.workflows.mdmx.run_make_waterbox")
    @patch("md_workflows.workflows.mdmx.run_make_crystal")
    @patch("md_workflows.workflows.mdmx.run_param_prot")
    def test_main_with_custom_parameters(
        self,
        mock_param_prot: Mock,
        mock_make_crystal: Mock,
        mock_make_waterbox: Mock,
        mock_solvate: Mock,
        mock_minimize: Mock,
        mock_equilibrate: Mock,
        mock_resolvate: Mock,
    ):
        """Test main() correctly passes custom parameters to workflow steps."""
        mdmx.main(
            ntomp=32,
            param_pdb_id="1ABC",
            crystal_ix=2,
            crystal_iy=3,
            crystal_iz=4,
            chimerax_exec="/custom/path/chimerax",
            resolv_ntmpi=16,
            resolv_ntomp=2,
            min_mdp="/path/to/min.mdp",
            min_water_mdp="/path/to/min_water.mdp",
            equil_mdp="/path/to/equil.mdp",
            equil_water_mdp="/path/to/equil_water.mdp",
        )

        mock_param_prot.assert_called_once_with(pdb_id="1ABC")
        mock_make_crystal.assert_called_once_with(
            ix=2, iy=3, iz=4, chimerax_exec="/custom/path/chimerax"
        )
        mock_make_waterbox.assert_called_once_with(
            ntomp=32, min_water_mdp="/path/to/min_water.mdp", equil_water_mdp="/path/to/equil_water.mdp"
        )
        mock_solvate.assert_called_once_with()
        mock_minimize.assert_called_once_with(ntomp=32, min_mdp="/path/to/min.mdp")
        mock_equilibrate.assert_called_once_with(ntomp=32, equil_mdp="/path/to/equil.mdp")
        mock_resolvate.assert_called_once_with(ntmpi=16, ntomp=2)

    @patch("md_workflows.workflows.mdmx.run_resolvate")
    @patch("md_workflows.workflows.mdmx.run_equilibrate")
    @patch("md_workflows.workflows.mdmx.run_minimize")
    @patch("md_workflows.workflows.mdmx.run_solvate")
    @patch("md_workflows.workflows.mdmx.run_make_waterbox")
    @patch("md_workflows.workflows.mdmx.run_make_crystal")
    @patch("md_workflows.workflows.mdmx.run_param_prot")
    def test_main_with_partial_crystal_dimensions(
        self,
        mock_param_prot: Mock,
        mock_make_crystal: Mock,
        mock_make_waterbox: Mock,
        mock_solvate: Mock,
        mock_minimize: Mock,
        mock_equilibrate: Mock,
        mock_resolvate: Mock,
    ):
        """Test main() handles partial crystal dimension specifications."""
        # Only specify ix and iy, iz remains None
        mdmx.main(
            crystal_ix=3,
            crystal_iy=5,
            min_mdp="tests/test_artifacts/min.mdp",
            min_water_mdp="tests/test_artifacts/min_water.mdp",
            equil_mdp="tests/test_artifacts/equil.mdp",
            equil_water_mdp="tests/test_artifacts/equil_water.mdp",
        )

        mock_make_crystal.assert_called_once_with(
            ix=3, iy=5, iz=None, chimerax_exec="/usr/bin/chimerax-daily"
        )


class TestMainWithMDPFiles:
    """Test MDP file parameter handling."""

    @patch("md_workflows.workflows.mdmx.run_resolvate")
    @patch("md_workflows.workflows.mdmx.run_equilibrate")
    @patch("md_workflows.workflows.mdmx.run_minimize")
    @patch("md_workflows.workflows.mdmx.run_solvate")
    @patch("md_workflows.workflows.mdmx.run_make_waterbox")
    @patch("md_workflows.workflows.mdmx.run_make_crystal")
    @patch("md_workflows.workflows.mdmx.run_param_prot")
    def test_main_with_all_mdp_files(
        self,
        mock_param_prot: Mock,
        mock_make_crystal: Mock,
        mock_make_waterbox: Mock,
        mock_solvate: Mock,
        mock_minimize: Mock,
        mock_equilibrate: Mock,
        mock_resolvate: Mock,
    ):
        """Test main() correctly passes all MDP file paths to workflow steps."""
        mdmx.main(
            min_mdp="custom/min.mdp",
            min_water_mdp="custom/min_water.mdp",
            equil_mdp="custom/equil.mdp",
            equil_water_mdp="custom/equil_water.mdp",
        )

        mock_minimize.assert_called_once_with(ntomp=26, min_mdp="custom/min.mdp")
        mock_equilibrate.assert_called_once_with(ntomp=26, equil_mdp="custom/equil.mdp")
        mock_make_waterbox.assert_called_once_with(
            ntomp=26,
            min_water_mdp="custom/min_water.mdp",
            equil_water_mdp="custom/equil_water.mdp",
        )

    @patch("md_workflows.workflows.mdmx.run_resolvate")
    @patch("md_workflows.workflows.mdmx.run_equilibrate")
    @patch("md_workflows.workflows.mdmx.run_minimize")
    @patch("md_workflows.workflows.mdmx.run_solvate")
    @patch("md_workflows.workflows.mdmx.run_make_waterbox")
    @patch("md_workflows.workflows.mdmx.run_make_crystal")
    @patch("md_workflows.workflows.mdmx.run_param_prot")
    def test_main_with_partial_mdp_files(
        self,
        mock_make_waterbox: Mock,
        mock_minimize: Mock,
        mock_equilibrate: Mock,
    ):
        """Test main() handles partial MDP file specifications."""
        mdmx.main(
            min_mdp="custom/min.mdp",
            equil_water_mdp="custom/equil_water.mdp",
        )

        mock_minimize.assert_called_once_with(ntomp=26, min_mdp="tests/test_artifacts/min.mdp")
        mock_equilibrate.assert_called_once_with(ntomp=26, equil_mdp="tests/test_artifacts/equil.mdp")
        mock_make_waterbox.assert_called_once_with(
            ntomp=26,
            min_water_mdp="tests/test_artifacts/min_water.mdp",
            equil_water_mdp="tests/test_artifacts/equil_water.mdp",
        )


class TestCLI:
    """Test the CLI argument parsing."""

    @patch("md_workflows.workflows.mdmx.main")
    @patch("sys.argv", ["mdmx"])
    def test_cli_with_no_arguments(self, mock_main: Mock):
        """Test CLI uses default values when no arguments provided."""
        mdmx._cli()

        mock_main.assert_called_once_with(
            ntomp=26,
            param_pdb_id="6B8X",
            crystal_ix=1,
            crystal_iy=None,
            crystal_iz=None,
            chimerax_exec="/usr/bin/chimerax-daily",
            resolv_ntmpi=8,
            resolv_ntomp=1,
            min_mdp="artifacts/min.mdp",
            min_water_mdp="artifacts/min_water.mdp",
            equil_mdp="artifacts/equil.mdp",
            equil_water_mdp="artifacts/equil_water.mdp",
        )

    @patch("md_workflows.workflows.mdmx.main")
    @patch(
        "sys.argv",
        [
            "mdmx",
            "--ntomp",
            "32",
            "--param-pdb-id",
            "7XYZ",
            "--ix",
            "2",
            "--iy",
            "3",
            "--iz",
            "4",
            "--chimerax-exec",
            "/opt/chimerax",
            "--resolv-ntmpi",
            "12",
            "--resolv-ntomp",
            "4",
            "--min",
            "/custom/min.mdp",
            "--min-water",
            "/custom/min_water.mdp",
            "--equil",
            "/custom/equil.mdp",
            "--equil-water",
            "/custom/equil_water.mdp",
        ],
    )
    def test_cli_with_all_arguments(self, mock_main: Mock):
        """Test CLI correctly parses all command-line arguments."""
        mdmx._cli()

        mock_main.assert_called_once_with(
            ntomp=32,
            param_pdb_id="7XYZ",
            crystal_ix=2,
            crystal_iy=3,
            crystal_iz=4,
            chimerax_exec="/opt/chimerax",
            resolv_ntmpi=12,
            resolv_ntomp=4,
            min_mdp="/custom/min.mdp",
            min_water_mdp="/custom/min_water.mdp",
            equil_mdp="/custom/equil.mdp",
            equil_water_mdp="/custom/equil_water.mdp",
        )

    @patch("md_workflows.workflows.mdmx.main")
    @patch(
        "sys.argv",
        ["mdmx", "--ntomp", "48", "--param-pdb-id", "3DEF"],
    )
    def test_cli_with_partial_arguments(self, mock_main: Mock):
        """Test CLI correctly handles mix of provided and default arguments."""
        mdmx._cli()

        mock_main.assert_called_once_with(
            ntomp=48,
            param_pdb_id="3DEF",
            crystal_ix=1,  # default
            crystal_iy=None,  # default
            crystal_iz=None,  # default
            chimerax_exec="/usr/bin/chimerax-daily",  # default
            resolv_ntmpi=8,  # default
            resolv_ntomp=1,  # default
            min_mdp="artifacts/min.mdp",  # default
            min_water_mdp="artifacts/min_water.mdp",  # default
            equil_mdp="artifacts/equil.mdp",  # default
            equil_water_mdp="artifacts/equil_water.mdp",  # default
        )

    @patch("md_workflows.workflows.mdmx.main")
    @patch(
        "sys.argv",
        [
            "mdmx",
            "--min",
            "artifacts/min.mdp",
            "--equil-water",
            "artifacts/equil_water.mdp",
        ],
    )
    def test_cli_with_only_mdp_arguments(self, mock_main: Mock):
        """Test CLI correctly handles only MDP file arguments."""
        mdmx._cli()

        mock_main.assert_called_once_with(
            ntomp=26,  # default
            param_pdb_id="6B8X",  # default
            crystal_ix=1,  # default
            crystal_iy=None,  # default
            crystal_iz=None,  # default
            chimerax_exec="/usr/bin/chimerax-daily",  # default
            resolv_ntmpi=8,  # default
            resolv_ntomp=1,  # default
            min_mdp="artifacts/min.mdp",
            min_water_mdp="artifacts/min_water.mdp",  # default
            equil_mdp="artifacts/equil.mdp",  # default
            equil_water_mdp="artifacts/equil_water.mdp",
        )
