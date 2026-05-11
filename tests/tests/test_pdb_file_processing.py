"""Integration tests for md_workflows.pdb_file_processing module."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from md_workflows import pdb_file_processing


# Sample PDB file content for testing
SAMPLE_PDB_CONTENT = """HEADER    OXIDOREDUCTASE                          07-OCT-17   6B8X
TITLE     CRYSTAL STRUCTURE OF MDM2 IN COMPLEX WITH SMALL MOLECULE
CRYST1   44.970   56.220   69.870  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 50.00           N
ATOM      2  CA  MET A   1      11.000  21.000  31.000  1.00 50.00           C
HETATM  100  C1  GOL A 201      15.000  25.000  35.000  1.00 30.00           C
HETATM  101  C2  GOL A 201      16.000  26.000  36.000  1.00 30.00           C
HETATM  102  C3  GOL A 201      17.000  27.000  37.000  1.00 30.00           C
HETATM  200  O   HOH A 301      20.000  30.000  40.000  1.00 20.00           O
HETATM  201  O   WAT B 401      21.000  31.000  41.000  1.00 20.00           O
HETATM  300  C1  LIG B 501      25.000  35.000  45.000  1.00 40.00           C
HETATM  301  C2  LIG B 501      26.000  36.000  46.000  1.00 40.00           C
END
"""

SAMPLE_PDB_WITH_WATER_ONLY = """HEADER    TEST STRUCTURE
ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 50.00           N
HETATM  100  O   HOH A 301      20.000  30.000  40.000  1.00 20.00           O
HETATM  101  O   WAT A 302      21.000  31.000  41.000  1.00 20.00           O
END
"""

SAMPLE_PDB_MULTIPLE_LIGANDS = """HEADER    TEST STRUCTURE
ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 50.00           N
HETATM  100  C1  GOL A 201      15.000  25.000  35.000  1.00 30.00           C
HETATM  101  C2  GOL A 201      16.000  26.000  36.000  1.00 30.00           C
HETATM  200  C1  NAG B 301      20.000  30.000  40.000  1.00 35.00           C
HETATM  201  C2  NAG B 301      21.000  31.000  41.000  1.00 35.00           C
END
"""

SAMPLE_PDB_SINGLE_LIGAND = """HEADER    TEST STRUCTURE
ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 50.00           N
HETATM  100  C1  GOL A 201      15.000  25.000  35.000  1.00 30.00           C
HETATM  101  C2  GOL A 201      16.000  26.000  36.000  1.00 30.00           C
HETATM  102  C3  GOL A 201      17.000  27.000  37.000  1.00 30.00           C
HETATM  200  O   HOH A 301      20.000  30.000  40.000  1.00 20.00           O
END
"""


class TestFindLigandsInLegacyPdbText(unittest.TestCase):
    """Test the find_ligands_in_legacy_pdb_text function."""

    def test_finds_hetero_residues_excluding_water(self):
        """Test that it finds hetero residues and excludes water by default."""
        ligands = pdb_file_processing.find_ligands_in_legacy_pdb_text(SAMPLE_PDB_CONTENT)

        self.assertEqual(len(ligands), 2)
        
        # First ligand: GOL chain A
        self.assertEqual(ligands[0]["resname"], "GOL")
        self.assertEqual(ligands[0]["chain_id"], "A")
        self.assertEqual(ligands[0]["residue_seq"], "201")
        self.assertEqual(ligands[0]["n_atoms"], 3)
        
        # Second ligand: LIG chain B
        self.assertEqual(ligands[1]["resname"], "LIG")
        self.assertEqual(ligands[1]["chain_id"], "B")
        self.assertEqual(ligands[1]["residue_seq"], "501")
        self.assertEqual(ligands[1]["n_atoms"], 2)

    def test_includes_water_when_requested(self):
        """Test that water residues are included when exclude_water_solvent=False."""
        ligands = pdb_file_processing.find_ligands_in_legacy_pdb_text(
            SAMPLE_PDB_CONTENT, exclude_water_solvent=False
        )

        self.assertEqual(len(ligands), 4)  # GOL + 2 waters + LIG
        
        # Should include HOH and WAT
        resnames = {lig["resname"] for lig in ligands}
        self.assertIn("HOH", resnames)
        self.assertIn("WAT", resnames)
        self.assertIn("GOL", resnames)
        self.assertIn("LIG", resnames)

    def test_handles_empty_pdb(self):
        """Test that empty PDB returns empty list."""
        ligands = pdb_file_processing.find_ligands_in_legacy_pdb_text("")
        self.assertEqual(ligands, [])

    def test_handles_pdb_with_only_atoms(self):
        """Test PDB with only ATOM records and no HETATM."""
        pdb_content = """HEADER    TEST
ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 50.00           N
END
"""
        ligands = pdb_file_processing.find_ligands_in_legacy_pdb_text(pdb_content)
        self.assertEqual(ligands, [])

    def test_ligand_sorting_by_chain_and_residue(self):
        """Test that ligands are sorted by chain ID and residue sequence."""
        pdb_content = """HEADER    TEST
HETATM  100  C1  LIG C 300      15.000  25.000  35.000  1.00 30.00           C
HETATM  200  C1  GOL A 100      15.000  25.000  35.000  1.00 30.00           C
HETATM  300  C1  NAG B 200      15.000  25.000  35.000  1.00 30.00           C
END
"""
        ligands = pdb_file_processing.find_ligands_in_legacy_pdb_text(pdb_content)
        
        # Should be sorted: A < B < C
        self.assertEqual(ligands[0]["chain_id"], "A")
        self.assertEqual(ligands[1]["chain_id"], "B")
        self.assertEqual(ligands[2]["chain_id"], "C")

    def test_counts_atoms_per_residue(self):
        """Test that atom counts are correct for each residue."""
        pdb_content = """HEADER    TEST
HETATM  100  C1  GOL A 201      15.000  25.000  35.000  1.00 30.00           C
HETATM  101  C2  GOL A 201      16.000  26.000  36.000  1.00 30.00           C
HETATM  102  C3  GOL A 201      17.000  27.000  37.000  1.00 30.00           C
HETATM  103  O1  GOL A 201      18.000  28.000  38.000  1.00 30.00           O
HETATM  104  O2  GOL A 201      19.000  29.000  39.000  1.00 30.00           O
END
"""
        ligands = pdb_file_processing.find_ligands_in_legacy_pdb_text(pdb_content)
        
        self.assertEqual(len(ligands), 1)
        self.assertEqual(ligands[0]["n_atoms"], 5)


class TestFindLigandsInLegacyPdbFile(unittest.TestCase):
    """Test the find_ligands_in_legacy_pdb_file function."""

    def test_reads_and_finds_ligands_from_file(self):
        """Test reading PDB file from disk and finding ligands."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_file = Path(tmpdir) / "test.pdb"
            pdb_file.write_text(SAMPLE_PDB_CONTENT)
            
            ligands = pdb_file_processing.find_ligands_in_legacy_pdb_file(pdb_file)
            
            self.assertEqual(len(ligands), 2)
            self.assertEqual(ligands[0]["resname"], "GOL")
            self.assertEqual(ligands[1]["resname"], "LIG")

    def test_handles_nonexistent_file(self):
        """Test that reading nonexistent file raises error."""
        with self.assertRaises(FileNotFoundError):
            pdb_file_processing.find_ligands_in_legacy_pdb_file("nonexistent.pdb")


class TestEnsureEntryPdbFile(unittest.TestCase):
    """Test the ensure_entry_pdb_file function."""

    @patch("md_workflows.pdb_file_processing.rcsb_api_request")
    def test_downloads_pdb_if_missing(self, mock_api_request: Mock):
        """Test that PDB is downloaded if file doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir)
            
            # Mock the download to create the file
            def mock_download(endpoint: str, out_path: Path) -> Path:
                out_path.write_text(SAMPLE_PDB_CONTENT)
                return out_path
            
            mock_api_request.side_effect = mock_download
            
            result = pdb_file_processing.ensure_entry_pdb_file("6B8X", out_dir)
            
            self.assertEqual(result, out_dir / "6B8X.pdb")
            self.assertTrue(result.exists())
            mock_api_request.assert_called_once()
            self.assertIn("6b8x.pdb", mock_api_request.call_args[0][0])

    def test_returns_existing_file_without_download(self):
        """Test that existing PDB file is returned without downloading."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir)
            existing_file = out_dir / "6B8X.pdb"
            existing_file.write_text(SAMPLE_PDB_CONTENT)
            
            with patch("md_workflows.pdb_file_processing.rcsb_api_request") as mock_api:
                result = pdb_file_processing.ensure_entry_pdb_file("6B8X", out_dir)
                
                self.assertEqual(result, existing_file)
                mock_api.assert_not_called()

    def test_validates_pdb_id_format(self):
        """Test that invalid PDB IDs raise ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir)
            
            # Too short
            with self.assertRaisesRegex(ValueError, "valid 4-letter PDB ID"):
                pdb_file_processing.ensure_entry_pdb_file("6B8", out_dir)
            
            # Too long
            with self.assertRaisesRegex(ValueError, "valid 4-letter PDB ID"):
                pdb_file_processing.ensure_entry_pdb_file("6B8XX", out_dir)
            
            # Invalid characters
            with self.assertRaisesRegex(ValueError, "valid 4-letter PDB ID"):
                pdb_file_processing.ensure_entry_pdb_file("6B8!", out_dir)

    def test_creates_output_directory_if_missing(self):
        """Test that output directory is created if it doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir) / "subdir" / "nested"
            
            with patch("md_workflows.pdb_file_processing.rcsb_api_request") as mock_api:
                def mock_download(endpoint: str, out_path: Path) -> Path:
                    out_path.write_text(SAMPLE_PDB_CONTENT)
                    return out_path
                
                mock_api.side_effect = mock_download
                
                result = pdb_file_processing.ensure_entry_pdb_file("6B8X", out_dir)
                
                self.assertTrue(out_dir.exists())
                self.assertTrue(result.exists())


class TestDownloadRcsbLegacyPdbAndFindLigands(unittest.TestCase):
    """Test the download_rcsb_legacy_pdb_and_find_ligands function."""

    @patch("md_workflows.pdb_file_processing.rcsb_api_request")
    def test_downloads_and_finds_ligands(self, mock_api_request: Mock):
        """Test downloading PDB and finding ligands."""
        with tempfile.TemporaryDirectory() as tmpdir:
            save_path = Path(tmpdir) / "test.pdb"
            
            # Mock the download
            def mock_download(endpoint: str, out_path: Path) -> Path:
                out_path.write_text(SAMPLE_PDB_CONTENT)
                return out_path
            
            mock_api_request.side_effect = mock_download
            
            pdb_path, ligands = pdb_file_processing.download_rcsb_legacy_pdb_and_find_ligands(
                "6B8X", save_path=save_path
            )
            
            self.assertEqual(pdb_path, save_path)
            self.assertTrue(pdb_path.exists())
            self.assertEqual(len(ligands), 2)
            self.assertEqual(ligands[0]["resname"], "GOL")
            self.assertEqual(ligands[1]["resname"], "LIG")

    @patch("md_workflows.pdb_file_processing.rcsb_api_request")
    def test_default_save_path(self, mock_api_request: Mock):
        """Test that default save path is used when not specified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Change to temp directory
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            
            try:
                def mock_download(endpoint: str, out_path: Path) -> Path:
                    out_path.write_text(SAMPLE_PDB_CONTENT)
                    return out_path
                
                mock_api_request.side_effect = mock_download
                
                pdb_path, ligands = pdb_file_processing.download_rcsb_legacy_pdb_and_find_ligands("6B8X")
                
                self.assertEqual(pdb_path.name, "6b8x_rcsb_legacy.pdb")
                self.assertTrue(pdb_path.exists())
            finally:
                os.chdir(original_dir)

    @patch("md_workflows.pdb_file_processing.rcsb_api_request")
    def test_validates_pdb_id(self, mock_api_request: Mock):
        """Test that invalid PDB IDs raise ValueError."""
        with self.assertRaisesRegex(ValueError, "valid 4-letter PDB ID"):
            pdb_file_processing.download_rcsb_legacy_pdb_and_find_ligands("INVALID123")


class TestPreparePdbAndResnFiles(unittest.TestCase):
    """Test the prepare_pdb_and_resn_files function."""

    @patch("md_workflows.pdb_file_processing.download_rcsb_legacy_pdb_and_find_ligands")
    def test_downloads_if_file_missing(self, mock_download: Mock):
        """Test that PDB is downloaded if file doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            lig_dir = Path(tmpdir)
            
            # Mock download to return single ligand
            mock_pdb_path = lig_dir / "6b8x_rcsb_legacy.pdb"
            mock_ligands = [{"resname": "GOL", "chain_id": "A", "residue_seq": "201", "insertion_code": "", "n_atoms": 3}]
            
            # Mock download creates the file
            def mock_download_side_effect(pdb_id, save_path):
                save_path.write_text(SAMPLE_PDB_SINGLE_LIGAND)
                return (save_path, mock_ligands)
            
            mock_download.side_effect = mock_download_side_effect
            
            # Don't create the file beforehand - let the mock create it
            resn, pdb_path = pdb_file_processing.prepare_pdb_and_resn_files(lig_dir, "6B8X")
            
            self.assertEqual(resn, "GOL")
            self.assertEqual(pdb_path, mock_pdb_path)
            mock_download.assert_called_once()

    def test_uses_existing_file(self):
        """Test that existing PDB file is used without downloading."""
        with tempfile.TemporaryDirectory() as tmpdir:
            lig_dir = Path(tmpdir)
            pdb_path = lig_dir / "6b8x_rcsb_legacy.pdb"
            # Use single ligand PDB content
            pdb_path.write_text(SAMPLE_PDB_SINGLE_LIGAND)
            
            with patch("md_workflows.pdb_file_processing.download_rcsb_legacy_pdb_and_find_ligands") as mock_download:
                resn, returned_path = pdb_file_processing.prepare_pdb_and_resn_files(lig_dir, "6B8X")
                
                self.assertEqual(resn, "GOL")
                self.assertEqual(returned_path, pdb_path)
                mock_download.assert_not_called()

    def test_raises_error_if_no_ligands(self):
        """Test that error is raised if no ligands are found."""
        with tempfile.TemporaryDirectory() as tmpdir:
            lig_dir = Path(tmpdir)
            pdb_path = lig_dir / "test_rcsb_legacy.pdb"
            pdb_path.write_text(SAMPLE_PDB_WITH_WATER_ONLY)
            
            with self.assertRaisesRegex(RuntimeError, "No non-solvent HETATM ligands found"):
                pdb_file_processing.prepare_pdb_and_resn_files(lig_dir, "test")

    def test_raises_error_if_multiple_ligand_types(self):
        """Test that error is raised if multiple ligand types are found."""
        with tempfile.TemporaryDirectory() as tmpdir:
            lig_dir = Path(tmpdir)
            pdb_path = lig_dir / "test_rcsb_legacy.pdb"
            pdb_path.write_text(SAMPLE_PDB_MULTIPLE_LIGANDS)
            
            with self.assertRaisesRegex(RuntimeError, "Expected exactly one ligand residue name"):
                pdb_file_processing.prepare_pdb_and_resn_files(lig_dir, "test")

    def test_accepts_single_ligand_type(self):
        """Test successful case with exactly one ligand type."""
        with tempfile.TemporaryDirectory() as tmpdir:
            lig_dir = Path(tmpdir)
            pdb_path = lig_dir / "test_rcsb_legacy.pdb"
            
            # Create PDB with only GOL ligands (multiple instances of same type)
            pdb_content = """HEADER    TEST
HETATM  100  C1  GOL A 201      15.000  25.000  35.000  1.00 30.00           C
HETATM  101  C2  GOL A 201      16.000  26.000  36.000  1.00 30.00           C
HETATM  200  C1  GOL B 301      15.000  25.000  35.000  1.00 30.00           C
HETATM  201  C2  GOL B 301      16.000  26.000  36.000  1.00 30.00           C
END
"""
            pdb_path.write_text(pdb_content)
            
            resn, returned_path = pdb_file_processing.prepare_pdb_and_resn_files(lig_dir, "test")
            
            self.assertEqual(resn, "GOL")
            self.assertEqual(returned_path, pdb_path)


if __name__ == "__main__":
    unittest.main()
