"""PDB ID download and ligand residue file preparation helpers."""

from __future__ import annotations

import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import TypedDict
import platform


class HeteroLigandHit(TypedDict):
    """One unique hetero residue from a legacy PDB (HETATM block)."""

    resname: str
    chain_id: str
    residue_seq: str
    insertion_code: str
    n_atoms: int


_WATER_SOLVENT_RESNAMES = frozenset({
    "HOH", "WAT", "SOL", "H2O", "DOD", "TIP", "TIP3", "SPC", "PE4", "P7G",
})

#Makes a GET request to the RCSB API and saves the response to a file.
def rcsb_api_request(endpoint: str, out_path: Path) -> Path:
    base_url = "https://files.rcsb.org/"
    req = urllib.request.Request(
        url=f"{base_url}{endpoint}",
        data=None,
        headers={"User-Agent": f"{platform.node()} {platform.system()}"},
        method="GET"
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        out_path.write_bytes(resp.read())
    return out_path

def download_rcsb_legacy_pdb_and_find_ligands(
    pdb_id: str,
    *,
    save_path: Path | str | None = None,
    exclude_water_solvent: bool = True,
) -> tuple[Path, list[HeteroLigandHit]]:
    """Fetch legacy-format PDB from RCSB, save locally, list hetero residues."""
    nid = pdb_id.strip().lower()
    if len(nid) != 4 or not nid.isalnum():
        raise ValueError(f"Expected a valid 4-letter PDB ID, got {pdb_id!r}")

    out = Path(save_path) if save_path is not None else Path(f"{nid}_rcsb_legacy.pdb")
    out.parent.mkdir(parents=True, exist_ok=True)

    out = rcsb_api_request(f"download/{nid}.pdb", out)
    return out, find_ligands_in_legacy_pdb_file(out, exclude_water_solvent=exclude_water_solvent)


def ensure_entry_pdb_file(pdb_id: str, out_dir: Path) -> Path:
    """Ensure <PDBID>.pdb exists in out_dir; download from RCSB if missing."""
    nid = pdb_id.strip()
    if len(nid) != 4 or not nid.isalnum():
        raise ValueError(f"Expected a valid 4-letter PDB ID, got {pdb_id!r}")
    out = out_dir / f"{nid.upper()}.pdb"
    if not out.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        rcsb_api_request(f"download/{nid.lower()}.pdb", out)
        print(f"Downloaded entry PDB from RCSB -> {out}")
    return out


def find_ligands_in_legacy_pdb_text(
    pdb_text: str,
    *,
    exclude_water_solvent: bool = True,
) -> list[HeteroLigandHit]:
    """Scan legacy PDB contents for hetero residues from HETATM records."""
    counts: dict[tuple[str, str, str, str], int] = defaultdict(int)

    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        if len(line) < 27:
            continue
        resname = line[17:20].strip()
        chain_id = line[21:22].strip() or "_"
        residue_seq = line[22:26].strip()
        icode = line[26:27].strip() if len(line) > 26 else ""

        if exclude_water_solvent and resname.upper() in _WATER_SOLVENT_RESNAMES:
            continue

        key = (resname, chain_id, residue_seq, icode)
        counts[key] += 1

    hits: list[HeteroLigandHit] = []
    for (resname, chain_id, residue_seq, icode), n_atoms in counts.items():
        hits.append(
            {
                "resname": resname,
                "chain_id": chain_id,
                "residue_seq": residue_seq,
                "insertion_code": icode,
                "n_atoms": n_atoms,
            }
        )

    def _sort_key(h: HeteroLigandHit) -> tuple[str, int, str, str]:
        try:
            snum = int(h["residue_seq"])
        except ValueError:
            snum = 0
        return (h["chain_id"], snum, h["insertion_code"], h["resname"])

    hits.sort(key=_sort_key)
    return hits


def find_ligands_in_legacy_pdb_file(
    path: Path | str, *, exclude_water_solvent: bool = True,
) -> list[HeteroLigandHit]:
    """Read a legacy PDB file from disk and list hetero residues."""
    text = Path(path).read_text(encoding="ascii", errors="replace")
    return find_ligands_in_legacy_pdb_text(text, exclude_water_solvent=exclude_water_solvent)


def prepare_pdb_and_resn_files(lig_dir: Path, pdb_id: str = "6B8X") -> tuple[str, Path]:
    """Ensure entry PDB exists and detect residue name from it."""
    pdb_path = lig_dir / f"{pdb_id.strip().lower()}_rcsb_legacy.pdb"
    ligands: list[HeteroLigandHit]
    if pdb_path.exists():
        ligands = find_ligands_in_legacy_pdb_file(pdb_path)
    else:
        pdb_path, ligands = download_rcsb_legacy_pdb_and_find_ligands(pdb_id, save_path=pdb_path)
        print(f"Downloaded legacy PDB from RCSB -> {pdb_path}")

    if not ligands:
        raise RuntimeError("No non-solvent HETATM ligands found in downloaded PDB")

    unique = sorted({h["resname"] for h in ligands})
    if len(unique) != 1:
        raise RuntimeError(f"Expected exactly one ligand residue name, found: {unique}")
    resn = unique[0]
    return resn, pdb_path
