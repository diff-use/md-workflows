"""Ligand parameterization using Gaussian and AmberTools.

Corresponds to run_all.sh lines:
  mkdir -p ligand
  cd ligand
  bash ../run_params_gaussian.sh
  cd -
"""

from __future__ import annotations

import os
import subprocess
import textwrap
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import TypedDict


class HeteroLigandHit(TypedDict):
    """One unique hetero residue from a legacy PDB (HETATM block)."""

    resname: str
    chain_id: str
    residue_seq: str
    insertion_code: str
    n_atoms: int


# Water and common solvent residue names in legacy PDBs (exclude from ligand hits).
_WATER_SOLVENT_RESNAMES = frozenset({
    "HOH", "WAT", "SOL", "H2O", "DOD", "TIP", "TIP3", "SPC", "PE4", "P7G",
})

# Base URL for legacy-format PDB coordinate files from RCSB PDB.
_RCSB_LEGACY_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def download_rcsb_legacy_pdb_and_find_ligands(
    pdb_id: str,
    *,
    save_path: Path | str | None = None,
    exclude_water_solvent: bool = True,
) -> tuple[Path, list[HeteroLigandHit]]:
    """Fetch legacy PDB text from RCSB, save locally, and list hetero small molecules.

    Uses the same legacy PDB endpoint as browsing ``RCSB.org``: ``download/<ID>.pdb``.
    Parses ``HETATM`` records and groups atoms by residue (chain + seq + insertion +
    residue name). Optionally drops standard water/solvent residue codes.

    Args:
        pdb_id: Four-character PDB identifier (e.g. ``"6b8x"``, case-insensitive).
        save_path: Where to write the file. Defaults to ``<pdb_id>_rcsb_legacy.pdb``
            in the current directory.
        exclude_water_solvent: If True, omit residues in ``_WATER_SOLVENT_RESNAMES``.

    Returns:
        Path to the downloaded legacy PDB file, and a sorted list of hetero summaries.

    Raises:
        urllib.error.HTTPError: If download fails (e.g. unknown PDB ID).
        ValueError: If ``pdb_id`` is not a plausible PDB identifier.
    """
    nid = pdb_id.strip().lower()
    if len(nid) != 4 or not nid.isalnum():
        raise ValueError(f"Expected a 4-letter PDB ID, got {pdb_id!r}")

    out = Path(save_path) if save_path is not None else Path(f"{nid}_rcsb_legacy.pdb")
    url = _RCSB_LEGACY_PDB_URL.format(pdb_id=nid)

    req = urllib.request.Request(url, headers={"User-Agent": "md-workflows/1.0"})
    out.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(req, timeout=120) as resp:
        body = resp.read().decode("ascii", errors="replace")

    out.write_text(body, encoding="ascii")

    ligands = find_ligands_in_legacy_pdb_text(
        body,
        exclude_water_solvent=exclude_water_solvent,
    )
    return out, ligands


def find_ligands_in_legacy_pdb_text(
    pdb_text: str,
    *,
    exclude_water_solvent: bool = True,
) -> list[HeteroLigandHit]:
    """Scan legacy PDB contents for hetero residues from ``HETATM`` records."""

    counts: dict[tuple[str, str, str, str], int] = defaultdict(int)

    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        if len(line) < 27:
            continue
        # PDB fixed columns (legacy): resName 18–20, chain 22, seq 23–26, iCode 27.
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


def find_ligands_in_legacy_pdb_file(path: Path | str) -> list[HeteroLigandHit]:
    """Read a legacy PDB file from disk and list hetero residues (same rules as above)."""

    text = Path(path).read_text(encoding="ascii", errors="replace")
    return find_ligands_in_legacy_pdb_text(text)


def _patch_gaussian_input(resn: str, nproc: int):
    """Insert NprocShared and switch HF -> B3LYP/6-31+G(d,p)."""
    with open(f"{resn}.gau") as fh:
        lines = fh.readlines()

    patched = []
    for line in lines:
        if "Link" in line:
            patched.append(line)
            patched.append(f"%NprocShared={nproc}\n")
        else:
            patched.append(line)

    text = "".join(patched)
    text = text.replace("#HF", "#B3LYP/6-31+G(d,p)")

    with open("tmp", "w") as fh:
        fh.write(text)


def _run_gaussian_opt(resn: str, nproc: int):
    """Run Gaussian geometry optimization if needed."""
    gau_file = f"{resn}.gau"
    log_file = f"{resn}.log"
    env = {**os.environ, "OMP_NUM_THREADS": str(nproc)}

    should_run = False
    if not os.path.exists(log_file):
        print("Running Gaussian optimization")
        should_run = True
    else:
        with open(gau_file) as a, open("tmp") as b:
            if a.read() != b.read():
                print("Gaussian input file is different. Running Gaussian optimization.")
                should_run = True
            else:
                print("Gaussian input file is the same. Skipping run.")

    if should_run:
        os.replace("tmp", gau_file)
        subprocess.run(["g16", gau_file], env=env, check=True)


def _run_gaussian_esp(resn: str, nproc: int):
    """Run Gaussian ESP / CHELPG calculation if needed."""
    with open(f"{resn}.gau") as fh:
        text = fh.read()

    text = text.replace("opt", "pop(chelpg,regular)")
    text = text.replace("molecule", "grid")

    resp_gau = f"{resn}_resp.gau"
    resp_log = f"{resn}_resp.log"
    env = {**os.environ, "OMP_NUM_THREADS": str(nproc)}

    should_run = False
    if not os.path.exists(resp_log):
        print("Running Gaussian ESP calculation.")
        should_run = True
    else:
        existing = ""
        if os.path.exists(resp_gau):
            with open(resp_gau) as fh:
                existing = fh.read()
        if existing != text:
            print("Gaussian input file is different. Running Gaussian ESP calculation.")
            should_run = True
        else:
            print("Gaussian input file is the same. Skipping run.")

    if should_run:
        with open(resp_gau, "w") as fh:
            fh.write(text)
        subprocess.run(["g16", resp_gau], env=env, check=True)


def _process_resp_charges(resn: str):
    """Derive RESP charges from Gaussian output and correct net charge."""
    subprocess.run([
        "antechamber", "-fi", "gout",
        "-i", f"{resn}_resp.log",
        "-cf", f"{resn}_resp.crg",
        "-c", "resp",
        "-o", f"{resn}_gauss.ac", "-fo", "ac", "-rn", resn,
    ], check=True)

    subprocess.run([
        "antechamber", "-fi", "gout",
        "-i", f"{resn}_resp.log",
        "-o", f"{resn}_gauss.pdb", "-fo", "pdb", "-rn", resn,
    ], check=True)

    orig_coords = _extract_coords_from_pdb(f"{resn}.pdb")
    _graft_coords_to_ac(f"{resn}_gauss.ac", orig_coords, f"{resn}_resp.ac")
    _correct_charge(f"{resn}_resp.ac")

    subprocess.run([
        "antechamber", "-fi", "ac", "-i", f"{resn}_resp.ac",
        "-fo", "mol2", "-o", f"{resn}_resp.mol2", "-rn", resn,
    ], check=True)

    subprocess.run([
        "atomtype", "-i", f"{resn}_resp.ac",
        "-o", f"{resn}_resp_gaff.ac", "-p", "gaff",
    ], check=True)

    subprocess.run([
        "prepgen", "-i", f"{resn}_resp_gaff.ac",
        "-o", f"{resn}_resp_gaff.prepc", "-f", "car",
    ], check=True)

    subprocess.run([
        "parmchk2", "-i", f"{resn}_resp_gaff.prepc",
        "-o", f"{resn}_resp.frcmod", "-f", "prepc",
    ], check=True)


def _extract_coords_from_pdb(pdb_file: str) -> list[str]:
    coords = []
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                coords.append(line[30:53])
    return coords


def _graft_coords_to_ac(ac_file: str, coords: list[str], out_file: str):
    anum = 0
    with open(ac_file) as fh, open(out_file, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                out.write(line[:30] + coords[anum] + line[53:])
                anum += 1
            else:
                out.write(line)


def _correct_charge(ac_file: str):
    """Adjust the largest-magnitude charge to make the total an exact integer."""
    charges = []
    with open(ac_file) as fh:
        lines = fh.readlines()
    for line in lines:
        if line.startswith("ATOM"):
            charges.append(float(line[54:63]))

    total = sum(charges)
    nearest_int = round(total)
    correction = nearest_int - total

    max_idx = max(range(len(charges)), key=lambda i: abs(charges[i]))
    new_charge = charges[max_idx] + correction

    anum = 0
    new_lines = []
    for line in lines:
        if line.startswith("ATOM"):
            if anum == max_idx:
                new_lines.append(f"{line[:54]}{new_charge:9.6f}{line[63:]}")
            else:
                new_lines.append(line)
            anum += 1
        else:
            new_lines.append(line)

    with open(ac_file, "w") as fh:
        fh.writelines(new_lines)


def _build_amber_lib(resn: str):
    """Run tleap to create Amber parameter/topology files."""
    tleap_input = textwrap.dedent(f"""\
        source leaprc.protein.ff19SB
        source leaprc.gaff
        loadamberparams {resn}_resp.frcmod
        loadamberprep {resn}_resp_gaff.prepc
        lig = loadmol2 {resn}_resp.mol2
        savepdb lig {resn}_resp.pdb
        saveamberparm lig {resn}_resp.parm7 {resn}_resp.rst7
        quit
    """)
    with open("tleap_lig.in", "w") as fh:
        fh.write(tleap_input)
    subprocess.run(["tleap", "-f", "tleap_lig.in"], check=True)


def _run_parameterization(resn: str, g16root: str, nproc: int):
    os.environ["g16root"] = g16root
    g16_profile = os.path.join(g16root, "g16", "bsd", "g16.profile")
    if os.path.exists(g16_profile):
        subprocess.run(["bash", "-c", f"source {g16_profile}"], check=True)

    subprocess.run([
        "antechamber", "-fi", "pdb", "-fo", "gcrt",
        "-i", f"{resn}.pdb", "-o", f"{resn}.gau",
        "-nc", "-2", "-m", "1",
    ], check=True)

    _patch_gaussian_input(resn, nproc)
    _run_gaussian_opt(resn, nproc)
    _run_gaussian_esp(resn, nproc)
    _process_resp_charges(resn)
    _build_amber_lib(resn)


def run(
    resn: str = "AR6",
    g16root: str = "/Users/mewall/packages",
    nproc: int = 8,
    pdb_id: str | None = None,
):
    base_dir = Path.cwd()
    lig_dir = base_dir / "ligand"
    lig_dir.mkdir(parents=True, exist_ok=True)

    if pdb_id:
        pdb_path, ligands = download_rcsb_legacy_pdb_and_find_ligands(
            pdb_id,
            save_path=lig_dir / f"{pdb_id.strip().lower()}_rcsb_legacy.pdb",
        )
        print(f"Downloaded legacy PDB from RCSB -> {pdb_path}")
        if ligands:
            print(
                "Hetero residues from HETATM records "
                "(common water/solvent residue names excluded): "
                f"{len(ligands)}"
            )
            for h in ligands:
                ins = h["insertion_code"] or ""
                print(
                    f"  {h['resname']:>4}  chain={h['chain_id']!r}  "
                    f"resSeq={h['residue_seq']}{ins}  atoms={h['n_atoms']}"
                )
        else:
            print(
                "No HETATM residues found after filtering (structure may be polymer-only "
                "or use only solvent codes)."
            )

    os.chdir(lig_dir)

    try:
        _run_parameterization(resn, g16root, nproc)
    finally:
        os.chdir(base_dir)



if __name__ == "__main__":
    run()
