"""Parameterize the protein: clean PDB, build Amber topology, convert to GROMACS.

Corresponds to run_all.sh line:
  bash scripts/param_prot.sh 6B8X
"""

import subprocess
import textwrap
from pathlib import Path


def run(pdb_id: str = "6B8X"):
    _clean_pdb(pdb_id)
    _initial_solvation()
    _amber_to_gromacs()
    _extract_solvent_pdbs()
    _final_protein_topology()
    _amber_to_gromacs()


def _clean_pdb(pdb_id: str):
    """Strip REMARK/KEYWDS, keep CRYST1/ATOM/HETATM/TER/END lines."""
    pdb_file = f"{pdb_id}.pdb"
    with open(pdb_file) as fh:
        lines = fh.readlines()

    kept = []
    for line in lines:
        if line.startswith(("REMARK", "KEYWDS")):
            continue
        if line.startswith(("CRYST1", "ATOM", "HETATM", "TER", "END")):
            kept.append(line)

    with open("pdb_clean.pdb", "w") as fh:
        fh.writelines(kept)

    subprocess.run([
        "pdb4amber", "-i", "pdb_clean.pdb", "--prot", "-o", "pdb_clean_amber.pdb",
    ], check=True)


def _initial_solvation():
    """Run tleap with minimal solvation to obtain solvent PDB templates."""
    tleap_input = textwrap.dedent("""\
        source leaprc.protein.ff19SB
        source leaprc.DNA.OL15
        source leaprc.RNA.OL3
        source leaprc.water.spceb
        source leaprc.gaff2
        p = loadpdb pdb_clean_amber.pdb
        x = combine{p}
        addions2 x Cl- 0
        addions2 x Na+ 0
        addions2 x Na+ 1
        addions2 x Cl- 1
        solvateBox x SPCBOX 1.
        set default PBradii mbondi3
        set default nocenter on
        saveAmberParm x prot.parm7 prot.rst7
        quit
    """)
    with open("tleap_temp.in", "w") as fh:
        fh.write(tleap_input)
    subprocess.run(["tleap", "-f", "tleap_temp.in"], check=True)

    for f in ["prot.top", "prot.pdb"]:
        Path(f).unlink(missing_ok=True)


def _amber_to_gromacs():
    """Convert Amber parm7/rst7 to GROMACS topology and PDB."""
    script = textwrap.dedent("""\
        import parmed as pmd
        parm = pmd.load_file("prot.parm7", "prot.rst7")
        parm.save("prot.top")
        parm.save("prot.pdb")
    """)
    with open("amber_to_gromacs.py", "w") as fh:
        fh.write(script)
    subprocess.run(["python", "amber_to_gromacs.py"], check=True)


def _extract_solvent_pdbs():
    """Pull single Na+, Cl-, and WAT coordinate templates from the solvated PDB."""
    with open("prot.pdb") as fh:
        lines = fh.readlines()

    hetatm_lines = [l for l in lines if l.startswith("HETATM")]

    na_lines = [l for l in hetatm_lines if "Na+" in l]
    cl_lines = [l for l in hetatm_lines if "Cl-" in l]
    wat_lines = [l for l in hetatm_lines if "WAT" in l]

    if na_lines:
        with open("Na+.pdb", "w") as fh:
            fh.write(na_lines[0])
    if cl_lines:
        with open("Cl-.pdb", "w") as fh:
            fh.write(cl_lines[0])
    if wat_lines:
        with open("WAT.pdb", "w") as fh:
            fh.writelines(wat_lines[:3])


def _final_protein_topology():
    """Build the final protein topology with one water + ions."""
    tleap_input = textwrap.dedent("""\
        source leaprc.protein.ff19SB
        source leaprc.DNA.OL15
        source leaprc.RNA.OL3
        source leaprc.water.spceb
        source leaprc.gaff2
        p = loadpdb pdb_clean_amber.pdb
        w = loadpdb WAT.pdb
        x = combine{p w}
        addions2 x Na+ 0
        addions2 x Cl- 0
        addions2 x Na+ 1
        addions2 x Cl- 1
        set default PBradii mbondi3
        set default nocenter on
        saveAmberParm x prot.parm7 prot.rst7
        quit
    """)
    with open("tleap_prot.in", "w") as fh:
        fh.write(tleap_input)
    subprocess.run(["tleap", "-f", "tleap_prot.in"], check=True)

    for f in ["prot.top", "prot.pdb"]:
        Path(f).unlink(missing_ok=True)


if __name__ == "__main__":
    run()
