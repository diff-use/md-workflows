"""Energy-minimize the solvated MD model.

Corresponds to run_all.sh line:
  bash scripts/minimize.sh
"""

import subprocess


def run(ntomp: int = 26):
    subprocess.run([
        "gmx", "grompp",
        "-f", "min.mdp",
        "-c", "md_model.pdb",
        "-o", "md_min.tpr",
        "-p", "md_model.top",
    ], capture_output=True, text=True, check=True)

    subprocess.run([
        "gmx", "mdrun",
        "-ntmpi", "1",
        "-ntomp", str(ntomp),
        "-deffnm", "md_min",
        "-v",
    ], check=True)


if __name__ == "__main__":
    run()
