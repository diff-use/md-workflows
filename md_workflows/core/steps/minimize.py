"""Energy-minimize the solvated MD model.

Reference implementation of the step contract (Inputs / resolve / guard / typed
Result) that the other steps follow. Corresponds to ``minimize.sh``:
``gmx grompp -f min.mdp -c md_model.pdb -o md_min.tpr -p md_model.top`` then
``gmx mdrun -deffnm md_min``.
"""

from __future__ import annotations

from pathlib import Path

from ..config import GromacsRunProfile, SystemConfig
from ..gmx import checksums, gmx_version, run_min
from ..paths import under
from ..results import StepInputs, StepResult, StepStatus

STEP = "minimize"
DEFFNM = "md_min"


class MinimizeInputs(StepInputs):
    model_pdb: Path
    model_top: Path
    min_mdp: Path

    def consumed_paths(self) -> list[Path]:
        return [self.model_pdb, self.model_top, self.min_mdp]


class MinimizeResult(StepResult):
    @property
    def gro(self) -> Path:
        return self.output("gro")

    @property
    def tpr(self) -> Path:
        return self.output("tpr")


def resolve_inputs(workdir: Path, cfg: SystemConfig) -> MinimizeInputs:
    """Map the conventional filenames in ``workdir`` to explicit input paths."""
    workdir = Path(workdir)
    mdp_dir = under(workdir, cfg.mdp_dir)
    return MinimizeInputs(
        workdir=workdir,
        model_pdb=workdir / "md_model.pdb",
        model_top=workdir / "md_model.top",
        min_mdp=mdp_dir / cfg.minimize.min_mdp,
    )


def check_inputs(inputs: MinimizeInputs) -> None:
    inputs.check_exists(STEP)


def minimize(
    inputs: MinimizeInputs,
    profile: GromacsRunProfile,
    *,
    resume: bool = False,
) -> MinimizeResult:
    """Run energy minimization and return a typed, metadata-carrying result."""
    check_inputs(inputs)

    consumed = {
        "model_pdb": inputs.model_pdb,
        "model_top": inputs.model_top,
        "min_mdp": inputs.min_mdp,
    }
    tpr = inputs.workdir / f"{DEFFNM}.tpr"
    gro = inputs.workdir / f"{DEFFNM}.gro"

    if resume and gro.exists():
        return MinimizeResult(
            step=STEP,
            status=StepStatus.SKIPPED,
            workdir=inputs.workdir,
            outputs={"tpr": tpr, "gro": gro},
            run_profiles_used={"min": profile},
            input_checksums=checksums(consumed),
        )

    run = run_min(
        inputs.model_pdb, inputs.model_top, inputs.min_mdp, DEFFNM, profile, inputs.workdir
    )
    return MinimizeResult(
        step=STEP,
        status=StepStatus.COMPLETED,
        workdir=inputs.workdir,
        outputs={"tpr": run.tpr, "gro": run.gro},
        run_profiles_used={"min": profile},
        input_checksums=checksums(consumed),
        tool_versions={"gmx": gmx_version(profile.gmx_bin)},
        log_paths=[run.log],
    )
