"""GROMACS-level primitives shared across steps.

These wrap one external operation each (grompp+mdrun, solvate, energy, topology edit)
so the higher-level steps and the resolvation flow reuse a single implementation
instead of duplicating the min/equil/solvate motif. All failures raise structured
exceptions; nothing here calls ``sys.exit``.
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

from pydantic import BaseModel

from .config import GromacsRunProfile
from .exceptions import GmxOutputParseError, StepToolError

_WATER_RESNAMES = frozenset({"WAT", "SOL", "HOH", "H2O"})
_WATER_OXYGENS = frozenset({"OW", "O", "OW1"})


# --------------------------------------------------------------------------- #
# Process / hashing helpers
# --------------------------------------------------------------------------- #
def run_tool(
    cmd: list[str],
    *,
    tool: str,
    cwd: str | Path | None = None,
    input: str | None = None,
    log_path: str | Path | None = None,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run an external tool, capturing output; raise :class:`StepToolError` on failure.

    When ``log_path`` is given, the combined stdout+stderr is written there (matching
    the shell scripts' ``>& tool.log`` redirects).
    """
    proc = subprocess.run(
        [str(c) for c in cmd],
        cwd=str(cwd) if cwd is not None else None,
        input=input,
        text=True,
        capture_output=True,
        env=env,
    )
    if log_path is not None:
        Path(log_path).write_text((proc.stdout or "") + (proc.stderr or ""))
    if proc.returncode != 0:
        raise StepToolError(
            tool,
            proc.returncode,
            cmd=[str(c) for c in cmd],
            stderr=proc.stderr,
            log_path=Path(log_path) if log_path is not None else None,
        )
    return proc


def sha256_file(path: str | Path) -> str:
    """Return the hex SHA-256 of a file (streamed)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def checksums(paths: dict[str, Path]) -> dict[str, str]:
    """SHA-256 each existing path in a name->path mapping; skip missing."""
    return {name: sha256_file(p) for name, p in paths.items() if Path(p).exists()}


def gmx_version(gmx_bin: str = "gmx") -> str:
    """Best-effort GROMACS version string; returns 'unknown' if it can't be read."""
    try:
        proc = subprocess.run([gmx_bin, "--version"], text=True, capture_output=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    for line in (proc.stdout or "").splitlines():
        if "version" in line.lower():
            return line.split(":", 1)[-1].strip() or line.strip()
    return "unknown"


# --------------------------------------------------------------------------- #
# Structure helpers
# --------------------------------------------------------------------------- #
def count_waters(structure: str | Path) -> int:
    """Count water molecules in a ``.gro`` or ``.pdb`` file (one oxygen per molecule).

    Unifies the two heuristics the scripts used (WAT-atoms/3 in make_waterbox, OW-count
    in pressure_interpolate): counting water oxygens gives the exact molecule count.
    """
    path = Path(structure)
    is_gro = path.suffix.lower() == ".gro"
    count = 0
    with open(path) as fh:
        lines = fh.readlines()
    if is_gro:
        # GRO fixed columns: resname [5:10], atom name [10:15]. Skip title + count +
        # the trailing box-vector line.
        for line in lines[2:-1]:
            resname = line[5:10].strip()
            atom = line[10:15].strip()
            if resname in _WATER_RESNAMES and atom in _WATER_OXYGENS:
                count += 1
    else:
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip()
            atom = line[12:16].strip()
            if resname in _WATER_RESNAMES and atom in _WATER_OXYGENS:
                count += 1
    return count


# --------------------------------------------------------------------------- #
# mdrun primitives
# --------------------------------------------------------------------------- #
class GmxRunResult(BaseModel):
    """Files produced by a grompp+mdrun invocation (paths; may not all exist)."""

    deffnm: str
    tpr: Path
    gro: Path
    edr: Path
    log: Path


def _grompp(
    profile: GromacsRunProfile,
    *,
    mdp: Path,
    structure: Path,
    top: Path,
    tpr: Path,
    ref: Path | None,
    maxwarn: int,
    workdir: Path,
    grompp_log: Path,
) -> None:
    cmd = [
        profile.gmx_bin,
        "grompp",
        "-f",
        str(mdp),
        "-c",
        str(structure),
        "-o",
        str(tpr),
        "-p",
        str(top),
    ]
    if ref is not None:
        cmd += ["-r", str(ref)]
    if maxwarn:
        cmd += ["-maxwarn", str(maxwarn)]
    run_tool(cmd, tool="gmx grompp", cwd=workdir, log_path=grompp_log)


def _mdrun(profile: GromacsRunProfile, *, tpr: Path, deffnm: str, workdir: Path) -> None:
    # -deffnm is relative, so mdrun must run with cwd=workdir; outputs land there.
    cmd = [
        profile.gmx_bin,
        "mdrun",
        "-s",
        str(tpr),
        *profile.to_mdrun_flags(),
        "-deffnm",
        deffnm,
        "-v",
    ]
    run_tool(cmd, tool="gmx mdrun", cwd=workdir, log_path=workdir / f"mdrun_{deffnm}.log")


def run_min(
    structure: Path,
    top: Path,
    mdp: Path,
    deffnm: str,
    profile: GromacsRunProfile,
    workdir: Path,
) -> GmxRunResult:
    """grompp + mdrun for an energy minimization; outputs ``<deffnm>.{tpr,gro,edr,log}``."""
    workdir = Path(workdir)
    tpr = workdir / f"{deffnm}.tpr"
    _grompp(
        profile,
        mdp=mdp,
        structure=structure,
        top=top,
        tpr=tpr,
        ref=None,
        maxwarn=0,
        workdir=workdir,
        grompp_log=workdir / f"grompp_{deffnm}.log",
    )
    _mdrun(profile, tpr=tpr, deffnm=deffnm, workdir=workdir)
    return GmxRunResult(
        deffnm=deffnm,
        tpr=tpr,
        gro=workdir / f"{deffnm}.gro",
        edr=workdir / f"{deffnm}.edr",
        log=workdir / f"{deffnm}.log",
    )


def run_equil(
    structure: Path,
    top: Path,
    mdp: Path,
    ref: Path,
    deffnm: str,
    profile: GromacsRunProfile,
    workdir: Path,
    *,
    maxwarn: int = 0,
) -> GmxRunResult:
    """grompp (with restraint ref ``-r``) + mdrun for an equilibration."""
    workdir = Path(workdir)
    tpr = workdir / f"{deffnm}.tpr"
    _grompp(
        profile,
        mdp=mdp,
        structure=structure,
        top=top,
        tpr=tpr,
        ref=ref,
        maxwarn=maxwarn,
        workdir=workdir,
        grompp_log=workdir / f"grompp_{deffnm}.log",
    )
    _mdrun(profile, tpr=tpr, deffnm=deffnm, workdir=workdir)
    return GmxRunResult(
        deffnm=deffnm,
        tpr=tpr,
        gro=workdir / f"{deffnm}.gro",
        edr=workdir / f"{deffnm}.edr",
        log=workdir / f"{deffnm}.log",
    )


# --------------------------------------------------------------------------- #
# solvate / energy primitives
# --------------------------------------------------------------------------- #
class AddWaterResult(BaseModel):
    out_pdb: Path
    nwat_added: int
    log: Path


def add_water(
    input_gro: Path,
    reservoir_gro: Path,
    maxsol: int | None,
    out_pdb: Path,
    workdir: Path,
    *,
    gmx_bin: str = "gmx",
    log_path: Path | None = None,
) -> AddWaterResult:
    """``gmx solvate`` to add water from a reservoir; parse the count actually added.

    ``maxsol=None`` runs an unbounded trial solvation (used to size the fill).
    """
    workdir = Path(workdir)
    log_path = Path(log_path) if log_path is not None else workdir / "gmx_solvate.log"
    cmd = [gmx_bin, "solvate", "-cp", str(input_gro), "-cs", str(reservoir_gro), "-o", str(out_pdb)]
    if maxsol is not None:
        cmd += ["-maxsol", str(maxsol)]
    run_tool(cmd, tool="gmx solvate", cwd=workdir, log_path=log_path)
    added = _parse_solvent_count(log_path)
    return AddWaterResult(out_pdb=out_pdb, nwat_added=added, log=log_path)


def _parse_solvent_count(log_path: Path) -> int:
    for line in Path(log_path).read_text().splitlines():
        if "Number of solvent molecules:" in line:
            for token in reversed(line.split()):
                if token.isdigit():
                    return int(token)
    raise GmxOutputParseError("number of solvent molecules", log_path)


def read_mean_pressure(edr: Path, workdir: Path, *, gmx_bin: str = "gmx") -> float:
    """Mean pressure (bar) from an ``.edr`` via ``gmx energy``."""
    workdir = Path(workdir)
    tmp = workdir / f".pressure_{Path(edr).stem}.xvg"
    try:
        proc = run_tool(
            [gmx_bin, "energy", "-f", str(edr), "-o", str(tmp)],
            tool="gmx energy",
            cwd=workdir,
            input="Pressure\n0\n",
        )
    finally:
        tmp.unlink(missing_ok=True)
    output = (proc.stdout or "") + (proc.stderr or "")
    for line in output.splitlines():
        if line.strip().startswith("Pressure"):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    return float(parts[1])
                except ValueError:
                    continue
    raise GmxOutputParseError("mean Pressure", edr, gmx_version=gmx_version(gmx_bin))


# --------------------------------------------------------------------------- #
# Topology WAT-block management
# --------------------------------------------------------------------------- #
def manage_wat_topology(top: Path, count: int, *, mode: str = "append") -> None:
    """Maintain the ordered ``WAT`` blocks in a topology's ``[ molecules ]`` section.

    Waters are added at several stages and sit in non-contiguous coordinate blocks
    (initial waters -> Cl-/Na+ -> each resolvation stage's waters), so the topology
    must carry a matching *sequence* of ``WAT`` lines, not one summed line.

    - ``mode="append"``: append ``WAT <count>`` (idempotent — skips if the file already
      ends with that exact block, so re-running a stage does not double-add).
    - ``mode="drop_last_then_append"``: drop the trailing ``WAT`` block (the discarded
      run-2 probe) then append ``WAT <count>`` — the final resolvation stage.
    """
    top = Path(top)
    lines = top.read_text().splitlines()
    if "[ molecules ]" not in "\n".join(lines):
        raise ValueError(f"{top}: no '[ molecules ]' section to manage WAT blocks in")

    while lines and not lines[-1].strip():
        lines.pop()

    if mode == "drop_last_then_append":
        if lines and lines[-1].split()[:1] == ["WAT"]:
            lines.pop()
    elif mode != "append":
        raise ValueError(f"unknown manage_wat_topology mode {mode!r}")

    target = f"WAT {count}"
    if not (lines and lines[-1].strip() == target):
        lines.append(target)

    top.write_text("\n".join(lines) + "\n")
