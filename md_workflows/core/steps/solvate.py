"""Solvate the crystal with water and add neutralizing/ionic-strength ions.

Corresponds to ``solvate.sh``: count symmetry copies, fill the crystal voids from the
equilibrated water reservoir, compute Na+/Cl- counts for the target ionic strength plus
charge neutralization, insert the ions, and finalize the model topology.

The target ionic strength is a tunable parameter (``solvate.ionic_strength``, default
0.1 M) rather than a hard-coded constant, so runs can request a different salt
concentration without editing code.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

from ..config import SolvateParams, SystemConfig
from ..exceptions import GmxOutputParseError
from ..gmx import checksums, count_waters, run_tool
from ..results import StepInputs, StepResult, StepStatus

STEP = "solvate"


class SolvateInputs(StepInputs):
    prot_dry: Path
    xtal: Path
    water_equil: Path
    prot_pdb: Path
    prot_top: Path
    cl_pdb: Path
    na_pdb: Path

    def consumed_paths(self) -> list[Path]:
        return [
            self.prot_dry,
            self.xtal,
            self.water_equil,
            self.prot_pdb,
            self.prot_top,
            self.cl_pdb,
            self.na_pdb,
        ]


class SolvateResult(StepResult):
    @property
    def md_model_pdb(self) -> Path:
        return self.output("md_model_pdb")

    @property
    def md_model_top(self) -> Path:
        return self.output("md_model_top")


def resolve_inputs(workdir: Path, cfg: SystemConfig) -> SolvateInputs:
    workdir = Path(workdir)
    return SolvateInputs(
        workdir=workdir,
        prot_dry=workdir / "prot_dry.pdb",
        xtal=workdir / "xtal.pdb",
        water_equil=workdir / "waterbox" / "water_equil.gro",
        prot_pdb=workdir / "prot.pdb",
        prot_top=workdir / "prot.top",
        cl_pdb=workdir / "Cl-.pdb",
        na_pdb=workdir / "Na+.pdb",
    )


def check_inputs(inputs: SolvateInputs) -> None:
    inputs.check_exists(STEP)


def solvate(
    inputs: SolvateInputs,
    params: SolvateParams,
    *,
    resume: bool = False,
) -> SolvateResult:
    check_inputs(inputs)
    wd = inputs.workdir

    xtal_solv = wd / "xtal_solv.pdb"
    xtal_solv_cl_na = wd / "xtal_solv_cl_na.pdb"
    md_model_pdb = wd / "md_model.pdb"
    md_model_top = wd / "md_model.top"
    outputs = {
        "xtal_solv": xtal_solv,
        "xtal_solv_cl_na": xtal_solv_cl_na,
        "md_model_pdb": md_model_pdb,
        "md_model_top": md_model_top,
    }
    consumed = {
        "prot_dry": inputs.prot_dry,
        "xtal": inputs.xtal,
        "water_equil": inputs.water_equil,
        "prot_pdb": inputs.prot_pdb,
        "prot_top": inputs.prot_top,
        "cl_pdb": inputs.cl_pdb,
        "na_pdb": inputs.na_pdb,
    }

    if resume and md_model_pdb.exists() and md_model_top.exists():
        return SolvateResult(
            step=STEP,
            status=StepStatus.SKIPPED,
            workdir=wd,
            outputs=outputs,
            params=params.model_dump(),
            input_checksums=checksums(consumed),
        )

    ncopies = _count_copies(inputs.prot_dry, inputs.xtal)

    # Fill crystal voids with the equilibrated water reservoir.
    solvate_log = wd / "gmx_solvate.log"
    run_tool(
        [
            "gmx",
            "solvate",
            "-cp",
            str(inputs.xtal),
            "-cs",
            str(inputs.water_equil),
            "-o",
            str(xtal_solv),
        ],
        tool="gmx solvate",
        cwd=wd,
        log_path=solvate_log,
    )
    nwat_initial = _read_solvate_nwat(solvate_log)

    _write_topology_header(inputs.prot_top, md_model_top)

    ions_pos, ions_neg = _count_ions(inputs.prot_pdb)
    net_ion_charge = (ions_pos - ions_neg) * ncopies
    nna, ncl = _compute_ion_counts(
        nwat_initial, net_ion_charge, params.ionic_strength, params.water_molarity
    )

    _insert_ions(wd, xtal_solv, inputs.cl_pdb, inputs.na_pdb, xtal_solv_cl_na, ncl, nna)

    nwat = count_waters(xtal_solv_cl_na)
    _finalize_topology(md_model_top, ncopies, nwat, ncl, nna)
    shutil.copy(xtal_solv_cl_na, md_model_pdb)

    return SolvateResult(
        step=STEP,
        status=StepStatus.COMPLETED,
        workdir=wd,
        outputs=outputs,
        params=params.model_dump(),
        input_checksums=checksums(consumed),
        metrics={
            "ncopies": ncopies,
            "nwat_initial": nwat_initial,
            "ions_pos": ions_pos,
            "ions_neg": ions_neg,
            "net_ion_charge": net_ion_charge,
            "nna": nna,
            "ncl": ncl,
            "nwat": nwat,
        },
        log_paths=[solvate_log],
    )


def _count_copies(prot_dry: Path, xtal: Path) -> int:
    def _natoms(path: Path) -> int:
        n = 0
        with open(path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    n += 1
        return n

    nats_one = _natoms(prot_dry)
    if nats_one == 0:
        raise GmxOutputParseError("asymmetric-unit atom count", prot_dry)
    return _natoms(xtal) // nats_one


def _read_solvate_nwat(log_path: Path) -> int:
    for line in Path(log_path).read_text().splitlines():
        m = re.search(r"Output configuration contains\s+(\d+)", line)
        if m:
            return int(m.group(1))
    raise GmxOutputParseError("output water count", log_path)


def _write_topology_header(prot_top: Path, md_model_top: Path) -> None:
    header: list[str] = []
    with open(prot_top) as fh:
        for line in fh:
            if "molecules" in line.lower():
                break
            header.append(line)
    with open(md_model_top, "w") as fh:
        fh.writelines(header)


def _count_ions(prot_pdb: Path) -> tuple[int, int]:
    pos = neg = 0
    with open(prot_pdb) as fh:
        for line in fh:
            if line.startswith("HETATM"):
                if "Na+" in line:
                    pos += 1
                elif "Cl-" in line:
                    neg += 1
    return pos, neg


def _compute_ion_counts(
    nwat: int, net_ion_charge: int, ionic_strength: float, water_molarity: float
) -> tuple[int, int]:
    """Na+/Cl- counts for the requested ionic strength plus charge neutralization.

    ``base`` is the salt-pair count implied by ``ionic_strength`` (mol/L) relative to
    water's molarity; the net charge is then neutralized by adding to the counter-ion.
    """
    base = nwat * ionic_strength // water_molarity
    if net_ion_charge >= 0:
        ncl = base
        nna = ncl + net_ion_charge
    else:
        nna = base
        ncl = nna - net_ion_charge
    return int(nna), int(ncl)


def _insert_ions(
    workdir: Path,
    solv_pdb: Path,
    cl_pdb: Path,
    na_pdb: Path,
    out_pdb: Path,
    ncl: int,
    nna: int,
) -> None:
    cl_out = workdir / "xtal_solv_cl.pdb"
    run_tool(
        [
            "gmx",
            "insert-molecules",
            "-f",
            str(solv_pdb),
            "-ci",
            str(cl_pdb),
            "-o",
            str(cl_out),
            "-replace",
            "SOL",
            "-nmol",
            str(ncl),
        ],
        tool="gmx insert-molecules",
        cwd=workdir,
        log_path=workdir / "insert_cl.log",
    )
    run_tool(
        [
            "gmx",
            "insert-molecules",
            "-f",
            str(cl_out),
            "-ci",
            str(na_pdb),
            "-o",
            str(out_pdb),
            "-replace",
            "SOL",
            "-nmol",
            str(nna),
        ],
        tool="gmx insert-molecules",
        cwd=workdir,
        log_path=workdir / "insert_na.log",
    )


def _finalize_topology(md_model_top: Path, ncopies: int, nwat: int, ncl: int, nna: int) -> None:
    with open(md_model_top, "a") as fh:
        fh.write("[ molecules ]\n")
        fh.write("; Compound       #mols\n")
        fh.write("system1              1\n" * ncopies)
        fh.write(f"WAT              {nwat}\n")
        fh.write(f"Cl-              {ncl}\n")
        fh.write(f"Na+              {nna}\n")
