"""Amber/OpenMM solvated systems: tleap (GAFF2 + OPC, optional ff19SB) and staged OpenMM MD.

The stage lengths, temperatures and restraints follow the GROMACS references:
POLYMER_IN_WATER_PROTOCOL the polymer benchmark (notebook 06A) and
ENZYME_POLYMER_IN_WATER_PROTOCOL the enzyme–polymer production runs (notebook 06B).
The force fields, thermostat, barostat and non-bonded settings differ (Amber/OpenMM),
so 06C results are not directly comparable with 06B. ``tleap`` and OpenMM run only
when a ``build_``/``run_`` function or the generated ``run_openmm_md.py`` is called;
everything else writes text files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from dataclasses import asdict
import importlib.util
import json
import re
import shutil
import subprocess


WATER_MOLARITY = 55.5
SYSTEM_PRMTOP = "system.prmtop"
SYSTEM_INPCRD = "system.inpcrd"
SYSTEM_PDB = "system.pdb"
TLEAP_INPUT = "tleap.in"
TLEAP_LOG = "tleap.log"
COUNT_INPUT = "tleap_count_waters.in"
COUNT_LOG = "tleap_count_waters.log"
RUN_SCRIPT = "run_openmm_md.py"
PROTOCOL_FILE = "protocol.json"
LOCAL_SCRIPT = "run_step6_local.sh"
HPC_SCRIPT = "run_hpc_equilibration_production.slurm"
RUN_SCRIPT_SOURCE = Path(__file__).resolve().parent / "data" / "amber_openmm" / RUN_SCRIPT

ADDED_RESIDUES_PATTERN = re.compile(r"Added\s+(\d+)\s+residues")
ION_PARAMETERS_PATTERN = re.compile(r"Loading parameters:\s*(\S*frcmod\.ion\S*)")
LEAP_ERRORS_PATTERN = re.compile(r"Exiting LEaP:\s*Errors\s*=\s*(\d+)")


@dataclass(frozen=True)
class AmberSolvationSettings:
    """Box and ion choices for the tleap build (lengths in nm, salt in mol/L)."""

    padding_nm: float = 1.2
    box_shape: str = "box"  # "box" (solvateBox) or "oct" (solvateOct)
    cubic: bool = True  # solvateBox ... iso: cubic box, like gmx editconf -bt cubic
    salt_molar: float = 0.15
    cation: str = "Na+"
    anion: str = "Cl-"

    def __post_init__(self) -> None:
        if self.box_shape not in ("box", "oct"):
            raise ValueError("box_shape must be 'box' or 'oct'")
        if self.padding_nm <= 0 or self.salt_molar < 0:
            raise ValueError("padding_nm must be > 0 and salt_molar >= 0")


@dataclass(frozen=True)
class AmberSystemOutputs:
    output_dir: Path
    tleap_input_path: Path
    tleap_log_path: Path
    prmtop_path: Path
    inpcrd_path: Path
    pdb_path: Path
    water_count: int
    salt_ion_pairs: int
    ion_parameter_files: tuple[str, ...]


@dataclass(frozen=True)
class TleapLogSummary:
    water_count: int | None
    ion_parameter_files: tuple[str, ...]
    errors: int | None


@dataclass(frozen=True)
class OpenMMProtocol:
    """Staged MD settings (times in ps/ns, time steps in fs, restraints in kJ mol^-1 nm^-2)."""

    name: str
    temperature_kelvin: float
    minimization_max_iterations: int
    nvt_ps: float
    nvt_timestep_fs: float
    npt_ps: float
    npt_timestep_fs: float
    production_ns: float
    production_frame_ps: float
    production_timestep_fs: float = 2.0
    equilibration_frame_ps: float = 0.0
    report_ps: float = 10.0
    pressure_bar: float = 1.0
    cutoff_nm: float = 1.0
    friction_per_ps: float = 1.0
    barostat_frequency: int = 25
    minimization_tolerance: float = 1000.0
    backbone_restraint: float = 0.0
    sidechain_restraint: float = 0.0


# Polymer benchmark (06A): 300 K, no restraints, 100 ps NVT, 500 ps NPT, 100 ns, frames every 2 ps.
POLYMER_IN_WATER_PROTOCOL = OpenMMProtocol(
    name="polymer_in_water",
    temperature_kelvin=300.0,
    minimization_max_iterations=50000,
    nvt_ps=100.0,
    nvt_timestep_fs=2.0,
    npt_ps=500.0,
    npt_timestep_fs=2.0,
    production_ns=100.0,
    production_frame_ps=2.0,
)

# Enzyme–polymer production runs / 06B: 303.15 K, restraints (N/CA/C/O and polymer
# heavy atoms 400, other protein heavy atoms 40) during minimisation, NVT and NPT,
# 125 ps NVT at 1 fs, 500 ps NPT, 200 ns, frames every 100 ps.
ENZYME_POLYMER_IN_WATER_PROTOCOL = OpenMMProtocol(
    name="enzyme_polymer_in_water",
    temperature_kelvin=303.15,
    minimization_max_iterations=5000,
    nvt_ps=125.0,
    nvt_timestep_fs=1.0,
    npt_ps=500.0,
    npt_timestep_fs=2.0,
    production_ns=200.0,
    production_frame_ps=100.0,
    equilibration_frame_ps=5.0,
    backbone_restraint=400.0,
    sidechain_restraint=40.0,
)


def salt_ion_pairs(water_count: int, salt_molar: float) -> int:
    """Ion pairs for ``salt_molar`` from the number of waters (55.5 mol/L water)."""

    return round(salt_molar * water_count / WATER_MOLARITY)


def density_g_per_ml(total_mass_dalton: float, volume_nm3: float) -> float:
    """Convert a system mass in Da and a box volume in nm^3 to g/mL (1 Da/nm^3 = 1.66054e-3 g/mL)."""

    return total_mass_dalton * 1.66053906660e-3 / volume_nm3


def write_tleap_input(
    path: str | Path,
    *,
    mol2_name: str,
    frcmod_name: str,
    settings: AmberSolvationSettings,
    protein_pdb_name: str | None = None,
    salt_pairs: int | None = None,
) -> str:
    """Write a tleap input and return its text.

    With ``salt_pairs=None`` the input only solvates (used to count waters);
    otherwise it neutralises, adds ``salt_pairs`` ion pairs and saves the system.
    """

    lines = []
    if protein_pdb_name:
        lines.append("source leaprc.protein.ff19SB")
    lines += [
        "source leaprc.gaff2",
        "source leaprc.water.opc",
        f"loadamberparams {frcmod_name}",
        f"POL = loadmol2 {mol2_name}",
    ]
    if protein_pdb_name:
        lines += [f"PROT = loadpdb {protein_pdb_name}", "SYS = combine { PROT POL }"]
    else:
        lines.append("SYS = copy POL")
    command = "solvateBox" if settings.box_shape == "box" else "solvateOct"
    iso = " iso" if settings.box_shape == "box" and settings.cubic else ""
    lines.append(f"{command} SYS OPCBOX {settings.padding_nm * 10:.1f}{iso}")
    if salt_pairs is not None:
        lines += [
            f"addIonsRand SYS {settings.cation} 0",
            f"addIonsRand SYS {settings.anion} 0",
        ]
        if salt_pairs > 0:
            lines.append(
                f"addIonsRand SYS {settings.cation} {salt_pairs} {settings.anion} {salt_pairs}"
            )
        lines += [
            "charge SYS",
            "check SYS",
            f"saveamberparm SYS {SYSTEM_PRMTOP} {SYSTEM_INPCRD}",
            f"savepdb SYS {SYSTEM_PDB}",
        ]
    lines += ["quit", ""]
    text = "\n".join(lines)
    Path(path).write_text(text)
    return text


def parse_tleap_log(text: str) -> TleapLogSummary:
    """Waters added by the (last) solvate command, ion frcmod files loaded, error count."""

    added = ADDED_RESIDUES_PATTERN.findall(text)
    errors = LEAP_ERRORS_PATTERN.findall(text)
    ion_files = tuple(dict.fromkeys(Path(name).name for name in ION_PARAMETERS_PATTERN.findall(text)))
    return TleapLogSummary(
        water_count=int(added[-1]) if added else None,
        ion_parameter_files=ion_files,
        errors=int(errors[-1]) if errors else None,
    )


def build_solvated_amber_system(
    mol2_path: str | Path,
    frcmod_path: str | Path,
    output_dir: str | Path,
    *,
    protein_pdb_path: str | Path | None = None,
    settings: AmberSolvationSettings = AmberSolvationSettings(),
    tleap: str = "tleap",
    runner=subprocess.run,
) -> AmberSystemOutputs:
    """Solvate the 05A polymer (optionally with a posed protein) in OPC and add ions.

    Copies the inputs into a new ``output_dir``, runs tleap once to count the
    waters for the requested padding, then builds the neutralised, salted system
    and writes ``system.prmtop`` / ``system.inpcrd`` / ``system.pdb`` with the
    tleap input and log next to them.
    """

    mol2 = Path(mol2_path).expanduser().resolve()
    frcmod = Path(frcmod_path).expanduser().resolve()
    protein = Path(protein_pdb_path).expanduser().resolve() if protein_pdb_path else None
    for required in (mol2, frcmod, *( [protein] if protein else [])):
        if not required.is_file():
            raise FileNotFoundError(f"Input not found: {required}")
    output = Path(output_dir).expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Output folder already exists; choose a new one: {output}")
    output.mkdir(parents=True)
    for source in (mol2, frcmod, *([protein] if protein else [])):
        shutil.copyfile(source, output / source.name)
    names = dict(
        mol2_name=mol2.name,
        frcmod_name=frcmod.name,
        settings=settings,
        protein_pdb_name=protein.name if protein else None,
    )

    write_tleap_input(output / COUNT_INPUT, **names)
    count = parse_tleap_log(_run_tleap(output, COUNT_INPUT, COUNT_LOG, tleap, runner))
    if count.water_count is None:
        raise RuntimeError(f"Could not read the number of added waters from {output / COUNT_LOG}")
    pairs = salt_ion_pairs(count.water_count, settings.salt_molar)

    write_tleap_input(output / TLEAP_INPUT, salt_pairs=pairs, **names)
    summary = parse_tleap_log(_run_tleap(output, TLEAP_INPUT, TLEAP_LOG, tleap, runner))
    if summary.errors:
        raise RuntimeError(f"tleap reported {summary.errors} error(s); see {output / TLEAP_LOG}")
    for name in (SYSTEM_PRMTOP, SYSTEM_INPCRD):
        if not (output / name).is_file():
            raise RuntimeError(f"tleap did not write {name}; see {output / TLEAP_LOG}")
    return AmberSystemOutputs(
        output_dir=output,
        tleap_input_path=output / TLEAP_INPUT,
        tleap_log_path=output / TLEAP_LOG,
        prmtop_path=output / SYSTEM_PRMTOP,
        inpcrd_path=output / SYSTEM_INPCRD,
        pdb_path=output / SYSTEM_PDB,
        water_count=count.water_count,
        salt_ion_pairs=pairs,
        ion_parameter_files=summary.ion_parameter_files,
    )


def write_openmm_run_files(
    system_dir: str | Path,
    protocol: OpenMMProtocol,
    *,
    job_name: str,
    polymer_residues: tuple[str, ...] = ("PHA",),
) -> tuple[Path, ...]:
    """Write ``run_openmm_md.py``, ``protocol.json``, ``run_step6_local.sh`` and the SLURM file.

    The layout mirrors the GROMACS run folders: step 6.0 locally, then steps 6.1,
    6.2 and 7 on the cluster. ``polymer_residues`` names the polymer residue(s)
    restrained with the backbone force constant.
    """

    folder = Path(system_dir)
    for name in (SYSTEM_PRMTOP, SYSTEM_INPCRD):
        if not (folder / name).is_file():
            raise FileNotFoundError(f"{name} not found in {folder}; build the system first")
    script = folder / RUN_SCRIPT
    shutil.copyfile(RUN_SCRIPT_SOURCE, script)
    protocol_path = folder / PROTOCOL_FILE
    protocol_path.write_text(
        json.dumps({**asdict(protocol), "polymer_residues": list(polymer_residues)}, indent=2) + "\n"
    )
    local = folder / LOCAL_SCRIPT
    local.write_text(LOCAL_SCRIPT_TEXT)
    local.chmod(0o755)
    hpc = folder / HPC_SCRIPT
    hpc.write_text(HPC_SCRIPT_TEMPLATE.replace("{JOB_NAME}", job_name))
    hpc.chmod(0o755)
    return script, protocol_path, local, hpc


def load_openmm_runner():
    """Import the packaged ``run_openmm_md.py`` (the same file the run folders use)."""

    spec = importlib.util.spec_from_file_location("iphasimulator_run_openmm_md", RUN_SCRIPT_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_openmm_stages(system_dir: str | Path, stages, *, platform_name: str | None = None) -> list[dict]:
    """Run stages (e.g. ``["step6.0_minimization"]``) in ``system_dir`` with its protocol.json."""

    runner = load_openmm_runner()
    return [runner.run_stage(stage, system_dir, platform_name=platform_name) for stage in stages]


def run_openmm_short_test(system_dir: str | Path, *, platform_name: str | None = None) -> list[dict]:
    """Minimise, then 10 ps NVT and 10 ps NPT in ``system_dir/short_test`` (protocol otherwise)."""

    return load_openmm_runner().run_short_test(system_dir, platform_name=platform_name)


def _run_tleap(folder: Path, input_name: str, log_name: str, tleap: str, runner) -> str:
    result = runner(
        [tleap, "-f", input_name],
        cwd=folder,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output = result.stdout or ""
    (folder / log_name).write_text(output)
    if result.returncode != 0:
        raise RuntimeError(f"tleap failed with return code {result.returncode}; see {folder / log_name}")
    return output


LOCAL_SCRIPT_TEXT = """#!/usr/bin/env bash
set -euo pipefail

python run_openmm_md.py step6.0_minimization
"""

HPC_SCRIPT_TEMPLATE = """#!/bin/bash -l
#SBATCH --job-name={JOB_NAME}
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=2-00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Activate a Python environment with OpenMM (CUDA build), for example:
# module load anaconda3 && conda activate ipha_clean

mkdir -p logs

# Steps 6.1 (NVT), 6.2 (NPT) and 7 (production); settings in protocol.json.
# Production continues from step7_production.chk when the job is resubmitted.
python run_openmm_md.py step6.1_nvt step6.2_npt step7_production --platform CUDA
"""
