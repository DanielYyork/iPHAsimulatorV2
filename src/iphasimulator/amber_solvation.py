"""Amber/OpenMM solvated systems: tleap (GAFF2 + OPC, optional ff19SB) and short OpenMM runs.

``tleap`` and OpenMM run only when :func:`build_solvated_amber_system`,
:func:`run_openmm_short_test` or a generated production script is called.
Everything else here writes text files into a new output folder.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
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
PRODUCTION_SCRIPT = "run_openmm_production.py"
PRODUCTION_SLURM = "run_openmm_production.slurm"

ADDED_RESIDUES_PATTERN = re.compile(r"Added\s+(\d+)\s+residues")
ION_PARAMETERS_PATTERN = re.compile(r"Loading parameters:\s*(\S*frcmod\.ion\S*)")
LEAP_ERRORS_PATTERN = re.compile(r"Exiting LEaP:\s*Errors\s*=\s*(\d+)")


@dataclass(frozen=True)
class AmberSolvationSettings:
    """Box and ion choices for the tleap build (lengths in nm, salt in mol/L)."""

    padding_nm: float = 1.2
    box_shape: str = "box"  # "box" (solvateBox) or "oct" (solvateOct)
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
class OpenMMTestSettings:
    """PME / HBonds / LangevinMiddle / Monte Carlo barostat settings for the Amber/OpenMM route."""

    cutoff_nm: float = 1.0
    temperature_kelvin: float = 300.0
    friction_per_ps: float = 1.0
    timestep_fs: float = 2.0
    pressure_bar: float = 1.0
    nvt_ps: float = 10.0
    npt_ps: float = 10.0
    report_interval_steps: int = 500
    platform_name: str | None = None

    def steps(self, picoseconds: float) -> int:
        return round(picoseconds * 1000 / self.timestep_fs)


@dataclass(frozen=True)
class OpenMMTestResult:
    output_dir: Path
    minimized_energy_kj_mol: float
    nvt_energy_kj_mol: float
    npt_energy_kj_mol: float
    density_g_ml: float
    report_path: Path


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
    lines.append(f"{command} SYS OPCBOX {settings.padding_nm * 10:.1f}")
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


def run_openmm_short_test(
    prmtop_path: str | Path,
    inpcrd_path: str | Path,
    output_dir: str | Path,
    *,
    settings: OpenMMTestSettings = OpenMMTestSettings(),
) -> OpenMMTestResult:
    """Minimise, then run short NVT and NPT stages and report energies and density."""

    import openmm
    from openmm import app, unit

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    prmtop = app.AmberPrmtopFile(str(prmtop_path))
    inpcrd = app.AmberInpcrdFile(str(inpcrd_path))
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=settings.cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
    )
    barostat = openmm.MonteCarloBarostat(
        settings.pressure_bar * unit.bar, settings.temperature_kelvin * unit.kelvin
    )
    barostat_index = system.addForce(barostat)
    barostat.setFrequency(0)  # off during NVT
    integrator = openmm.LangevinMiddleIntegrator(
        settings.temperature_kelvin * unit.kelvin,
        settings.friction_per_ps / unit.picosecond,
        settings.timestep_fs * unit.femtosecond,
    )
    platform = openmm.Platform.getPlatformByName(settings.platform_name) if settings.platform_name else None
    simulation = app.Simulation(prmtop.topology, system, integrator, *([platform] if platform else []))
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    def energy() -> float:
        state = simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    simulation.minimizeEnergy()
    minimized = energy()
    report_path = output / "short_test.csv"
    simulation.reporters.append(
        app.StateDataReporter(
            str(report_path),
            settings.report_interval_steps,
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            volume=True,
        )
    )
    simulation.context.setVelocitiesToTemperature(settings.temperature_kelvin * unit.kelvin)
    simulation.step(settings.steps(settings.nvt_ps))
    nvt = energy()
    system.getForce(barostat_index).setFrequency(25)
    simulation.context.reinitialize(preserveState=True)
    simulation.step(settings.steps(settings.npt_ps))
    npt = energy()

    state = simulation.context.getState(getPositions=True)
    volume = state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
    mass = sum(
        system.getParticleMass(index).value_in_unit(unit.dalton)
        for index in range(system.getNumParticles())
    )
    with open(output / "short_test_final.pdb", "w") as handle:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle)
    with open(output / "short_test_state.xml", "w") as handle:
        handle.write(
            openmm.XmlSerializer.serialize(
                simulation.context.getState(getPositions=True, getVelocities=True)
            )
        )
    return OpenMMTestResult(
        output_dir=output,
        minimized_energy_kj_mol=minimized,
        nvt_energy_kj_mol=nvt,
        npt_energy_kj_mol=npt,
        density_g_ml=density_g_per_ml(mass, volume),
        report_path=report_path,
    )


def write_openmm_production_files(
    system_dir: str | Path,
    *,
    production_ns: float = 100.0,
    job_name: str = "amber_openmm",
    settings: OpenMMTestSettings = OpenMMTestSettings(),
    report_interval_ps: float = 100.0,
) -> tuple[Path, Path]:
    """Write a restartable OpenMM production script and a SLURM file into ``system_dir``."""

    folder = Path(system_dir)
    for name in (SYSTEM_PRMTOP, SYSTEM_INPCRD):
        if not (folder / name).is_file():
            raise FileNotFoundError(f"{name} not found in {folder}; build the system first")
    report_steps = settings.steps(report_interval_ps)
    script = PRODUCTION_SCRIPT_TEMPLATE.format(
        prmtop=SYSTEM_PRMTOP,
        inpcrd=SYSTEM_INPCRD,
        cutoff_nm=settings.cutoff_nm,
        temperature=settings.temperature_kelvin,
        friction=settings.friction_per_ps,
        timestep_fs=settings.timestep_fs,
        pressure=settings.pressure_bar,
        total_steps=settings.steps(production_ns * 1000),
        report_steps=report_steps,
    )
    script_path = folder / PRODUCTION_SCRIPT
    script_path.write_text(script)
    slurm_path = folder / PRODUCTION_SLURM
    slurm_path.write_text(PRODUCTION_SLURM_TEMPLATE.replace("{JOB_NAME}", job_name))
    return script_path, slurm_path


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


PRODUCTION_SCRIPT_TEMPLATE = '''#!/usr/bin/env python
"""Amber/OpenMM production run written by iPHASimulator notebook 06B.

Restarts from production.chk when it exists. Settings: PME {cutoff_nm} nm,
HBonds constraints, LangevinMiddle {temperature} K, {timestep_fs} fs,
Monte Carlo barostat {pressure} bar.
"""

from pathlib import Path

import openmm
from openmm import app, unit

TOTAL_STEPS = {total_steps}
REPORT_STEPS = {report_steps}
CHECKPOINT = Path("production.chk")

prmtop = app.AmberPrmtopFile("{prmtop}")
inpcrd = app.AmberInpcrdFile("{inpcrd}")
system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff={cutoff_nm} * unit.nanometer,
    constraints=app.HBonds,
)
system.addForce(openmm.MonteCarloBarostat({pressure} * unit.bar, {temperature} * unit.kelvin))
integrator = openmm.LangevinMiddleIntegrator(
    {temperature} * unit.kelvin, {friction} / unit.picosecond, {timestep_fs} * unit.femtosecond
)

platform = None
for name in ("CUDA", "OpenCL", "CPU"):
    try:
        platform = openmm.Platform.getPlatformByName(name)
        break
    except Exception:
        continue
print("Platform:", platform.getName())
simulation = app.Simulation(prmtop.topology, system, integrator, platform)

if CHECKPOINT.exists():
    simulation.loadCheckpoint(str(CHECKPOINT))
    print("Restarted at step", simulation.currentStep)
else:
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature({temperature} * unit.kelvin)

append = CHECKPOINT.exists()
simulation.reporters += [
    app.DCDReporter("production.dcd", REPORT_STEPS, append=append),
    app.StateDataReporter(
        "production.csv", REPORT_STEPS, step=True, time=True, potentialEnergy=True,
        temperature=True, density=True, speed=True, append=append,
    ),
    app.CheckpointReporter(str(CHECKPOINT), REPORT_STEPS),
]
remaining = TOTAL_STEPS - simulation.currentStep
if remaining > 0:
    simulation.step(remaining)
simulation.saveCheckpoint(str(CHECKPOINT))
with open("production_final.pdb", "w") as handle:
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle)
print("Finished at step", simulation.currentStep)
'''

PRODUCTION_SLURM_TEMPLATE = """#!/bin/bash -l
#SBATCH --job-name={JOB_NAME}
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=2-00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Activate a Python environment with OpenMM (CUDA build) here, for example:
# module load anaconda3 && conda activate ipha_clean

mkdir -p logs
python run_openmm_production.py
"""
