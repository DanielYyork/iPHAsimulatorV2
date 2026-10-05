#!/usr/bin/env python
"""Staged OpenMM MD for an Amber system (system.prmtop / system.inpcrd).

Written by iPHASimulator (notebook 06C). The stages and file names mirror
the GROMACS run folders: step6.0_minimization, step6.1_nvt, step6.2_npt and
step7_production. All settings come from protocol.json next to this script.

    python run_openmm_md.py step6.0_minimization                        # local
    python run_openmm_md.py step6.1_nvt step6.2_npt step7_production    # HPC
    python run_openmm_md.py --short-test                                # quick check

Each stage starts from the previous stage's state (<stage>.xml) and writes
<stage>.xml, <stage>.pdb and <stage>.csv; NVT/NPT/production also write
<stage>.dcd frames. Every MD stage resumes from its own <stage>.chk and validates
its input fingerprint from <stage>.restart.json. Completed stages take no new steps.
Keep system.prmtop, system.inpcrd, this file and protocol.json together; OpenMM is
required, but an iPHASimulator installation is not.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

STAGES = ("step6.0_minimization", "step6.1_nvt", "step6.2_npt", "step7_production")
PREVIOUS_STAGE = dict(zip(STAGES[1:], STAGES[:-1]))
PROTEIN_RESIDUES = frozenset(
    "ALA ARG ASN ASP ASH CYS CYX CYM GLN GLU GLH GLY HIS HID HIE HIP ILE LEU LYS LYN "
    "MET PHE PRO SER THR TRP TYR VAL ACE NME".split()
)
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})
SHORT_TEST_PS = 10.0
SHORT_TEST_MINIMIZATION_ITERATIONS = 1000


def load_protocol(folder: str | Path = ".") -> dict:
    return json.loads((Path(folder) / "protocol.json").read_text())


def stage_steps(protocol: dict, stage: str, picoseconds: float | None = None) -> int:
    """Number of MD steps for a stage (0 for minimisation)."""

    if stage == "step6.0_minimization":
        return 0
    if stage == "step6.1_nvt":
        length, timestep = protocol["nvt_ps"], protocol["nvt_timestep_fs"]
    elif stage == "step6.2_npt":
        length, timestep = protocol["npt_ps"], protocol["npt_timestep_fs"]
    elif stage == "step7_production":
        length, timestep = protocol["production_ns"] * 1000, protocol["production_timestep_fs"]
    else:
        raise ValueError(f"Unknown stage: {stage}")
    if picoseconds is not None:
        length = picoseconds
    return round(length * 1000 / timestep)


def stage_intervals(protocol: dict, stage: str) -> tuple[float, float]:
    """CSV and trajectory intervals in ps, with defaults for older protocol files."""

    report = protocol.get({"step6.1_nvt": "nvt_report_ps", "step6.2_npt": "npt_report_ps"}.get(stage, "report_ps"))
    report = protocol["report_ps"] if report is None else report
    frames = protocol["production_frame_ps"] if stage == "step7_production" else protocol["equilibration_frame_ps"]
    if stage == "step6.2_npt" and protocol.get("npt_frame_ps") is not None:
        frames = protocol["npt_frame_ps"]
    return report, frames


def input_signature(folder: Path, protocol: dict) -> str:
    """Bind a restart to its topology, original coordinates and physical settings."""

    digest = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode())
    for name in ("system.prmtop", "system.inpcrd"):
        digest.update((folder / name).read_bytes())
    return digest.hexdigest()


def validate_inputs(prmtop, inpcrd, cutoff_nm: float) -> dict:
    """Check particle counts, finite coordinates and a PME-compatible periodic box."""

    from openmm import unit

    atoms = prmtop.topology.getNumAtoms()
    if len(inpcrd.positions) != atoms:
        raise ValueError(f"Topology/coordinate mismatch: {atoms} vs {len(inpcrd.positions)} particles")
    if not all(math.isfinite(float(x)) for p in inpcrd.positions.value_in_unit(unit.nanometer) for x in p):
        raise ValueError("Coordinates contain non-finite values")
    if inpcrd.boxVectors is None or prmtop.topology.getPeriodicBoxVectors() is None:
        raise ValueError("PME requires a periodic solvated topology and coordinate box")
    box = inpcrd.boxVectors.value_in_unit(unit.nanometer)
    if any(float(box[i][i]) <= 2 * cutoff_nm for i in range(3)):
        raise ValueError("Periodic box dimensions must exceed twice the nonbonded cutoff")
    return {"particles": atoms, "residues": prmtop.topology.getNumResidues()}


def restraint_groups(atoms, polymer_residues) -> tuple[list[int], list[int]]:
    """Split heavy atoms into backbone-type and side-chain restraint groups.

    ``atoms`` yields (index, residue_name, atom_name, element_symbol). As in the
    CHARMM-GUI packages: protein N/CA/C/O and polymer heavy atoms use the backbone
    force constant, other protein heavy atoms the side-chain one; water, ions and
    hydrogens are never restrained.
    """

    polymer = set(polymer_residues)
    backbone: list[int] = []
    sidechain: list[int] = []
    for index, residue, name, element in atoms:
        if element == "H":
            continue
        if residue in PROTEIN_RESIDUES:
            (backbone if name in BACKBONE_ATOMS else sidechain).append(index)
        elif residue in polymer:
            backbone.append(index)
    return backbone, sidechain


def run_stage(
    stage: str,
    folder: str | Path = ".",
    *,
    output_dir: str | Path | None = None,
    picoseconds: float | None = None,
    max_iterations: int | None = None,
    platform_name: str | None = None,
) -> dict:
    """Run one stage and return a summary (energy in kJ/mol, density in g/mL)."""

    import openmm
    from openmm import app, unit

    if stage not in STAGES:
        raise ValueError(f"Unknown stage: {stage}")
    folder = Path(folder)
    output = Path(output_dir) if output_dir else folder
    output.mkdir(parents=True, exist_ok=True)
    protocol = load_protocol(folder)
    temperature = protocol["temperature_kelvin"] * unit.kelvin
    checkpoint = output / f"{stage}.chk"
    provenance = output / f"{stage}.restart.json"
    signature = input_signature(folder, protocol)
    resumed = stage != "step6.0_minimization" and checkpoint.exists()
    if resumed and (not provenance.exists() or json.loads(provenance.read_text()).get("input_signature") != signature):
        raise ValueError(f"Restart inputs/settings changed or provenance missing: {checkpoint}")

    prmtop = app.AmberPrmtopFile(str(folder / "system.prmtop"))
    inpcrd = app.AmberInpcrdFile(str(folder / "system.inpcrd"))
    validation = validate_inputs(prmtop, inpcrd, protocol["cutoff_nm"])
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=protocol["cutoff_nm"] * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=protocol.get("ewald_error_tolerance", 5e-4),
    )

    restrained = stage != "step7_production" and (
        protocol["backbone_restraint"] > 0 or protocol["sidechain_restraint"] > 0
    )
    if restrained:
        atoms = (
            (atom.index, atom.residue.name, atom.name, atom.element.symbol if atom.element else "")
            for atom in prmtop.topology.atoms()
        )
        backbone, sidechain = restraint_groups(atoms, protocol["polymer_residues"])
        residue_names = {a.residue.name for a in prmtop.topology.atoms()}
        if not residue_names.intersection(protocol["polymer_residues"]):
            raise ValueError("polymer_residues does not match the topology; PHA restraints would be missing")
        force = openmm.CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
        for name in ("k", "x0", "y0", "z0"):
            force.addPerParticleParameter(name)
        reference = inpcrd.positions.value_in_unit(unit.nanometer)
        for indices, constant in ((backbone, protocol["backbone_restraint"]), (sidechain, protocol["sidechain_restraint"])):
            for index in indices:
                force.addParticle(index, [constant, *reference[index]])
        system.addForce(force)
    if stage in ("step6.2_npt", "step7_production"):
        system.addForce(
            openmm.MonteCarloBarostat(
                protocol["pressure_bar"] * unit.bar, temperature, protocol["barostat_frequency"]
            )
        )

    timestep = {
        "step6.1_nvt": protocol["nvt_timestep_fs"],
        "step6.2_npt": protocol["npt_timestep_fs"],
        "step7_production": protocol["production_timestep_fs"],
    }.get(stage, protocol["nvt_timestep_fs"])
    integrator = openmm.LangevinMiddleIntegrator(
        temperature, protocol["friction_per_ps"] / unit.picosecond, timestep * unit.femtosecond
    )
    integrator.setConstraintTolerance(protocol.get("constraint_tolerance", 1e-5))
    simulation = None
    platform_error = None
    for name in ([platform_name] if platform_name else ["CUDA", "OpenCL", "CPU"]):
        try:
            platform = openmm.Platform.getPlatformByName(name)
            simulation = app.Simulation(prmtop.topology, system, integrator, platform)
            break
        except openmm.OpenMMException as exc:
            platform_error = exc
            if platform_name:
                raise
    if simulation is None:
        raise RuntimeError("No usable OpenMM platform; choose an installed CPU/GPU platform") from platform_error

    if resumed:
        simulation.loadCheckpoint(str(checkpoint))
    elif stage == "step6.0_minimization":
        if (output / f"{stage}.xml").exists():
            raise FileExistsError("Minimisation output already exists; choose a new output directory")
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    else:
        if any((output / f"{stage}.{suffix}").exists() for suffix in ("xml", "csv", "dcd")):
            raise FileExistsError(f"Existing {stage} outputs without a checkpoint; choose a new directory")
        previous = output / f"{PREVIOUS_STAGE[stage]}.xml"
        if not previous.exists():
            raise FileNotFoundError(f"{previous} not found; run {PREVIOUS_STAGE[stage]} first")
        state = openmm.XmlSerializer.deserialize(previous.read_text())
        simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        simulation.context.setPositions(state.getPositions())
        if stage == "step6.1_nvt":
            simulation.context.setVelocitiesToTemperature(temperature)
        else:
            simulation.context.setVelocities(state.getVelocities())
    # Each stage has its own step/time origin; checkpoint continuation retains both.

    steps = stage_steps(protocol, stage, picoseconds)
    executed = 0
    if stage == "step6.0_minimization":
        simulation.minimizeEnergy(
            tolerance=protocol["minimization_tolerance"] * unit.kilojoule_per_mole / unit.nanometer,
            maxIterations=max_iterations or protocol["minimization_max_iterations"],
        )
    else:
        report_ps, frame_ps = stage_intervals(protocol, stage)
        report = max(1, round(report_ps * 1000 / timestep))
        remaining = steps - (simulation.currentStep if resumed else 0)
        csv_path = output / f"{stage}.csv"
        dcd_path = output / f"{stage}.dcd"
        if remaining > 0:
            simulation.reporters.append(
                app.StateDataReporter(
                    str(csv_path), report, step=True, time=True, potentialEnergy=True,
                    temperature=True, density=True, volume=True, speed=True,
                    append=resumed and csv_path.exists() and csv_path.stat().st_size > 0,
                )
            )
        if remaining > 0 and frame_ps > 0:
            frames = max(1, round(frame_ps * 1000 / timestep))
            simulation.reporters.append(app.DCDReporter(str(dcd_path), frames,
                                       append=resumed and dcd_path.exists() and dcd_path.stat().st_size > 0))
        provenance.write_text(json.dumps({"input_signature": signature, "platform": platform.getName(),
                                         "openmm_version": openmm.__version__}, indent=2) + "\n")
        simulation.reporters.append(app.CheckpointReporter(str(checkpoint), 10 * report))
        if remaining > 0:
            simulation.step(remaining)
            executed = remaining
        simulation.saveCheckpoint(str(checkpoint))

    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
    (output / f"{stage}.xml").write_text(openmm.XmlSerializer.serialize(state))
    with open(output / f"{stage}.pdb", "w") as handle:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle)
    mass = sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))
    volume = state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy) or not math.isfinite(volume) or volume <= 0:
        raise RuntimeError(f"Non-finite energy or invalid volume after {stage}")
    return {
        "stage": stage,
        "steps": steps,
        "steps_executed": executed,
        "current_step": simulation.currentStep,
        "resumed": resumed,
        **validation,
        "potential_energy_kj_mol": energy,
        "density_g_ml": mass * 1.66053906660e-3 / volume,
        "platform": simulation.context.getPlatform().getName(),
    }


def run_short_test(folder: str | Path = ".", *, platform_name: str | None = None) -> list[dict]:
    """Minimise, then 10 ps NVT and 10 ps NPT, into short_test/ (protocol settings otherwise)."""

    output = Path(folder) / "short_test"
    results = [run_stage("step6.0_minimization", folder, output_dir=output,
                         max_iterations=SHORT_TEST_MINIMIZATION_ITERATIONS, platform_name=platform_name)]
    for stage in ("step6.1_nvt", "step6.2_npt"):
        results.append(run_stage(stage, folder, output_dir=output, picoseconds=SHORT_TEST_PS, platform_name=platform_name))
    return results


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # argparse checks an empty nargs='*' list against choices, breaking --short-test.
    parser.add_argument("stages", nargs="*", metavar="STAGE", help="stages to run, in order: " + ", ".join(STAGES))
    parser.add_argument("--short-test", action="store_true", help="minimise + 10 ps NVT + 10 ps NPT in short_test/")
    parser.add_argument("--platform", default=None, help="CUDA, OpenCL or CPU (default: first available)")
    args = parser.parse_args(argv)
    if any(stage not in STAGES for stage in args.stages):
        parser.error("unknown stage; choose from " + ", ".join(STAGES))
    if not args.stages and not args.short_test:
        parser.error("give one or more stages, or --short-test")
    results = run_short_test(platform_name=args.platform) if args.short_test else [
        run_stage(stage, platform_name=args.platform) for stage in args.stages
    ]
    for result in results:
        print(
            f"{result['stage']}: {result['steps']} steps, Epot {result['potential_energy_kj_mol']:.1f} kJ/mol, "
            f"density {result['density_g_ml']:.3f} g/mL ({result['platform']})"
        )


if __name__ == "__main__":
    main()
