#!/usr/bin/env python
"""Staged OpenMM MD for an Amber system (system.prmtop / system.inpcrd).

Written by iPHASimulator (notebooks 06A1/06A2). The stages and file names mirror
the GROMACS run folders: step6.0_minimization, step6.1_nvt, step6.2_npt and
step7_production. All settings come from protocol.json next to this script.

    python run_openmm_md.py step6.0_minimization                        # local
    python run_openmm_md.py step6.1_nvt step6.2_npt step7_production    # HPC
    python run_openmm_md.py --short-test                                # quick check

Each stage starts from the previous stage's state (<stage>.xml) and writes
<stage>.xml, <stage>.pdb and <stage>.csv; NVT/NPT/production also write
<stage>.dcd frames. Production restarts from step7_production.chk when it exists.
Only this file and protocol.json are needed (no iPHASimulator install).
"""

from __future__ import annotations

import argparse
import json
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
    else:
        length, timestep = protocol["production_ns"] * 1000, protocol["production_timestep_fs"]
    if picoseconds is not None:
        length = picoseconds
    return round(length * 1000 / timestep)


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

    folder = Path(folder)
    output = Path(output_dir) if output_dir else folder
    output.mkdir(parents=True, exist_ok=True)
    protocol = load_protocol(folder)
    temperature = protocol["temperature_kelvin"] * unit.kelvin

    prmtop = app.AmberPrmtopFile(str(folder / "system.prmtop"))
    inpcrd = app.AmberInpcrdFile(str(folder / "system.inpcrd"))
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=protocol["cutoff_nm"] * unit.nanometer,
        constraints=app.HBonds,
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
    platform = None
    for name in ([platform_name] if platform_name else ["CUDA", "OpenCL", "CPU"]):
        try:
            platform = openmm.Platform.getPlatformByName(name)
            break
        except Exception:
            continue
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)

    checkpoint = output / "step7_production.chk"
    resumed = stage == "step7_production" and checkpoint.exists()
    if resumed:
        simulation.loadCheckpoint(str(checkpoint))
    elif stage == "step6.0_minimization":
        simulation.context.setPositions(inpcrd.positions)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    else:
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

    steps = stage_steps(protocol, stage, picoseconds)
    if stage == "step6.0_minimization":
        simulation.minimizeEnergy(
            tolerance=protocol["minimization_tolerance"] * unit.kilojoule_per_mole / unit.nanometer,
            maxIterations=max_iterations or protocol["minimization_max_iterations"],
        )
    else:
        report = max(1, round(protocol["report_ps"] * 1000 / timestep))
        frame_ps = protocol["production_frame_ps"] if stage == "step7_production" else protocol["equilibration_frame_ps"]
        simulation.reporters.append(
            app.StateDataReporter(
                str(output / f"{stage}.csv"), report, step=True, time=True, potentialEnergy=True,
                temperature=True, density=True, volume=True, speed=True, append=resumed,
            )
        )
        if frame_ps > 0:
            frames = max(1, round(frame_ps * 1000 / timestep))
            simulation.reporters.append(app.DCDReporter(str(output / f"{stage}.dcd"), frames, append=resumed))
        if stage == "step7_production":
            simulation.reporters.append(app.CheckpointReporter(str(checkpoint), max(report, 10 * report)))
        remaining = steps - (simulation.currentStep if resumed else 0)
        if remaining > 0:
            simulation.step(remaining)
        if stage == "step7_production":
            simulation.saveCheckpoint(str(checkpoint))

    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
    (output / f"{stage}.xml").write_text(openmm.XmlSerializer.serialize(state))
    with open(output / f"{stage}.pdb", "w") as handle:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle)
    mass = sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))
    volume = state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
    return {
        "stage": stage,
        "steps": steps,
        "potential_energy_kj_mol": state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
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
    parser.add_argument("stages", nargs="*", choices=STAGES, help="stages to run, in order")
    parser.add_argument("--short-test", action="store_true", help="minimise + 10 ps NVT + 10 ps NPT in short_test/")
    parser.add_argument("--platform", default=None, help="CUDA, OpenCL or CPU (default: first available)")
    args = parser.parse_args(argv)
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
