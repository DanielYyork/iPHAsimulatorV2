"""Amber/OpenMM tleap helpers; tleap is mocked and OpenMM is never run."""

import ast
import subprocess

import pytest

import json

from iphasimulator.amber_solvation import (
    ENZYME_POLYMER_IN_WATER_PROTOCOL,
    POLYMER_IN_WATER_PROTOCOL,
    AmberSolvationSettings,
    build_solvated_amber_system,
    density_g_per_ml,
    load_openmm_runner,
    parse_tleap_log,
    salt_ion_pairs,
    write_openmm_run_files,
    write_tleap_input,
)


COUNT_LOG = """Loading parameters: /amber/dat/leap/parm/frcmod.opc
Loading parameters: /amber/dat/leap/parm/frcmod.ionslm_126_opc
  Solute vdw bounding box:              12.1 9.8 7.5
  Added 3700 residues.
Exiting LEaP: Errors = 0; Warnings = 0; Notes = 0.
"""
BUILD_LOG = COUNT_LOG + "Adding 11 counter ions to \"SYS\".\n"


@pytest.fixture
def inputs(tmp_path):
    mol2 = tmp_path / "P3HB_4.gaff2.mol2"
    frcmod = tmp_path / "P3HB_4.gaff2.frcmod"
    mol2.write_text("@<TRIPOS>MOLECULE\nPHA\n")
    frcmod.write_text("MASS\n")
    return mol2, frcmod


def fake_tleap(calls, log=BUILD_LOG, errors=0):
    def runner(command, cwd, **kwargs):
        calls.append((command, (cwd / command[2]).read_text()))
        if command[2] == "tleap.in":
            (cwd / "system.prmtop").write_text("%VERSION\n")
            (cwd / "system.inpcrd").write_text("PHA\n")
            (cwd / "system.pdb").write_text("END\n")
            text = log.replace("Errors = 0", f"Errors = {errors}")
        else:
            text = COUNT_LOG
        return subprocess.CompletedProcess(command, 0, text, None)

    return runner


def test_salt_ion_pairs_and_density():
    assert salt_ion_pairs(3700, 0.15) == 10
    assert salt_ion_pairs(3700, 0.0) == 0
    assert density_g_per_ml(602.214076, 1.0) == pytest.approx(1.0, rel=1e-6)


def test_write_tleap_input_for_polymer_in_opc(tmp_path):
    text = write_tleap_input(
        tmp_path / "tleap.in",
        mol2_name="p.mol2",
        frcmod_name="p.frcmod",
        settings=AmberSolvationSettings(padding_nm=1.2, box_shape="oct"),
        salt_pairs=10,
    )

    lines = text.splitlines()
    assert lines[:2] == ["source leaprc.gaff2", "source leaprc.water.opc"]
    assert "leaprc.protein.ff19SB" not in text
    assert "solvateOct SYS OPCBOX 12.0" in lines
    assert lines.index("addIonsRand SYS Na+ 0") < lines.index("addIonsRand SYS Na+ 10 Cl- 10")
    assert "saveamberparm SYS system.prmtop system.inpcrd" in lines


def test_write_tleap_input_with_protein_and_count_only(tmp_path):
    text = write_tleap_input(
        tmp_path / "count.in",
        mol2_name="p.mol2",
        frcmod_name="p.frcmod",
        settings=AmberSolvationSettings(),
        protein_pdb_name="enzyme.pdb",
    )
    cubic = write_tleap_input(tmp_path / "c.in", mol2_name="p.mol2", frcmod_name="p.frcmod",
                              settings=AmberSolvationSettings(cubic=False))

    assert text.splitlines()[0] == "source leaprc.protein.ff19SB"
    assert "PROT = loadpdb enzyme.pdb" in text
    assert "SYS = combine { PROT POL }" in text
    assert "solvateBox SYS OPCBOX 12.0 iso" in text  # cubic by default, like editconf -bt cubic
    assert "solvateBox SYS OPCBOX 12.0\n" in cubic
    assert "addIonsRand" not in text and "saveamberparm" not in text


def test_settings_reject_unknown_box_shape():
    with pytest.raises(ValueError, match="box_shape"):
        AmberSolvationSettings(box_shape="sphere")


def test_parse_tleap_log_reports_waters_ion_files_and_errors():
    summary = parse_tleap_log(COUNT_LOG)

    assert summary.water_count == 3700
    assert summary.ion_parameter_files == ("frcmod.ionslm_126_opc",)
    assert summary.errors == 0


def test_build_solvated_amber_system_counts_waters_then_builds(inputs, tmp_path):
    mol2, frcmod = inputs
    calls = []

    outputs = build_solvated_amber_system(
        mol2, frcmod, tmp_path / "solvated", runner=fake_tleap(calls)
    )

    assert [call[0] for call in calls] == [
        ["tleap", "-f", "tleap_count_waters.in"],
        ["tleap", "-f", "tleap.in"],
    ]
    assert "addIonsRand SYS Na+ 10 Cl- 10" in calls[1][1]
    assert outputs.water_count == 3700
    assert outputs.salt_ion_pairs == 10
    assert outputs.ion_parameter_files == ("frcmod.ionslm_126_opc",)
    assert outputs.prmtop_path.is_file()
    assert (outputs.output_dir / mol2.name).is_file()
    assert outputs.tleap_log_path.read_text() == BUILD_LOG


def test_build_solvated_amber_system_raises_on_tleap_errors(inputs, tmp_path):
    mol2, frcmod = inputs

    with pytest.raises(RuntimeError, match="2 error"):
        build_solvated_amber_system(
            mol2, frcmod, tmp_path / "solvated", runner=fake_tleap([], errors=2)
        )


def test_build_solvated_amber_system_refuses_existing_output(inputs, tmp_path):
    mol2, frcmod = inputs
    (tmp_path / "solvated").mkdir()

    with pytest.raises(FileExistsError):
        build_solvated_amber_system(mol2, frcmod, tmp_path / "solvated", runner=fake_tleap([]))


def test_protocols_mirror_the_gromacs_references():
    polymer, enzyme = POLYMER_IN_WATER_PROTOCOL, ENZYME_POLYMER_IN_WATER_PROTOCOL
    # polymer benchmark (gromacs_mdp / charmm_gromacs_polymer): 300 K, 100 ps NVT, 500 ps NPT, 100 ns, 2 ps frames
    assert (polymer.temperature_kelvin, polymer.nvt_ps, polymer.nvt_timestep_fs) == (300.0, 100.0, 2.0)
    assert (polymer.npt_ps, polymer.production_ns, polymer.production_frame_ps) == (500.0, 100.0, 2.0)
    assert (polymer.backbone_restraint, polymer.sidechain_restraint, polymer.minimization_max_iterations) == (0.0, 0.0, 50000)
    # enzyme runs (charmm_gromacs): 303.15 K, 125 ps NVT at 1 fs, 500 ps NPT, 200 ns, 100 ps frames, restraints 400/40
    assert (enzyme.temperature_kelvin, enzyme.nvt_ps, enzyme.nvt_timestep_fs) == (303.15, 125.0, 1.0)
    assert (enzyme.npt_ps, enzyme.production_ns, enzyme.production_frame_ps) == (500.0, 200.0, 100.0)
    assert (enzyme.backbone_restraint, enzyme.sidechain_restraint, enzyme.minimization_max_iterations) == (400.0, 40.0, 5000)
    for protocol in (polymer, enzyme):
        assert (protocol.minimization_tolerance, protocol.pressure_bar, protocol.cutoff_nm) == (1000.0, 1.0, 1.0)


def test_runner_stage_steps_match_protocol_lengths():
    runner = load_openmm_runner()
    protocol = {**ENZYME_POLYMER_IN_WATER_PROTOCOL.__dict__}

    assert runner.stage_steps(protocol, "step6.0_minimization") == 0
    assert runner.stage_steps(protocol, "step6.1_nvt") == 125000  # GROMACS step6.1_nvt nsteps
    assert runner.stage_steps(protocol, "step6.2_npt") == 250000
    assert runner.stage_steps(protocol, "step7_production") == 100000000
    assert runner.stage_steps(protocol, "step6.1_nvt", picoseconds=10) == 10000


def test_runner_restraint_groups_follow_charmm_gui_selection():
    runner = load_openmm_runner()
    atoms = [
        (0, "LEU", "N", "N"), (1, "LEU", "H", "H"), (2, "LEU", "CA", "C"), (3, "LEU", "CB", "C"),
        (4, "HID", "ND1", "N"), (5, "LEU", "C", "C"), (6, "LEU", "O", "O"),
        (7, "PHA", "C1", "C"), (8, "PHA", "H1", "H"), (9, "WAT", "O", "O"), (10, "Na+", "Na+", "Na"),
    ]

    backbone, sidechain = runner.restraint_groups(atoms, ["PHA"])

    assert backbone == [0, 2, 5, 6, 7]
    assert sidechain == [3, 4]


def test_write_openmm_run_files(tmp_path):
    (tmp_path / "system.prmtop").write_text("")
    (tmp_path / "system.inpcrd").write_text("")

    script, protocol_path, local, hpc = write_openmm_run_files(
        tmp_path, ENZYME_POLYMER_IN_WATER_PROTOCOL, job_name="GK13_PHO4_amber", polymer_residues=("PHA",)
    )

    ast.parse(script.read_text())
    protocol = json.loads(protocol_path.read_text())
    assert protocol["name"] == "enzyme_polymer_in_water" and protocol["polymer_residues"] == ["PHA"]
    assert protocol["backbone_restraint"] == 400.0
    assert "python run_openmm_md.py step6.0_minimization" in local.read_text()
    slurm = hpc.read_text()
    assert "#SBATCH --job-name=GK13_PHO4_amber" in slurm
    assert "python run_openmm_md.py step6.1_nvt step6.2_npt step7_production" in slurm
    assert load_openmm_runner().load_protocol(tmp_path)["nvt_timestep_fs"] == 1.0


def test_write_openmm_run_files_needs_system(tmp_path):
    with pytest.raises(FileNotFoundError, match="system.prmtop"):
        write_openmm_run_files(tmp_path, POLYMER_IN_WATER_PROTOCOL, job_name="x")
