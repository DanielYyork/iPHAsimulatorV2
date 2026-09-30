"""Route A tleap/OpenMM helpers; tleap is mocked and OpenMM is never run."""

import ast
import subprocess

import pytest

from iphasimulator.amber_solvation import (
    AmberSolvationSettings,
    OpenMMTestSettings,
    build_solvated_amber_system,
    density_g_per_ml,
    parse_tleap_log,
    salt_ion_pairs,
    write_openmm_production_files,
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

    assert text.splitlines()[0] == "source leaprc.protein.ff19SB"
    assert "PROT = loadpdb enzyme.pdb" in text
    assert "SYS = combine { PROT POL }" in text
    assert "solvateBox SYS OPCBOX 12.0" in text
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


def test_openmm_settings_convert_picoseconds_to_steps():
    settings = OpenMMTestSettings()

    assert settings.steps(10) == 5000
    assert settings.timestep_fs == 2.0 and settings.cutoff_nm == 1.0


def test_write_openmm_production_files(tmp_path):
    (tmp_path / "system.prmtop").write_text("")
    (tmp_path / "system.inpcrd").write_text("")

    script, slurm = write_openmm_production_files(tmp_path, production_ns=100, job_name="P3HB4_OPC")

    source = script.read_text()
    ast.parse(source)  # valid Python
    assert "TOTAL_STEPS = 50000000" in source
    assert "app.PME" in source and "app.HBonds" in source
    assert "LangevinMiddleIntegrator" in source and "MonteCarloBarostat(1.0 * unit.bar" in source
    assert "loadCheckpoint" in source
    assert "#SBATCH --job-name=P3HB4_OPC" in slurm.read_text()
    assert "python run_openmm_production.py" in slurm.read_text()


def test_write_openmm_production_files_needs_system(tmp_path):
    with pytest.raises(FileNotFoundError, match="system.prmtop"):
        write_openmm_production_files(tmp_path)
