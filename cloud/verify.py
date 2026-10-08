"""Exercise a disposable checkout in the full scientific/GUI environment.

Run from the repository root: python cloud/verify.py --output /tmp/ipha-results
Never run against a working database: the GUI can create registry files.
This checks representative workflows, not every GUI action or Cloud limits.
"""

import argparse
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
os.environ["IPHASIMULATOR_PYTHON"] = sys.executable
os.environ["IPHA_OPENMM_PLATFORM"] = "CPU"
os.environ["OPENMM_CPU_THREADS"] = "1"
os.environ["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]


def imports():
    modules = ("parmed", "rdkit", "openmm", "MDAnalysis", "mdtraj", "py3Dmol",
               "stmol", "openbabel.openbabel", "acpype", "pysmiles", "cgsmiles",
               "polyply", "kneed", "streamlit", "pdb2pqr", "Bio")
    result = {}
    for name in modules:
        try:
            module = importlib.import_module(name)
            result[name] = getattr(module, "__version__", "imported")
        except Exception as exc:
            result[name] = f"FAILED: {exc}"
    for name in ("tleap", "antechamber", "parmchk2", "prepgen", "obabel", "acpype", "polyply"):
        result[name + " executable"] = shutil.which(name) or "FAILED: absent from PATH"
    if any(str(value).startswith("FAILED:") for value in result.values()):
        raise RuntimeError(json.dumps(result, indent=2))
    return result


def build(output):
    from iphasimulator.build_pha import PHAPolymerBuilder
    import parmed

    database = output / "structure_database"
    database.mkdir()
    for name in ("residue_codes.csv", "polymer_smiles.csv"):
        shutil.copy2(ROOT / "structure_database" / name, database / name)
    shutil.copytree(ROOT / "structure_database/PHA_types/3HB", database / "PHA_types/3HB")
    builder = PHAPolymerBuilder(database)
    # Exercise the full AmberTools distribution, including the PREPGEN binary
    # omitted by some smaller pip bundles. Only the disposable fixture changes.
    prepins = builder.generate_polymer_prepins("3HB")
    for name in ("head_prepin", "mainchain_prepin", "tail_prepin"):
        assert Path(prepins[name]).stat().st_size > 0, prepins[name]
    result = builder.build_PHA_polymer("3HB", 4)
    for name in ("prmtop_file", "rst7_file", "pdb_file"):
        path = Path(result[name])
        assert path.is_file() and path.stat().st_size, f"Missing build output: {path}"
    structure = parmed.load_file(str(result["prmtop_file"]), xyz=str(result["rst7_file"]))
    assert len(structure.atoms) > 0
    return {key: str(value) for key, value in result.items() if key.endswith("_file")}


def simulate(built, output):
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    import MDAnalysis as mda

    topology = app.AmberPrmtopFile(built["prmtop_file"])
    coordinates = app.AmberInpcrdFile(built["rst7_file"])
    system = topology.createSystem(nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                            0.001 * unit.picoseconds)
    simulation = app.Simulation(topology.topology, system, integrator,
                                mm.Platform.getPlatformByName("CPU"))
    simulation.context.setPositions(coordinates.positions)
    simulation.minimizeEnergy(maxIterations=100)
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin, 42)
    trajectory = output / "short.dcd"
    reporter = app.DCDReporter(str(trajectory), 5)
    simulation.reporters.append(reporter)
    simulation.step(20)
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    assert np.isfinite(energy)
    simulation.reporters.clear()
    del reporter
    universe = mda.Universe(built["prmtop_file"], str(trajectory))
    assert len(universe.trajectory) == 4
    return {"platform": "CPU", "steps": 20, "frames": 4, "energy_kj_mol": energy}


def analysis():
    import numpy as np
    from iphasimulator.analysis.descriptors import ChainDistanceDescriptors
    from iphasimulator.analysis.pca import calculate_chain_pca
    from iphasimulator.analysis.tg_analysis.glass_transition import TemperatureResponse, fit_historical_tg

    distances = np.random.default_rng(42).normal(size=(30, 6))
    descriptors = ChainDistanceDescriptors("test", 0, "PHA", np.arange(30), distances, 4)
    pca = calculate_chain_pca(descriptors, 3)
    assert pca.transformed_data.shape == (30, 3)
    temperatures = np.linspace(150, 450, 61)
    response = 5 * (1 - np.tanh(0.025 * temperatures - 7.5)) - 1
    fit = fit_historical_tg(TemperatureResponse(temperatures, response, np.ones(61)))
    assert abs(fit.tg - 300) < 1, fit.tg
    return {"pca_shape": list(pca.transformed_data.shape), "synthetic_tg_kelvin": fit.tg}


def script_generation(output):
    from iphasimulator.openmmscript_builder import OpenMMScriptBuilder

    builder = OpenMMScriptBuilder("P3HB_4_dry", "dry", run_name="CloudSmoke")
    builder.add_minimization()
    builder.add_basic_NVT(total_steps=20, temp=300, filename="short_nvt")
    path = builder.write_script(output / "generated_openmm.py")
    compile(Path(path).read_text(), str(path), "exec")
    return {"generated_script": str(path), "validation": "generated and compiled; not executed"}


def gui():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_file(str(ROOT / "pha_gui.py"), default_timeout=90).run()
    assert not app.exception, [item.message for item in app.exception]
    assert len(app.tabs) == 7, [item.value for item in app.error]
    return {"tabs": [item.label for item in app.tabs], "displayed_errors": [item.value for item in app.error]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    # Keep third-party caches out of the installed environment and home folder.
    for variable, folder in (("NUMBA_CACHE_DIR", "numba"), ("MPLCONFIGDIR", "matplotlib"),
                             ("XDG_CACHE_HOME", "cache")):
        os.environ.setdefault(variable, str(output / folder))
    os.chdir(ROOT)
    report = {"python": sys.version, "platform": platform.platform(), "checks": {}}
    built = {}

    def run(name, callback):
        print(f"CHECK {name}", flush=True)
        try:
            details = callback()
            report["checks"][name] = {"status": "passed", "details": details}
            return details
        except Exception:
            error = traceback.format_exc()
            print(error, flush=True)
            report["checks"][name] = {"status": "failed", "error": error}
        finally:
            (output / "report.json").write_text(json.dumps(report, indent=2))

    run("dependencies", imports)
    built = run("polymer_build", lambda: build(output))
    if built:
        run("cpu_simulation_and_trajectory", lambda: simulate(built, output))
    else:
        report["checks"]["cpu_simulation_and_trajectory"] = {"status": "blocked", "reason": "polymer build failed"}
    run("analysis", analysis)
    run("script_generation", lambda: script_generation(output))
    run("gui", gui)
    passed = all(check["status"] == "passed" for check in report["checks"].values())
    report["passed"] = passed
    (output / "report.json").write_text(json.dumps(report, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
