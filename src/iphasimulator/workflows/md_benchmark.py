"""Batch MD benchmark workflow for PHA validation systems."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse

from iphasimulator.parameterization_gaff2 import (
    AmberToolsError,
    parameterize_gaff2,
)
from iphasimulator.simulation_gromacs_runner import (
    prepare_gromacs_run_folder,
    write_gromacs_solvation_files,
)
from iphasimulator.simulation_openmm_amber_runner import (
    OpenMMRunnerError,
    openmm_available,
    run_openmm_with_amber_topology,
)


SYSTEMS = ("P3HB_4", "P3HB_8", "P3HO_4", "P3HO_8", "P3HDD_4", "P3HDD_8")
CHARGE_METHODS = ("abcg2", "bcc", "gas", "mul")
DEFAULT_REPO_ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class BenchmarkTarget:
    """Input and output paths for one benchmark system."""

    name: str
    sdf_path: Path
    output_dir: Path

    @property
    def gaff2_dir(self) -> Path:
        return self.output_dir / "gaff2"

    @property
    def prmtop_path(self) -> Path:
        return self.gaff2_dir / f"{self.name}.prmtop"

    @property
    def inpcrd_path(self) -> Path:
        return self.gaff2_dir / f"{self.name}.inpcrd"

    @property
    def openmm_dir(self) -> Path:
        return self.output_dir / "openmm" / "dry_polymer"

    @property
    def gromacs_dir(self) -> Path:
        return self.output_dir / "gromacs"


def _find_sdf(output_root: Path, system_name: str) -> Path | None:
    structures_root = output_root / "polymer_structures"
    candidates = (
        structures_root / f"{system_name}.sdf",
        structures_root / system_name / f"{system_name}.sdf",
        output_root / f"{system_name}.sdf",
        output_root / system_name / f"{system_name}.sdf",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def discover_targets(
    system_names: tuple[str, ...],
    *,
    repo_root: Path = DEFAULT_REPO_ROOT,
) -> list[BenchmarkTarget]:
    output_root = repo_root / "examples" / "output"
    md_root = output_root / "benchmark"
    targets: list[BenchmarkTarget] = []

    for system_name in system_names:
        sdf_path = _find_sdf(output_root, system_name)
        if sdf_path is None:
            print(
                f"[{system_name}] missing SDF input under "
                f"{output_root / 'polymer_structures'}; skipping"
            )
            continue
        targets.append(
            BenchmarkTarget(
                name=system_name,
                sdf_path=sdf_path,
                output_dir=md_root / system_name,
            )
        )

    return targets


def selected_system_names(args: argparse.Namespace) -> tuple[str, ...]:
    names: list[str] = []
    for option in (args.system, args.systems, args.target):
        if option:
            names.extend(option)

    if not names:
        return SYSTEMS

    return tuple(dict.fromkeys(names))


def _run_gaff2(
    target: BenchmarkTarget,
    *,
    charge_method: str,
    atom_count_warning_threshold: int | None,
    verbose: bool,
    reuse_existing_gaff2: bool,
) -> None:
    if (
        reuse_existing_gaff2
        and target.prmtop_path.exists()
        and target.inpcrd_path.exists()
    ):
        print(f"[{target.name}] reusing {target.prmtop_path}")
        print(f"[{target.name}] reusing {target.inpcrd_path}")
        return

    print(f"[{target.name}] GAFF2 parameterisation from {target.sdf_path}")
    outputs = parameterize_gaff2(
        target.sdf_path,
        target.gaff2_dir,
        name=target.name,
        net_charge=0,
        residue_name="PHA",
        charge_method=charge_method,
        atom_count_warning_threshold=atom_count_warning_threshold,
        verbose=verbose,
    )
    print(f"[{target.name}] wrote {outputs.prmtop_path}")
    print(f"[{target.name}] wrote {outputs.inpcrd_path}")
    print(f"[{target.name}] timings {outputs.timing_log}")


def _run_openmm(
    target: BenchmarkTarget,
    *,
    minimization_iterations: int,
    nvt_steps: int,
    npt_steps: int,
    production_steps: int,
) -> None:
    print(f"[{target.name}] OpenMM dry-polymer MD")
    outputs = run_openmm_with_amber_topology(
        target.prmtop_path,
        target.inpcrd_path,
        target.openmm_dir,
        minimization_max_iterations=minimization_iterations,
        nvt_steps=nvt_steps,
        npt_steps=npt_steps,
        production_steps=production_steps,
    )
    print(f"[{target.name}] wrote {outputs.minimized_pdb_path}")
    print(f"[{target.name}] wrote {outputs.final_pdb_path}")
    print(f"[{target.name}] summary {outputs.summary_log_path}")


def _prepare_gromacs(target: BenchmarkTarget) -> None:
    print(f"[{target.name}] GROMACS workflow folder preparation")
    outputs = prepare_gromacs_run_folder(
        target.prmtop_path,
        target.inpcrd_path,
        target.gromacs_dir,
        target.name,
    )
    print(f"[{target.name}] wrote {outputs.dry_polymer_dir}")
    print(f"[{target.name}] wrote {outputs.solvated_polymer_dir}")
    print(f"[{target.name}] wrote {outputs.charmm_gui_membrane_dir}")

    solvation_outputs = write_gromacs_solvation_files(
        outputs.solvated_polymer_dir,
        workflow_type="polymer",
        box_padding_nm=1.2,
        ion_concentration_molar=0.15,
        clean=True,
    )
    print(f"[{target.name}] wrote {solvation_outputs.solvate_script_path}")
    print(f"[{target.name}] wrote {solvation_outputs.ions_mdp_path}")
    for solvation_itp_path in solvation_outputs.solvation_itp_paths:
        print(f"[{target.name}] wrote {solvation_itp_path}")


def run_benchmark(args: argparse.Namespace) -> int:
    system_names = selected_system_names(args)
    charge_method = "gas" if args.skip_am1_bcc else args.charge_method
    atom_count_warning_threshold = (
        None
        if args.atom_count_warning_threshold <= 0
        else args.atom_count_warning_threshold
    )
    targets = discover_targets(system_names)

    if not targets:
        print("No benchmark inputs were found.")
        return 1

    can_run_openmm = openmm_available()
    if args.skip_openmm:
        can_run_openmm = False
    elif not can_run_openmm:
        print("OpenMM is not installed; dry-polymer MD will be skipped.")

    failures: list[tuple[str, str]] = []
    for target in targets:
        try:
            target.output_dir.mkdir(parents=True, exist_ok=True)
            _run_gaff2(
                target,
                charge_method=charge_method,
                atom_count_warning_threshold=atom_count_warning_threshold,
                verbose=not args.quiet_ambertools,
                reuse_existing_gaff2=args.reuse_existing_gaff2,
            )
            if can_run_openmm:
                _run_openmm(
                    target,
                    minimization_iterations=args.minimization_iterations,
                    nvt_steps=args.nvt_steps,
                    npt_steps=args.npt_steps,
                    production_steps=args.production_steps,
                )
            if args.prepare_gromacs:
                _prepare_gromacs(target)
        except (
            AmberToolsError,
            OpenMMRunnerError,
            ImportError,
            FileNotFoundError,
            ValueError,
        ) as exc:
            failures.append((target.name, str(exc)))
            print(f"[{target.name}] FAILED: {exc}")
            if args.stop_on_error:
                break

    if failures:
        print("\nFailures:")
        for system_name, message in failures:
            print(f"- {system_name}: {message}")
        return 1

    print("\nBenchmark workflow completed for:")
    for target in targets:
        print(f"- {target.name}: {target.output_dir}")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--system",
        action="append",
        choices=SYSTEMS,
        help=(
            "Run one system. May be passed more than once. Defaults to all six "
            "if no system option is provided."
        ),
    )
    parser.add_argument(
        "--systems",
        nargs="+",
        choices=SYSTEMS,
        help=(
            "Run one or more systems, for example --systems P3HB_4 P3HB_8. "
            "Defaults to all six if no system option is provided."
        ),
    )
    parser.add_argument(
        "--target",
        action="append",
        choices=SYSTEMS,
        help="Deprecated alias for --system.",
    )
    parser.add_argument(
        "--skip-openmm",
        action="store_true",
        help="Run GAFF2 parameterisation only.",
    )
    parser.add_argument(
        "--prepare-gromacs",
        action="store_true",
        help=(
            "Also prepare GROMACS dry, solvated, and CHARMM-GUI-style folders, "
            "including solvated_polymer solvation files."
        ),
    )
    parser.add_argument(
        "--reuse-existing-gaff2",
        action="store_true",
        help="Reuse existing prmtop/inpcrd files instead of rerunning AmberTools.",
    )
    parser.add_argument(
        "--charge-method",
        choices=CHARGE_METHODS,
        default="abcg2",
        help="antechamber charge method (default abcg2). Existing project benchmark/enzyme runs used bcc.",
    )
    parser.add_argument(
        "--skip-am1-bcc",
        action="store_true",
        help="Shortcut for --charge-method gas; bypasses the AM1 step used by both abcg2 and bcc.",
    )
    parser.add_argument(
        "--atom-count-warning-threshold",
        type=int,
        default=120,
        help="Warn when prepared molecules exceed this atom count. Use 0 to disable.",
    )
    parser.add_argument(
        "--quiet-ambertools",
        action="store_true",
        help="Do not print each AmberTools command before it runs.",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop after the first system failure instead of reporting all failures.",
    )
    parser.add_argument(
        "--minimization-iterations",
        type=int,
        default=200,
        help="Maximum OpenMM minimization iterations.",
    )
    parser.add_argument(
        "--nvt-steps",
        type=int,
        default=100,
        help="OpenMM NVT equilibration steps.",
    )
    parser.add_argument(
        "--npt-steps",
        type=int,
        default=100,
        help="OpenMM NPT equilibration steps when box vectors are present.",
    )
    parser.add_argument(
        "--production-steps",
        type=int,
        default=100,
        help="OpenMM production steps.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    raise SystemExit(run_benchmark(parse_args(argv)))


if __name__ == "__main__":
    main()
