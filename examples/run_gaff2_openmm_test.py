"""Run the staged GAFF2 -> OpenMM dry-polymer test workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse

from iphasimulator.parameterization_gaff2 import (
    AmberToolsError,
    ambertools_available,
    parameterize_gaff2,
)
from iphasimulator.simulation_openmm_amber_runner import (
    OpenMMRunnerError,
    openmm_available,
    run_openmm_with_amber_topology,
)


TARGETS = ("P3HB_4", "P3HB_8", "P3HO_4", "P3HO_8", "P3HDD_4", "P3HDD_8")
DEFAULT_TARGETS = ("P3HB_4", "P3HO_4", "P3HDD_4")
CHARGE_METHODS = ("abcg2", "bcc", "gas", "mul")


@dataclass(frozen=True)
class WorkflowTarget:
    name: str
    sdf_path: Path
    output_dir: Path

    @property
    def gaff2_dir(self) -> Path:
        return self.output_dir / "gaff2"

    @property
    def openmm_dir(self) -> Path:
        return self.output_dir / "openmm" / "dry_polymer"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _find_sdf(output_root: Path, name: str) -> Path | None:
    structures_root = output_root / "polymer_structures"
    candidates = [
        structures_root / name / f"{name}.sdf",
        structures_root / f"{name}.sdf",
        output_root / name / f"{name}.sdf",
        output_root / f"{name}.sdf",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def discover_targets(
    output_root: Path,
    md_root: Path,
    *,
    target_names: tuple[str, ...] = DEFAULT_TARGETS,
) -> list[WorkflowTarget]:
    targets: list[WorkflowTarget] = []
    for name in target_names:
        sdf_path = _find_sdf(output_root, name)
        if sdf_path is None:
            print(f"Skipping {name}: no SDF found under {output_root}")
            continue
        targets.append(
            WorkflowTarget(
                name=name,
                sdf_path=sdf_path,
                output_dir=md_root / name,
            )
        )
    return targets


def run_workflow(
    targets: list[WorkflowTarget],
    *,
    run_openmm: bool,
    charge_method: str,
    atom_count_warning_threshold: int | None,
    verbose: bool,
    minimization_max_iterations: int,
    nvt_steps: int,
    npt_steps: int,
    production_steps: int,
) -> None:
    if not ambertools_available():
        raise AmberToolsError(
            "AmberTools is required for Stage 1. Install AmberTools and ensure "
            "antechamber, parmchk2, and tleap are on PATH."
        )

    can_run_openmm = openmm_available()
    if run_openmm and not can_run_openmm:
        print("OpenMM is not installed; Stage 2 will be skipped.")

    for target in targets:
        print(f"[{target.name}] Stage 1 GAFF2 from {target.sdf_path}")
        gaff2_outputs = parameterize_gaff2(
            target.sdf_path,
            target.gaff2_dir,
            name=target.name,
            net_charge=0,
            residue_name="PHA",
            charge_method=charge_method,
            atom_count_warning_threshold=atom_count_warning_threshold,
            verbose=verbose,
        )
        print(f"[{target.name}] Wrote {gaff2_outputs.prmtop_path}")
        print(f"[{target.name}] Wrote {gaff2_outputs.inpcrd_path}")
        print(f"[{target.name}] Timings: {gaff2_outputs.timing_log}")

        if not run_openmm or not can_run_openmm:
            continue

        print(f"[{target.name}] Stage 2 OpenMM")
        openmm_outputs = run_openmm_with_amber_topology(
            gaff2_outputs.prmtop_path,
            gaff2_outputs.inpcrd_path,
            target.openmm_dir,
            minimization_max_iterations=minimization_max_iterations,
            nvt_steps=nvt_steps,
            npt_steps=npt_steps,
            production_steps=production_steps,
        )
        print(f"[{target.name}] Wrote {openmm_outputs.minimized_pdb_path}")
        print(f"[{target.name}] Wrote {openmm_outputs.nvt_trajectory_path}")
        print(f"[{target.name}] Wrote {openmm_outputs.production_trajectory_path}")
        print(f"[{target.name}] Summary: {openmm_outputs.summary_log_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GAFF2 parameterisation and optional OpenMM tests."
    )
    parser.add_argument(
        "--skip-openmm",
        action="store_true",
        help="Run only Stage 1 GAFF2 parameterisation.",
    )
    parser.add_argument(
        "--target",
        action="append",
        choices=TARGETS,
        help=(
            "Run one target. May be passed more than once. Defaults to P3HB_4, "
            "P3HO_4, and P3HDD_4."
        ),
    )
    parser.add_argument(
        "--charge-method",
        choices=CHARGE_METHODS,
        default="abcg2",
        help="antechamber charge method (default abcg2). The existing polymer-only benchmark systems used bcc.",
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
        help="OpenMM NPT equilibration steps. Requires periodic AMBER coordinates.",
    )
    parser.add_argument(
        "--production-steps",
        type=int,
        default=100,
        help="OpenMM production steps after NVT/NPT.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = _repo_root()
    output_root = repo_root / "examples" / "output"
    md_root = output_root / "md_tests"
    target_names = tuple(args.target) if args.target else DEFAULT_TARGETS
    charge_method = "gas" if args.skip_am1_bcc else args.charge_method
    atom_count_warning_threshold = (
        None
        if args.atom_count_warning_threshold <= 0
        else args.atom_count_warning_threshold
    )
    targets = discover_targets(output_root, md_root, target_names=target_names)

    if not targets:
        raise SystemExit(f"No target SDF files found under {output_root}")

    try:
        run_workflow(
            targets,
            run_openmm=not args.skip_openmm,
            charge_method=charge_method,
            atom_count_warning_threshold=atom_count_warning_threshold,
            verbose=not args.quiet_ambertools,
            minimization_max_iterations=args.minimization_iterations,
            nvt_steps=args.nvt_steps,
            npt_steps=args.npt_steps,
            production_steps=args.production_steps,
        )
    except (AmberToolsError, OpenMMRunnerError) as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    main()
