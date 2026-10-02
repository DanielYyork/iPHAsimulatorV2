#!/usr/bin/env python
"""Turn a CHARMM-GUI Solution Builder download (enzyme + PHA) into a GROMACS run folder.

Route 05B -> 06B: protein CHARMM36m, PHA CGenFF (CHARMM-GUI Ligand Reader & Modeler),
CHARMM TIP3P water, SOD/CLA ions, GROMACS.

What it does
  1. Reports the CGenFF penalty scores of the PHA (lig/lig.rtf, lig/lig.prm).
  2. Builds a NEW run folder from the download's gromacs/ folder:
     step3_input.gro -> step5_input.gro, CHARMM-GUI's step4/step5 mdp files are left out,
     and the project's run files are added:
       step6.0_minimization.mdp   (minimisation, restraints)
       step6.1_nvt.mdp            (NVT 125 ps, 303.15 K, restraints)
       step6.2_npt.mdp            (NPT 500 ps, 1 bar, restraints)
       step7_production.mdp       (NPT production, 200 ns)
       run_step6_local.sh         (local minimisation)
       run_hpc_equilibration_production.slurm (HPC: 6.1 -> 6.2 -> 7)
  3. Checks the folder: CHARMM [ defaults ], includes, atom counts, molecules,
     net charge, CGenFF atom types, PHA charges equal lig.rtf, PHA stereocentres
     still R (needs RDKit and the SDF that was uploaded to CHARMM-GUI).
  4. Optionally runs the local minimisation (needs gmx on PATH).

The download is never modified. The run folder must not exist yet.

Example (ANC55 + P3HB4):
  python examples/charmm_gui_enzyme_pha/prepare_enzyme_pha_run.py \
      --download ~/Downloads/Anc55_ANC_56_PHB4/Anc55_PHB4-9000069110 \
      --sdf      ~/Data/MD_projects/PHA/PHA_Simulation_set_up/polymer_structures/PHB4_R.sdf \
      --out      ~/Data/Input/PHA_new_Enzyme/ANC55_P3HB4_gromacs_new
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from iphasimulator.charmmgui_import import (  # noqa: E402
    check_cgenff_penalties,
    check_ligand_charges_match_rtf,
    check_ligand_signed_volumes,
    describe_run_folder_protocol,
    find_lig_topology,
    format_cgenff_penalty_report,
    format_checks,
    prepare_charmm_gui_run_folder,
    read_cgenff_penalties,
    read_minimization_result,
    validate_charmm_gui_run_folder,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--download", required=True, type=Path,
                        help="Unpacked CHARMM-GUI Solution Builder download (contains gromacs/ and lig/).")
    parser.add_argument("--out", required=True, type=Path,
                        help="New run folder to create (must not exist).")
    parser.add_argument("--sdf", type=Path, default=None,
                        help="The R-configured PHA SDF uploaded to CHARMM-GUI (for the stereo check).")
    parser.add_argument("--job-name", default=None,
                        help="SLURM job name (default: the run folder's name).")
    parser.add_argument("--minimise", action="store_true",
                        help="Run step 6.0 minimisation locally with run_step6_local.sh (needs gmx).")
    args = parser.parse_args()

    download = args.download.expanduser().resolve()
    lig_dir = download / "lig"
    sdf = args.sdf.expanduser().resolve() if args.sdf else None

    print("=" * 70)
    print("1. CGenFF penalties of the PHA (<10 fair, 10-50 check, >50 validate)")
    print("=" * 70)
    if lig_dir.is_dir():
        print(format_cgenff_penalty_report(read_cgenff_penalties(lig_dir)))
    else:
        print(f"No lig/ folder in {download}; skipping.")

    print("\n" + "=" * 70)
    print("2. Build the run folder")
    print("=" * 70)
    run_dir = prepare_charmm_gui_run_folder(download, args.out, job_name=args.job_name)
    print(f"Created {run_dir}")
    for stage in describe_run_folder_protocol(run_dir):
        print(f"  {stage.name:28s} {stage.length_text()}")

    print("\n" + "=" * 70)
    print("3. Checks")
    print("=" * 70)
    checks = validate_charmm_gui_run_folder(run_dir, sdf_path=sdf)
    checks.append(check_cgenff_penalties(lig_dir))
    checks.append(check_ligand_charges_match_rtf(run_dir / "topol.top", find_lig_topology(lig_dir)))
    checks.append(check_ligand_signed_volumes(run_dir / "step5_input.gro", run_dir / "topol.top", sdf))
    print(format_checks(checks))
    failed = [c for c in checks if c.status == "FAIL"]
    if failed:
        print(f"\n{len(failed)} check(s) FAILED - inspect before running.")
        return 1

    if args.minimise:
        print("\n" + "=" * 70)
        print("4. Local minimisation (step 6.0)")
        print("=" * 70)
        subprocess.run(["bash", "run_step6_local.sh"], cwd=run_dir, check=True)
        result = read_minimization_result(run_dir / "step6.0_minimization.log")
        print(f"Potential energy: {result.potential_energy}  Max force: {result.maximum_force}")

    print("\n" + "=" * 70)
    print("Next steps")
    print("=" * 70)
    if not args.minimise:
        print(f"  cd {run_dir} && bash run_step6_local.sh      # local minimisation")
    print(f"  copy {run_dir.name}/ to the cluster, then:")
    print("  sbatch run_hpc_equilibration_production.slurm  # NVT -> NPT -> production")
    return 0


if __name__ == "__main__":
    sys.exit(main())
