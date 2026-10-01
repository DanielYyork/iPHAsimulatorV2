# CHARMM/GROMACS run templates

Run files for a CHARMM-GUI Solution Builder GROMACS package (CHARMM36m +
CGenFF + CHARMM TIP3P + SOD/CLA). Notebook 06B copies them next to the package's
`gromacs/` files after renaming `step3_input.gro` to `step5_input.gro`.

Source: the GK13_P3HO_4 enzyme–PHA production run (CHARMM-GUI job 8214536317).

| File | Stage | Notes |
| --- | --- | --- |
| `step6.0_minimization.mdp` | steepest descent | identical to CHARMM-GUI's `step4.0_minimization.mdp`; backbone/side-chain restraints |
| `step6.1_nvt.mdp` | NVT, 303.15 K | dt 1 fs, 125000 steps (125 ps); restraints |
| `step6.2_npt.mdp` | NPT, 1 bar, C-rescale | dt 2 fs, 250000 steps (500 ps); restraints |
| `step7_production.mdp` | NPT production | dt 2 fs, 100000000 steps (200 ns) |
| `run_step6_local.sh` | local minimisation | `grompp` + `mdrun` for step 6.0 |
| `run_hpc_equilibration_production.slurm` | HPC steps 6.1–7 | `{JOB_NAME}` is filled in by 06B; module line is for the KCL GPU cluster |

The mdp files use the `SOLU`/`SOLV` index groups that CHARMM-GUI writes to `index.ndx`.
