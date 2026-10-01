# CHARMM/GROMACS polymer-in-water run templates

mdp files for a CGenFF (CHARMM) PHA oligomer solvated locally in CHARMM TIP3P water
with SOD/CLA ions (notebook 06B1). They follow the GAFF2 polymer benchmark protocol
(`../gromacs_mdp/`, used for `examples/output/benchmark/`) so results are comparable;
only the non-bonded block differs, because CGenFF/CHARMM36 require CHARMM settings.

| File | Stage | Protocol |
| --- | --- | --- |
| `step6.0_minimization.mdp` | steepest descent | emtol 1000, up to 50000 steps, no restraints |
| `step6.1_nvt.mdp` | NVT | 100 ps, dt 2 fs, 300 K, V-rescale (`tc-grps = System`) |
| `step6.2_npt.mdp` | NPT | 500 ps, dt 2 fs, 300 K, 1 bar, C-rescale (τp 5 ps) |
| `step7_production.mdp` | NPT production | 100 ns, dt 2 fs, frames every 2 ps |

Non-bonded (all stages): PME, `rlist = rcoulomb = rvdw = 1.2` nm, `vdw-modifier = Force-switch`,
`rvdw_switch = 1.0` nm. The benchmark used 1.0 nm plain cut-offs, which suit GAFF2.
