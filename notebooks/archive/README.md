# Archived notebooks

Hybrid GAFF2 → GROMACS route with CHARMM-style TIP3P/SOD/CLA. Used to generate the
six polymer-only benchmark systems (notebook 10, examples/output/benchmark/, AM1-BCC
charges) and the P3HB_4_01 example trajectory analysed in 08/09. Kept for provenance.
Superseded by Route A (05A→06A) and Route C (05B→06B); the dry GROMACS check is
superseded by 06B §1.

- `06B_gromacs_dry_polymer.ipynb`: AMBER → GROMACS conversion and dry GROMACS minimisation.
- `06C_gromacs_solvated_system.ipynb`: explicit-solvent GROMACS preparation (TIP3P, SOD/CLA) and NVT/NPT/production scripts.
- `06D_openmm_solvated_system.ipynb`: disabled explicit-solvent OpenMM template.
