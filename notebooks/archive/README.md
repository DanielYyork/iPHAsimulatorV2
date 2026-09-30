# Archived notebooks

Hybrid GAFF2 → GROMACS route with CHARMM-style TIP3P/SOD/CLA. Used to generate the
six polymer-only benchmark systems (notebook 10, examples/output/benchmark/, AM1-BCC
charges) and the P3HB_4_01 example trajectory analysed in 08/09. Kept for provenance.
Superseded by Route A (05A → 06A/06B) and Route C (05B → 06C/06D); the dry GROMACS
check is superseded by 06C.

The `hybrid_` prefix separates these from the current 06B/06C/06D notebooks, which
belong to Routes A and C.

- `hybrid_06B_gromacs_dry_polymer.ipynb`: AMBER → GROMACS conversion and dry GROMACS minimisation.
- `hybrid_06C_gromacs_solvated_system.ipynb`: explicit-solvent GROMACS preparation (TIP3P, SOD/CLA) and NVT/NPT/production scripts.
- `hybrid_06D_openmm_solvated_system.ipynb`: disabled explicit-solvent OpenMM template.
