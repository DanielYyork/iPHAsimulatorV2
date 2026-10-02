# Archived notebooks

Hybrid GAFF2 → GROMACS route with CHARMM-style TIP3P/SOD/CLA. Used to generate the
six polymer-only benchmark systems (notebook 10, examples/output/benchmark/, AM1-BCC
charges) and the P3HB_4_01 example trajectory analysed in 08/09. Kept for provenance.
Superseded by 06A (GAFF2 PHA in water, GROMACS), 06B (CGenFF PHA + enzyme, GROMACS) and
06C (GAFF2 PHA + enzyme, OpenMM).

The `hybrid_` prefix separates these from the current 06 notebooks.

- `hybrid_06B_gromacs_dry_polymer.ipynb`: AMBER → GROMACS conversion and dry GROMACS minimisation.
- `hybrid_06C_gromacs_solvated_system.ipynb`: explicit-solvent GROMACS preparation (TIP3P, SOD/CLA) and NVT/NPT/production scripts.
- `hybrid_06D_openmm_solvated_system.ipynb`: disabled explicit-solvent OpenMM template.
- `06A1_amber_openmm_polymer_in_water.ipynb`: GAFF2 polymer in OPC water with staged OpenMM run files (replaced by `06A_gaff2_gromacs_pha_in_water.ipynb`, the benchmark method). Its code (`amber_solvation`, `run_openmm_md.py`) is still used by 06C.
