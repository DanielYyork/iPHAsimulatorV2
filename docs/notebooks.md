# Notebook catalogue

The original notebooks remain in the repository. Download one for local use, or
view its existing contents on GitHub. **Documentation builds never execute notebook
cells, simulations, shell commands in notebooks, or stored outputs.** The site
links/downloads the source files rather than rendering them through an execution engine.

Install the relevant [dependencies](installation.md), start `jupyter lab` from the
repository root, and choose your environment's kernel. Review paths, selections,
run switches and outputs before executing cells. [Known import blockers](capabilities.md)
apply even if a notebook contains historical successful output.

## Construction and design

| Notebook | Download | Status / prerequisites |
| --- | --- | --- |
| [Build a PHA oligomer](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/01_build_pha_oligomer.ipynb) | {download}`01_build_pha_oligomer.ipynb <../notebooks/01_build_pha_oligomer.ipynb>` | RDKit route; import blockers |
| [Design a custom PHA](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/02_design_custom_pha.ipynb) | {download}`02_design_custom_pha.ipynb <../notebooks/02_design_custom_pha.ipynb>` | RDKit route; import blockers |
| [Validate structures](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/03_validate_structures.ipynb) | {download}`03_validate_structures.ipynb <../notebooks/03_validate_structures.ipynb>` | RDKit route; import blockers |
| [Export structures](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/04_export_structures.ipynb) | {download}`04_export_structures.ipynb <../notebooks/04_export_structures.ipynb>` | RDKit route; import blockers |

## Parameterisation and simulation

| Notebook | Download | Status / prerequisites |
| --- | --- | --- |
| [Amber/GAFF2 parameters](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/05A_amber_gaff2_parameters.ipynb) | {download}`05A_amber_gaff2_parameters.ipynb <../notebooks/05A_amber_gaff2_parameters.ipynb>` | AmberTools; import blocker |
| [CHARMM/CGenFF parameters](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/05B_charmm_cgenff_parameters.ipynb) | {download}`05B_charmm_cgenff_parameters.ipynb <../notebooks/05B_charmm_cgenff_parameters.ipynb>` | Manual CHARMM-GUI step; checks run on the example dataset by default |
| [PHA (GAFF2) in water, GROMACS](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/06A_gaff2_gromacs_pha_in_water.ipynb) | {download}`06A_gaff2_gromacs_pha_in_water.ipynb <../notebooks/06A_gaff2_gromacs_pha_in_water.ipynb>` | Polymer benchmark method; ParmEd + GROMACS; run flags off by default |
| [PHA (CGenFF) + enzyme in water, GROMACS](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/06B_cgenff_gromacs_pha_enzyme_in_water.ipynb) | {download}`06B_cgenff_gromacs_pha_enzyme_in_water.ipynb <../notebooks/06B_cgenff_gromacs_pha_enzyme_in_water.ipynb>` | Docked complex PDB (notebook 11); CHARMM-GUI Solution Builder download (default: the example dataset `examples/data/charmm_gui_ANC55_P3HB4/`); GROMACS optional |
| [Optional: PHA (GAFF2) + enzyme in water, OpenMM](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/06C_optional_gaff2_openmm_pha_enzyme_in_water.ipynb) | {download}`06C_optional_gaff2_openmm_pha_enzyme_in_water.ipynb <../notebooks/06C_optional_gaff2_openmm_pha_enzyme_in_water.ipynb>` | Docked complex PDB (notebook 11); tleap + OpenMM; a different force-field setup from 06B; run flags off by default |
| [Optional: quick check of 05A parameters, OpenMM](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/05A_quick_check_openmm.ipynb) | {download}`05A_quick_check_openmm.ipynb <../notebooks/05A_quick_check_openmm.ipynb>` | 05A GAFF2 files in OpenMM without water; not a physical result; run flag off by default |

## Execution and analysis

For enzyme–PHA trajectory preparation, use the [standalone Bash templates](enzyme_trajectory_processing.md); the superseded GK13 preparation notebook is no longer needed.

| Notebook | Download | Status / prerequisites |
| --- | --- | --- |
| [HPC execution](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/07_hpc_execution.ipynb) | {download}`07_hpc_execution.ipynb <../notebooks/07_hpc_execution.ipynb>` | Cluster-specific; configured runner imports blocked |
| [Trajectory preprocessing](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/08_trajectory_preprocessing.ipynb) | {download}`08_trajectory_preprocessing.ipynb <../notebooks/08_trajectory_preprocessing.ipynb>` | Requires GROMACS and MD inputs |
| [Solvated polymer analysis](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/09_solvated_polymer_analysis.ipynb) | {download}`09_solvated_polymer_analysis.ipynb <../notebooks/09_solvated_polymer_analysis.ipynb>` | MDTraj notebook workflow |
| [Polymer benchmark batch](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/10_polymer_benchmark_batch.ipynb) | {download}`10_polymer_benchmark_batch.ipynb <../notebooks/10_polymer_benchmark_batch.ipynb>` | Fixed workflow assumptions; import blockers |
| [Enzyme docking setup](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/11_enzyme_docking_setup.ipynb) | {download}`11_enzyme_docking_setup.ipynb <../notebooks/11_enzyme_docking_setup.ipynb>` | Manual; exports the polymer only from `examples/output/benchmark/<SYSTEM>/` after a whole-molecule check; benchmark import blocker |
| [Enzyme–polymer analysis](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/12_enzyme_polymer_analysis.ipynb) | {download}`12_enzyme_polymer_analysis.ipynb <../notebooks/12_enzyme_polymer_analysis.ipynb>` | Energy and RMSD diagnostics; system-specific inputs |
| [PHA–enzyme contacts](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/src/md_simulation_scripts/enzyme_contacts/enzyme_contacts.ipynb) | {download}`enzyme_contacts.ipynb <../src/md_simulation_scripts/enzyme_contacts/enzyme_contacts.ipynb>` | Reusable, validated preview; matching TPR/XTC required |

The two `01_APO_*` files in `src/md_simulation_scripts/` are currently empty
placeholders. They are preserved but are not runnable tutorials.

## Workflow order

The notebooks are workflow modules, not a 01→12 sequence. Notebooks 01–04 build the PHA;
then pick one workflow; 07 runs it on HPC; 08–12 analyse:

- PHA alone in water: 05A → 06A
- Enzyme–PHA in water (the method used for the project's research simulations): 05B → 06B
- Optional enzyme–PHA alternative with OpenMM: 05A → 06C
- Optional quick check of 05A parameters: 05A → 05A_quick_check_openmm (outside the main workflows)

Both enzyme workflows need a docked enzyme–PHA complex PDB first (docking preparation:
notebook 11, then manual HADDOCK). For enzyme contacts, use the original periodic coordinates rather than
an independently fitted trajectory with an untransformed box.

The earlier [notebook guide](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/README.md)
is retained in the repository.
