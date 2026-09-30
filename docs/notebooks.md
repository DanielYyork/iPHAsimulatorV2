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
| [Build a PHA oligomer](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/01_examples_pha_oligomers.ipynb) | {download}`01_examples_pha_oligomers.ipynb <../notebooks/01_examples_pha_oligomers.ipynb>` | RDKit route; import blockers |
| [Choose a polymer design](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/02_design_polymer_for_user_request.ipynb) | {download}`02_design_polymer_for_user_request.ipynb <../notebooks/02_design_polymer_for_user_request.ipynb>` | RDKit route; import blockers |
| [Validate and visualise](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/03_validate_and_visualize.ipynb) | {download}`03_validate_and_visualize.ipynb <../notebooks/03_validate_and_visualize.ipynb>` | RDKit route; import blockers |
| [Export structures](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/04_export_structures.ipynb) | {download}`04_export_structures.ipynb <../notebooks/04_export_structures.ipynb>` | RDKit route; import blockers |

## Parameterisation and simulation

| Notebook | Download | Status / prerequisites |
| --- | --- | --- |
| [GAFF2 parameterisation](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/05A_amber_gaff2_parameterisation.ipynb) | {download}`05A_amber_gaff2_parameterisation.ipynb <../notebooks/05A_amber_gaff2_parameterisation.ipynb>` | AmberTools; import blocker |
| [CHARMM/CGenFF handoff](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/05B_charmm_cgenff_parameterisation.ipynb) | {download}`05B_charmm_cgenff_parameterisation.ipynb <../notebooks/05B_charmm_cgenff_parameterisation.ipynb>` | Incomplete manual workflow |
| [OpenMM dry polymer (Route A)](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/06A_openmm_dry_polymer.ipynb) | {download}`06A_openmm_dry_polymer.ipynb <../notebooks/06A_openmm_dry_polymer.ipynb>` | AMBER runner; import blocker |
| [GROMACS dry preparation](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/archive/hybrid_06B_gromacs_dry_polymer.ipynb) | {download}`hybrid_06B_gromacs_dry_polymer.ipynb <../notebooks/archive/hybrid_06B_gromacs_dry_polymer.ipynb>` | Archived (hybrid route, provenance) |
| [GROMACS solvation preparation](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/archive/hybrid_06C_gromacs_solvated_system.ipynb) | {download}`hybrid_06C_gromacs_solvated_system.ipynb <../notebooks/archive/hybrid_06C_gromacs_solvated_system.ipynb>` | Archived (hybrid route, provenance) |
| [OpenMM solvation](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/archive/hybrid_06D_openmm_solvated_system.ipynb) | {download}`hybrid_06D_openmm_solvated_system.ipynb <../notebooks/archive/hybrid_06D_openmm_solvated_system.ipynb>` | Archived (hybrid route, provenance) |

## Execution and analysis

For enzyme–PHA trajectory preparation, use the [standalone Bash templates](enzyme_trajectory_processing.md); the superseded GK13 preparation notebook is no longer needed.

| Notebook | Download | Status / prerequisites |
| --- | --- | --- |
| [HPC execution and restart](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/07_hpc_workflows.ipynb) | {download}`07_hpc_workflows.ipynb <../notebooks/07_hpc_workflows.ipynb>` | Cluster-specific; configured runner imports blocked |
| [Trajectory preprocessing](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/08_trajectory_preprocessing.ipynb) | {download}`08_trajectory_preprocessing.ipynb <../notebooks/08_trajectory_preprocessing.ipynb>` | Requires GROMACS and MD inputs |
| [Basic polymer analysis](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/09_basic_polymer_analysis.ipynb) | {download}`09_basic_polymer_analysis.ipynb <../notebooks/09_basic_polymer_analysis.ipynb>` | MDTraj notebook workflow |
| [Six-system MD benchmark](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/10_batch_md_benchmark.ipynb) | {download}`10_batch_md_benchmark.ipynb <../notebooks/10_batch_md_benchmark.ipynb>` | Fixed workflow assumptions; import blockers |
| [Docking preparation](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/11_PHA_Enzyme_Docking.ipynb) | {download}`11_PHA_Enzyme_Docking.ipynb <../notebooks/11_PHA_Enzyme_Docking.ipynb>` | Manual; verify ligand-only input; benchmark import blocker |
| [Enzyme/polymer stability](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/12_enzyme_polymer_stable_analysis.ipynb) | {download}`12_enzyme_polymer_stable_analysis.ipynb <../notebooks/12_enzyme_polymer_stable_analysis.ipynb>` | Energy and RMSD diagnostics; system-specific inputs |
| [PHA–enzyme contacts](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/src/md_simulation_scripts/enzyme_contacts/enzyme_contacts.ipynb) | {download}`enzyme_contacts.ipynb <../src/md_simulation_scripts/enzyme_contacts/enzyme_contacts.ipynb>` | Reusable, validated preview; matching TPR/XTC required |

The two `01_APO_*` files in `src/md_simulation_scripts/` are currently empty
placeholders. They are preserved but are not runnable tutorials.

## Workflow order

Construction → validation/export → parameterisation → one engine branch →
execution → preprocessing/analysis. Docking preparation is a separate manual
handoff. For enzyme contacts, use the original periodic coordinates rather than
an independently fitted trajectory with an untransformed box.

The earlier [notebook guide](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/notebooks/README.md)
is retained in the repository.
