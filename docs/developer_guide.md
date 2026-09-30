# iPHASimulator v2 Design

> **Historical design notes:** this document is preserved from the earlier
> implementation. Some scope/status statements below predate the current
> GROMACS, packing, GUI and docking-preparation code. Use
> [current capabilities and blockers](capabilities.md) for the inspected status,
> and [documentation maintenance](contributing_docs.md) to edit or build the site.

## Scope

This package builds PHA oligomers and prepares them for molecular simulation.

## Architecture

The builder has three levels:

1. Built-in curated PHA library for canonical monomer/residue codes: `3HB`, `3HV`, `3HHx`, `3HHep`, `3HO`, `3HN`, `3HD`, and `3HDD`.
2. Generic alkyl-side-chain PHA generation with `build_pha_by_sidechain(side_chain_carbons, degree)`.
3. Custom PHA generation with `build_custom_pha(monomer_smiles, degree, name)` for user-defined 3-hydroxy acid monomers with branched, unsaturated, or functionalised side chains.

The design goal is curated common PHA chemistry plus generalisable custom generation, not manual enumeration of every possible PHA.

## Naming convention

Name generation and validation live in `src/iphasimulator/naming.py`.

| Meaning | Examples |
|---|---|
| Monomer / residue code | `3HB`, `3HO`, `3HDD` |
| Polymer code | `P3HB`, `P3HO`, `P3HDD` |
| Single oligomer chain | `P3HB_4`, `P3HO_8`, `P3HDD_4` |
| Multi-chain system | `25_P3HB_3`, `10_P3HO_8` |
| Head/main/tail residue entries | `3HB_H`, `3HB_M`, `3HB_T` |

Notebook and script outputs should use polymer-level names such as `P3HB_4`.
Residue-level database entries should use the head/main/tail suffixes.

## Project organisation

- `src/iphasimulator/` contains reusable Python modules.
- `src/iphasimulator/parameterization_gaff2.py` contains AmberTools/GAFF2 code.
- `src/iphasimulator/simulation_openmm_amber_runner.py` contains OpenMM execution code.
- `src/iphasimulator/workflows/` contains reusable workflow helpers for tutorials and examples.
- `notebooks/` contains separate step-by-step Jupyter tutorials written for
  users who may not be computational specialists.
- `examples/` contains thin command-line scripts only.
- `examples/output/` contains generated structures and MD outputs.
- `tests/` contains automated checks for package behaviour.

Tutorial notebooks should call package functions from `src/iphasimulator` instead of
owning reusable chemistry, export, parameterization, or simulation logic.

User-directed polymer design should be exposed through small workflow helpers
such as `PolymerDesign` and `design_polymer`, so notebooks can offer plain input
choices without duplicating builder logic.

Scripted and HPC workflows should use YAML configuration files plus command-line
runners. Notebooks may explain the workflow, but repeatable production-style runs
should be launched through scripts such as `examples/run_configured_workflow.py`.

## Debug validation systems

- Default workflow: `P3HB_4`, `P3HO_4`, and `P3HDD_4`.
- Larger opt-in systems: `P3HB_8`, `P3HO_8`, and `P3HDD_8`.

The n = 8 systems are intentionally excluded from the default AmberTools/OpenMM
workflow until the smaller systems are validated. AM1-BCC charge generation in
`antechamber` can be very slow for larger `P3HDD` oligomers.

## Core design

RDKit-first polymer construction.

The polymer molecule is built chemically first, then exported to PDB/SDF.
The initial MD workflow is split into workflow-oriented notebooks:

1. GAFF2 parameterisation with AmberTools (`antechamber`, `parmchk2`, `tleap`).
2. OpenMM dry validation from generated AMBER files.
3. GROMACS dry conversion and validation.
4. Engine-specific explicit-solvent preparation.
5. HPC execution and restart/benchmark management.

An alternative GROMACS route is also supported after GAFF2 parameterisation:

1. Convert AMBER `prmtop`/`inpcrd` to GROMACS `top`/`gro` with ParmEd.
2. Write separate workflow folders under
   `examples/output/md_tests/<SYSTEM>/gromacs/`: `dry_polymer/`,
   `solvated_polymer/`, and `charmm_gui_membrane/`.
3. Keep each workflow folder self-contained with its own `topol.top`,
   `step5_input.gro`, `index.ndx`, MDP files, and run scripts. The dry route
   writes a boxed polymer to `dry_polymer/step5_input.gro`; the solvation route
   resets from `dry_polymer/`, adds water/ions only in `solvated_polymer/`, and
   writes final ionised coordinates back to `solvated_polymer/step5_input.gro`.

The GROMACS route is an alternative to OpenMM, not a replacement.

OpenMM outputs are also split by workflow under
`examples/output/md_tests/<SYSTEM>/openmm/`: `dry_polymer/`,
`solvated_polymer/`, `advanced_workflows/`, and `debug/`. Debug output is not
production workflow state.

For solvated GROMACS systems, the simple GROMACS solvate/genion-style route
remains the default quick-validation path. Intermediate files such as
`step5_input_box.gro`, `step5_solvated.gro`, `genion.tpr`, and
`system_neutralized.gro` stay inside `solvated_polymer/`. The generated
solvation script deletes stale generated files before rerunning and starts from
the fresh polymer-only dry topology, which prevents duplicated solvent or ion
molecule counts in `topol.top`. A separate advanced Packmol builder in
`src/iphasimulator/system_builder_packmol.py` can prepare future
realistic PHA oligomer systems with a defined cubic box, TIP3P water, NaCl, and
multiple polymer copies. It is intentionally not a membrane builder yet, but the
API leaves room for future enzyme/polymer systems.

Notebook `05A_amber_gaff2_parameterisation.ipynb` is the roadmap before the workflow
branches into the 06-series notebooks. It explains when to use OpenMM vs
GROMACS, when dry validation is sufficient, when explicit solvent is required,
and why the project keeps both simplified polymer and CHARMM-GUI-style staged
GROMACS workflows.

Two GROMACS equilibration workflows are intentionally kept available:

- `solvated_polymer/` is the simplified polymer workflow. It is the default for
  polymer benchmarking and rapid iteration: `step6.0_minimization`,
  `step6.1_nvt`, `step6.2_npt`, then `step7_production`.
- `charmm_gui_membrane/` is the CHARMM-GUI-style staged workflow. It is an
  advanced template for membrane proteins, enzyme/polymer systems, and sensitive
  complexes: `step6.0_minimization`, `step6.1_equilibration` through
  `step6.6_equilibration`, then `step7_production`.

The simplified workflow is shorter and easier to debug, but it is less
conservative for heterogeneous systems. The CHARMM-GUI-style workflow uses more
equilibration stages so temperature, density, interfaces, and pressure coupling
are relaxed gradually; it costs more wall-clock time and has more files to
inspect. The detailed conceptual comparison belongs in notebook 05 so the 06
notebooks can stay focused on reproducible execution.

The GAFF2 workflow records the exact AmberTools command for each stage, separate
raw logs for `antechamber` and `sqm`, and per-stage timings. It defaults to
`abcg2` (ABCG2) charges. For debugging large
systems, the example runner can use faster temporary charge methods such as
`--charge-method gas` or `--skip-am1-bcc`.

OpenMM AMBER execution lives in `src/iphasimulator/simulation_openmm_amber_runner.py`
and uses existing `prmtop`/`inpcrd` files. The v2 runner performs minimisation,
NVT equilibration, NPT equilibration when periodic box vectors are present, and a
short production stage. Single-molecule GAFF2 outputs generated without a box are
nonperiodic; for those inputs the runner writes an explicit skipped-NPT log.

## Main requirements

1. Generate chemically correct PHA oligomers.
2. Preserve R stereochemistry.
3. Validate repeat number.
4. Validate side-chain length for curated and generic linear alkyl PHAs.
5. Export PDB and SDF.
6. Keep AMBER separate from polymer construction.
7. Allow custom 3-hydroxy acid monomers without requiring every possible PHA to be registered manually.
8. Keep reusable MD workflow code in root-level `src/iphasimulator/` modules with descriptive prefixes; examples only contain runnable scripts.

## MD dependencies

Use conda-forge for the full MD stack:

```bash
conda install -c conda-forge ambertools openmm parmed mdtraj -y
```

`ambertools` provides the Stage 1 command-line tools. `openmm`, `parmed`, and `mdtraj` provide the Stage 2 Python MD ecosystem.

## Out of scope for v2

- docking
- enzyme complex building
- Packmol
- GROMACS
- GUI
- general polymer chemistry outside 3-hydroxyalkanoate PHA construction
