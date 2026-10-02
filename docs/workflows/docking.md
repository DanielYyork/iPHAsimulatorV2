# Docking preparation

**Status: notebook-based, manual preparation; incomplete as an end-to-end pipeline.**
There is no automated HADDOCK submission, validated enzyme–polymer complex
builder, or catalytic-activity predictor in this workflow.

## 1. Choose the polymer structure

Both enzyme workflows (05B → 06B and 05A → 06C) start from a docked enzyme–PHA complex
PDB. Notebook **11** reads the polymer benchmark folder that notebook 10 writes,
`examples/output/benchmark/<system>/gromacs/solvated_polymer/`, by default its final
structure `step7_production.gro` (`GRO_NAME`; `step5_input.gro` is the structure before MD).
Review that path and confirm the intended system and structure. Its
benchmark import also depends on modules with [runtime blockers](../capabilities.md).

## 2. Prepare and inspect the PDB

The notebook writes the polymer only (residue `POLYMER_RESNAME`: `PHA` for 06A
systems, `LIG` for CHARMM-GUI systems) to `examples/output/docking_inputs/<system>/`,
with {py:func}`iphasimulator.docking_inputs.export_polymer_pdb`. Before writing it
checks that the polymer is whole: every bonded heavy-atom pair (bonds from the
topology) must be shorter than 0.3 nm. If not, it stops; use a whole-molecule frame
instead, such as the representative frame from notebook 08 or `gmx trjconv -pbc mol`.

Inspect chain/residue identifiers, atom names, elements and terminal groups in the
exported PDB. Elements come from the topology masses (or the atom names); atoms and
residues are renumbered from 1 in chain A, and no box record is written.

## 3. Prepare the enzyme separately

Supply the reviewed enzyme structure and decide the docking restraints for the
specific research question. The enzyme annotations and binding statements in
the historical notebook are user-supplied context, not independently validated
results of this package. Confirm them against your own evidence before use.

## 4. Record the manual docking job

Upload inputs through your chosen docking service when ready, and retain the
input files, restraints, job identifiers and returned results. The notebook
provides preparation notes and records; it does not submit the job.

For an existing MD trajectory of a complex, use [contact analysis](analysis.md)
to examine proximity. Contact occupancy cannot establish affinity or catalytic
activity.

Notebook: [11 docking preparation](../notebooks.md#execution-and-analysis).
