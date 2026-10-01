# OpenMM and GROMACS preparation

Choose one engine route after [parameterisation](gaff2.md). The 06-series
notebooks are alternatives; you do not need to run all four in sequence.

```{warning}
The AMBER OpenMM runner, AMBER-to-GROMACS converter and GROMACS preparation module
currently have [import blockers](../capabilities.md#import-blockers). This page
describes their implemented interfaces and required handoffs. No simulation is
run during a documentation build.
```

## OpenMM: AMBER topology route

1. Obtain matching `P3HB_4.prmtop` and `P3HB_4.inpcrd` from the GAFF2 step.
2. Review §0 of notebook **06A1** and choose a fresh output folder.
3. Start with a deliberately short run after the source and engine checks pass:

```python
from pathlib import Path
from iphasimulator.simulation_openmm_amber_runner import run_openmm_with_amber_topology

base = Path("examples/output/md_tests/P3HB_4")
outputs = run_openmm_with_amber_topology(
    base / "gaff2/P3HB_4.prmtop",
    base / "gaff2/P3HB_4.inpcrd",
    base / "openmm/dry_polymer",
    minimization_max_iterations=200,
    nvt_steps=100, npt_steps=100, production_steps=100,
    report_interval=10, platform_name="CPU",
)
```

That call executes minimisation and dynamics; it is not just a file writer.
The runner uses periodic handling when box vectors are present and records NPT
as skipped for nonperiodic input. Review stage logs, the final structure and
`openmm_summary.log` before increasing the run length. These short settings are
for a smoke test, not equilibration or convergence evidence.

**Solvated OpenMM with Amber force fields is notebooks 06A1 (polymer in water) and
06A2 (enzyme + polymer in water).** They build the system in OPC water with tleap
(ff19SB for the protein) and write staged OpenMM run files whose protocols match
the GROMACS notebooks 06B1 and 06B2.

## GROMACS with CHARMM-GUI inputs

For CGenFF/CHARMM36m systems prepared in CHARMM-GUI, use notebooks **06B1** (vacuum
check from a Ligand Reader & Modeler download) and **06B2** (Solution Builder package).

## GROMACS from GAFF2/AMBER files (conversion helpers)

The package can also convert GAFF2/AMBER files to GROMACS; the polymer-only benchmark
systems were prepared this way.

1. Use the converter/preparation helpers to convert AMBER parameters
   with ParmEd and create a self-contained dry run folder.
2. Validate topology includes, atom counts and the box relative to nonbonded
   cutoffs before any `grompp`/minimisation run.
3. Generate the explicit-water/ion preparation scripts from the dry
   inputs. Review solvent group names and force-field/water compatibility.
4. Run the generated solvation script deliberately, then check coordinate and
   topology atom counts and solvent/ion molecule counts.
5. Carry the prepared folder to [HPC execution](hpc.md).

The standard folders are under `examples/output/md_tests/P3HB_4/gromacs/`, with
`dry_polymer/` and `solvated_polymer/` branches. The preparation module also
contains a CHARMM-GUI-style staged workflow; that is distinct from an automated
enzyme-complex builder. Solvation scripts reset intermediate files in their
target folder, so keep original inputs and use a dedicated workflow folder.

Notebook links: [06A1–06B2](../notebooks.md#parameterisation-and-simulation).
API: {py:func}`iphasimulator.conversion_amber_to_gromacs.convert_amber_to_gromacs`,
{py:func}`iphasimulator.simulation_gromacs_runner.prepare_gromacs_run_folder`,
{py:func}`iphasimulator.simulation_gromacs_runner.validate_gromacs_coordinate_topology_counts`.
