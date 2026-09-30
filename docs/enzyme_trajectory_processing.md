# Standalone enzyme–PHA trajectory preparation

The preparation workflow now uses **Bash and GROMACS**, with two templates in
`src/md_simulation_scripts/trajectory_preparation/`. Copy them into a simulation's
`analysis/` folder, edit the plain `KEY=value` instructions, then run
`bash process_trajectory.sh`. No Python package or preparation notebook is required.

- {download}`process_trajectory.sh <../src/md_simulation_scripts/trajectory_preparation/process_trajectory.sh>`
- {download}`instrcution.txt <../src/md_simulation_scripts/trajectory_preparation/instrcution.txt>`

The complete setup, settings, atom groups, error handling, VMD inspection and
contact-analysis instructions are maintained in the README below.

```{include} ../src/md_simulation_scripts/README.md
:start-after: "# Simulation analysis scripts"
:end-before: "## 5. Run enzyme-contact analysis"
:relative-docs: docs/
:relative-images:
```

For Python analysis, continue with the [enzyme-contact guide](enzyme_contacts.md).
The complete [script-folder README](https://github.com/MMLabCodes/iPHAsimulatorV2/blob/main/src/md_simulation_scripts/README.md) also covers the moved entry points.
