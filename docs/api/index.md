# Python API reference

These pages are generated from the package's Python signatures and docstrings
using Sphinx AutoAPI. Static parsing does not import the package, load molecular
data, execute notebook cells, or start simulation tools.

Start with `naming`, `build`, `export` and `workflows.design` for the RDKit route;
`analysis_contacts` and `trajectory_*` for completed simulations. The
`build_pha`, `pha_filepath_manager`, `build_single_PHA_systems`, `pha_melt_builder`,
`openmmscript_builder` and `sw_openmm` modules expose additional database and MD
interfaces. Their presence in this reference does not remove external-tool or
input-data prerequisites. See [current limits](../capabilities.md).

```{toctree}
:maxdepth: 3

../autoapi/iphasimulator/index
../autoapi/iphasimulator/visualisation/visualiser/index
../autoapi/iphasimulator/data/amber_openmm/run_openmm_md/index
```
