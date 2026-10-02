# GAFF2 parameterisation

**Input:** an exported SDF. **Outputs:** MOL2, FRCMOD, PRMTOP, INPCRD, PDB and
command/timing logs. The existing helper runs AmberTools; it does not perform MD.

```{warning}
`parameterization_gaff2.py` currently has an [import-order blocker](../capabilities.md#import-blockers).
Resolve that source issue before attempting the commands below. This documentation
does not run AmberTools.
```

## 1. Check prerequisites and choose paths

Use the [installation guide](../installation.md) to install AmberTools. The helper
needs `antechamber`, `parmchk2` and `tleap` on `PATH`.

```python
from pathlib import Path
from iphasimulator.parameterization_gaff2 import ambertools_available, parameterize_gaff2

print("AmberTools found:", ambertools_available())
input_sdf = Path("examples/output/quickstart/P3HB_4/P3HB_4.sdf")
output_dir = Path("examples/output/md_tests/P3HB_4/gaff2")
```

## 2. Select the charge model deliberately

The default is `charge_method="abcg2"` (ABCG2). AM1-BCC (`charge_method="bcc"`)
is also supported; the project's polymer-only benchmark systems used `bcc`.
Faster/debug charge choices are not automatically interchangeable with a
validated production model.

When ready to run the external tools:

```python
outputs = parameterize_gaff2(
    input_sdf, output_dir,
    name="P3HB_4", residue_name="PHA",
    net_charge=0, charge_method="abcg2", verbose=True,
)
print(outputs.prmtop_path)
print(outputs.inpcrd_path)
```

## 3. Review the outputs

Inspect the antechamber, parmchk2, tleap and timing logs, atom types and charges.
The helper includes checks for failed commands and MOL2 charge-rounding residue;
successful file creation alone does not validate a force field for every PHA.
Continue with [engine preparation](simulation.md) using the matching topology
and coordinates.

The [05A notebook](../notebooks.md#parameterisation-and-simulation) has a disabled
execution switch and explains the branch into OpenMM and GROMACS. The 05B
CHARMM/CGenFF notebook is an **incomplete manual handoff**, not an automated
replacement for this route.

API: {py:func}`iphasimulator.parameterization_gaff2.parameterize_gaff2`.
