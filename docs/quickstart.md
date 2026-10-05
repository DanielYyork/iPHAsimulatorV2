# Quick-start: construct, validate and export P3HB_4

**Goal:** create one linear oligomer with four 3HB repeat units and export
`P3HB_4.sdf` and `P3HB_4.pdb`. This is structure preparation, not an MD simulation.



This tutorial combines the existing construction, validation and export notebooks
([01, 03 and 04](notebooks.md#construction-and-design)). Run these steps in order from the repository root using the
[core environment](installation.md).

## 1. Construct the oligomer

```python
from rdkit import Chem
from iphasimulator.naming import oligomer_name
from iphasimulator.build import build_pha_chain

monomer = "3HB"
repeat_units = 4
name = oligomer_name(monomer, repeat_units)
pha = build_pha_chain(monomer, repeat_units, stereochemistry="R")
pha.SetProp("_Name", name)
print(name, pha.GetNumAtoms())
```

Expected chemistry: `P3HB_4`, 25 atoms with hydrogens implicit. The molecule is a
hydroxyl/carboxyl-terminated oligomer. `build_pha_chain` already checks its core
chemistry while constructing it.

## 2. Validate the structure

```python
from iphasimulator.build import validate_pha_chain
from iphasimulator.monomers import get_monomer
from iphasimulator.stereochemistry import validate_r_chiral_centres

validate_pha_chain(pha, get_monomer(monomer), repeat_units, "R")
validate_r_chiral_centres(pha, expected_count=repeat_units)
print(Chem.FindMolChiralCenters(pha, includeUnassigned=True))
print(Chem.MolToSmiles(pha, isomericSmiles=True))
```

These functions raise an exception on failure. For this example the invariants
are four R stereocentres, three linking ester bonds and four side-chain carbon
atoms. Atom indices in the printed stereocentre list are RDKit indices, not
residue numbers.

## 3. Export PDB and SDF

```python
from pathlib import Path
from iphasimulator.export import to_pdb, to_sdf

output_dir = Path("examples/output/quickstart/P3HB_4")
output_dir.mkdir(parents=True, exist_ok=False)
to_sdf(pha, output_dir / f"{name}.sdf")
to_pdb(pha, output_dir / f"{name}.pdb")
print(output_dir.resolve())
```

Use a new output directory for each attempt: the export helpers themselves do
not protect existing files against replacement. They add explicit hydrogens,
generate a 3D conformer with RDKit ETKDG, and run a short UFF optimisation before
writing. This is a starting structure, not an equilibrated conformation or a
parameterised force field.

## 4. Inspect the output

Open the PDB in a molecular viewer and confirm the chain and terminal groups.
Use the SDF as the input to [GAFF2 parameterisation](workflows/gaff2.md).
For a different PHA, continue with [polymer design](workflows/design.md).

API: {py:func}`iphasimulator.build.build_pha_chain`,
{py:func}`iphasimulator.build.validate_pha_chain`,
{py:func}`iphasimulator.export.to_sdf`, {py:func}`iphasimulator.export.to_pdb`.
