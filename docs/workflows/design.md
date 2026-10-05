# Polymer design

**Inputs:** one monomer code, alkyl side-chain length, or custom R 3-hydroxy-acid
SMILES, plus a repeat count. **Output:** a named RDKit molecule.



## 1. Choose one design mode

```python
from iphasimulator.workflows import PolymerDesign, design_polymer

request = PolymerDesign(common_name="3HB", degree=4)
name, molecule = design_polymer(request)
```

To specify a linear alkyl side chain instead, replace the request with
`PolymerDesign(side_chain_carbons=5, degree=4)`. For a custom monomer, use
`PolymerDesign(custom_monomer_smiles="COC[C@H](O)CC(=O)O", degree=4, name="ether_pha")`.
Supply exactly one design mode.

The curated set includes 3HB, 3HV, 3HHx, 3HHep, 3HO, 3HN, 3HD and 3HDD.
`supported_polymer_table()` exposes the registry and naming fields. The common
monomer route accepts R or S, while the side-chain and custom routes construct
R-PHA. Custom input is restricted to an R 3-hydroxy-acid backbone; it is not a
general arbitrary-polymer builder.

## 2. Validate and export

Follow the [quick-start validation and export steps](../quickstart.md). Confirm
the repeat count and stereochemistry before creating MD parameters. Use the
[design notebook](../notebooks.md#construction-and-design) to compare requests
interactively.

## 3. If you use the AmberTools database route

`build_pha.PHAPolymerBuilder` offers a separate workflow:

1. `parameterise_trimer(...)` generates trimer parameter files and registers
   residue codes.
2. Supply the manual prepgen definition files under
   `structure_database/PHA_types/<PHA_type>/input/`.
3. `generate_polymer_prepins(...)` prepares head/mainchain/tail residue files.
4. `build_PHA_polymer(PHA_type, length)` assembles the desired polymer in tleap.

The class also implements `build_PHA_copolymer(...)` with patterned or random
sequences from **existing** monomer-unit prepins. Each requested PHA type needs
its residue/parameter files first. This does not mean the RDKit `build_pha_chain`
interface supports mixed sequences, or that arbitrary comonomers are already
parameterised. The database route remains a separate advanced interface.

This route requires Open Babel, AmberTools and a repository-root environment
because some imports are written as `src.iphasimulator`. It is not substituted
silently for the RDKit tutorial. See the generated
{py:class}`iphasimulator.build_pha.PHAPolymerBuilder` reference and the
`src/dan_example_scripts` directory in the repository for its developing examples.
