# Current capabilities and limitations

This page is based on the code and notebook sources in this checkout. **Implemented**
means a concrete helper exists; it is not a claim that a production MD protocol
has been validated. Historical notebook outputs may predate current source errors.

| Area | Current implementation | Limitations |
| --- | --- | --- |
| Naming | Canonical monomer, polymer, oligomer and multi-chain names | A naming helper does not build a mixed polymer |
| RDKit construction / export | Curated monomers, alkyl side-chain length, custom R 3-hydroxy acid; PDB/SDF export | Currently blocked by import-order errors |
| Database construction | `PHAPolymerBuilder`, residue/path managers; trimer parameterisation, manual prepgen definitions, tleap assembly; patterned/random copolymers from existing prepins | Requires Open Babel and AmberTools; uses repository-root `src.iphasimulator` imports; separate from the RDKit notebooks |
| GAFF2 parameterisation | SDF → MOL2/FRCMOD/PRMTOP/INPCRD, charge handling and logs | Import blocker; external AmberTools required; default charge method is currently `abcg2` |
| OpenMM AMBER runner | Minimisation, NVT, NPT when periodic, production and logs | Import blocker; requires valid input parameters and OpenMM |
| GROMACS preparation | ParmEd conversion; dry/solvated input folders, MDP/scripts and validation helpers | Import blockers; production settings and prepared systems still need validation |
| HPC | YAML planning, SLURM generation, execution and restart guidance | Workflow imports depend on blocked modules; cluster settings must be edited |
| Trajectories | GROMACS merge/check, centering, wrapping, optional fitting, representative frames | Requires GROMACS and matching simulation inputs |
| Enzyme–PHA contacts | Reusable minimum-image distances, both heavy-atom modes, occupancies, plots and provenance | Proximity only; no affinity/catalysis/convergence inference |
| Basic polymer analysis | Notebook Rg, end-to-end distance and SASA using MDTraj | Notebook-level workflow; requires meaningful selections and preprocessed coordinates |
| Enzyme stability | Notebook energy and RMSD diagnostics | System-specific paths/selections; manual interpretation |
| Docking preparation | Notebook GRO-to-PDB conversion and manual job records | No HADDOCK submission or validated complex builder; inspect/extract the polymer first |
| CHARMM/CGenFF with GROMACS | 05B documents the CHARMM-GUI route with penalty and stereocentre checks; 06B1/06D prepare and check GROMACS folders from CHARMM-GUI downloads | CGenFF and CHARMM-GUI steps are manual; minimisation runs only on request |
| Solvated OpenMM (Amber force fields) | Notebook 06A2: tleap GAFF2 + OPC (optional ff19SB protein), short OpenMM test, production script | tleap/OpenMM are mocked in tests; GAFF2 + OPC is not validated for PHA |
| Packing / script builders | Additional single-chain, melt and OpenMM script classes | Advanced, separate interfaces; not a fully unified end-to-end tutorial |
| APO comparisons, catalytic geometry, repeat-unit contacts, ML/DFT | Not part of the reusable contact workflow | Planned or outside the documented v2 scope |

## Import blockers

The following files have an additional standalone comment string before
`from __future__ import annotations`, which Python rejects at import time:

- `build.py`, `monomers.py`, `stereochemistry.py`, `export.py`
- `parameterization_gaff2.py`
- `conversion_amber_to_gromacs.py`, `simulation_gromacs_runner.py`
- `simulation_openmm_amber_runner.py`, `system_builder_packmol.py`

Other modules may be affected transitively, including `iphasimulator.workflows`
and the configured execution runner. Even a runner's `--dry-run` needs its imports
to succeed. These source files are intentionally unchanged by the documentation
implementation. The [API reference](api/index.md) uses static parsing and can
still show their signatures and docstrings.

## Choosing a construction route

The [P3HB_4 tutorial](quickstart.md) follows notebooks 01–04 and the existing
RDKit helper interface. `build_pha.PHAPolymerBuilder` is an additional AmberTools
database route; it requires prepared trimer data and manually supplied head,
mainchain and tail definitions. Do not treat the two routes' parameters or files
as interchangeable. See [polymer design](workflows/design.md).

## Reading older material

The [repository README](project_readme.md), [developer guide](developer_guide.md)
and [contact validation record](enzyme_contacts.md) are preserved. Their historical
status statements should be read alongside this page. The documentation does not
promote notebook claims about enzyme binding or catalytic residues into verified
package capabilities.
