# Current capabilities and limitations

The teaching notebooks are complete workflow guides. Their reusable helpers import
in the installed `ipha_clean` environment; the old import-order blocker list no
longer describes this checkout. Completion of a tutorial does not validate a
force field, a production trajectory or a binding prediction.

| Area | Current implementation | Prerequisites / limits |
| --- | --- | --- |
| Construction and export | 01–04: curated monomers, side-chain/custom design, R checks, SDF/PDB export | RDKit and the installed package; inspect geometry before parameterisation |
| GAFF2 | 05A: MOL2/FRCMOD/PRMTOP/INPCRD and command/timing logs | AmberTools; default ABCG2 requires version ≥23; benchmark used AM1-BCC |
| Quick dry test | 05A_quick_check_openmm: load, minimise and short dynamics | 05A files + OpenMM; no solution-physics interpretation |
| Polymer in water | 06A: ParmEd conversion to GROMACS, CHARMM-style TIP3P/SOD/CLA | 05A inputs, ParmEd, GROMACS; explicit run/write flags |
| Enzyme–PHA reference | 05B → 06B: CGenFF checks and CHARMM-GUI Solution Builder preparation | Website steps manual; packaged ANC55/P3HB4 inputs run by default; reviewed docking and protonation for new systems |
| Optional enzyme–PHA | 05A → 06C: ff19SB/GAFF2/OPC, tleap + staged OpenMM | Docked complex and same-oligomer 05A files; [matched conditions and differences](workflows/optional_openmm.md); GAFF2/OPC is not validated for PHA |
| HPC | 07: reviewed SLURM execution, restart and benchmarking; YAML helper interface | Prepared run folder; cluster environment/directives must be edited |
| Polymer benchmark | 10: six-system GAFF2/GROMACS batch workflow | External engines and storage; its execution cell runs preparation and may submit with sbatch; separate from a single enzyme simulation |
| Docking inputs | 11: whole-polymer GRO-to-PDB export and manual HADDOCK records | Benchmark frame/topology or other reviewed docking source; no automated submission |
| Teaching analysis | 08/09: GROMACS preprocessing and Rg/end-to-end/SASA; 12: enzyme stability | Existing matching trajectories/topologies; machine-specific paths and selections must be replaced |
| Research contacts | Separate contact notebook/CLI under `src/md_simulation_scripts/enzyme_contacts/` | Matching periodic trajectory/topology; sampled proximity does not establish affinity or catalysis |
| Database construction | `PHAPolymerBuilder` and prepin/trimer route | Open Babel/AmberTools and manually supplied head/main/tail definitions; distinct from RDKit tutorials |
| Packing, GUI, script builders | Additional advanced interfaces | Not a unified validated production pipeline |

## Validation and remaining limits

The maintained `tests/` suite is the test entry point. Bare repository-wide pytest
also discovers a legacy `src/dan_example_scripts/test_openmm_script_builder.py`
whose constructor call uses an obsolete `root_dir` argument; it is outside these
teaching routes. Research datasets are not distributed as tutorial trajectories.
A notebook's saved outputs do not validate new inputs.

Schema/syntax checks and short execution can establish runnable preparation and
restart behaviour. They do not establish full equilibration, convergence or
PHA-specific accuracy. [06C](workflows/optional_openmm.md) retains deliberate
force-field and engine differences from 06B. GROMACS and OpenMM output formats
also differ; GROMACS analysis cells cannot consume DCD/CSV by changing filenames.

## Teaching and research material

Use the [notebook catalogue](notebooks.md) for names, prerequisites and run order.
The [contact validation record](enzyme_contacts.md) records a particular research
preview and its historical environment; it is separate from the current teaching
validation. Notebook docking tables and catalytic-residue notes are project
context, not validated package predictions.
