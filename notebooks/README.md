# Tutorial Notebooks

## Design goal

iPHASimulator v2 serves two purposes:

1. **A package** that builds PHA oligomers with RDKit and prepares molecular dynamics (MD) systems from them.
2. **The author's research**: MD of enzyme + PHA systems.

The notebooks follow one shared build stage (01–04), then split into two force-field
routes. The routes are kept separate: each uses one consistent force-field family from
polymer to protein to water and ions, so parameters from the two families are never mixed.

## The two routes

| | Amber/OpenMM | CHARMM/GROMACS |
| --- | --- | --- |
| Polymer | GAFF2, ABCG2 charges by default | CGenFF via CHARMM-GUI Ligand Reader & Modeler |
| Protein | ff19SB | CHARMM36m |
| Water / ions | OPC (with the ion parameters loaded by `leaprc.water.opc`) | CHARMM TIP3P + SOD/CLA |
| Engine | OpenMM | GROMACS |
| Parameters | `05A_amber_gaff2_parameterisation.ipynb` | `05B_charmm_cgenff_parameterisation.ipynb` |
| Dry polymer check | `06A_openmm_dry_polymer.ipynb` | `06C_gromacs_dry_polymer.ipynb` |
| Solvated system (polymer in water; polymer + protein) | `06B_openmm_solvated_system.ipynb` | `06D_gromacs_solvated_system.ipynb` |
| HPC | `07_hpc_workflows.ipynb`, Amber/OpenMM section | `07_hpc_workflows.ipynb`, CHARMM/GROMACS section |

The enzyme–PHA production simulations (GK13/ANC45 × P3HO_4/P3HB_4) used the **CHARMM/GROMACS** route.
Their input structures were iPHASimulator's R-configured SDF/PDB files (for example
`PHB4_R.sdf`); CGenFF assigned all atom types, charges and parameters, so no GAFF2
charges enter them.

```mermaid
flowchart LR
    S["01–04 shared<br/>RDKit build / validate / export"]
    S --> A5["05A GAFF2 (ABCG2)"]
    S --> C5["05B CGenFF (CHARMM-GUI)"]
    subgraph RA["Amber / OpenMM"]
        A5 --> A6["06A OpenMM dry polymer"]
        A5 --> A7["06B OpenMM solvated system"]
    end
    subgraph RC["CHARMM / GROMACS"]
        C5 --> C6["06C GROMACS dry polymer"]
        C5 --> C7["06D GROMACS solvated system<br/>(CHARMM-GUI Solution Builder)"]
    end
    A7 --> H["07 HPC workflows"]
    C7 --> H
    H --> X["08–12 analysis and tools"]
```

## Notebooks

Shared build stage:

1. `01_examples_pha_oligomers.ipynb`: built-in PHA oligomer generation examples.
2. `02_design_polymer_for_user_request.ipynb`: select or define a PHA target.
3. `03_validate_and_visualize.ipynb`: validate generated oligomers and inspect structures.
4. `04_export_structures.ipynb`: export validated oligomers to SDF/PDB.

Amber/OpenMM:

5. `05A_amber_gaff2_parameterisation.ipynb`: AmberTools GAFF2 parameterisation (ABCG2 charges by default).
6. `06A_openmm_dry_polymer.ipynb`: dry (vacuum) OpenMM check of the GAFF2 polymer. A sanity check, not a physical result.
7. `06B_openmm_solvated_system.ipynb`: polymer in OPC water via tleap, optionally with a posed ff19SB protein; short OpenMM test (minimise, 10 ps NVT, 10 ps NPT) and restartable production files for HPC.

CHARMM/GROMACS:

8. `05B_charmm_cgenff_parameterisation.ipynb`: CGenFF parameters through CHARMM-GUI Ligand Reader & Modeler and the Solution Builder handoff.
9. `06C_gromacs_dry_polymer.ipynb`: dry (vacuum) GROMACS check of the CGenFF polymer, from the Ligand Reader's native `gromacs/` files.
10. `06D_gromacs_solvated_system.ipynb`: prepare and validate a CHARMM-GUI Solution Builder GROMACS package (polymer in water; polymer + protein), as used for the enzyme–PHA runs.

Execution and analysis:

11. `07_hpc_workflows.ipynb`: HPC execution, SLURM submission, restart continuation, benchmarking and performance tuning; Amber/OpenMM and CHARMM/GROMACS sections.
12. `08_trajectory_preprocessing.ipynb`: GROMACS trajectory analysis (CHARMM/GROMACS; example data from the archived hybrid route). PBC reconstruction, centering with reusable `[ center ]` index groups, compact wrapping, optional fitting and representative frames.
13. `09_basic_polymer_analysis.ipynb`: GROMACS trajectory analysis (CHARMM/GROMACS; example data from the archived hybrid route). Rg, end-to-end distance and SASA from the centered trajectory.
14. `10_batch_md_benchmark.ipynb`: launcher and progress checker for the six-system polymer-only benchmark, which uses the archived hybrid route.
15. `11_PHA_Enzyme_Docking.ipynb`: prepares PHA oligomer PDB inputs and job notes for manual HADDOCK docking.
16. `12_enzyme_polymer_stable_analysis.ipynb`: stability diagnostics for one enzyme–polymer GROMACS production run (total energy, protein backbone RMSD, polymer RMSD relative to the protein).

## Archive

[`archive/`](archive/README.md) keeps the old 06B, 06C and 06D notebooks, renamed `hybrid_06B_…`, `hybrid_06C_…` and `hybrid_06D_…` so they are not confused with the current Amber/OpenMM and CHARMM/GROMACS notebooks. They form a
hybrid GAFF2 → GROMACS route with CHARMM-style TIP3P/SOD/CLA, which was used for the
six polymer-only benchmark systems (notebook 10, AM1-BCC charges) and the P3HB_4_01
example trajectory analysed in 08/09. They are kept for provenance and superseded by
the Amber/OpenMM and CHARMM/GROMACS routes.

## Conventions

These tutorials are written for users who may not be computational specialists.
They explain what each step means and keep editable input cells small.

Reusable Python logic belongs in `src/iphasimulator`. Notebook cells should guide
the workflow, not duplicate chemistry or simulation code.
