# Optional enzyme–PHA in water with OpenMM (06C)

Use `06C_optional_gaff2_openmm_pha_enzyme_in_water.ipynb` after 05A and a reviewed
**docked enzyme–PHA complex PDB**. It preserves the GAFF2/Amber route: ff19SB
protein, GAFF2 PHA (ABCG2 by default in 05A), OPC water and Na⁺/Cl⁻. 06B uses
CHARMM36m/CGenFF/CHARMM TIP3P in GROMACS. Matching macroscopic conditions does not
make these physical models identical or establish comparable binding predictions.

## Inputs, preparation and hand-off

1. Install the package and activate one environment containing RDKit, AmberTools
   (`antechamber`, `parmchk2`, `tleap`) and OpenMM. Start Jupyter in that environment.
2. Run 04 → 05A for the same oligomer as the complex. 06C reads
   `<SYSTEM>.gaff2.mol2`, `<SYSTEM>.gaff2.frcmod` and `<SYSTEM>.antechamber.sdf`
   from `examples/output/md_tests/<SYSTEM>/gaff2/`. Its default system is `P3HB_4`,
   matching 05A. Dry-polymer PRMTOP/INPCRD are not solvated complex inputs.
3. Set `COMPLEX_PDB`, `LIGAND_RESNAME` (normally `LIG`) and a fresh `OUTPUT_DIR`.
   Use notebook 11 plus manual HADDOCK when starting from the polymer benchmark;
   11 exports a whole polymer from a completed 10 run or a chosen starting frame.
   Another reviewed docking source is also acceptable.
4. Repair the protein first: resolve missing atoms/residues, alternate locations,
   termini, disulfides and any cofactors deliberately. Review protonation against
   the pH 7.0 CHARMM-GUI setup in 06B. The helper preserves named HID/HIE/HIP and
   converts HSD/HSE/HSP; generic HIS uses tleap's HIE. Optional `CHARMM_PROTEIN_ITP`
   imports histidine states by residue number for one matching protein chain.
   Other protonation choices are not automatically copied from CHARMM-GUI.
5. Enable `RUN_PREPARE_INPUTS`. The helper extracts protein heavy atoms, preserves
   chain boundaries and transfers the ligand pose into the 05A MOL2. Element/order,
   bond connectivity and all-R stereochemistry are checked; parameters and charges
   remain from 05A. The preparation directory must be new.
6. Enable `RUN_TLEAP`. Two tleap calls count waters, then neutralise and add the
   requested NaCl pairs. Inspect both logs. The new system directory contains
   `system.prmtop`, `system.inpcrd`, `system.pdb`, the copied inputs and build logs.
7. Keep the generated `run_openmm_md.py`, `protocol.json`, `run_step6_local.sh` and
   `run_hpc_equilibration_production.slurm` with those inputs. The runner checks
   particle counts, finite coordinates and periodic box dimensions before PME.
   OPC's fourth sites remain part of the topology/coordinates. The polymer residue
   name in `polymer_residues` must match the 05A MOL2 (normally `PHA`).
8. Use `RUN_OPENMM_TEST` for minimisation + 10 ps NVT + 10 ps NPT in `short_test/`.
   It tests numerical execution; it is not equilibration or scientific validation.
   Run full step 6.0 separately, inspect it, then use 07's OpenMM section to review
   and submit the cluster script. Edit environment activation and SLURM directives.

All notebook execution flags default to false. Run All checks paths and previews
text without writing a prepared complex. Existing input/output folders are protected;
keep new outputs under ignored `examples/output/` or external scratch storage.

## Final 06B–06C conditions

06B values below come from its Solution Builder instructions, the actual packaged
`src/iphasimulator/data/charmm_gromacs/step6.*.mdp` / `step7_production.mdp` files
and the example's `PROA.itp`/`LIG.itp` restraint definitions. 06C values come from
`ENZYME_SOLVATION_SETTINGS`, `ENZYME_POLYMER_IN_WATER_PROTOCOL` and its runner.

| Condition | 06B (GROMACS reference) | 06C (Amber/OpenMM) |
| --- | --- | --- |
| Protein / PHA / water | CHARMM36m / CGenFF / CHARMM TIP3P | ff19SB / GAFF2 / OPC; retained physical-model difference |
| Temperature / pressure | 303.15 K / 1 bar | 303.15 K / 1 bar |
| Salt | Neutralised, 0.05 M NaCl, SOD/CLA | Neutralised, 0.05 M NaCl, Na⁺/Cl⁻; pair count rounded from waters/55.5 M |
| Box | Rectangular, 30 Å edge distance fitted to protein; example 11.4 nm cube | Rectangular, 3.0 nm padding around complete complex; no `iso` option; dimensions/water counts need not match |
| Minimisation | Steepest descent; emtol 1000 kJ mol⁻¹ nm⁻¹; max 5000 steps | OpenMM L-BFGS; RMS force tolerance 1000 kJ mol⁻¹ nm⁻¹; max 5000 iterations; stopping criterion differs |
| NVT | 125 ps, 1 fs | 125 ps, 1 fs |
| NPT | 500 ps, 2 fs | 500 ps, 2 fs |
| Production | 200 ns, 2 fs | 200 ns, 2 fs |
| Restraints in minimisation/NVT/NPT | Protein N/CA/C/O + PHA heavy atoms: 400; other protein heavy atoms: 40 kJ mol⁻¹ nm⁻²; none in production | Same selections/strengths, periodic harmonic force; fixed original reference positions |
| Restraint reference under NPT | `refcoord_scaling = com` | Fixed original coordinates; MC volume move/reference treatment differs |
| Constraints | H bonds, LINCS; water SETTLE | H bonds, rigid OPC water, OpenMM constraint solver, tolerance 10⁻⁵; no hydrogen mass repartitioning |
| Electrostatics | PME, 1.2 nm real-space cutoff | PME, 1.0 nm cutoff, error tolerance 5×10⁻⁴ |
| Lennard-Jones | CHARMM force switch 1.0–1.2 nm | Amber 1.0 nm cutoff with default dispersion correction; no CHARMM force switch |
| Thermostat | v-rescale, SOLU/SOLV groups, tau 1 ps | LangevinMiddle, whole system, friction 1 ps⁻¹; stochastic dynamics differ |
| Barostat | Isotropic C-rescale, tau 5 ps, compressibility 4.5×10⁻⁵ bar⁻¹ | Isotropic Monte Carlo, attempt every 25 steps; no equivalent tau/compressibility |
| Centre-of-mass motion removal | Separate SOLU/SOLV groups every 100 steps | Whole-system OpenMM CMMotionRemover, every step |
| Trajectory interval: NVT / NPT / production | 5 / 100 / 100 ps, XTC | 5 / 100 / 100 ps, DCD |
| Energy/log interval: NVT / NPT / production | 1 / 2 / 2 ps | 1 / 2 / 2 ps, CSV |

Salt and padding were changed from the older 06C defaults (0.15 M and cubic
1.2 nm padding). NPT frames and CSV intervals were aligned to 06B. Existing run
folders retain their copied script and JSON; create a new run folder to use these
settings. Do not silently replace the protocol of an existing simulation.

OPC and Amber nonbonded treatment are retained deliberately. Applying CHARMM's
water parameters or force switch would be a consequential force-field change,
not an engine translation. GAFF2/OPC hydration-dependent properties remain
unvalidated for PHA. See the [OpenMM AMBER example](https://docs.openmm.org/latest/userguide/application/02_running_sims.html#using-amber-files)
for the preserved PME/HBonds route.

## Restart and analysis formats

Each MD stage saves `<stage>.chk` and `<stage>.restart.json`, along with final XML,
PDB, CSV and DCD. Resubmitting the same script loads each stage's checkpoint,
appends output and runs only its remaining steps; completed stages take zero new
steps. Input/physical-setting fingerprints reject changed topology, original
coordinates or protocol. A missing checkpoint with existing stage outputs fails
instead of overwriting them. Binary checkpoints require compatible hardware,
platform and OpenMM version; see [OpenMM checkpoint compatibility](https://docs.openmm.org/latest/api-python/generated/openmm.app.simulation.Simulation.html).

DCD append assumes outputs correspond to the checkpoint. After an interrupted
job, inspect and reconcile any frames/CSV rows written beyond the checkpoint
before resubmission; this runner does not trim them automatically. Keep the
original folder and checkpoint together. A portable XML state is a stage hand-off,
not an exact replacement for an interrupted checkpoint's random state.

06C produces DCD/PDB/CSV. Notebooks 08 and 12 expect GROMACS XTC/TPR/EDR and are
not direct consumers of these files. Use an appropriate OpenMM/MDTraj analysis
with matching topology and periodic coordinates; review virtual-site selection.
Research enzyme-contact analysis is described [separately](analysis.md).

## Validation scope

The validation record for this update distinguishes schema/syntax checks, default
notebook execution, actual tleap preparation and shortened OpenMM execution.
Short runs use temporary directories and shortened stage durations; the published
production protocol is never executed by documentation builds. Finite energy and
successful restart demonstrate execution, not equilibrium, convergence, binding
or force-field accuracy.

Validation on 5 October 2026 used Python 3.11, RDKit 2025.09.5, OpenMM 8.5.1
and AmberTools from `ipha_clean` (the tools were placed on PATH before execution):

- All 06C default code cells executed with preparation/MD flags false.
- Actual pose/protein preparation and two tleap calls built a synthetic Ala–Ala
  protein plus the existing all-R P3HB tetramer in OPC: **55,740 particles**,
  13,936 pre-ion-insertion waters and 13 salt pairs at the requested 0.05 M;
  OPC virtual sites loaded correctly. This is a mechanics fixture, not an enzyme
  docking validation. Original GAFF2 input hashes were unchanged.
- Actual CPU execution used at most 50 minimisation iterations, **0.01 ps NVT**,
  **0.10 ps NPT** (including Monte Carlo volume attempts) and **0.02 ps production**
  split across a checkpoint restart. Energies stayed finite. DCD frames appended
  from 5 to 10; completed NVT/production stages took zero further steps. A changed
  temperature was rejected for checkpoint continuation.
- Durations and reporting intervals were shortened in a temporary protocol for
  that check; the 200 ns production protocol was inspected, not executed. The
  notebook’s separate 10 ps + 10 ps short-test preset was not run at full length.
- The maintained test suite passed; notebook schemas/code syntax, documentation
  quickstart exports, warning-free Sphinx build, and built-site link/download checks
  were also checked. External research trajectories were not executed.
