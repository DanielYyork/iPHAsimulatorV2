# Simulation analysis scripts

Trajectory preparation needs only **Bash and GROMACS**. The research analysis
and enzyme-contact notebooks use Python analysis dependencies.

## Parameterisation provenance

The GK13/ANC45 × P3HO4/P3HB4 production simulations analysed here used the CHARMM/GROMACS route:
CHARMM36m + CGenFF (via CHARMM-GUI Ligand Reader & Modeler, from iPHASimulator's
R-configured SDF/PDB) + CHARMM TIP3P + SOD/CLA in GROMACS. No GAFF2 charges enter
them.

```text
src/md_simulation_scripts/
├── README.md
├── trajectory_preparation/
│   ├── process_trajectory.sh
│   └── instrcution.txt
├── enzyme_contacts/
│   ├── enzyme_contacts.ipynb
│   ├── run_enzyme_contacts.py
│   └── GK13_P3HO_4.yaml
└── enzyme_pha_analysis/
    └── enzyme_pha_analysis.ipynb
```

## 1. Prepare a new simulation's analysis folder

Keep the combined XTC in the simulation folder. Copy the two templates and
matching production TPR, viewing GRO and index into a **new** `analysis/` folder.
For example, after setting your repository and simulation paths:

```bash
REPO=/path/to/iPHASimulator_v2
SIM=/path/to/GK13_P3HO_4_gromacs
mkdir "$SIM/analysis"
cp "$REPO/md_simulation_scripts/trajectory_preparation/process_trajectory.sh" "$SIM/analysis/"
cp "$REPO/md_simulation_scripts/trajectory_preparation/instrcution.txt" "$SIM/analysis/"
cp "$SIM/step7_production.tpr" "$SIM/analysis/"
cp "$SIM/step6.2_npt.gro" "$SIM/analysis/"
cp "$SIM/index.ndx" "$SIM/analysis/"
cd "$SIM/analysis"
```

Use filenames from **your** system. The production TPR must describe the same
atoms in the same order as the combined XTC; the GRO is for viewing, and must also
match that order. Matching filenames or atom counts alone do not prove this.
If `analysis/` already exists, inspect its contents and preserve previous files
before copying templates into it. These examples are setup instructions, not an
instruction to replace your existing settings.

```text
simulation/
├── production_combined_1us.xtc       # raw combined trajectory stays here
└── analysis/
    ├── process_trajectory.sh
    ├── instrcution.txt
    ├── step7_production.tpr
    ├── step6.2_npt.gro
    └── index.ndx
```

## 2. Edit instrcution.txt

This spelling is retained for compatibility. It is **plain KEY=value data, not
YAML**. Do not use quotes, spaces around `=`, or inline comments; paths containing
spaces are otherwise supported. Settings are read as data, never sourced as code.
Paths resolve from the script/instruction folder, including when the script is
invoked from another working directory.

```text
INPUT=../production_combined_1us.xtc
TPR=step7_production.tpr
INDEX=index.ndx
CENTER_GROUP=SOLU
OUTPUT_GROUP=SYSTEM
STRIDE=100
PROCESSED=processed.xtc
PREVIEW=processed_every100.xtc
GMX=gmx
```

Change the input filenames and groups for a new system. `GMX` is an executable
name or path (for example `gmx_mpi`), not a command with extra arguments. Output
settings must be distinct `.xtc` filenames inside `analysis/`.

In GK13–P3HO₄, the inspected index contains:

| Group | Atoms | Use |
| --- | ---: | --- |
| `SOLU` | 3,800 | Enzyme (3,701 atoms) + `LIG` (99 atoms); default centre |
| `SOLV` | 133,560 | Solvent and ions |
| `SYSTEM` | 137,360 | Entire simulation; default output |

Use actual names from your index, not guessed group numbers. **Retain the full
system in `OUTPUT_GROUP`**, so the same TPR/index also match the processed XTC in
step 2. A subset needs its own matching topology/index and is outside this simple
two-pass template. There is no automatic choice of centring group or fitting.
For enzyme-only centring, explicitly create an enzyme group in your copied index
first; `Protein` and protein Cα groups are not present in this example's index.
The previously inspected protein has 245 Cα atoms, but this script does not
create selections from Python expressions.

## 3. Check and run

From the analysis folder:

```bash
bash process_trajectory.sh --dry-run
bash process_trajectory.sh
```

The dry run checks settings, executable availability, input file readability and
output conflicts, then prints the two commands and selections. It executes no
GROMACS command and writes no files. It cannot verify trajectory integrity,
topology compatibility or index-group contents; check those for your system.

The normal run performs, in this order:

1. Read the combined XTC from the parent folder; use the production TPR and index
   with `-pbc mol -ur compact -center`. Select `CENTER_GROUP`, then `OUTPUT_GROUP`.
   Write **every processed frame** to `processed.xtc`, with no frame reduction or
   fitting in this step.
2. Read `processed.xtc` and write every 100th frame to `processed_every100.xtc`
   using `-skip 100` (or the configured `STRIDE`). This keeps zero-based frames
   0, 100, 200, …; it is a frame stride, not a time interval.

The script stops on errors and refuses existing output files, including symlinks.
Rename earlier outputs or choose new `PROCESSED`/`PREVIEW` names to rerun. A failed
run can leave a partial file: inspect it before deciding what to keep. Run only
one preparation job per analysis folder at a time. Raw inputs are not modified.
The full processed XTC may be roughly as large as the raw file; allow space for
both full and reduced outputs.

## 4. Analyse the full trajectory; view the reduced one

Use **`processed.xtc` for analysis**. In MDAnalysis you can sample without creating
another stored trajectory:

```python
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

u = mda.Universe("step7_production.tpr", "processed.xtc")
enzyme = u.select_atoms("protein and not name H*")
pha = u.select_atoms("resname LIG and not name H*")
for ts in u.trajectory[::100]:
    distances = distance_array(enzyme.positions, pha.positions, box=ts.dimensions)
    # Consume the sampled distances here (MDAnalysis length units are Angstrom).
```

This is a simple sampling/PBC illustration; the contact workflow below provides
more complete heavy-atom selection, summaries and plots. **Enzyme–PHA distances
must account for periodic boundaries using each frame's box**, even after
centring. Centring does not bind a separated ligand to the enzyme or guarantee
that the displayed pair shares the nearest periodic image. Do not apply rotational
fitting and then calculate periodic distances using an unrotated box.

In VMD, load the matching `step6.2_npt.gro`, then add
**`processed_every100.xtc` to the same molecule**. The GRO's initial coordinates
need not be the first processed frame; start viewing the added XTC frames (or
remove the initial GRO frame). Display `protein or resname LIG` for the pair, or
`resname LIG` for the PHA. Scrub through the beginning, middle and end, checking
molecular continuity, periodic images and the chosen centre. Visual inspection
is a manual step; running this script does not perform it.

Do not infer time coverage from `1us` in a filename. Inspect the actual times.
Do not concatenate the combined XTC with continuation files it already contains.

## 5. Run enzyme-contact analysis

Install the repository's analysis dependencies in your Python environment:

```bash
python -m pip install -e '.[analysis]'
```

Open [enzyme_contacts/enzyme_contacts.ipynb](enzyme_contacts/enzyme_contacts.ipynb),
set `PROJECT_DIR`, `TPR_PATH`, `TRAJECTORY_PATH` and `OUTPUT_DIR`, then run its steps.
For the workflow above, point `TPR_PATH` at `analysis/step7_production.tpr` and
`TRAJECTORY_PATH` at **`analysis/processed.xtc`**. The moved notebook retains its
previous raw-input defaults and saved outputs; those historical outputs do not
represent a new run on your processed trajectory.

Alternatively, edit [enzyme_contacts/GK13_P3HO_4.yaml](enzyme_contacts/GK13_P3HO_4.yaml).
Set `topology` and `trajectory` to the same TPR and full processed XTC. Absolute
paths are simplest; relative paths resolve against **this YAML's directory**.
From the repository root:

```bash
python src/md_simulation_scripts/enzyme_contacts/run_enzyme_contacts.py src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml
```

Or, from `src/md_simulation_scripts/enzyme_contacts/`:

```bash
python run_enzyme_contacts.py GK13_P3HO_4.yaml
```

The default is a 31-frame preview of the selected window. Use `--full` to analyse
all configured samples; `--sample-interval-ns` controls time sampling. This Python
interface is separate from the Bash `STRIDE`, which only controls the VMD file.
The moved YAML preserves its original raw inputs and the ignored
`examples/output/enzyme_contacts/` result location. Each contact run creates a
new result directory. See [the contact guide](../../docs/enzyme_contacts.md).

## 6. Compare the four enzyme–PHA research systems

Open [enzyme_pha_analysis/enzyme_pha_analysis.ipynb](enzyme_pha_analysis/enzyme_pha_analysis.ipynb)
using the `ipha_clean` environment. It analyses GK13/ANC45 with P3HO₄/P3HB₄, using
each system's raw-folder `production_combined_1us.edr` and `step7_production.tpr`
and the verified **`analysis/processed_protein_centered.xtc`**. Preparation has
already been completed for these research systems; the notebook does not repeat it.

Set `STRIDE` in the first code cell (initially 100) and run all cells. `START_NS=0` and
`END_NS=1000` restrict every figure and CSV to the same 0–1000 ns window. The first
sampled frame is the backbone RMSD reference. Every EDR sample within this window
is retained independently of trajectory stride. The plots follow notebook 12: a
red dashed cumulative energy mean and a thicker RMSD rolling mean.
`ROLLING_WINDOW_FRAMES=10` controls the trend over analysed samples for both RMSD
and minimum distance; their rolling means are also exported in the CSVs.

Each simulation's own `analysis/` folder receives three PNGs and matching CSVs:

- `<system>_research_0_1000ns_total_energy`: total energy and cumulative sample mean, kJ/mol.
- `<system>_research_0_1000ns_protein_backbone_rmsd`: N/CA/C backbone RMSD after alignment
  to the first analysed frame, Å, plus its rolling mean.
- `<system>_research_0_1000ns_enzyme_pha_min_distance`: minimum enzyme–PHA heavy-atom
  distance with PBC, Å, calculated from unaligned coordinates and same-frame boxes,
  plus a rolling mean of the per-frame minima.

Time is in ns. The notebook defaults to `RUN_ANALYSIS=False`: Run All redraws plots
from existing CSVs. For plot edits after restarting the kernel, run **1. Paths**,
then **5. Plot settings**, then the separate **Figure 1**, **2**, or **3** cell.
Edit plot colours/labels/axes directly in that cell. Plots keep left/bottom borders
and hide top/right borders. `PLOT_ROLLING_WINDOW` controls the displayed mean;
plotting replaces PNGs only, leaving CSVs unchanged. Set `SAVE_FIGURES=False` for
preview only. For fresh calculations, set `RUN_ANALYSIS=True` and choose a new
`OUTPUT_TAG` (or explicitly enable `OVERWRITE`).

An energy extraction log and a JSON summary record inputs,
selections, sampling and independent numerical checks. Existing results are
protected by default; change `OUTPUT_TAG` for a new run, or explicitly enable
`OVERWRITE`. Check the figures before interpreting drift, stability or contacts;
the notebook does not establish binding or convergence. It uses teaching notebook
12 only as a reference and does not modify teaching notebooks 01–12.

## What this replaces

The GK13 preparation notebook and `examples/enzyme_trajectory_GK13_P3HO_4.yaml`
are superseded by the two standalone templates. Use this script with the plain
instruction file; do not pass it to the older Python YAML workflow. Existing
shared package helpers and `notebooks/08_trajectory_preprocessing.ipynb` remain
available for their separate workflows. The templates do not depend on them.

## CGenFF / CHARMM-GUI provenance notes

Moved here from `notebooks/05B_charmm_cgenff_parameters.ipynb` (2026-10-01), when 05B was
simplified for non-computational users.

### Production runs and route choice

The enzyme–PHA production simulations (GK13/ANC45 × P3HO_4/P3HB_4) used the CHARMM/GROMACS route:
polymer **CGenFF** via CHARMM-GUI Ligand Reader & Modeler · protein **CHARMM36m** · water/ions
**CHARMM TIP3P + SOD/CLA** · engine **GROMACS**. CGenFF assigned all atom types, charges and
parameters, so no GAFF2 charges enter them.

Use the CHARMM/GROMACS route for protein–polymer systems that should stay in one consistent
CHARMM force-field family. Use the Amber/OpenMM route (`05A_amber_gaff2_parameters.ipynb` →
`06_optional_gaff2_openmm_pha_enzyme_in_water.ipynb`) for the Amber family (GAFF2 / ff19SB / OPC)
in OpenMM.

### CHARMM-GUI tools

- [CHARMM-GUI](https://www.charmm-gui.org)
- [Ligand Reader & Modeler](https://www.charmm-gui.org/?doc=input/ligandrm): CGenFF atom types,
  charges and parameters for a non-standard molecule.
- [Solution Builder](https://www.charmm-gui.org/?doc=input/solution): solvated box, ions and
  engine-specific inputs (GROMACS here).
- [Polymer Builder](https://www.charmm-gui.org/?doc=input/polymer): Choi et al., *J. Chem. Theory
  Comput.* 2021, [doi:10.1021/acs.jctc.1c00169](https://doi.org/10.1021/acs.jctc.1c00169). PHA
  monomer support not verified.

### Ligand Reader & Modeler input and output

Input: the R-configured SDF exported by `04_export_structures.ipynb`, one per oligomer:
`PHB4_R.sdf`, `PHB8_R.sdf`, `PHO4_R.sdf`, `PHO8_R.sdf`, `PHDD4_R.sdf`, `PHDD8_R.sdf`. Check the
structure shown in Marvin JS: CHARMM-GUI builds the topology from that drawing, so the chirality
(all backbone stereocentres R), the bond orders and the explicit hydrogens must be correct.

| File (in `lig/`) | Content |
| --- | --- |
| `lig.rtf` | CGenFF topology, with a charge penalty on every atom |
| `lig_g.rtf` | the same topology with charge groups |
| `lig.prm` | CGenFF parameters by analogy, with a penalty on every parameter |
| `lig.log`, `ndihe.str` | CGenFF log and extra dihedrals |

`lig_g.rtf` and `lig.prm` are what Solution Builder asks for in 06B ("Upload CHARMM top & par for
hetero chain").

### Open items

- TODO(Zhiwen): CGenFF version used for the production systems (read it from the `lig.rtf` /
  `lig.prm` header).
  *Evidence from the PHO4 Ligand Reader download (Apr 2026):* CGenFF program version 4.0
  (released June 2024), for CGenFF topology and parameter files version 5.0.
- TODO(Zhiwen): which rtf (`lig.rtf` or `lig_g.rtf`) the GROMACS `LIG.itp` charges came from.
  *Evidence from the PHO4 Ligand Reader download (Apr 2026):* `lig.rtf` and `lig_g.rtf` give
  identical per-atom charges (99/99 atoms; `lig_g.rtf` only regroups them), and `gromacs/LIG.itp`
  matches both (net charge 0). Confirm for the production systems.
- TODO(Zhiwen): how the enzyme–polymer starting pose was obtained (docking; the complexes are
  named `..._pose2_complex.pdb`).
  *Evidence from the Solution Builder downloads:* the GK13–PHO4 builds start from
  `GK13_PHO4_pose2_complex.pdb` (docking pose 2).
- *Production provenance (checked 2026-09-30):* the GK13_P3HO_4 production input is
  byte-identical to CHARMM-GUI job **8214536317** (22 June 2026, 137,360 atoms), not to the later
  rebuild 8221445616 (23 June, 137,363 atoms). Both builds have the same LIG types and charges,
  and all four LIG stereocentres are R in both.
- *Caution:* the April PHO4 Ligand Reader download (job 7687539766) used for the version and
  charge evidence above contains the **S** enantiomer (all four stereocentres S in its
  `drawing_3D.mol` and `ligandrm.pdb`). Its CGenFF version, penalties and per-atom charges are
  unaffected by chirality, but do not use its coordinates for an R system.

### Quality checklist

- **CGenFF penalties:** below 10 the analogy is fair; 10–50 basic validation is recommended;
  above 50 the parameters need extensive validation or optimisation.
- **Stereochemistry:** every backbone stereocentre must be R in the SDF sent to CHARMM-GUI
  (RDKit CIP labels from the 3D coordinates) and still R after CHARMM-GUI (signed-volume check in
  06B).
- **`[ defaults ]`** in the GROMACS topology is CHARMM-style: `1 2 yes 1.0 1.0`.
- **Net charge and ion counts** match the intended system (checked in 06B).

### Files to keep

Keep the whole Ligand Reader & Modeler download (`lig/`, `toppar/`, `gromacs/`, `ligandrm.*`)
together with the exact SDF you uploaded. Treat them as one matched set: atom names, charges and
parameters in these files belong together.

### Limitations

- CGenFF targets drug-like small molecules; PHA oligomer parameters are assigned by analogy.
- There is no PHA-specific validation of these parameters. Review penalties and, for production
  work, consider targeted validation of the highest-penalty terms.
