# Simulation analysis scripts

Trajectory preparation needs only **Bash and GROMACS**. The research analysis
and enzyme-contact notebooks use Python analysis dependencies.

```text
md_simulation_scripts/
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
python md_simulation_scripts/enzyme_contacts/run_enzyme_contacts.py md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml
```

Or, from `md_simulation_scripts/enzyme_contacts/`:

```bash
python run_enzyme_contacts.py GK13_P3HO_4.yaml
```

The default is a 31-frame preview of the selected window. Use `--full` to analyse
all configured samples; `--sample-interval-ns` controls time sampling. This Python
interface is separate from the Bash `STRIDE`, which only controls the VMD file.
The moved YAML preserves its original raw inputs and the ignored
`examples/output/enzyme_contacts/` result location. Each contact run creates a
new result directory. See [the contact guide](../docs/enzyme_contacts.md).

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
