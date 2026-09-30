# PHA–enzyme distance and contact workflow

The workflow answers which enzyme residues stay close to the **whole selected
PHA** in an existing trajectory. It measures proximity and contacts, without
estimating affinity, assigning catalytic activity, or claiming convergence.

## Run it

Use the existing `ipha_clean` environment, which has the required dependencies.
The system Homebrew `python3` lacks MDAnalysis and matplotlib. From the project
root, run the preview:

```bash
cd /Users/k20098771/opt/iPHASimulator_v2
/opt/homebrew/Caskroom/miniconda/base/envs/ipha_clean/bin/python src/md_simulation_scripts/enzyme_contacts/run_enzyme_contacts.py src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml
```

Run the **full analysis** with:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/ipha_clean/bin/python src/md_simulation_scripts/enzyme_contacts/run_enzyme_contacts.py src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml --full
```

With the supplied YAML, full means both modes on all 14,806 stored frames over
0–1480.5 ns, including the initial frame. It does not discard equilibration.
`python -m iphasimulator.analysis_contacts CONFIG --full` is equivalent.

Alternatively open
[`src/md_simulation_scripts/enzyme_contacts/enzyme_contacts.ipynb`](../src/md_simulation_scripts/enzyme_contacts/enzyme_contacts.ipynb),
select a kernel using the `ipha_clean` interpreter, set `PREVIEW_FRAMES = None`
in its configuration cell, and run that cell and all following cells. The
notebook defaults to a 31-frame preview. Both interfaces use the same package
functions and create a new timestamped directory on every run.

In another environment, install from the repository with
`python -m pip install -e ".[analysis]"`. A notebook also needs Jupyter/IPython
and an appropriate kernel. No package installation was needed for validation.

## Configuration and sampling

Edit [`src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml`](../src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml).
Relative paths are resolved against the YAML file's directory. The reusable
module contains no simulation-specific paths or ligand-residue assumption.

After using the [standalone preparation workflow](enzyme_trajectory_processing.md),
set `topology` to the matching production TPR and `trajectory` to the full
`analysis/processed.xtc`, not its reduced VMD preview. The moved configuration
and notebook retain their previous raw-input defaults and ignored output location.
Centring/wrapping does not remove the need for the per-frame periodic-distance
calculation below; do not supply rotationally fitted coordinates with an unrotated box.

| Setting | Meaning |
| --- | --- |
| `topology`, `trajectory` | One matching topology and one trajectory; no discovery, concatenation or fallback to another topology. |
| `output_root` | Parent for new timestamp/UUID run directories. |
| `protein_selection`, `pha_selection` | MDAnalysis topology selections; both must be nonempty and disjoint. |
| `start_ns`, `end_ns` | Inclusive time window in ns; null uses the corresponding available boundary. Windows outside the data are rejected. |
| `sample_interval_ns` | Null retains every stored frame. A positive value must be an integer multiple of the regular native interval, anchored at the first frame within the window. |
| `modes` | `all_heavy`, `sidechain_heavy`, or both. |
| `primary_cutoff_A` | Primary strict contact cutoff, initially 4.5 Å. |
| `extra_cutoffs_A` | Additional cutoffs; 4.0 and 5.0 Å are always included. |
| `preview_frames` | At least 2: spread up to this many samples across eligible frame indices, including both ends. Null: use all eligible frames. |
| `distance_block_size` | Atom blocks bound temporary pair-distance memory; default 1024. |

CLI overrides include `--full`, `--preview-frames 31`, `--start-ns`, `--end-ns`,
`--sample-interval-ns`, `--mode both`, `--cutoff-A`, and `--output-root`.
For example, add `--start-ns 100 --end-ns 1000 --sample-interval-ns 1.0` to
analyse a deliberately chosen window at 1 ns intervals. This is an illustration
of configuration, not evidence for equilibration at 100 ns.

Every stored frame is decoded during inspection, even for a preview. The
workflow verifies topology/trajectory atom counts, strictly increasing finite
times and valid periodic dimensions, and saves all timestamps and boxes. A
duplicate/reversed time or invalid box fails the run. Irregular native intervals
are reported and can be retained by all-frame sampling; specifying a fixed
interval for irregular input is rejected instead of silently rounding.

The calculation makes a second streaming pass over selected frames. Coordinates
are not retained across frames; only residue × frame distance matrices and
metadata accumulate. Full calculation therefore reads the XTC twice. The preview
reads the full XTC once plus its sampled frames. Temporary symlinks direct
MDAnalysis's XTC/TRR offset caches to scratch, preserving the input directory.

## Definitions and implementation

For residue `r` and sampled frame `t`, the stored value is the minimum over
selected atom pairs between `r` and the whole PHA, using the current frame's
periodic cell. No protein fitting or coordinate transformations are applied.
An externally fitted file cannot be detected from its filename; always supply
the original trajectory or a trajectory whose coordinates and box remain in
the same frame of reference.

Hydrogen exclusion uses topology elements (H/D/T), falls back to mass greater
than 3.5 u if elements are absent, then atom names with leading digits stripped.
Zero-mass virtual sites are excluded. For unusual or coarse-grained topologies,
review the saved atom inventory before interpreting a heavy-atom analysis.

Side-chain mode further removes N, CA, C, O and common terminal oxygen aliases:
OXT, OT1, OT2, O1, O2, OC1, OC2, OCT1 and OCT2. Glycine always has an unavailable
side-chain row. The initial protein selection defines which residues are
included; rows are retained even if heavy/side-chain filtering leaves no atoms.
Narrowing the initial protein selection itself also narrows the residue set.

The code retains filtered atom indices, calculates atom-to-PHA minima in
vectorised distance-array blocks, then reduces those minima by topology residue
index. Metadata access through residues never expands the calculation atoms.
Both modes share the atom-distance calculation. This fixes the reference
notebook's pattern of selecting atoms, traversing their residues, and then
reintroducing excluded atoms with `residue.atoms`.

The numerical routine follows the documented [MDAnalysis minimum-image
distance-array API](https://docs.mdanalysis.org/2.8.0/documentation_pages/lib/distances.html),
which accepts orthorhombic and triclinic boxes. Reader-cache handling follows
the [XDR reader documentation](https://docs.mdanalysis.org/2.9.0/documentation_pages/coordinates/XDR.html).

Contact occupancy is `100 × count(distance < cutoff) / analysed_frame_count`.
The inequality is strict. Each residue contributes at most one contact per
sampled frame, irrespective of how many atom pairs fall inside the cutoff.
The denominator is the number of **analysed** frames, including the initial
sample if it belongs to the window. Any unavailable values leave that residue's
occupancy unavailable; finite-frame counts are never silently substituted as
the denominator. Mean and median also remain unavailable for incomplete rows.

These percentages are sampled-frame occupancies, not time-weighted occupancies
or independent observations. Temporal correlation and deliberately sparse
preview sampling limit interpretation; no standard errors or convergence claim
are attached. The preview's 31 samples are 49.3–49.4 ns apart. Its endpoint-aware
sampling is therefore slightly irregular, which is recorded explicitly.

## Outputs and reuse

All outputs are inside a fresh directory under `examples/output/enzyme_contacts/`.
These generated files are ignored by Git.

| File | Contents |
| --- | --- |
| `distances.npz` | Float64 residue × sampled-frame matrices for each mode, times in ns, original frame indices, same-frame boxes, residue metadata, selected atom counts and provenance. JSON string metadata; no pickle needed. |
| `residue_summary.csv` | One row per residue/mode: resindex, resid, resname, segment index/ID, chain IDs, insertion code, label, selected atom count, analysed/valid frame counts, availability, mean/median distance and occupancy at every cutoff. Missing values are explicit `NaN`. |
| `distance_heatmap_<mode>.png` | Distance heatmap for each mode, time in ns, distance colour bar in Å, shared colour scale, grey unavailable rows. |
| `contact_occupancy_<mode>.png` | Each cutoff's occupancy profile with a 0–100% axis and unavailable rows marked grey. |
| `config.json`, `provenance.json` | Resolved paths/settings, requested and actual windows, actual sampling, units, cutoffs, selections, versions, input fingerprints and implementation SHA-256. |
| `inspection.json`, `trajectory_frames.csv` | Complete trajectory validation summary plus every stored frame's time and periodic dimensions. |
| `selections.json`, `selected_atoms.csv`, `protein_residues.csv` | Selection policy, unavailable residue identities, exact selected atom indices and topology residue identities. |
| `status.json` | Running, complete or failed status; failed runs retain a diagnostic and their original directory. |

NPZ keys are `times_ns`, `frame_indices`, `boxes_A_degrees`, `residues_json`,
`provenance_json`, `distance_<mode>_A`, and `selected_atom_count_<mode>`.
The topology index `resindex` is zero-based and is the unambiguous row key.
`resid` is the original residue number, never a fabricated matrix index plus
one. Labels include the residue name/number, chain (or segment if absent), and
topology index. All chain and segment fields remain available in the metadata.

Plots use selected identity ticks for readability. Heatmap time bins are bounded
by sample midpoints; no distance values are interpolated or smoothed. Preview
titles explicitly give the label, number of samples and time range. Wide bins
in a sparse preview do not imply measurements between sampled frames.

Change cutoffs and plots without rereading the trajectory:

```python
from pathlib import Path
from iphasimulator.analysis_contacts import load_distances, summarise, plot_contacts

run = Path("examples/output/enzyme_contacts/<run-directory>")
saved = load_distances(run / "distances.npz")
derived = run / "replot_6A"
derived.mkdir(exist_ok=False)
cutoffs = (4.0, 4.5, 5.0, 6.0)
summarise(saved, cutoffs).to_csv(derived / "residue_summary.csv", index=False, na_rep="NaN")
plot_contacts(saved, derived, cutoffs, primary_cutoff_A=4.5, distance_vmax_A=15)
```

The optional 15 Å colour maximum changes display saturation only. Existing
NPZ files and figures are protected against replacement. Use a new derived
directory for each revised plot set.

## Verified GK13–P3HO_4 preview, 7 September 2026

Run directory:
`examples/output/enzyme_contacts/20260907T132220_010222Z_preview_895d7808/`.

The matching inputs were inspected directly in
`/Users/k20098771/Data/MD_projects/PHA/MD_data/gromacs/PHA_Enzyme/GK13_P3HO_4_gromacs/`:
`step7_production.tpr` and **only** `production_combined_1us.xtc`. Continuation
files and other merged trajectories were not concatenated or analysed.

| Verified property | Result |
| --- | --- |
| TPR / XTC atom counts | 137,360 / 137,360 |
| Stored frames and time range | 14,806; 0–1480.5 ns |
| Native interval | Uniform 0.1 ns across every stored frame |
| Periodic cell | Cubic, 90° angles; lengths 110.6547–111.0368 Å; valid at every frame |
| Protein | `protein`: 245 residues, 3,701 atoms; segment `seg_0_PROA`, chain `PROA`; LEU1–PRO245 |
| Protein heavy atoms | 1,898 all-heavy; 917 side-chain-heavy |
| PHA | `resname LIG and segid seg_1_LIG`: 99 total atoms; 41 heavy atoms (32 C, 9 O) |
| PHA identity | resindex 245, resid 246, resname LIG, segment `seg_1_LIG` |
| Missing side chains | 25 glycines, retained as NaN |
| Preview | 31 frames, including 0 and 1480.5 ns; actual spacing 49.3–49.4 ns |

`topol.top` declares one PROA and one LIG molecule and includes `toppar/LIG.itp`.
The ligand composition is consistent with the P3HO tetramer. No repeat-unit
mapping or catalytic annotation is inferred. XTC does not encode atom identities,
so atom count alone cannot prove atom ordering. As a separate check, all 3,869
protein/PHA topology bonds were evaluated in five widely spaced frames; lengths
were 0.9496–1.9008 Å, supporting compatibility of the supplied atom order.

Validation evidence in the preview directory includes `validation.json`,
`independent_distance_checks.csv`, `topology_bond_checks.csv`, and
`validate_preview.py`. Five frames (indices 0, 3454, 7402, 11350, 14805) and six
residues (resindices 0, 48, 91, 98, 139, 244) were cross-checked in both modes.
An independent float64 NumPy implementation wrapped displacements with
`delta -= box * rint(delta / box)` and took Euclidean atom-pair minima.
The **55 finite distance checks agreed within 5.07 × 10⁻⁶ Å**; the five missing
glycine checks remained NaN. This independent real-data check uses the observed
orthorhombic cell; triclinic behavior is separately covered by a known-coordinate
test against explicit lattice-image enumeration.

All summary occupancy values were independently recomputed from the saved
matrices. Every available side-chain minimum was at least its corresponding
all-heavy minimum, as required by atom-set inclusion. All four PNGs were visually
inspected for axes, units, readable identities, shared distance scale, missing
rows and occupancy limits.

All seven notebook code cells were also executed sequentially and headlessly
in `ipha_clean`, starting from `src/md_simulation_scripts/`, without substituting
their source. The notebook produced a second fresh preview at
`examples/output/enzyme_contacts/20260907T133200_068659Z_preview_154678c6/`;
its frame indices and both distance matrices match the terminal preview exactly,
including NaNs. The primary preview's `notebook_execution.json` and
`notebook_execution.log` record this check. The source notebook retains clean,
unexecuted cells. `test_results.json` records the automated-test results.

In this **sparse preview**, PHE49 and PHE92 have the highest 4.5 Å occupancy:
12/31 frames (38.71%) in both modes. Other prominent sampled contacts include
SER48, ASN96, ILE97 and VAL197. These are preliminary sampled proximities;
use the full analysis before describing persistence or ranking regions.

The environment used Python 3.11.15, MDAnalysis 2.10.0, NumPy 2.4.3,
SciPy 1.17.1, pandas 2.3.3 and matplotlib 3.10.9; provenance records these.

Automated checks:

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/ipha_clean/bin/python -m pytest tests/test_analysis_contacts.py -q
```

**17 focused tests passed**, including known distances, orthorhombic/triclinic
periodicity, changed/invalid boxes, hydrogen/terminal exclusions, glycine,
duplicate residue numbers across chains, original-selection preservation,
strict cutoff/denominator behavior, irregular sampling, serialization, plotting,
input preservation, atom-count mismatch and fresh output directories. The 34
existing trajectory preprocessing/merging and naming tests also passed.

The full `pytest tests -q` attempt stopped at **10 collection errors** caused by
pre-existing misplaced `from __future__ import annotations` statements in other
modules (`build`, `monomers`, `stereochemistry`, `parameterization_gaff2`,
`conversion_amber_to_gromacs`, `simulation_gromacs_runner`, and
`system_builder_packmol`). Those files have no changes in this implementation;
these errors do not affect the new module, example or notebook.

The two earlier reference notebooks and original simulation files were preserved.
This first workflow intentionally stops at residue-to-whole-PHA distances and
contacts. Repeat-unit mapping, catalytic geometry, APO RMSF and batch comparisons
remain subsequent extensions.
