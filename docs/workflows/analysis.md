# Analysis

Choose the analysis that matches your research question. Structure fluctuations
and minimum-image contacts have different coordinate-preparation requirements.

## Enzyme–PHA distances and contacts

This reusable workflow measures each enzyme residue's minimum distance to all
selected PHA heavy atoms. It supports all-heavy and side-chain-heavy modes,
preserves residue identities and unavailable glycine side chains, and reports
the percentage of sampled frames below each contact cutoff.

1. Open [the contact notebook](../notebooks.md#execution-and-analysis).
2. Set `TPR_PATH`, `TRAJECTORY_PATH` and `OUTPUT_DIR` in step 1. Use a matching
   TPR and one original production XTC, or the full `processed.xtc` from the
   [standalone preparation workflow](../enzyme_trajectory_processing.md). The
   notebook retains its previous `step6.2_npt.tpr`/raw-XTC defaults; update both
   paths for your prepared system.
3. Review the protein/PHA selections, time window, sampling and cutoffs. The
   notebook starts with a 31-frame preview; `PREVIEW_FRAMES = None` requests the
   full selected analysis.
4. Run the input checks and analysis. Four PNG plots are saved and displayed.
5. Inspect `residue_summary.csv`, `distances.npz` and the saved provenance.

For terminal runs, edit a copy of the
{download}`contact configuration <../../src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml>`:

```bash
python src/md_simulation_scripts/enzyme_contacts/run_enzyme_contacts.py src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml
```

Add `--full` for all configured sampled frames. The original YAML uses
`step7_production.tpr`; notebook overrides are local to the notebook and do not
edit that YAML. Configure your own paths explicitly in whichever interface you use.

Distances use Å, time uses ns, occupancy uses 0–100%. A residue contributes once
per sampled frame when `distance < cutoff`; default comparisons are 4.0, 4.5 and
5.0 Å. These are correlated sampled-frame proportions, not affinity, catalytic
activity or independent observations. No equilibration period is discarded
automatically. See the [full contact method and validation record](../enzyme_contacts.md).

## Basic polymer structure

Notebook **09** uses MDTraj and the centered trajectory from notebook 08 to
calculate radius of gyration, an example end-to-end distance, and SASA. Review
the atom selection and units in that notebook: its MDTraj calculations are not
the contact module's Å/ns interface. Its end-to-end example uses the first and
last selected atoms, which need not be the most meaningful chemical endpoints.
Its SASA is computed on the polymer-only trajectory.

## Enzyme/polymer stability diagnostics

Notebook **12** examines total energy, protein-backbone RMSD and polymer RMSD
after protein alignment. Its inputs and selections are system-specific. These
diagnostics support manual review; the notebook does not establish convergence.
Do not reuse fitted coordinates with an untransformed periodic box for contacts.

APO RMSF comparisons, catalytic geometry, PHA repeat-unit mapping and batch
contact comparisons are **future extensions**. Existing empty APO notebook
placeholders do not implement them.

API: {py:func}`iphasimulator.analysis_contacts.run_contacts`,
{py:func}`iphasimulator.analysis_contacts.load_distances`,
{py:func}`iphasimulator.analysis_contacts.plot_contacts`.
