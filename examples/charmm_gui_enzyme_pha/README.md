# Example: enzyme + PHA in water with CHARMM-GUI (05B → 06B)

Force fields: protein **CHARMM36m**, PHA **CGenFF**, water **CHARMM TIP3P**, ions **SOD/CLA**. Engine: **GROMACS**.

Worked example: **ANC55 + P3HB4** (enzyme ANC55, PHA = 3-hydroxybutyrate tetramer, all R).

```
PHB4_R.sdf ──(A) Ligand Reader──► lig/ ──(B) Solution Builder──► download/ ──(C) this script──► run folder ──(D) HPC
```

---

## A. PHA parameters: CHARMM-GUI Ligand Reader & Modeler (notebook 05B)

1. Go to https://www.charmm-gui.org → *Input Generator* → **Ligand Reader & Modeler**.
2. **Upload MOL/MOL2/SDF**: `PHB4_R.sdf` (from notebook 04). Use the SDF, not the PDB: it keeps bond orders, hydrogens and the R stereocentres.
3. In Marvin JS, check that all stereocentres are R, then continue.
4. Download the `lig/` folder. It contains `lig.rtf` and/or `lig_g.rtf`, plus `lig.prm` (CGenFF topology and parameters).

## B. Build the system: CHARMM-GUI Solution Builder (notebook 06B)

| Page | Setting |
| --- | --- |
| PDB | Upload the docked enzyme–PHA complex, e.g. `Anc55_PHB4_pose2_complex.pdb` |
| Chains | Tick **Protein** (PROA) and **Hetero** (LIG); tick **Check pKa** |
| PDB manipulation | System pH **7.0**. *Reading hetero chain residues*: **Upload CHARMM top & par** → topology `lig_g.rtf` (or `lig.rtf` if absent), parameter `lig.prm` |
| Water box | **Fit to protein size**, rectangular, edge distance **30 Å** |
| Ions | **NaCl**, **0.05 M**, Monte-Carlo placement, neutralising |
| Temperature | **303.15 K** |
| Input Generator | **GROMACS** |

Download and unpack the package, e.g. `~/Downloads/Anc55_ANC_56_PHB4/Anc55_PHB4-9000069110/`.
It contains `gromacs/` (topology, `toppar/`, `step3_input.gro`, `index.ndx`) and `lig/`.

## C. Make the run folder and check it (this script)

From the repository root, in the `ipha_clean` environment:

```bash
python examples/charmm_gui_enzyme_pha/prepare_enzyme_pha_run.py \
    --download ~/Downloads/Anc55_ANC_56_PHB4/Anc55_PHB4-9000069110 \
    --sdf      ~/Data/MD_projects/PHA/PHA_Simulation_set_up/polymer_structures/PHB4_R.sdf \
    --out      ~/Data/Input/PHA_new_Enzyme/ANC55_P3HB4_gromacs_new
```

Add `--minimise` to also run the step 6.0 minimisation locally (needs `gmx`).
The download is never changed, and the `--out` folder must not exist yet.
If any check fails, the script exits before minimisation. Without `--sdf`, the
stereochemistry checks are skipped; provide it for a complete check.

To try the packaged input dataset without a download, run from the repository root:

```bash
python examples/charmm_gui_enzyme_pha/prepare_enzyme_pha_run.py \
    --download examples/data/charmm_gui_ANC55_P3HB4 \
    --sdf examples/data/charmm_gui_ANC55_P3HB4/PHB4_R.sdf \
    --out examples/output/ANC55_P3HB4_gromacs
```

Keep generated run folders under `examples/output/` (ignored by Git) or outside the repository.

Expected output for ANC55 + P3HB4 (abridged):

```
1. CGenFF penalties
Max charge penalty: 3.539 (fair analogy)
Max parameter penalty: 4.5 (fair analogy)

2. Build the run folder
  step6.0_minimization   up to 5000 steps
  step6.1_nvt            125 ps (125000 steps x 0.001 ps)
  step6.2_npt            500 ps (250000 steps x 0.002 ps)
  step7_production       200 ns (100000000 steps x 0.002 ps)

3. Checks
[PASS] [ defaults ] CHARMM-style: 1 2 yes 1.0 1.0
[PASS] #include files resolve
[PASS] step5_input.gro atoms == topology atoms: 141760 vs 141760
[PASS] molecule counts: PROA 1, LIG 1, SOD 50, CLA 43, TIP3 45977
[PASS] net charge ~ 0
[PASS] LIG atom types CGenFF-style
[PASS] LIG atom names/order match LIG.itp: 51 atoms
[PASS] LIG stereocentres all R
[PASS] run files and index groups: SOLU, SOLV present
[PASS] LIG charges equal lig.rtf: 51/51 atoms
[PASS] LIG stereocentres R (signed volume vs SDF): 4/4 R
```

The run folder contains:

| File | What it is |
| --- | --- |
| `step5_input.gro`, `topol.top`, `toppar/`, `index.ndx` | system from CHARMM-GUI (`step3_input.gro` renamed) |
| `step6.0_minimization.mdp` | minimisation, protein restraints |
| `step6.1_nvt.mdp` | NVT 125 ps, 303.15 K, restraints |
| `step6.2_npt.mdp` | NPT 500 ps, 1 bar, restraints |
| `step7_production.mdp` | production 200 ns (edit `nsteps` for longer runs) |
| `run_step6_local.sh` | local minimisation |
| `run_hpc_equilibration_production.slurm` | HPC job: NVT → NPT → production |

## D. Run

```bash
cd ~/Data/Input/PHA_new_Enzyme/ANC55_P3HB4_gromacs_new
bash run_step6_local.sh                           # minimisation (if not done with --minimise)
# copy the folder to the cluster, then:
sbatch run_hpc_equilibration_production.slurm     # NVT -> NPT -> production
```

The SLURM script is set up for the KCL GPU cluster (`module load gromacs/2021.5-gcc-11.4.0-cuda-11.8.0`).
Edit the `#SBATCH` lines and the `module load` line for other clusters.
Then analyse with notebooks 08 and 12.
