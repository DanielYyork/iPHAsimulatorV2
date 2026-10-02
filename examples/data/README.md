# Example data

## `charmm_gui_ANC55_P3HB4/`: enzyme + PHA in water, CHARMM-GUI Solution Builder

A small dataset for the CHARMM/GROMACS method (notebooks 05B and 06B), taken from one
CHARMM-GUI download (job 9000069110). It is the default input of
`notebooks/06B_cgenff_gromacs_pha_enzyme_in_water.ipynb`.

**The system**

- Enzyme **ANC55** with one PHA oligomer, **P3HB4** (3-hydroxybutyrate tetramer, all four
  stereocentres R), from the docked complex `Anc55_PHB4_pose2_complex.pdb`.
- Force fields: protein **CHARMM36m**; PHA **CGenFF** (CHARMM-GUI Ligand Reader & Modeler:
  CGenFF program 4.0, parameter files 5.0, penalties param 4.5 / charge 3.539); water
  **CHARMM TIP3P**; ions **SOD/CLA** (NaCl).
- Built with CHARMM-GUI **Solution Builder**, GROMACS output: rectangular box 11.4 nm,
  141 760 atoms; molecules PROA 1, LIG 1, SOD 50, CLA 43, TIP3 45 977.

**Files** (same layout as the unpacked download, about 9 MB)

| File | What it is |
| --- | --- |
| `gromacs/step3_input.gro` | Solvated, ionised coordinates (06B renames it `step5_input.gro` in the run folder) |
| `gromacs/topol.top` | GROMACS topology |
| `gromacs/toppar/*.itp` | `forcefield.itp`, `PROA.itp`, `LIG.itp`, `TIP3.itp`, `SOD.itp`, `CLA.itp` |
| `gromacs/index.ndx` | Index groups, including `SOLU`/`SOLV` used by the mdp files |
| `lig/lig.rtf`, `lig/lig.prm` | CGenFF topology and parameters for the PHA (`LIG`), with penalty scores |
| `PHB4_R.sdf` | The R-configured PHA structure uploaded to Ligand Reader & Modeler (from iPHASimulator) |

The files are byte-identical copies of the download (and of `PHB4_R.sdf`). This download
has no `lig/lig_g.rtf`. In Solution Builder: Upload `lig_g.rtf` as topology and `lig.prm` as
parameters. If your download has no `lig_g.rtf`, use `lig.rtf`; the charges are identical.

**Left out**, to keep the dataset small and limited to what 06B reads:

- `.psf`, `.pdb` and `.crd` files (`gromacs/step3_input.psf`/`.pdb` alone are about 29 MB):
  GROMACS and 06B use the `.gro` and `.top` files.
- The CHARMM `step1_*`, `step2_*` and `step3_pbcsetup.*` inputs, outputs and structures: the
  intermediate CHARMM-GUI build steps.
- The CHARMM `toppar/` folder and `.str` files (about 35 MB): the force field is already in
  `gromacs/toppar/*.itp`.
- CHARMM-GUI's `step4.0`, `step4.1` and `step5` mdp files: 06B leaves them out of the run folder
  and adds the production mdp files and scripts from `src/iphasimulator/data/charmm_gromacs/`.

Do not write into this folder. 06B writes its run folder to a temporary folder by default
and refuses a `RUN_DIR` inside `examples/data/`.

**Citations**

If you use this dataset or the method, cite CHARMM-GUI and CGenFF:

- S. Jo, T. Kim, V. G. Iyer, W. Im. CHARMM-GUI: A web-based graphical user interface for
  CHARMM. *J. Comput. Chem.* **29**, 1859–1865 (2008). doi:10.1002/jcc.20945
- J. Lee *et al.* CHARMM-GUI Input Generator for NAMD, GROMACS, AMBER, OpenMM, and
  CHARMM/OpenMM simulations using the CHARMM36 additive force field. *J. Chem. Theory Comput.*
  **12**, 405–413 (2016). doi:10.1021/acs.jctc.5b00935
- S. Kim *et al.* CHARMM-GUI Ligand Reader and Modeler for CHARMM force field generation of
  small molecules. *J. Comput. Chem.* **38**, 1879–1886 (2017). doi:10.1002/jcc.24829
- K. Vanommeslaeghe *et al.* CHARMM General Force Field: A force field for drug-like molecules
  compatible with the CHARMM all-atom additive biological force fields. *J. Comput. Chem.*
  **31**, 671–690 (2010). doi:10.1002/jcc.21367
- K. Vanommeslaeghe, A. D. MacKerell Jr. Automation of the CHARMM General Force Field (CGenFF) I:
  bond perception and atom typing. *J. Chem. Inf. Model.* **52**, 3144–3154 (2012).
  doi:10.1021/ci300363c
- K. Vanommeslaeghe, E. P. Raman, A. D. MacKerell Jr. Automation of the CHARMM General Force
  Field (CGenFF) II: assignment of bonded parameters and partial atomic charges.
  *J. Chem. Inf. Model.* **52**, 3155–3168 (2012). doi:10.1021/ci3003649
- J. Huang *et al.* CHARMM36m: an improved force field for folded and intrinsically
  disordered proteins. *Nat. Methods* **14**, 71–73 (2017). doi:10.1038/nmeth.4067
