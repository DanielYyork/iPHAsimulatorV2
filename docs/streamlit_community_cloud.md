# Deploy the existing GUI to Streamlit Community Cloud

Use **`cloud/streamlit_cloud.py`** as the Cloud main file. The root launcher
still exists, but deploying it directly uses the broad development environment
that stalled during the reported Cloud builds.

This deployment retains the existing GUI, including polymer construction,
system building, viewing, script generation, local simulation execution and
analysis. The **Submit to Slurm** button remains present, but cannot submit
jobs without a configured Slurm installation. No Sunbird connection is used.

## What is already present

- `pha_gui.py`: the existing application.
- `gui/` and `src/iphasimulator/`: GUI and scientific code.
- `environment.yml`: the full Conda environment, including Streamlit, RDKit,
  AmberTools, ACPYPE, OpenMM, MDAnalysis, analysis packages and viewers. Its
  pip section installs the project with `-e .`.
- `pyproject.toml`: package installation metadata.
- `structure_database/residue_codes.csv`, `polymer_smiles.csv` and tracked
  `PHA_types/` files: the monomer database and available structure parameters.
- `streamlit_cloud.py`: a Cloud launcher that sets `IPHASIMULATOR_PYTHON` to
  the running interpreter, then runs the existing `pha_gui.py` unchanged.

The `cloud/streamlit_cloud.py` entry point selects the adjacent
`cloud/environment.yml` before any root dependency file. This Linux-only
snapshot pins all 374 Conda packages resolved for Python 3.12 and Streamlit
1.58.0, keeping every application/scientific dependency and omitting only
pytest and JupyterLab. Its pip requirements are retained from the original
environment and are not a complete pip lockfile. The original unpinned solve
also selected Streamlit 1.9.0, so an explicit modern Streamlit version is
important.

Cloud does not automatically select `environment-gui.yml`. Do not add a
competing requirements file. The root environment remains for local/Sunbird
use. This snapshot reduces Conda version choices; it does not bypass solving
or guarantee installation into Cloud’s existing base environment.
See [Streamlit dependency selection](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/app-dependencies).

## 1. Put the required files on GitHub

Community Cloud downloads the selected GitHub branch, not the files on your
computer. Review your changes before committing. This checkout uses branch
`dan`; deploy whichever branch you actually push the workshop files to.

From the project directory, add the new deployment files:

```bash
git add cloud streamlit_cloud.py docs/streamlit_community_cloud.md
git diff --cached --stat
git commit -m "Add Streamlit Community Cloud launcher and guide"
git push origin dan
```

These commands assume `dan` is still your branch and that the staged changes
are the ones you intend to publish. Review any previously staged changes too.
Commit other intended application changes separately before deployment.

At review time, `src/iphasimulator/analysis/tg_analysis/system_analysis.py`
was untracked. Include it if the deployed workflows need system-level analysis;
an untracked local file will not be available on Cloud.

The existing `.gitignore` excludes built polymers, MD systems and trajectories.
Consequently, a new Cloud instance does not initially have your locally built
systems or melt simulations. `md_systems.csv` is created empty if absent.
For workshop examples, distribute selected structures/results deliberately
with a matching registry; copying the registry alone does not copy its systems.
Do not upload the complete trajectory collection just to start the GUI.

## 2. Create the Cloud application

1. Open [Streamlit Community Cloud](https://share.streamlit.io) and sign in
   with the GitHub account that has access to the repository.
2. Select the workspace for the repository owner, then **Create app**.
3. Choose **Yup, I have an app** if prompted.
4. Enter:

   | Setting | Value |
   | --- | --- |
   | Repository | `MMLabCodes/iPHAsimulatorV2` |
   | Branch | `dan`, or the branch containing your deployment changes |
   | Main file path | `cloud/streamlit_cloud.py` |
   | App URL | An available workshop subdomain of your choice |

5. Open **Advanced settings**, select **Python 3.12** to match
   `cloud/environment.yml`, and save.
6. No secrets or Sunbird credentials are required for this setup. In
   particular, do not set `IPHASIMULATOR_PYTHON` to a Mac or `/scratch/` path.
7. Click **Deploy** and watch the build logs.

You can create a new app with this entry point and keep the failed app until
the replacement works. In the build log, confirm the selected dependency file
is `/mount/src/iphasimulatorv2/cloud/environment.yml`. If the root
`environment.yml` is selected, the entry point is still incorrect.

The full scientific environment has substantial dependencies; allow time for
installation. Successful local startup does not establish that the Linux
Conda solve will succeed. Inspect the actual Cloud build log before changing
dependency versions.

Official reference: [Deploy your app](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy).

## 3. Verify the deployed application

Check the following in order:

1. All seven main tabs appear and the monomer list loads.
2. The sidebar's Python executable belongs to the Cloud environment, rather
   than your local Miniconda directory or a Sunbird scratch directory.
3. Choose a monomer with complete parameter files, create a short polymer
   sequence, and inspect its molecular preview.
4. Build a small polymer and confirm that coordinates and topology are produced.
5. Build a small MD system and inspect it in the viewer.
6. Generate an OpenMM script and run a very short CPU simulation with
   **Run locally**. In this deployment, “locally” means on the Cloud server.
7. Check analysis using a suitably small compatible simulation. This does not
   establish that the full melt trajectories fit within Cloud resources.

The existing OpenMM code attempts CUDA, then OpenCL, then CPU when no platform
is explicitly requested. Do not force a GPU platform for the first Cloud test.

**Submit to Slurm** still calls the existing Sunbird-specific launcher and is
expected to fail here. Installing an `sbatch` executable alone would not supply
a working scheduler or a Sunbird connection.

## 4. Troubleshoot using the Cloud log

| Symptom | What to check |
| --- | --- |
| Conda solve or installation failure | The first failing package in the build log; check that the log names `cloud/environment.yml`, not the root development environment. |
| Missing `iphasimulator` module | The `-e .` installation in `environment.yml` completed, and `src/` plus `pyproject.toml` are in the deployed branch. |
| Missing scientific Python executable | Main file is `cloud/streamlit_cloud.py`; remove any obsolete `IPHASIMULATOR_PYTHON` override from Cloud secrets. |
| No existing systems or polymers | Generated directories are excluded from Git; build examples in the app or supply selected example files. |
| Missing monomer parameters | The selected monomer has complete parameter files in the deployed `PHA_types/` tree. A monomer appearing in the CSV alone is insufficient. |
| Missing `tleap` or `acpype` | Full environment installation succeeded, and the Cloud interpreter is selected. |
| Resource-limit error during building or analysis | Check the workload and Cloud's current resource allowance; changing the launcher does not provide more memory. |
| Slurm submission failure | Expected without a scheduler; use **Run locally** for Cloud execution. |

Use **Manage app** to view logs. Code changes pushed to the deployed branch
are picked up automatically; dependency changes trigger environment updates.

## Storage and runtime limits

The GUI's current files and registry are shared on the server; they are not
separate private workspaces for each browser session. Community Cloud does not
guarantee persistence of locally generated files. Download needed outputs and
keep durable copies elsewhere.

All application features remain in the code, but Cloud resource limits still
apply. The previously analysed melt runs used a descriptor matrix of
20,000 × 5,050 float64 values per chain for P3HO, around 20.2 GB across 25
chains before PCA overhead. Those exact analyses do not fit the approximately
2.7 GB maximum memory allowance currently documented for ordinary Community
Cloud apps; Streamlit labels its published limits approximate and subject to
change. This does not prevent testing GUI startup and small examples.

References:

- [Resource limits and app management](https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app)
- [Local storage persistence](https://docs.streamlit.io/develop/concepts/connections/connecting-to-data)

## Scope of verification

Local Streamlit AppTest checks passed for the launcher: all seven main tabs
rendered without unhandled exceptions, and the backend interpreter matched
the running Python. A second check used an isolated copy containing tracked
application files and tracked monomer data plus the new launcher. Existing
warnings about incomplete monomer parameters were still shown.

These checks used the installed GUI environment on macOS. They do not verify
the Linux Conda build, scientific execution or rendering in a real Cloud
browser. The Conda portion of the Cloud snapshot was separately resolved for Linux-64
with a glibc 2.31 compatibility constraint. This is a dry-run, not a Linux
package installation. Pip installation and execution on Cloud remain to be
verified. No Community Cloud deployment, Git commit or push has been performed.
