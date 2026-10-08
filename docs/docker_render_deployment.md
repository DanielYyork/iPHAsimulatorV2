# Deploy the workshop GUI with Docker on Render

Use **Render → Web Service → Docker**, with the repository-root **`Dockerfile`**.
This deployment does not use Streamlit Community Cloud. You do not need Docker
installed on your Mac: Render builds the container from your GitHub branch.

The existing GUI and its actions are retained. It runs the full scientific
environment, including AmberTools/PREPGEN, Open Babel, RDKit, OpenMM, ACPYPE,
Polyply and analysis tools. CPU execution is selected for the hosted workshop.
Submit to Slurm remains present and requires an actual Slurm cluster to work.

## What changes from the failed Cloud installation

`deploy/conda-linux-64.lock` contains explicit download URLs and checksums for
374 packages derived from the previously resolved Linux environment. Micromamba
installs those archives **without a dependency solve**. The Conda packages
alone total approximately 1.03 GB of compressed downloads; the first build
will still take time to download and extract them.

The image currently retains the Conda package cache. The first Render build
completed the Conda transaction, then failed in `micromamba clean --all` on
AmberTools 26's cached `bin/amber.conda` entry. This matches an
[upstream cleanup issue](https://github.com/conda-forge/ambertools-feedstock/issues/193).
The optional cleanup command has been removed; this increases image size but
keeps installation failures and the application checks fatal to the build.

The Docker lock uses MDTraj 1.10.3 because 1.11.1's Python metadata requires
NumPy 2, while this environment's PyDeck requires NumPy below 2. This avoids a
conflict that Conda's package metadata alone did not identify. The original
Community Cloud environment file is unchanged; the override is recorded in
`deploy/environment-provenance.json`.

`deploy/pip-requirements.txt` installs the remaining Python packages. Its
scientific-package constraints prevent pip from silently replacing the main
Conda scientific stack. Direct pip requirements are pinned; the entire pip
dependency tree and base container digest are not fully locked. This is not
a claim of a byte-for-byte reproducible image.

The Docker build stops on any failed installation, `pip check`, or functional
smoke check. It exercises PREPIN generation, construction of P3HB with four
units, 20 OpenMM CPU steps, MDAnalysis trajectory loading, synthetic PCA/Tg
analysis, workflow-script generation/compilation and all seven GUI tabs.
These checks are representative; they do not establish every scientific
workflow, browser interaction or workshop capacity.

## 1. Push the deployment files

Commit these exact additions on the branch you want to deploy (currently `dan`):

```bash
git add Dockerfile .dockerignore render.yaml deploy cloud/verify.py docs/docker_render_deployment.md
git diff --cached --stat
git commit -m "Add Docker workshop deployment with explicit scientific environment"
git push origin dan
```

Review staged changes before committing; do not use `git add .`, because this
checkout contains unrelated research data and scripts.

The existing `pha_gui.py`, `streamlit_cloud.py`, `cloud/streamlit_cloud.py`,
`gui/`, `src/`, `cluster/`, `pyproject.toml`, residue registry and PHA chemistry
inputs must already be on the branch. The `.dockerignore` deliberately excludes
`.git`, secrets, generated MD systems, trajectories, scratch folders and local
registries. A new deployment starts with the supplied chemistry and empty
polymer/system registries, rather than links to research files absent from it.

The checkout's `src/iphasimulator/analysis/tg_analysis/system_analysis.py` is
currently untracked. If you want that separate system-summary workflow deployed,
review and commit it too. It is not needed for the current GUI's replica-analysis
startup check. GitHub deployments cannot include untracked local files.

## 2. Create a Render Web Service

1. Sign in at [Render](https://dashboard.render.com).
2. Select **New → Web Service** and connect your existing GitHub repository.
3. Use these settings:

   | Setting | Value |
   | --- | --- |
   | Name | `iphasimulator-workshop`, or another available name |
   | Branch | `dan` |
   | Language / Runtime | **Docker** |
   | Root directory | Leave empty |
   | Dockerfile path | `./Dockerfile` |
   | Docker build context | Repository root (`.`), if shown |
   | Docker command | Leave empty; the image supplies its startup command |
   | Health-check path | `/_stcore/health` |
   | Auto-deploy | Off for the initial workshop setup |
   | Starting instance | `1c-2g`: 1 CPU, 2 GB RAM (formerly Standard) |

4. Under the service's disk settings, add a **5 GB** persistent disk mounted at
   **`/app/structure_database`**. It preserves generated polymers, systems,
   simulations, analysis results and saved workflow definitions.
5. Review the price and click **Deploy / Create Web Service** when ready.
6. Watch the build log, then open the supplied `https://...onrender.com` URL.

No build/start command, Python-version selector or Sunbird credentials are
needed. The Dockerfile supplies Python and the launch command. The launcher
binds to `0.0.0.0` and Render's `PORT` value.

The proposed compute plan is currently approximately **US$25/month**, with
persistent storage and any other applicable charges additional. Check Render's
displayed total before accepting. This is a starting point for GUI use and small
examples, not a capacity guarantee. The large research trajectories and full
melt analyses can require substantially more RAM and disk.

Alternatively, **New → Blueprint** can read the included `render.yaml` and
pre-fill the same paid plan, region and disk. Review the resources and prices
before approving creation. No account, service or paid resource was created
while preparing these files.

## 3. Confirm success

The build log should show actual archive downloads/extraction, pip installation,
and the six `CHECK ...` stages. This route should not show the earlier long
Conda `Solving environment` stage. Any failed check must fail the build.

After deployment:

1. Open the URL and confirm all seven tabs appear.
2. Build a short polymer from a chemistry with complete parameters, such as 3HB.
3. Build a small dry system, generate a short OpenMM workflow, and use **Run locally**.
   Here, “locally” means the Render server.
4. Download any outputs you need, and check the saved workflow/results after
   restarting the service.

The persistent disk stores the structure database, including workflow JSON.
Generated Python scripts in `/app/md_simulation_scripts` are outside that disk
and may disappear on redeployment. Regenerate them from the saved workflow.
Existing persistent chemistry is reused; startup never overwrites it with a
new image's chemistry library. A partially populated disk without the residue
registry causes an explicit error for review instead of silent replacement.

The current application shares its database across browser sessions. This
deployment preserves that behaviour. Avoid concurrent builds of the same system
name. Docker preserves software functionality, not unlimited compute capacity.

## 4. If deployment fails

Use Render's build/deploy logs. Copy the **first failed command and its error**:

- `Transaction finished` followed by a cleanup error mentioning `amber.conda`:
  deploy the updated Dockerfile, which omits `micromamba clean --all`. Push it
  to `dan`, then choose **Manual Deploy → Deploy latest commit** on Render.
- Archive download/checksum failure: the log identifies the exact package URL.
- Pip conflict: the resolver or `pip check` names the conflicting dependencies.
- `CHECK ...` failure: the traceback names the failing scientific or GUI step.
- Disk error: confirm the mount is `/app/structure_database` and inspect existing
  contents. Do not delete a populated disk just to restart.
- Runtime memory limit: reduce the example workload or choose a larger instance.

Do not return to editing Community Cloud's `environment.yml` for this deployment;
Render follows the Dockerfile. The old Cloud files remain available for reference.

## Optional local Docker commands

Only for a computer that already has Docker installed:

```bash
docker build --platform linux/amd64 -t iphasimulator-workshop .
docker run --rm --platform linux/amd64 -p 8501:8501 \
  -v ipha-workshop-data:/app/structure_database iphasimulator-workshop
```

Open `http://localhost:8501`. Linux x86_64 is required by the explicit package
lock; Apple Silicon uses Docker's amd64 emulation.

## Verification status

Prepared on 8 October 2026. The Conda lock's archive identities and checksums
were matched to conda-forge package metadata. Local checks passed for fresh
database creation, restart preservation and refusal to overwrite a partial
database. In a disposable copy of the committed application, PREPIN generation,
polymer construction, a 20-step CPU simulation, trajectory loading, PCA/Tg
analysis, workflow-script compilation and all seven GUI tabs passed. The GUI
still displays its existing notices for incompletely parameterised monomers.

Those application checks used the existing Mac scientific and GUI environments
separately; they do not validate the new Linux container. A Docker runtime is
not installed on this Mac, so **the complete Docker build and Render deployment
have not yet been executed here**.
The build includes blocking checks so a failed scientific installation cannot
quietly advance to deployment as a working app.

References:

- [Micromamba explicit locks bypass the solver](https://micromamba-docker.readthedocs.io/en/latest/advanced_usage.html#using-a-lockfile)
- [Docker on Render](https://render.com/docs/docker)
- [Render compute plans](https://render.com/docs/compute-plans)
- [Render pricing](https://render.com/pricing)
- [Persistent disks](https://render.com/docs/disks)
- [HTTP health checks](https://render.com/docs/health-checks)
