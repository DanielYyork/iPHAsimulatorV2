# Installation

## 1. Get the repository

Use a terminal on Linux, macOS, or Windows with WSL. For an existing checkout,
change into its directory; otherwise:

```bash
git clone https://github.com/MMLabCodes/iPHAsimulatorV2.git
cd iPHAsimulatorV2
```

## 2. Create a Python environment

The package requires Python 3.10 or newer; Python 3.11 is used for these guides.

```bash
conda create -n ipha python=3.11
conda activate ipha
python -m pip install -e .
```

This installs the core dependencies declared in `pyproject.toml`: NumPy, SciPy,
NetworkX, pandas, PyYAML and RDKit. Editable installation keeps imports connected
to this checkout. No simulation is started by these commands.

## 3. Add the tools for your workflow

| What you want to do | Install |
| --- | --- |
| Run notebooks and tests | `python -m pip install -e ".[dev]"` |
| Use OpenMM, ParmEd and MDTraj | `python -m pip install -e ".[md]"` |
| Analyse enzyme–PHA contacts | `python -m pip install -e ".[analysis]"` |
| Launch the Streamlit interface | `python -m pip install -e ".[gui]"` |

AmberTools and GROMACS executables are separate from these Python extras. Install
the tools your selected tutorial needs, for example through conda-forge:

```bash
conda install -c conda-forge ambertools gromacs openbabel
```

The database builder additionally uses Open Babel's Python `pybel` interface.
Packmol is needed only for the advanced packing route. GPU platforms and cluster
module names depend on your local installation.

## 4. Check the installation

```bash
python -c "import iphasimulator; from iphasimulator.naming import oligomer_name; print(oligomer_name('3HB', 4))"
```

Expected output: `P3HB_4`. This checks the package and naming helper; it does not
certify every module. For contact analysis, also check:

```bash
python -c "import MDAnalysis, matplotlib; from iphasimulator.analysis_contacts import ContactConfig; print('Contact imports OK')"
```

For 06C, verify `tleap` is on PATH and `import openmm` succeeds in the
same notebook kernel. For 06A/06B execution, verify `gmx --version`. Activate
the environment before starting Jupyter; a Python import check does not prove
that external command-line tools are available.


## 5. Open a notebook

```bash
jupyter lab
```

Choose a notebook from the [catalogue](notebooks.md), then select the kernel for
the environment you installed. Review its editable file paths and execution
switches first. Several notebooks contain machine-specific paths and historical
outputs; downloading or viewing them does not validate those inputs for your
system.

For building only this documentation, use the separate, lightweight
[documentation environment](contributing_docs.md). It needs no MD engines.
