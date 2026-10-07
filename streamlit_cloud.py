"""Run the existing GUI with Community Cloud's single Python environment."""

import os
import importlib.util
from pathlib import Path
import runpy
import sys


PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

# Community Cloud installs environment.yml into its own environment. It does
# not have the separate ~/miniconda3/envs/iphasimulator used by local installs.
os.environ.setdefault("IPHASIMULATOR_PYTHON", sys.executable)

# Record only runtime/package information, never environment variables or
# secrets. Cloud's public error page hides the missing module's name.
print(
    f"[iPHA Cloud] Python {sys.version.split()[0]}; executable={sys.executable}",
    flush=True,
)
package_status = {}
for module in ("parmed", "rdkit", "openmm", "MDAnalysis", "py3Dmol", "stmol"):
    status = "available" if importlib.util.find_spec(module) is not None else "MISSING"
    package_status[module] = status
    print(f"[iPHA Cloud] {module}: {status}", flush=True)

# Keep every tab and action in the existing application.
try:
    runpy.run_path(str(PROJECT_ROOT / "pha_gui.py"), run_name="__main__")
except ModuleNotFoundError as error:
    print(f"[iPHA Cloud] Import failed: missing module {error.name!r}", flush=True)
    # Cloud may omit stdout from downloadable logs and redact exceptions.
    # Show a bounded diagnostic on startup failure; never expose secrets or
    # the full process environment. Successful startup keeps the existing UI.
    import streamlit as st

    st.error("The Cloud Python environment is missing an application dependency.")
    st.write("Copy the diagnostic below when reporting this startup failure.")
    st.json({
        "missing_module": error.name,
        "python_version": sys.version.split()[0],
        "expected_python_version": "3.12 (environment.yml)",
        "python_executable": sys.executable,
        "packages": package_status,
    })
    st.info(
        "The repository already declares ParmEd in environment.yml. "
        "This diagnostic checks the running environment; it does not confirm "
        "that the Conda installation completed."
    )
    st.stop()
