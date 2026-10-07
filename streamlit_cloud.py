"""Run the existing GUI with Community Cloud's single Python environment."""

import os
from pathlib import Path
import runpy
import sys


PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

# Community Cloud installs environment.yml into its own environment. It does
# not have the separate ~/miniconda3/envs/iphasimulator used by local installs.
os.environ.setdefault("IPHASIMULATOR_PYTHON", sys.executable)

# Keep every tab and action in the existing application.
runpy.run_path(str(PROJECT_ROOT / "pha_gui.py"), run_name="__main__")
