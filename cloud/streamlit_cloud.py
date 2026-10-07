"""Select the Cloud-specific environment and run the unchanged GUI launcher."""

from pathlib import Path
import runpy


PROJECT_ROOT = Path(__file__).resolve().parents[1]
runpy.run_path(str(PROJECT_ROOT / "streamlit_cloud.py"), run_name="__main__")
