"""Run ``python run_enzyme_contacts.py GK13_P3HO_4.yaml`` from this folder."""

from pathlib import Path
import sys

# Support both repository-root and script-folder invocation in a source checkout.
source = Path(__file__).resolve().parents[2] / "src"
if str(source) not in sys.path:
    sys.path.insert(0, str(source))

from iphasimulator.analysis_contacts import main


if __name__ == "__main__":
    raise SystemExit(main())
