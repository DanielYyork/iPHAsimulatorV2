"""Docker health probe; Render uses the same path through its HTTP health check."""

import os
from urllib.request import urlopen

port = int(os.environ.get("PORT", "8501"))
with urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=4) as response:
    if response.status != 200 or response.read().strip() != b"ok":
        raise SystemExit("Streamlit health check failed")
