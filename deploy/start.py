"""Prepare container storage, drop privileges and launch the existing GUI."""

import argparse
import os
from pathlib import Path
import pwd
import shutil
import sys
import tempfile


def prepare_database(seed: Path, destination: Path) -> None:
    """Stage a new database; never overwrite workshop data on restart."""
    if destination.exists() and any(destination.iterdir()):
        if not (destination / "residue_codes.csv").is_file():
            raise RuntimeError(
                f"{destination} contains files but no residue_codes.csv. "
                "Inspect the existing disk contents before starting the app."
            )
        return
    if not (seed / "residue_codes.csv").is_file() or not (seed / "PHA_types").is_dir():
        raise RuntimeError(f"Incomplete chemistry seed: {seed}")
    # A mounted empty directory cannot be renamed, so stage the seed inside it.
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ipha-seed-", dir=destination) as stage:
        staged = Path(stage) / "database"
        shutil.copytree(seed, staged)
        # The GUI creates its own empty build/system registries as needed.
        (staged / "polymer_smiles.csv").write_text("polymer_name,smiles\n")
        # Publish the registry last. An interrupted copy cannot be mistaken for
        # a complete database on the next startup; partial data is never erased.
        items = sorted(staged.iterdir(), key=lambda item: item.name == "residue_codes.csv")
        for item in items:
            item.rename(destination / item.name)


def prepare_storage(app_root: Path, seed_root: Path) -> None:
    data = app_root / "structure_database"
    # The mount itself can be root-owned even though its previous files are not.
    # Only these top-level directories change owner, never existing science data.
    if os.geteuid() == 0:
        user = pwd.getpwnam(os.environ.get("MAMBA_USER", "mambauser"))
        for path in (data, app_root / "md_simulation_scripts"):
            path.mkdir(parents=True, exist_ok=True)
            os.chown(path, user.pw_uid, user.pw_gid)
        os.setgroups([])
        os.setgid(user.pw_gid)
        os.setuid(user.pw_uid)
        os.environ["HOME"] = user.pw_dir
    prepare_database(seed_root, data)
    (app_root / "md_simulation_scripts").mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    prepare_storage(root, Path("/opt/ipha-seed"))
    if args.prepare_only:
        return
    port = int(os.environ.get("PORT", "8501"))
    if not 1 <= port <= 65535:
        raise ValueError("PORT must be between 1 and 65535")
    os.chdir(root)
    print(f"Starting iPHAsimulator on port {port} with {sys.executable}", flush=True)
    os.execv(sys.executable, [
        sys.executable, "-m", "streamlit", "run", str(root / "pha_gui.py"),
        "--server.address=0.0.0.0", f"--server.port={port}",
        "--server.headless=true", "--browser.gatherUsageStats=false",
    ])


if __name__ == "__main__":
    main()
