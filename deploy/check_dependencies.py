"""Run pip check, verifying the one known Conda/PyPI Open Babel mismatch.

ACPYPE's Conda recipe depends on `openbabel`; its Python metadata still names
`openbabel-wheel`. Accept that exact mismatch only after testing the Conda
installation. Every other pip check error remains fatal.
"""

import json
from pathlib import Path
import subprocess
import sys


CONDA_OPENBABEL_MISMATCH = (
    "acpype 2026.9.4 requires openbabel-wheel, which is not installed."
)


def conda_record(prefix, name):
    records = [json.loads(path.read_text())
               for path in (prefix / "conda-meta").glob(f"{name}-*.json")]
    records = [record for record in records if record.get("name") == name]
    if len(records) != 1:
        raise RuntimeError(f"Expected one Conda record for {name} in {prefix}")
    return records[0]


def verify_conda_openbabel():
    from openbabel import openbabel as ob

    prefix = Path(sys.prefix).resolve()
    if conda_record(prefix, "acpype")["version"] != "2026.9.4":
        raise RuntimeError("The Open Babel exception only covers Conda ACPYPE 2026.9.4")
    record = conda_record(prefix, "openbabel")
    module = Path(ob.__file__).resolve().relative_to(prefix).as_posix()
    binary = prefix / "bin" / "obabel"
    if module not in record["files"] or "bin/obabel" not in record["files"]:
        raise RuntimeError("Open Babel bindings and executable must come from this Conda environment")

    conversion = ob.OBConversion()
    molecule = ob.OBMol()
    if not conversion.SetInFormat("smi") or not conversion.ReadString(molecule, "CCO"):
        raise RuntimeError("Open Babel Python SMILES conversion failed")
    if molecule.NumAtoms() != 3:
        raise RuntimeError("Open Babel Python conversion produced the wrong atom count")
    result = subprocess.run(
        [str(binary), "-ismi", "-osmi"], input="CCO\n", text=True,
        capture_output=True, check=True, timeout=30,
    )
    if not result.stdout.split() or result.stdout.split()[0] != "CCO":
        raise RuntimeError("Open Babel executable SMILES conversion failed")
    print(f"Verified Conda Open Babel {record['version']}: Python and executable conversions passed.")


def main():
    result = subprocess.run(
        [sys.executable, "-m", "pip", "check"], text=True, capture_output=True,
    )
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    if result.returncode == 0:
        return 0
    if (result.returncode != 1 or result.stderr.strip()
            or result.stdout.strip() != CONDA_OPENBABEL_MISMATCH):
        return result.returncode
    try:
        verify_conda_openbabel()
    except Exception as exc:
        print(f"Conda Open Babel verification failed: {exc}", file=sys.stderr)
        return 1
    print("Dependency check passed with the verified Conda Open Babel provider.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
