"""Read-only helpers for CHARMM-GUI / CGenFF outputs (Route C)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


RESI_PENALTY_PATTERN = re.compile(
    r"param penalty=\s*([-\d.]+)\s*;\s*charge penalty=\s*([-\d.]+)"
)
PARAMETER_PENALTY_PATTERN = re.compile(r"penalty=\s*([-\d.]+)")
PROGRAM_VERSION_PATTERN = re.compile(r"CGenFF\) program version\s+(\S+)")
FILES_VERSION_PATTERN = re.compile(r"parameter files version\s+(\S+)")

TOPOLOGY_FILES = ("lig.rtf", "lig_g.rtf")
PARAMETER_FILE = "lig.prm"


@dataclass(frozen=True)
class PenaltyEntry:
    """One penalised atom charge or bonded parameter."""

    source: str
    label: str
    penalty: float


@dataclass(frozen=True)
class CGenFFPenaltyReport:
    """Penalty scores parsed from a Ligand Reader & Modeler ``lig/`` folder."""

    lig_dir: Path
    files_read: tuple[str, ...]
    program_version: str | None
    parameter_files_version: str | None
    resi_penalties: dict[str, tuple[float, float]]
    charge_penalties: tuple[PenaltyEntry, ...]
    parameter_penalties: tuple[PenaltyEntry, ...]

    @property
    def max_charge_penalty(self) -> float | None:
        values = [entry.penalty for entry in self.charge_penalties]
        values += [charge for _, charge in self.resi_penalties.values()]
        return max(values, default=None)

    @property
    def max_parameter_penalty(self) -> float | None:
        values = [entry.penalty for entry in self.parameter_penalties]
        values += [param for param, _ in self.resi_penalties.values()]
        return max(values, default=None)

    def top_charge_penalties(self, n: int = 5) -> tuple[PenaltyEntry, ...]:
        return _top(self.charge_penalties, n)

    def top_parameter_penalties(self, n: int = 5) -> tuple[PenaltyEntry, ...]:
        return _top(self.parameter_penalties, n)


def penalty_category(penalty: float | None) -> str:
    """Classify a CGenFF penalty with the thresholds printed in CGenFF output."""

    if penalty is None:
        return "not reported"
    if penalty < 10:
        return "fair analogy"
    if penalty <= 50:
        return "basic validation recommended"
    return "extensive validation required"


def read_cgenff_penalties(lig_dir: str | Path) -> CGenFFPenaltyReport:
    """Parse penalty scores from ``lig.rtf``, ``lig_g.rtf`` and ``lig.prm``.

    Missing files are skipped; at least one of them must exist.
    """

    lig_dir = Path(lig_dir).expanduser()
    files_read: list[str] = []
    program_version = None
    parameter_files_version = None
    resi_penalties: dict[str, tuple[float, float]] = {}
    charge_penalties: list[PenaltyEntry] = []
    parameter_penalties: list[PenaltyEntry] = []

    for filename in (*TOPOLOGY_FILES, PARAMETER_FILE):
        path = lig_dir / filename
        if not path.is_file():
            continue
        files_read.append(filename)
        for line in path.read_text(errors="replace").splitlines():
            if line.startswith("*"):
                program_version = program_version or _search(PROGRAM_VERSION_PATTERN, line)
                parameter_files_version = parameter_files_version or _search(
                    FILES_VERSION_PATTERN, line
                )
                continue
            if line.lstrip().startswith("!") or "!" not in line:
                continue
            body, comment = line.split("!", 1)
            fields = body.split()
            if not fields:
                continue
            keyword = fields[0].upper()
            if filename in TOPOLOGY_FILES and keyword == "RESI":
                match = RESI_PENALTY_PATTERN.search(comment)
                if match:
                    resi_penalties[filename] = (float(match[1]), float(match[2]))
            elif filename in TOPOLOGY_FILES and keyword == "ATOM" and len(fields) >= 4:
                penalty = _leading_float(comment)
                if penalty is not None:
                    label = f"{fields[1]} ({fields[2]}, q={fields[3]})"
                    charge_penalties.append(PenaltyEntry(filename, label, penalty))
            elif filename == PARAMETER_FILE:
                match = PARAMETER_PENALTY_PATTERN.search(comment)
                if match:
                    label = " ".join(fields)
                    parameter_penalties.append(PenaltyEntry(filename, label, float(match[1])))

    if not files_read:
        expected = ", ".join((*TOPOLOGY_FILES, PARAMETER_FILE))
        raise FileNotFoundError(f"None of {expected} found in {lig_dir}")

    return CGenFFPenaltyReport(
        lig_dir=lig_dir,
        files_read=tuple(files_read),
        program_version=program_version,
        parameter_files_version=parameter_files_version,
        resi_penalties=resi_penalties,
        charge_penalties=tuple(charge_penalties),
        parameter_penalties=tuple(parameter_penalties),
    )


def format_cgenff_penalty_report(report: CGenFFPenaltyReport, top: int = 5) -> str:
    """Return a short plain-text summary suitable for notebook output."""

    lines = [
        f"Folder: {report.lig_dir}",
        f"Files read: {', '.join(report.files_read)}",
        f"CGenFF program {report.program_version or '?'}; "
        f"topology/parameter files {report.parameter_files_version or '?'}",
    ]
    for filename, (param, charge) in report.resi_penalties.items():
        lines.append(f"{filename} RESI: param penalty {param:g}; charge penalty {charge:g}")
    for title, value, entries in (
        ("charge", report.max_charge_penalty, report.top_charge_penalties(top)),
        ("parameter", report.max_parameter_penalty, report.top_parameter_penalties(top)),
    ):
        shown = "n/a" if value is None else f"{value:g}"
        lines.append(f"Max {title} penalty: {shown} ({penalty_category(value)})")
        lines.extend(f"  {entry.penalty:>7g}  {entry.source}: {entry.label}" for entry in entries)
    return "\n".join(lines)


def _top(entries: tuple[PenaltyEntry, ...], n: int) -> tuple[PenaltyEntry, ...]:
    return tuple(sorted(entries, key=lambda entry: entry.penalty, reverse=True)[:n])


def _search(pattern: re.Pattern[str], text: str) -> str | None:
    match = pattern.search(text)
    return match[1] if match else None


def _leading_float(text: str) -> float | None:
    fields = text.split()
    try:
        return float(fields[0]) if fields else None
    except ValueError:
        return None
