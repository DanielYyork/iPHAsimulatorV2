"""Helpers for CHARMM-GUI / CGenFF outputs (CHARMM/GROMACS route).

Nothing here modifies a CHARMM-GUI download. The ``prepare_*`` functions write
only into a new folder chosen by the caller, and ``gmx`` runs only when the
``run_*`` functions are called explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import subprocess

from iphasimulator.simulation_gromacs_runner import (
    INCLUDE_PATTERN,
    validate_gromacs_coordinate_topology_counts,
    validate_gromacs_run_folder,
)


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


# ---------------------------------------------------------------------------
# CHARMM-GUI GROMACS systems (notebooks 06C and 06D)
# ---------------------------------------------------------------------------

CHARMM_GROMACS_TEMPLATE_DIR = Path(__file__).resolve().parent / "data" / "charmm_gromacs"
CHARMM_GROMACS_MDP_FILES = (
    "step6.0_minimization.mdp",
    "step6.1_nvt.mdp",
    "step6.2_npt.mdp",
    "step7_production.mdp",
)
LOCAL_SCRIPT = "run_step6_local.sh"
HPC_SCRIPT = "run_hpc_equilibration_production.slurm"
HPC_SCRIPT_TEMPLATE = HPC_SCRIPT + ".template"
CHARMM_GUI_INDEX_GROUPS = ("SOLU", "SOLV")

CHARMM_DEFAULTS = ("1", "2", "yes", 1.0, 1.0)
CGENFF_TYPE_PATTERN = re.compile(r"^[A-Z]G[A-Z0-9]+$")
WATER_MOLECULES = ("TIP3", "SOL", "HOH", "WAT")
ION_MOLECULES = ("SOD", "CLA", "POT", "NA", "CL")

DRY_COORDINATES = "lig.gro"
DRY_MDP = "dry_minimization.mdp"
DRY_DEFFNM = "dry_minimization"
DRY_MINIMIZATION_MDP = """\
; Dry (no solvent) steepest-descent minimisation of the CGenFF polymer.
; Sanity check only: GROMACS' Verlet scheme needs a periodic box, so
; the molecule sits alone in a box larger than twice the cut-off.
integrator              = steep
emtol                   = 100.0
emstep                  = 0.01
nsteps                  = 5000
nstlist                 = 10
cutoff-scheme           = Verlet
pbc                     = xyz
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
constraints             = h-bonds
constraint_algorithm    = LINCS
"""

POTENTIAL_ENERGY_PATTERN = re.compile(r"Potential Energy\s*=\s*([-+\d.eE]+)")
MAXIMUM_FORCE_PATTERN = re.compile(r"Maximum force\s*=\s*([-+\d.eE]+)")
ELEMENT_BY_MASS = {1: "H", 12: "C", 14: "N", 16: "O", 19: "F", 31: "P", 32: "S", 35: "Cl"}


@dataclass(frozen=True)
class ItpAtom:
    name: str
    type: str
    charge: float
    mass: float


@dataclass(frozen=True)
class ItpMolecule:
    """One ``[ moleculetype ]`` with its atoms and bonds (0-based indices)."""

    name: str
    atoms: tuple[ItpAtom, ...]
    bonds: tuple[tuple[int, int], ...]

    @property
    def charge(self) -> float:
        return sum(atom.charge for atom in self.atoms)


@dataclass(frozen=True)
class GroAtom:
    residue_number: int
    residue_name: str
    atom_name: str
    position: tuple[float, float, float]


@dataclass(frozen=True)
class CheckResult:
    """One PASS/FAIL/SKIP line with the evidence behind it."""

    name: str
    status: str
    evidence: str

    @property
    def passed(self) -> bool:
        return self.status != "FAIL"


@dataclass(frozen=True)
class MinimizationResult:
    log_path: Path
    potential_energy: float | None
    maximum_force: float | None


def format_checks(results: list[CheckResult]) -> str:
    """Return one line per check, followed by an overall verdict."""

    lines = [f"[{result.status:<4}] {result.name}: {result.evidence}" for result in results]
    failed = [result.name for result in results if result.status == "FAIL"]
    lines.append("Overall: " + ("FAIL (" + ", ".join(failed) + ")" if failed else "PASS"))
    return "\n".join(lines)


def read_topology_molecules(topology_path: str | Path) -> dict[str, ItpMolecule]:
    """Read every ``[ moleculetype ]`` in a topology and its ``#include`` files."""

    molecules: dict[str, ItpMolecule] = {}
    for source in _topology_files(Path(topology_path)):
        name: str | None = None
        atoms: list[ItpAtom] = []
        bonds: list[tuple[int, int]] = []
        section = None
        for fields in _itp_lines(source):
            if fields[0] == "[":
                if fields[1] == "moleculetype" and name is not None:
                    molecules[name] = ItpMolecule(name, tuple(atoms), tuple(bonds))
                    name, atoms, bonds = None, [], []
                section = fields[1]
                continue
            if section == "moleculetype" and name is None:
                name = fields[0]
            elif section == "atoms" and len(fields) >= 8:
                atoms.append(ItpAtom(fields[4], fields[1], float(fields[6]), float(fields[7])))
            elif section == "bonds" and len(fields) >= 2:
                bonds.append((int(fields[0]) - 1, int(fields[1]) - 1))
        if name is not None:
            molecules[name] = ItpMolecule(name, tuple(atoms), tuple(bonds))
    return molecules


def read_topology_defaults(topology_path: str | Path) -> tuple[str, ...] | None:
    """Return the first ``[ defaults ]`` line found in the topology or its includes."""

    for source in _topology_files(Path(topology_path)):
        in_defaults = False
        for fields in _itp_lines(source):
            if fields[0] == "[":
                in_defaults = fields[1] == "defaults"
            elif in_defaults:
                return tuple(fields)
    return None


def read_topology_molecule_counts(topology_path: str | Path) -> dict[str, int]:
    """Read the ``[ molecules ]`` table in order."""

    counts: dict[str, int] = {}
    in_molecules = False
    for fields in _itp_lines(Path(topology_path)):
        if fields[0] == "[":
            in_molecules = fields[1] == "molecules"
        elif in_molecules and len(fields) >= 2:
            counts[fields[0]] = counts.get(fields[0], 0) + int(fields[1])
    return counts


def read_gro(gro_path: str | Path) -> tuple[tuple[GroAtom, ...], tuple[float, ...]]:
    """Read atoms (positions in nm) and box vectors from a ``.gro`` file."""

    lines = Path(gro_path).read_text().splitlines()
    count = int(lines[1].strip())
    atoms = tuple(
        GroAtom(
            int(line[0:5]),
            line[5:10].strip(),
            line[10:15].strip(),
            (float(line[20:28]), float(line[28:36]), float(line[36:44])),
        )
        for line in lines[2 : 2 + count]
    )
    box = tuple(float(value) for value in lines[2 + count].split())
    return atoms, box


def check_charmm_defaults(topology_path: str | Path) -> CheckResult:
    name = "[ defaults ] CHARMM-style"
    defaults = read_topology_defaults(topology_path)
    if defaults is None:
        return CheckResult(name, "FAIL", "no [ defaults ] section found")
    shown = " ".join(defaults)
    try:
        ok = (
            tuple(defaults[:3]) == CHARMM_DEFAULTS[:3]
            and float(defaults[3]) == CHARMM_DEFAULTS[3]
            and float(defaults[4]) == CHARMM_DEFAULTS[4]
        )
    except (IndexError, ValueError):
        ok = False
    return CheckResult(name, "PASS" if ok else "FAIL", f"{shown} (expected 1 2 yes 1.0 1.0)")


def check_includes(folder: str | Path, topology_name: str = "topol.top") -> CheckResult:
    validation = validate_gromacs_run_folder(folder, topology_name=topology_name)
    if validation.missing_files:
        missing = ", ".join(path.name for path in validation.missing_files)
        return CheckResult("#include files resolve", "FAIL", f"missing: {missing}")
    names = ", ".join(path.name for path in validation.included_files) or "none (standalone)"
    return CheckResult("#include files resolve", "PASS", names)


def check_atom_counts(
    folder: str | Path,
    coordinate_name: str,
    topology_name: str = "topol.top",
) -> CheckResult:
    name = f"{coordinate_name} atoms == topology atoms"
    validation = validate_gromacs_coordinate_topology_counts(
        folder, coordinate_name=coordinate_name, topology_name=topology_name
    )
    if not validation.can_compare:
        return CheckResult(name, "FAIL", "a molecule in [ molecules ] has no [ moleculetype ]")
    status = "PASS" if validation.valid else "FAIL"
    return CheckResult(
        name,
        status,
        f"{validation.coordinate_atom_count} vs {validation.expected_atom_count}",
    )


def check_molecule_counts(
    topology_path: str | Path,
    *,
    ligand: str = "LIG",
    require_water: bool,
) -> CheckResult:
    counts = read_topology_molecule_counts(topology_path)
    upper = {key.upper(): value for key, value in counts.items()}
    evidence = ", ".join(f"{key} {value}" for key, value in counts.items())
    problems = []
    if upper.get(ligand.upper(), 0) < 1:
        problems.append(f"no {ligand}")
    if require_water and not any(upper.get(water, 0) for water in WATER_MOLECULES):
        problems.append("no water")
    protein = [key for key in counts if key.upper().startswith("PRO")]
    evidence += f" | protein chains: {', '.join(protein) if protein else 'none'}"
    if problems:
        evidence += " | " + "; ".join(problems)
    return CheckResult("molecule counts", "FAIL" if problems else "PASS", evidence)


def check_net_charge(topology_path: str | Path, tolerance: float = 1e-3) -> CheckResult:
    molecules = read_topology_molecules(topology_path)
    counts = read_topology_molecule_counts(topology_path)
    missing = [name for name in counts if name not in molecules]
    if missing:
        return CheckResult("net charge ~ 0", "FAIL", f"no [ moleculetype ] for {', '.join(missing)}")
    parts = {name: count * molecules[name].charge for name, count in counts.items()}
    total = sum(parts.values())
    detail = ", ".join(f"{name} {value:+.3f}" for name, value in parts.items())
    status = "PASS" if abs(total) <= tolerance else "FAIL"
    return CheckResult("net charge ~ 0", status, f"{total:+.4f} e ({detail})")


def check_cgenff_types(topology_path: str | Path, ligand: str = "LIG") -> CheckResult:
    name = f"{ligand} atom types CGenFF-style"
    molecule = _find_molecule(read_topology_molecules(topology_path), ligand)
    if molecule is None:
        return CheckResult(name, "FAIL", f"no [ moleculetype ] {ligand}")
    types = sorted({atom.type for atom in molecule.atoms})
    offenders = [atom_type for atom_type in types if not CGENFF_TYPE_PATTERN.match(atom_type)]
    if offenders:
        return CheckResult(name, "FAIL", f"not CGenFF-style: {', '.join(offenders)}")
    return CheckResult(name, "PASS", f"{len(types)} types, e.g. {', '.join(types[:6])}")


def check_ligand_atom_order(
    gro_path: str | Path,
    topology_path: str | Path,
    ligand: str = "LIG",
) -> CheckResult:
    name = f"{ligand} atom names/order match {ligand}.itp"
    molecule = _find_molecule(read_topology_molecules(topology_path), ligand)
    if molecule is None:
        return CheckResult(name, "FAIL", f"no [ moleculetype ] {ligand}")
    atoms, _ = read_gro(gro_path)
    ligand_atoms = _first_residue(atoms, ligand)
    gro_names = [atom.atom_name for atom in ligand_atoms]
    itp_names = [atom.name for atom in molecule.atoms]
    if gro_names == itp_names:
        return CheckResult(name, "PASS", f"{len(itp_names)} atoms in the same order")
    mismatches = sum(a != b for a, b in zip(gro_names, itp_names))
    return CheckResult(
        name,
        "FAIL",
        f"{len(gro_names)} coordinates vs {len(itp_names)} topology atoms; {mismatches} name mismatches",
    )


def ligand_stereocentres(
    gro_path: str | Path,
    topology_path: str | Path,
    sdf_path: str | Path,
    ligand: str = "LIG",
) -> list[tuple[str, str]]:
    """CIP labels of the ligand in ``gro_path``, using the SDF only for bond orders.

    The ligand graph comes from the ``.itp`` bonds and its geometry from the
    ``.gro`` coordinates (made whole across the box), so the labels describe the
    structure GROMACS will simulate, not the uploaded SDF.
    """

    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    molecule = _find_molecule(read_topology_molecules(topology_path), ligand)
    if molecule is None:
        raise ValueError(f"No [ moleculetype ] {ligand} in {topology_path}")
    atoms, box = read_gro(gro_path)
    ligand_atoms = _first_residue(atoms, ligand)
    if len(ligand_atoms) != len(molecule.atoms):
        raise ValueError(
            f"{ligand}: {len(ligand_atoms)} coordinates vs {len(molecule.atoms)} topology atoms"
        )
    positions = _make_whole([atom.position for atom in ligand_atoms], molecule.bonds, box[:3])

    editable = Chem.RWMol()
    for atom in molecule.atoms:
        editable.AddAtom(Chem.Atom(_element(atom.mass)))
    for i, j in sorted({tuple(sorted(bond)) for bond in molecule.bonds}):
        editable.AddBond(i, j, Chem.BondType.SINGLE)
    mol = editable.GetMol()
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for index, (x, y, z) in enumerate(positions):
        conformer.SetAtomPosition(index, Point3D(x * 10, y * 10, z * 10))
    mol.AddConformer(conformer)

    template = Chem.MolFromMolFile(str(sdf_path), removeHs=False)
    if template is None:
        raise ValueError(f"RDKit could not read {sdf_path}")
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    centres = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    return [(molecule.atoms[index].name, label) for index, label in centres]


def check_ligand_stereocentres(
    gro_path: str | Path,
    topology_path: str | Path,
    sdf_path: str | Path | None,
    ligand: str = "LIG",
    expected: str = "R",
) -> CheckResult:
    name = f"{ligand} stereocentres all {expected}"
    if sdf_path is None:
        return CheckResult(name, "SKIP", "no SDF given (needed for bond orders)")
    try:
        centres = ligand_stereocentres(gro_path, topology_path, sdf_path, ligand)
    except Exception as exc:  # report, do not crash the notebook
        return CheckResult(name, "FAIL", f"could not map {ligand} to the SDF: {exc}")
    shown = ", ".join(f"{atom} {label}" for atom, label in centres) or "none found"
    ok = bool(centres) and all(label == expected for _, label in centres)
    return CheckResult(name, "PASS" if ok else "FAIL", shown)


def check_gromacs_run_files(folder: str | Path) -> CheckResult:
    folder = Path(folder)
    required = ("step5_input.gro", "topol.top", "index.ndx", *CHARMM_GROMACS_MDP_FILES, LOCAL_SCRIPT, HPC_SCRIPT)
    missing = [name for name in required if not (folder / name).is_file()]
    groups = _index_groups(folder / "index.ndx") if (folder / "index.ndx").is_file() else set()
    missing_groups = [group for group in CHARMM_GUI_INDEX_GROUPS if group not in groups]
    problems = []
    if missing:
        problems.append("missing " + ", ".join(missing))
    if missing_groups:
        problems.append("index.ndx lacks " + ", ".join(missing_groups))
    if problems:
        return CheckResult("run files and index groups", "FAIL", "; ".join(problems))
    return CheckResult(
        "run files and index groups",
        "PASS",
        f"{len(required)} files; index groups {', '.join(CHARMM_GUI_INDEX_GROUPS)} present",
    )


def validate_dry_ligand_folder(
    folder: str | Path,
    *,
    sdf_path: str | Path | None = None,
    ligand: str = "LIG",
) -> list[CheckResult]:
    """Read-only checks for a folder made by :func:`prepare_dry_ligand_folder`."""

    folder = Path(folder)
    topology = folder / "topol.top"
    gro = folder / DRY_COORDINATES
    return [
        check_charmm_defaults(topology),
        check_includes(folder),
        check_atom_counts(folder, DRY_COORDINATES),
        check_molecule_counts(topology, ligand=ligand, require_water=False),
        check_net_charge(topology),
        check_cgenff_types(topology, ligand),
        check_ligand_atom_order(gro, topology, ligand),
        check_ligand_stereocentres(gro, topology, sdf_path, ligand),
    ]


def validate_charmm_gui_run_folder(
    folder: str | Path,
    *,
    sdf_path: str | Path | None = None,
    ligand: str = "LIG",
) -> list[CheckResult]:
    """Read-only checks for a folder made by :func:`prepare_charmm_gui_run_folder`."""

    folder = Path(folder)
    topology = folder / "topol.top"
    gro = folder / "step5_input.gro"
    return [
        check_charmm_defaults(topology),
        check_includes(folder),
        check_atom_counts(folder, "step5_input.gro"),
        check_molecule_counts(topology, ligand=ligand, require_water=True),
        check_net_charge(topology),
        check_cgenff_types(topology, ligand),
        check_ligand_atom_order(gro, topology, ligand),
        check_ligand_stereocentres(gro, topology, sdf_path, ligand),
        check_gromacs_run_files(folder),
    ]


def prepare_dry_ligand_folder(
    ligand_reader_dir: str | Path,
    output_dir: str | Path,
    *,
    padding_nm: float = 1.7,
) -> Path:
    """Build a solvent-free GROMACS folder from a Ligand Reader & Modeler download.

    Copies ``gromacs/topol.top``, ``charmm36.itp`` and ``LIG.itp``, writes
    ``lig.gro`` from ``ligandrm.pdb`` centred in a cubic box of the molecule's
    extent plus ``2 * padding_nm``, and writes ``dry_minimization.mdp``.
    """

    source = Path(ligand_reader_dir).expanduser().resolve()
    gromacs_dir = source / "gromacs"
    pdb_path = source / "ligandrm.pdb"
    for required in (gromacs_dir / "topol.top", pdb_path):
        if not required.is_file():
            raise FileNotFoundError(f"Not a Ligand Reader & Modeler download: missing {required}")
    output = _new_output_folder(output_dir, source)

    topology = gromacs_dir / "topol.top"
    shutil.copyfile(topology, output / "topol.top")
    for line in topology.read_text().splitlines():
        match = INCLUDE_PATTERN.match(line)
        if match:
            target = output / match.group(1)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(gromacs_dir / match.group(1), target)

    records = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            records.append(
                (
                    line[17:21].strip(),
                    line[12:16].strip(),
                    tuple(float(line[start : start + 8]) / 10 for start in (30, 38, 46)),
                )
            )
    low = [min(position[axis] for _, _, position in records) for axis in range(3)]
    high = [max(position[axis] for _, _, position in records) for axis in range(3)]
    edge = max(h - l for h, l in zip(high, low)) + 2 * padding_nm
    shift = [edge / 2 - (h + l) / 2 for h, l in zip(high, low)]
    lines = ["Dry CGenFF polymer from Ligand Reader & Modeler", f"{len(records):5d}"]
    for index, (residue, atom_name, position) in enumerate(records, start=1):
        x, y, z = (value + offset for value, offset in zip(position, shift))
        lines.append(f"{1:5d}{residue:<5.5s}{atom_name:>5.5s}{index % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append(f"{edge:10.5f}{edge:10.5f}{edge:10.5f}")
    (output / DRY_COORDINATES).write_text("\n".join(lines) + "\n")
    (output / DRY_MDP).write_text(DRY_MINIMIZATION_MDP)
    return output


def prepare_charmm_gui_run_folder(
    solution_builder_dir: str | Path,
    run_dir: str | Path,
    *,
    template_dir: str | Path | None = None,
    job_name: str = "charmm_gromacs",
) -> Path:
    """Turn a Solution Builder download into a run folder, as in the production runs.

    Copies ``gromacs/`` (topology, ``toppar/``, ``index.ndx``, coordinates) into the
    new ``run_dir``, renames ``step3_input.gro`` to ``step5_input.gro``, leaves out
    CHARMM-GUI's ``step4.x``/``step5`` mdp files, and adds the CHARMM/GROMACS mdp files and
    scripts from ``template_dir`` (default: the packaged CHARMM/GROMACS templates).
    """

    source = Path(solution_builder_dir).expanduser().resolve()
    if (source / "gromacs").is_dir():
        source = source / "gromacs"
    if not (source / "step3_input.gro").is_file() or not (source / "topol.top").is_file():
        raise FileNotFoundError(
            f"Not a Solution Builder gromacs/ folder (needs step3_input.gro and topol.top): {source}"
        )
    templates = Path(template_dir).expanduser().resolve() if template_dir else CHARMM_GROMACS_TEMPLATE_DIR
    missing = [
        name
        for name in (*CHARMM_GROMACS_MDP_FILES, LOCAL_SCRIPT)
        if not (templates / name).is_file()
    ]
    if not ((templates / HPC_SCRIPT).is_file() or (templates / HPC_SCRIPT_TEMPLATE).is_file()):
        missing.append(HPC_SCRIPT)
    if missing:
        raise FileNotFoundError(f"Template folder {templates} lacks: {', '.join(missing)}")
    output = _new_output_folder(run_dir, source)

    for item in source.iterdir():
        if item.name == "step3_input.gro":
            shutil.copyfile(item, output / "step5_input.gro")
        elif item.is_dir():
            _copy_tree_contents(item, output / item.name)
        elif item.suffix == ".mdp" and item.name.startswith(("step4", "step5")):
            continue
        else:
            shutil.copyfile(item, output / item.name)

    for name in (*CHARMM_GROMACS_MDP_FILES, LOCAL_SCRIPT):
        shutil.copyfile(templates / name, output / name)
    if (templates / HPC_SCRIPT).is_file():
        shutil.copyfile(templates / HPC_SCRIPT, output / HPC_SCRIPT)
    else:
        text = (templates / HPC_SCRIPT_TEMPLATE).read_text().replace("{JOB_NAME}", job_name)
        (output / HPC_SCRIPT).write_text(text)
    (output / LOCAL_SCRIPT).chmod(0o755)
    return output


def run_dry_minimization(
    folder: str | Path,
    *,
    gmx: str = "gmx",
    runner=subprocess.run,
) -> MinimizationResult:
    """Run ``grompp`` + ``mdrun`` for the dry folder and return Epot and Fmax."""

    folder = Path(folder)
    commands = (
        [gmx, "grompp", "-f", DRY_MDP, "-c", DRY_COORDINATES, "-p", "topol.top", "-o", f"{DRY_DEFFNM}.tpr"],
        [gmx, "mdrun", "-deffnm", DRY_DEFFNM],
    )
    for command in commands:
        result = runner(command, cwd=folder, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            raise RuntimeError(
                f"{' '.join(command)} failed with return code {result.returncode}.\n"
                f"STDERR:\n{result.stderr}"
            )
    return read_minimization_result(folder / f"{DRY_DEFFNM}.log")


def read_minimization_result(log_path: str | Path) -> MinimizationResult:
    """Read the final potential energy and maximum force from an ``mdrun`` log."""

    log_path = Path(log_path)
    text = log_path.read_text(errors="replace")
    energies = POTENTIAL_ENERGY_PATTERN.findall(text)
    forces = MAXIMUM_FORCE_PATTERN.findall(text)
    return MinimizationResult(
        log_path=log_path,
        potential_energy=float(energies[-1]) if energies else None,
        maximum_force=float(forces[-1]) if forces else None,
    )


def _topology_files(topology_path: Path) -> list[Path]:
    files = [topology_path]
    for line in topology_path.read_text().splitlines():
        match = INCLUDE_PATTERN.match(line)
        if match:
            include = Path(match.group(1))
            include = include if include.is_absolute() else topology_path.parent / include
            if include.is_file():
                files.append(include)
    return files


def _itp_lines(path: Path):
    """Yield split, comment-free lines; section headers come back as ``["[", name]``."""

    for raw_line in path.read_text(errors="replace").splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            yield ["[", line.strip("[]").strip().lower()]
        else:
            yield line.split()


def _find_molecule(molecules: dict[str, ItpMolecule], name: str) -> ItpMolecule | None:
    for key, molecule in molecules.items():
        if key.upper() == name.upper():
            return molecule
    return None


def _first_residue(atoms: tuple[GroAtom, ...], residue_name: str) -> list[GroAtom]:
    selected: list[GroAtom] = []
    for atom in atoms:
        if atom.residue_name.upper() != residue_name.upper():
            if selected:
                break
            continue
        if selected and atom.residue_number != selected[0].residue_number:
            break
        selected.append(atom)
    return selected


def _make_whole(positions, bonds, box):
    """Unwrap a molecule split across a rectangular periodic box by walking its bonds."""

    positions = [list(position) for position in positions]
    neighbours: dict[int, list[int]] = {index: [] for index in range(len(positions))}
    for i, j in bonds:
        neighbours[i].append(j)
        neighbours[j].append(i)
    seen = {0}
    stack = [0]
    while stack:
        current = stack.pop()
        for neighbour in neighbours[current]:
            if neighbour in seen:
                continue
            for axis in range(3):
                length = box[axis] if axis < len(box) else 0.0
                if length > 0:
                    delta = positions[neighbour][axis] - positions[current][axis]
                    positions[neighbour][axis] -= length * round(delta / length)
            seen.add(neighbour)
            stack.append(neighbour)
    return [tuple(position) for position in positions]


def _element(mass: float) -> str:
    nearest = min(ELEMENT_BY_MASS, key=lambda value: abs(value - mass))
    return ELEMENT_BY_MASS[nearest]


def _index_groups(index_path: Path) -> set[str]:
    return {
        line.strip().strip("[]").strip()
        for line in index_path.read_text().splitlines()
        if line.strip().startswith("[")
    }


def _copy_tree_contents(source: Path, target: Path) -> None:
    """Copy file contents only, so read-only downloads give an editable copy."""

    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_dir():
            _copy_tree_contents(item, target / item.name)
        else:
            shutil.copyfile(item, target / item.name)


def _new_output_folder(output_dir: str | Path, source: Path) -> Path:
    output = Path(output_dir).expanduser().resolve()
    if output == source or source in output.parents:
        raise ValueError(f"Refusing to write inside the CHARMM-GUI download: {output}")
    if output.exists():
        raise FileExistsError(f"Output folder already exists; choose a new one: {output}")
    output.mkdir(parents=True)
    return output
