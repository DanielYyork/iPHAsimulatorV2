"""Inputs for an Amber/OpenMM enzyme + polymer system built from a docked complex.

The same complex PDB used for CHARMM-GUI (protein + polymer pose, e.g. from
docking) is split into a protein PDB for ``tleap``/ff19SB and a posed copy of the
05A GAFF2 polymer mol2, so both force-field routes start from the same structure.
Nothing here runs tleap or OpenMM.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


CHARMM_TO_AMBER_HISTIDINE = {"HSD": "HID", "HSE": "HIE", "HSP": "HIP"}
HISTIDINE_NAMES = ("HIS", "HSD", "HSE", "HSP", "HID", "HIE", "HIP")
SOLVENT_AND_IONS = ("HOH", "WAT", "TIP3", "SOL", "SOD", "CLA", "NA", "CL", "POT", "K")


@dataclass(frozen=True)
class ProteinPdbSummary:
    output_path: Path
    residues: int
    removed_hydrogens: int
    histidines: dict[int, str]


@dataclass(frozen=True)
class PoseTransferResult:
    output_mol2: Path
    heavy_atoms: int
    hydrogens_placed: int
    stereocentres: list[tuple[str, str]]


def read_charmm_histidine_states(protein_itp: str | Path) -> dict[int, str]:
    """Residue number -> Amber histidine name (HID/HIE/HIP) from a CHARMM-GUI ``PROA.itp``."""

    states: dict[int, str] = {}
    in_atoms = False
    for raw in Path(protein_itp).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            in_atoms = line.strip("[] ").lower() == "atoms"
            continue
        fields = line.split()
        if in_atoms and len(fields) >= 5 and fields[3] in CHARMM_TO_AMBER_HISTIDINE:
            states[int(fields[2])] = CHARMM_TO_AMBER_HISTIDINE[fields[3]]
    return states


def write_protein_pdb(
    complex_pdb: str | Path,
    output_pdb: str | Path,
    *,
    ligand_resname: str = "LIG",
    histidine_states: dict[int, str] | None = None,
) -> ProteinPdbSummary:
    """Write the protein part of a complex PDB for ``tleap`` (ff19SB).

    Drops the ligand, water, ions and hydrogens (tleap adds hydrogens), and renames
    histidines to HID/HIE/HIP from ``histidine_states`` (residue number -> name).
    Histidines without a state stay ``HIS``, which tleap treats as HIE.
    """

    histidine_states = histidine_states or {}
    lines: list[str] = []
    residues: set[tuple[str, int]] = set()
    histidines: dict[int, str] = {}
    removed = 0
    previous_chain = None
    for line in Path(complex_pdb).read_text().splitlines():
        if line.startswith("TER"):
            if lines and lines[-1] != "TER":
                lines.append("TER")
            previous_chain = None
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:21].strip()
        if resname.upper() in (ligand_resname.upper(), *SOLVENT_AND_IONS):
            continue
        if _pdb_element(line) == "H":
            removed += 1
            continue
        resnum = int(line[22:26])
        if resname in HISTIDINE_NAMES:
            resname = histidine_states.get(resnum, CHARMM_TO_AMBER_HISTIDINE.get(resname, resname))
            histidines[resnum] = resname
            line = f"{line[:17]}{resname:<3} {line[21:]}"
        if previous_chain is not None and line[21] != previous_chain:
            lines.append("TER")
        previous_chain = line[21]
        residues.add((line[21], resnum))
        lines.append(line)
    if not lines:
        raise ValueError(f"No protein atoms found in {complex_pdb}")
    if histidine_states and len({chain for chain, _ in residues}) > 1:
        raise ValueError("Residue-number histidine mapping requires a single protein chain")
    if lines[-1] != "TER":
        lines.append("TER")
    Path(output_pdb).write_text("\n".join(lines) + "\nEND\n")
    return ProteinPdbSummary(Path(output_pdb), len(residues), removed, histidines)


def pose_polymer_mol2(
    mol2_path: str | Path,
    template_sdf: str | Path,
    complex_pdb: str | Path,
    output_mol2: str | Path,
    *,
    ligand_resname: str = "LIG",
) -> PoseTransferResult:
    """Copy the docked polymer pose onto the 05A GAFF2 mol2.

    ``template_sdf`` is the antechamber input written by 05A (same atom order as the
    mol2; it carries bond orders and stereochemistry). Heavy atoms are matched to the
    pose by element and connectivity; hydrogens are placed from the posed heavy
    atoms. Atom names, types and charges in the mol2 are unchanged.
    """

    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Geometry import Point3D

    template = Chem.MolFromMolFile(str(template_sdf), removeHs=False)
    if template is None:
        raise ValueError(f"RDKit could not read {template_sdf}")
    mol2_lines, atom_rows = _read_mol2_atoms(Path(mol2_path))
    mol2_elements = [_mol2_element(row) for row in atom_rows]
    if mol2_elements != [atom.GetSymbol() for atom in template.GetAtoms()]:
        raise ValueError("The mol2 and template SDF atoms are not in the same order")
    expected_bonds = {tuple(sorted((b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1))) for b in template.GetBonds()}
    mol2_bonds = set()
    in_bonds = False
    for line in mol2_lines:
        if line.startswith("@<TRIPOS>"):
            in_bonds = line.strip() == "@<TRIPOS>BOND"
        elif in_bonds and line.strip():
            fields = line.split()
            mol2_bonds.add(tuple(sorted((int(fields[1]), int(fields[2])))))
    if mol2_bonds != expected_bonds:
        raise ValueError("The mol2 and template SDF bond connectivity/order do not agree")

    pose_lines = [
        line for line in Path(complex_pdb).read_text().splitlines()
        if line.startswith(("ATOM", "HETATM")) and line[17:21].strip().upper() == ligand_resname.upper()
        and _pdb_element(line) != "H"
    ]
    heavy_indices = [atom.GetIdx() for atom in template.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(pose_lines) != len(heavy_indices):
        raise ValueError(
            f"{ligand_resname}: {len(pose_lines)} heavy atoms in the complex vs "
            f"{len(heavy_indices)} in the 05A polymer; is it the same oligomer?"
        )

    pose = Chem.RWMol()
    conformer = Chem.Conformer(len(pose_lines))
    for index, line in enumerate(pose_lines):
        pose.AddAtom(Chem.Atom(_pdb_element(line)))
        conformer.SetAtomPosition(index, Point3D(float(line[30:38]), float(line[38:46]), float(line[46:54])))
    pose = pose.GetMol()
    pose.AddConformer(conformer)
    rdDetermineBonds.DetermineConnectivity(pose)
    pose.UpdatePropertyCache(strict=False)

    heavy_template = Chem.RemoveHs(template)
    graph = Chem.RWMol(heavy_template)
    for bond in graph.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in graph.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    graph = graph.GetMol()
    graph.UpdatePropertyCache(strict=False)
    match = pose.GetSubstructMatch(graph)
    if len(match) != graph.GetNumAtoms():
        raise ValueError(f"{ligand_resname} in the complex does not match the 05A polymer's connectivity")

    posed_heavy = Chem.Mol(heavy_template)
    posed_heavy.RemoveAllConformers()
    heavy_conformer = Chem.Conformer(posed_heavy.GetNumAtoms())
    for template_index, pose_index in enumerate(match):
        heavy_conformer.SetAtomPosition(template_index, pose.GetConformer().GetAtomPosition(pose_index))
    posed_heavy.AddConformer(heavy_conformer)
    with_h = Chem.AddHs(posed_heavy, addCoords=True)

    positions = {}
    for heavy_position, template_index in enumerate(heavy_indices):
        positions[template_index] = with_h.GetConformer().GetAtomPosition(heavy_position)
        template_h = sorted(n.GetIdx() for n in template.GetAtomWithIdx(template_index).GetNeighbors() if n.GetAtomicNum() == 1)
        placed_h = sorted(n.GetIdx() for n in with_h.GetAtomWithIdx(heavy_position).GetNeighbors() if n.GetAtomicNum() == 1)
        if len(template_h) != len(placed_h):
            raise ValueError(f"Hydrogen count differs on template atom {template_index}")
        for t_index, p_index in zip(template_h, placed_h):
            positions[t_index] = with_h.GetConformer().GetAtomPosition(p_index)

    posed = Chem.Mol(template)
    posed_conformer = posed.GetConformer()
    for index, point in positions.items():
        posed_conformer.SetAtomPosition(index, point)
    Chem.AssignStereochemistryFrom3D(posed)
    centres = Chem.FindMolChiralCenters(posed, includeUnassigned=True, useLegacyImplementation=False)
    stereocentres = [(atom_rows[index][1], label) for index, label in centres]

    template_centres = Chem.FindMolChiralCenters(template, includeUnassigned=True, useLegacyImplementation=False)
    if not centres or len(centres) != len(template_centres) or any(label != "R" for _, label in centres):
        raise ValueError("The docked PHA does not preserve all template R stereocentres")

    _write_mol2_with_coordinates(mol2_lines, atom_rows, [positions[i] for i in range(template.GetNumAtoms())], Path(output_mol2))
    return PoseTransferResult(Path(output_mol2), len(heavy_indices), template.GetNumAtoms() - len(heavy_indices), stereocentres)


def _pdb_element(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if not element:
        element = re.sub(r"[^A-Za-z]", "", line[12:16])[:1]
    return element.capitalize()


def _read_mol2_atoms(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text().splitlines()
    rows: list[list[str]] = []
    in_atoms = False
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.strip() == "@<TRIPOS>ATOM"
            continue
        if in_atoms and line.strip():
            rows.append(line.split())
    return lines, rows


def _mol2_element(row: list[str]) -> str:
    return re.match(r"[A-Za-z][a-z]?", row[1]).group(0).capitalize() if row[1][1:2].islower() else row[1][0].upper()


def _write_mol2_with_coordinates(lines, rows, positions, output: Path) -> None:
    out: list[str] = []
    in_atoms = False
    atom_index = 0
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.strip() == "@<TRIPOS>ATOM"
            out.append(line)
            continue
        if in_atoms and line.strip():
            fields = rows[atom_index]
            point = positions[atom_index]
            rest = " ".join(fields[5:])
            out.append(f"{int(fields[0]):>7d} {fields[1]:<8s}{point.x:>10.4f}{point.y:>10.4f}{point.z:>10.4f} {rest}")
            atom_index += 1
        else:
            out.append(line)
    output.write_text("\n".join(out) + "\n")
