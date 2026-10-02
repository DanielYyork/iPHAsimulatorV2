"""Polymer-only PDB for docking, from a GROMACS ``.gro`` and its topology (notebook 11).

Only the polymer is written: the atoms with one residue name, ``PHA`` for systems
built with 06A (and the polymer benchmark) or ``LIG`` for CHARMM-GUI systems. Before
writing, the polymer must be whole: a molecule split across the periodic box would
give the docking program a broken structure.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from iphasimulator.charmmgui_import import ELEMENT_BY_MASS, ItpMolecule, read_gro, read_topology_molecules


MAX_BOND_NM = 0.3
# Hydrogens weigh about 1 u, or about 3 u with hydrogen mass repartitioning.
HYDROGEN_MAX_MASS = 4.0
WHOLE_MOLECULE_ADVICE = (
    "Use a whole-molecule frame instead, e.g. the representative frame from notebook 08 "
    "or a frame written with `gmx trjconv -pbc mol`."
)


@dataclass(frozen=True)
class PolymerWholeCheck:
    """Bonded heavy-atom distances of the polymer (nm)."""

    residue_name: str
    atoms: int
    heavy_bonds: int
    longest_bond: tuple[str, str, float]
    broken_bonds: tuple[tuple[str, str, float], ...]
    max_bond_nm: float = MAX_BOND_NM

    @property
    def is_whole(self) -> bool:
        return not self.broken_bonds


@dataclass(frozen=True)
class PolymerPdbExport:
    pdb_path: Path
    atoms: int
    check: PolymerWholeCheck


def check_polymer_whole(
    gro_path: str | Path,
    topology_path: str | Path,
    residue_name: str = "PHA",
    max_bond_nm: float = MAX_BOND_NM,
) -> PolymerWholeCheck:
    """Measure every bonded heavy-atom pair of the polymer; whole if all are < ``max_bond_nm``.

    Bonds come from the polymer's ``[ moleculetype ]`` (named like its residue) in the
    topology or its ``#include`` files; coordinates come from the ``.gro`` file.
    """

    polymer, molecule = _polymer_atoms(gro_path, topology_path, residue_name)
    return _whole_check(polymer, molecule, residue_name, max_bond_nm, topology_path)


def _whole_check(polymer, molecule, residue_name, max_bond_nm, topology_path) -> PolymerWholeCheck:
    heavy = [atom.mass > HYDROGEN_MAX_MASS for atom in molecule.atoms]
    lengths = []
    for i, j in sorted({tuple(sorted(bond)) for bond in molecule.bonds}):
        if heavy[i] and heavy[j]:
            length = math.dist(polymer[i].position, polymer[j].position)
            lengths.append((polymer[i].atom_name, polymer[j].atom_name, length))
    if not lengths:
        raise ValueError(f"No bonded heavy-atom pairs for {residue_name} in {topology_path}")
    return PolymerWholeCheck(
        residue_name=residue_name,
        atoms=len(polymer),
        heavy_bonds=len(lengths),
        longest_bond=max(lengths, key=lambda bond: bond[2]),
        broken_bonds=tuple(bond for bond in lengths if bond[2] >= max_bond_nm),
        max_bond_nm=max_bond_nm,
    )


def export_polymer_pdb(
    gro_path: str | Path,
    topology_path: str | Path,
    pdb_path: str | Path,
    *,
    residue_name: str = "PHA",
    max_bond_nm: float = MAX_BOND_NM,
) -> PolymerPdbExport:
    """Write the polymer only (``residue_name``) as a PDB, after checking it is whole.

    Raises ``ValueError`` without writing anything when a bonded heavy-atom pair is
    ``max_bond_nm`` or longer. Atoms and residues are renumbered from 1, in chain A;
    coordinates are converted from nm to Å; no box (CRYST1) record is written.
    """

    polymer, molecule = _polymer_atoms(gro_path, topology_path, residue_name)
    check = _whole_check(polymer, molecule, residue_name, max_bond_nm, topology_path)
    if not check.is_whole:
        first, second, length = max(check.broken_bonds, key=lambda bond: bond[2])
        raise ValueError(
            f"The polymer ({residue_name}) is not whole in {Path(gro_path).name}: "
            f"{len(check.broken_bonds)} bonded heavy-atom pair(s) are {max_bond_nm} nm or longer "
            f"(longest {first}-{second}, {length:.2f} nm). {WHOLE_MOLECULE_ADVICE}"
        )

    residue_numbers: dict[int, int] = {}
    lines = []
    for serial, (atom, itp_atom) in enumerate(zip(polymer, molecule.atoms), start=1):
        residue = residue_numbers.setdefault(atom.residue_number, len(residue_numbers) + 1)
        name = atom.atom_name if len(atom.atom_name) == 4 else f" {atom.atom_name:<3s}"
        x, y, z = (10.0 * value for value in atom.position)
        lines.append(
            f"ATOM  {serial:5d} {name:4s} {residue_name[:3]:>3s} A{residue:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          "
            f"{_element(atom.atom_name, itp_atom.mass):>2s}"
        )
    pdb_path = Path(pdb_path)
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text("\n".join(lines) + "\nTER\nEND\n")
    return PolymerPdbExport(pdb_path, len(polymer), check)


def _polymer_atoms(gro_path, topology_path, residue_name):
    """The polymer's .gro atoms and its [ moleculetype ], checked against each other."""

    atoms, _ = read_gro(gro_path)
    polymer = [atom for atom in atoms if atom.residue_name.upper() == residue_name.upper()]
    if not polymer:
        names = sorted({atom.residue_name for atom in atoms})
        raise ValueError(
            f"No residue {residue_name} in {Path(gro_path).name} (residues: {', '.join(names)}). "
            "Set the polymer residue name: PHA for 06A systems, LIG for CHARMM-GUI systems."
        )
    molecule = _molecule(read_topology_molecules(topology_path), residue_name)
    if molecule is None:
        raise ValueError(f"No [ moleculetype ] {residue_name} in {topology_path} or its #include files")
    if len(polymer) != len(molecule.atoms):
        raise ValueError(
            f"{residue_name}: {len(polymer)} atoms in {Path(gro_path).name} vs {len(molecule.atoms)} "
            "in the topology; export works for one polymer molecule"
        )
    differing = [atom.atom_name for atom, itp in zip(polymer, molecule.atoms) if atom.atom_name != itp.name]
    if differing:
        raise ValueError(f"{residue_name}: atom names differ between the .gro and topology, e.g. {differing[0]}")
    return polymer, molecule


def _molecule(molecules: dict[str, ItpMolecule], name: str) -> ItpMolecule | None:
    return next((molecule for key, molecule in molecules.items() if key.upper() == name.upper()), None)


def _element(atom_name: str, mass: float) -> str:
    element = ELEMENT_BY_MASS.get(round(mass))
    if element:
        return element
    letters = "".join(character for character in atom_name if character.isalpha())
    return letters[:1].upper()
