"""Docked-complex inputs for the Amber/OpenMM enzyme + polymer system; RDKit only."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from iphasimulator.amber_complex import (
    pose_polymer_mol2,
    read_charmm_histidine_states,
    write_protein_pdb,
)


R_LIGAND = "C[C@@H](O)CC(=O)O"


def _template(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles(R_LIGAND))
    AllChem.EmbedMolecule(mol, randomSeed=3)
    sdf = tmp_path / "PHA.antechamber.sdf"
    Chem.MolToMolFile(mol, str(sdf))
    lines = ["@<TRIPOS>MOLECULE", "PHA", f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0", "SMALL", "bcc", "", "@<TRIPOS>ATOM"]
    for atom in mol.GetAtoms():
        p = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        name = f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
        lines.append(f"{atom.GetIdx() + 1:>7d} {name:<8s}{p.x:>10.4f}{p.y:>10.4f}{p.z:>10.4f} {atom.GetSymbol().lower()}3 1 PHA 0.100000")
    lines.append("@<TRIPOS>BOND")
    lines += [f"{i + 1} {b.GetBeginAtomIdx() + 1} {b.GetEndAtomIdx() + 1} 1" for i, b in enumerate(mol.GetBonds())]
    mol2 = tmp_path / "PHA.gaff2.mol2"
    mol2.write_text("\n".join(lines) + "\n")
    return mol, sdf, mol2


def _complex_pdb(tmp_path, ligand):
    """Protein (2 residues, one HIS, with an H) + heavy-atom-only ligand in a new pose."""

    pose = Chem.Mol(ligand)
    conf = pose.GetConformer()
    for i in range(pose.GetNumAtoms()):  # rigid shift = a different position
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 10.0, p.y - 4.0, p.z + 2.5))
    rows = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  H   ALA A   1       0.500   0.500   0.000  1.00  0.00           H",
        "ATOM      3  CA  ALA A   1       1.400   0.000   0.000  1.00  0.00           C",
        "ATOM      4  N   HIS A   2       3.000   0.000   0.000  1.00  0.00           N",
        "ATOM      5  CA  HIS A   2       4.400   0.000   0.000  1.00  0.00           C",
    ]
    heavy = [a for a in pose.GetAtoms() if a.GetAtomicNum() > 1]
    for k, atom in enumerate(reversed(heavy)):  # different atom order than the template
        p = conf.GetAtomPosition(atom.GetIdx())
        rows.append(f"HETATM{6 + k:>5d} {atom.GetSymbol() + str(k + 1):<4} LIG X   1    {p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {atom.GetSymbol():>2}")
    path = tmp_path / "complex.pdb"
    path.write_text("\n".join(rows) + "\nEND\n")
    return path, pose


def test_pose_polymer_mol2_moves_heavy_atoms_onto_pose_and_keeps_r(tmp_path):
    ligand, sdf, mol2 = _template(tmp_path)
    complex_pdb, pose = _complex_pdb(tmp_path, ligand)

    result = pose_polymer_mol2(mol2, sdf, complex_pdb, tmp_path / "posed.mol2")

    assert (result.heavy_atoms, result.hydrogens_placed) == (7, 8)  # C4H8O3
    assert [label for _, label in result.stereocentres] == ["R"]
    rows = [l.split() for l in result.output_mol2.read_text().split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().splitlines()]
    for atom in ligand.GetAtoms():
        row = rows[atom.GetIdx()]
        assert row[1] == f"{atom.GetSymbol()}{atom.GetIdx() + 1}" and row[5:] == [f"{atom.GetSymbol().lower()}3", "1", "PHA", "0.100000"]
        if atom.GetAtomicNum() > 1:
            p = pose.GetConformer().GetAtomPosition(atom.GetIdx())
            assert [float(v) for v in row[2:5]] == pytest.approx([p.x, p.y, p.z], abs=1e-3)
    assert "@<TRIPOS>BOND" in result.output_mol2.read_text()


def test_pose_polymer_mol2_rejects_a_different_oligomer(tmp_path):
    ligand, sdf, mol2 = _template(tmp_path)
    complex_pdb, _ = _complex_pdb(tmp_path, ligand)
    complex_pdb.write_text("\n".join(l for l in complex_pdb.read_text().splitlines() if "HETATM    6" not in l))

    with pytest.raises(ValueError, match="same oligomer"):
        pose_polymer_mol2(mol2, sdf, complex_pdb, tmp_path / "posed.mol2")


def test_write_protein_pdb_drops_ligand_and_hydrogens_and_sets_histidines(tmp_path):
    ligand, _, _ = _template(tmp_path)
    complex_pdb, _ = _complex_pdb(tmp_path, ligand)

    summary = write_protein_pdb(complex_pdb, tmp_path / "protein.pdb", histidine_states={2: "HID"})

    text = summary.output_path.read_text()
    assert "LIG" not in text and " H   ALA" not in text
    assert "HID A   2" in text
    assert (summary.residues, summary.removed_hydrogens, summary.histidines) == (2, 1, {2: "HID"})
    assert write_protein_pdb(complex_pdb, tmp_path / "p2.pdb").histidines == {2: "HIS"}


def test_read_charmm_histidine_states(tmp_path):
    itp = tmp_path / "PROA.itp"
    itp.write_text("[ moleculetype ]\nPROA 3\n[ atoms ]\n"
                   "1 NH1 58 HSD N 1 -0.47 14.007\n2 NH1 78 HSE N 2 -0.47 14.007\n3 NH1 80 HSP N 3 -0.47 14.007\n4 NH1 81 LEU N 4 -0.47 14.007\n")

    assert read_charmm_histidine_states(itp) == {58: "HID", 78: "HIE", 80: "HIP"}
