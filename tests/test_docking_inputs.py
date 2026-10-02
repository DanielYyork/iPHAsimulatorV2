"""Polymer-only PDB export for docking (notebook 11); no GROMACS is run."""

from pathlib import Path

import pytest

from iphasimulator.docking_inputs import check_polymer_whole, export_polymer_pdb


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_P3HB_4 = ROOT / "examples/output/benchmark/P3HB_4/gromacs/solvated_polymer"
EXAMPLE_CHARMM_GUI = ROOT / "examples/data/charmm_gui_ANC55_P3HB4/gromacs"

TOPOLOGY = """[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
PHA 3

[ atoms ]
1 os 1 PHA O1 1 -0.50 16.00
2 c3 1 PHA C1 2 0.20 12.01
3 c3 1 PHA C2 3 0.20 12.01
4 hc 1 PHA H1 4 0.10 1.008

[ bonds ]
1 2 1
2 3 1
3 4 1

[ system ]
test

[ molecules ]
PHA 1
SOL 1
"""


def write_gro(path, polymer_x, box=3.0):
    """Polymer O1, C1, C2, H1 at the given x (nm), then one water."""

    rows = [(1, "PHA", name, x, 1.0, 1.0) for name, x in zip(("O1", "C1", "C2", "H1"), polymer_x)]
    rows += [(2, "SOL", "OW", 2.0, 2.0, 2.0), (2, "SOL", "HW1", 2.1, 2.0, 2.0), (2, "SOL", "HW2", 1.97, 2.09, 2.0)]
    lines = ["test", f"{len(rows):5d}"]
    lines += [f"{res:5d}{resname:<5s}{name:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}" for i, (res, resname, name, x, y, z) in enumerate(rows, 1)]
    lines.append(f"{box:10.5f}{box:10.5f}{box:10.5f}")
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def system(tmp_path):
    (tmp_path / "topol.top").write_text(TOPOLOGY)
    return tmp_path


def test_export_writes_the_polymer_only(system):
    gro = write_gro(system / "step5_input.gro", (1.00, 1.14, 1.29, 1.39))

    result = export_polymer_pdb(gro, system / "topol.top", system / "out" / "pha.pdb")

    atoms = [line for line in result.pdb_path.read_text().splitlines() if line.startswith("ATOM")]
    assert result.atoms == 4 and len(atoms) == 4
    assert {line[17:20] for line in atoms} == {"PHA"}
    assert [line[76:78].strip() for line in atoms] == ["O", "C", "C", "H"]
    assert atoms[1][30:38].strip() == "11.400"  # nm -> Å
    assert result.check.heavy_bonds == 2 and result.check.is_whole


def test_split_polymer_stops_with_advice_and_writes_nothing(system):
    gro = write_gro(system / "frame.gro", (2.80, 2.94, 0.09, 0.19))  # C1-C2 across the box edge
    pdb = system / "out" / "pha.pdb"

    with pytest.raises(ValueError, match="not whole.*C1-C2.*gmx trjconv -pbc mol"):
        export_polymer_pdb(gro, system / "topol.top", pdb)
    assert not pdb.exists()


def test_only_heavy_atom_bonds_count(system):
    gro = write_gro(system / "frame.gro", (1.00, 1.14, 1.29, 2.50))  # C2-H1 stretched, heavy atoms fine

    assert check_polymer_whole(gro, system / "topol.top").is_whole


def test_wrong_residue_name_is_explained(system):
    gro = write_gro(system / "frame.gro", (1.00, 1.14, 1.29, 1.39))

    with pytest.raises(ValueError, match="No residue LIG.*PHA for 06A systems, LIG for CHARMM-GUI"):
        export_polymer_pdb(gro, system / "topol.top", system / "lig.pdb", residue_name="LIG")


@pytest.mark.skipif(not (BENCHMARK_P3HB_4 / "step5_input.gro").is_file(), reason="benchmark output not present")
def test_benchmark_p3hb_4_exports_51_pha_atoms(tmp_path):
    gro = BENCHMARK_P3HB_4 / "step5_input.gro"
    before = gro.read_bytes()

    result = export_polymer_pdb(gro, BENCHMARK_P3HB_4 / "topol.top", tmp_path / "P3HB_4.pdb", residue_name="PHA")

    assert result.atoms == 51 and result.check.is_whole
    assert sum(line.startswith("ATOM") for line in result.pdb_path.read_text().splitlines()) == 51
    assert gro.read_bytes() == before


@pytest.mark.skipif(not EXAMPLE_CHARMM_GUI.is_dir(), reason="example dataset not present")
def test_charmm_gui_example_exports_51_lig_atoms(tmp_path):
    result = export_polymer_pdb(
        EXAMPLE_CHARMM_GUI / "step3_input.gro", EXAMPLE_CHARMM_GUI / "topol.top", tmp_path / "lig.pdb", residue_name="LIG"
    )

    assert result.atoms == 51 and result.check.is_whole
