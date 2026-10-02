"""06B steps on the packaged CHARMM-GUI example (examples/data/charmm_gui_ANC55_P3HB4); gmx is never run."""

from pathlib import Path

import pytest

from iphasimulator.charmmgui_import import (
    check_atom_counts,
    check_cgenff_penalties,
    check_cgenff_types,
    check_charmm_defaults,
    check_gromacs_run_files,
    check_includes,
    check_ligand_atom_order,
    check_ligand_charges_match_rtf,
    check_ligand_signed_volumes,
    check_molecule_counts,
    check_net_charge,
    prepare_charmm_gui_run_folder,
    read_topology_molecule_counts,
)


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "data" / "charmm_gui_ANC55_P3HB4"

pytestmark = pytest.mark.skipif(not EXAMPLE.is_dir(), reason="example dataset not present")


def _files(folder: Path) -> dict[str, bytes]:
    return {str(path.relative_to(folder)): path.read_bytes() for path in sorted(folder.rglob("*")) if path.is_file()}


def test_example_dataset_passes_every_06b_check(tmp_path):
    before = _files(EXAMPLE)
    run_dir = prepare_charmm_gui_run_folder(EXAMPLE, tmp_path / "ANC55_P3HB4_gromacs")
    topology, gro = run_dir / "topol.top", run_dir / "step5_input.gro"

    results = [
        check_charmm_defaults(topology),
        check_includes(run_dir),
        check_molecule_counts(topology, require_water=True),
        check_net_charge(topology),
        check_cgenff_penalties(EXAMPLE / "lig"),
        check_ligand_charges_match_rtf(topology, EXAMPLE / "lig" / "lig.rtf"),
        check_cgenff_types(topology),
        check_atom_counts(run_dir, "step5_input.gro"),
        check_ligand_atom_order(gro, topology),
        check_ligand_signed_volumes(gro, topology, EXAMPLE / "PHB4_R.sdf"),
        check_gromacs_run_files(run_dir),
    ]

    assert [result.status for result in results] == ["PASS"] * len(results), results
    assert read_topology_molecule_counts(topology) == {"PROA": 1, "LIG": 1, "SOD": 50, "CLA": 43, "TIP3": 45977}
    assert results[7].evidence == "141760 vs 141760"
    assert results[4].evidence.startswith("4.5 / 3.539")
    assert results[9].evidence.startswith("4/4 R")
    assert "#SBATCH --job-name=ANC55_P3HB4_gromacs" in (run_dir / "run_hpc_equilibration_production.slurm").read_text()
    assert _files(EXAMPLE) == before


def test_charge_check_accepts_lig_g_rtf(tmp_path):
    """lig_g.rtf has the same ATOM lines as lig.rtf without the penalty comments."""

    lig_g = tmp_path / "lig_g.rtf"
    lig_g.write_text("\n".join(line.split("!", 1)[0].rstrip() for line in (EXAMPLE / "lig" / "lig.rtf").read_text().splitlines()))

    result = check_ligand_charges_match_rtf(EXAMPLE / "gromacs" / "topol.top", lig_g)

    assert result.status == "PASS" and result.name == "LIG charges equal lig_g.rtf"
