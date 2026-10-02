"""The worked-example CLI must validate inputs before starting GROMACS."""

import importlib.util
from pathlib import Path
import sys

import pytest

from iphasimulator.charmmgui_import import CheckResult


SCRIPT = Path(__file__).resolve().parents[1] / "examples/charmm_gui_enzyme_pha/prepare_enzyme_pha_run.py"


@pytest.mark.parametrize("failed_check", ["structure", "penalties", "charges", "stereochemistry"])
def test_failed_checks_prevent_minimisation(tmp_path, monkeypatch, failed_check):
    spec = importlib.util.spec_from_file_location("prepare_enzyme_pha_run", SCRIPT)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    lig = tmp_path / "download" / "lig"
    lig.mkdir(parents=True)
    topology = lig / "lig_g.rtf"
    topology.write_text("fixture")
    monkeypatch.setattr(example, "read_cgenff_penalties", lambda *_: None)
    monkeypatch.setattr(example, "format_cgenff_penalty_report", lambda *_: "fixture penalties")
    monkeypatch.setattr(example, "prepare_charmm_gui_run_folder", lambda *_args, **_kwargs: tmp_path / "run")
    monkeypatch.setattr(example, "describe_run_folder_protocol", lambda *_: [])

    def result(name):
        return CheckResult(name, "FAIL" if name == failed_check else "PASS", "fixture")

    monkeypatch.setattr(example, "validate_charmm_gui_run_folder", lambda *_args, **_kwargs: [result("structure")])
    monkeypatch.setattr(example, "check_cgenff_penalties", lambda *_: result("penalties"))
    monkeypatch.setattr(example, "check_ligand_signed_volumes", lambda *_: result("stereochemistry"))

    def check_charges(_topology, rtf):
        assert rtf == topology  # Accept downloads with only lig_g.rtf.
        return result("charges")

    monkeypatch.setattr(example, "check_ligand_charges_match_rtf", check_charges)

    def unexpected_run(*_args, **_kwargs):
        pytest.fail("Minimisation must not run after a failed check")

    monkeypatch.setattr(example.subprocess, "run", unexpected_run)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--download", str(lig.parent),
                                    "--out", str(tmp_path / "run"), "--minimise"])
    assert example.main() == 1
