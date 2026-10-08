"""Regression checks for the Render dependency-check failures."""

import contextlib
import importlib.util
import io
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location(
    "deployment_dependencies", Path(__file__).resolve().parents[1] / "deploy/check_dependencies.py",
)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class DeploymentDependencyChecks(unittest.TestCase):
    def run_check(self, stdout, returncode=1, stderr="", provider_error=None):
        result = subprocess.CompletedProcess([], returncode, stdout, stderr)
        with (patch.object(checker.subprocess, "run", return_value=result),
              patch.object(checker, "verify_conda_openbabel", side_effect=provider_error) as provider,
              contextlib.redirect_stdout(io.StringIO()),
              contextlib.redirect_stderr(io.StringIO())):
            status = checker.main()
        return status, provider.call_count

    def test_clean_environment_passes_without_exception(self):
        self.assertEqual(self.run_check("No broken requirements found.\n", 0), (0, 0))

    def test_reported_missing_packages_and_version_conflict_still_fail(self):
        # The original Render failure must never be accepted wholesale.
        output = "\n".join([
            "packmol-memgen 2026.3.25 requires pdb2pqr, which is not installed.",
            "proprep 1.0.0 requires pdb2pqr, which is not installed.",
            checker.CONDA_OPENBABEL_MISMATCH,
            "proprep 1.0.0 has requirement biopython<1.86,>=1.83, but you have biopython 1.88.",
        ])
        self.assertEqual(self.run_check(output), (1, 0))

    def test_known_mismatch_requires_working_conda_provider(self):
        self.assertEqual(self.run_check(checker.CONDA_OPENBABEL_MISMATCH + "\n"), (0, 1))

    def test_missing_or_broken_conda_provider_fails(self):
        status, calls = self.run_check(
            checker.CONDA_OPENBABEL_MISMATCH, provider_error=RuntimeError("missing library"),
        )
        self.assertEqual((status, calls), (1, 1))

    def test_other_acpype_versions_are_not_silently_exempted(self):
        output = checker.CONDA_OPENBABEL_MISMATCH.replace("2026.9.4", "2027.1.0")
        self.assertEqual(self.run_check(output), (1, 0))

    def test_pip_execution_failure_still_fails(self):
        self.assertEqual(self.run_check("", returncode=2, stderr="pip failed"), (2, 0))

    def test_additional_warning_does_not_hide_an_error(self):
        self.assertEqual(self.run_check(
            checker.CONDA_OPENBABEL_MISMATCH, stderr="unexpected failure",
        ), (1, 0))


if __name__ == "__main__":
    unittest.main()
