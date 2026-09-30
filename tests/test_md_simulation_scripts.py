"""Standalone-template and relocated contact entry-point checks; no MD data needed."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
PREPARATION = ROOT / 'src/md_simulation_scripts/trajectory_preparation'
CONTACTS = ROOT / 'src/md_simulation_scripts/enzyme_contacts'


@pytest.fixture
def analysis(tmp_path):
    folder = tmp_path / 'simulation with spaces' / 'analysis'
    folder.mkdir(parents=True)
    for filename in ('process_trajectory.sh', 'instrcution.txt'):
        shutil.copyfile(PREPARATION / filename, folder / filename)
    (folder.parent / 'production_combined_1us.xtc').write_text('raw fixture')
    (folder / 'step7_production.tpr').write_text('topology fixture')
    (folder / 'index.ndx').write_text('[ SOLU ]\n1 2\n[ SYSTEM ]\n1 2 3\n')
    stub = folder / 'mock gmx'
    stub.write_text('''#!/usr/bin/env bash
set -eu
printf '%s\\n' "$@" >> commands.log
input="$(cat)"
printf '%s\\n' "$input" >> selections.log
[[ "${FAIL_GMX:-false}" == false ]] || exit 9
while [[ $# -gt 0 ]]; do
    if [[ "$1" == -o ]]; then printf 'processed fixture' > "$2"; exit 0; fi
    shift
done
exit 8
''')
    stub.chmod(0o755)
    change(folder, 'GMX', str(stub))
    return folder


def change(folder, key, value):
    path = folder / 'instrcution.txt'
    path.write_text('\n'.join(f'{key}={value}' if line.startswith(f'{key}=') else line
                              for line in path.read_text().splitlines()) + '\n')


def run(folder, *args, **environment):
    return subprocess.run(['bash', str(folder / 'process_trajectory.sh'), *args],
                          cwd=folder.parent, env={**os.environ, **environment},
                          text=True, capture_output=True)


def snapshot(folder):
    return {str(p.relative_to(folder)): p.read_bytes() for p in folder.rglob('*') if p.is_file()}


def test_bash_syntax():
    subprocess.run(['bash', '-n', str(PREPARATION / 'process_trajectory.sh')], check=True)


def test_dry_run_is_read_only_and_preserves_sequence(analysis):
    before = snapshot(analysis.parent)
    result = run(analysis, '--dry-run')
    assert result.returncode == 0, result.stderr
    assert snapshot(analysis.parent) == before
    first, second = result.stdout.split('Step 2 selection:')
    assert '-pbc mol -ur compact -center' in first
    assert '-skip' not in first
    assert 'SOLU, then SYSTEM' in first
    assert '-f processed.xtc' in second and '-skip 100' in second
    assert '-o processed_every100.xtc' in second
    assert 'SYSTEM' in second
    assert not (analysis / 'commands.log').exists()


def test_two_commands_and_stdin_with_stub_only(analysis):
    raw = (analysis.parent / 'production_combined_1us.xtc').read_bytes()
    result = run(analysis)
    assert result.returncode == 0, result.stderr
    commands = (analysis / 'commands.log').read_text().splitlines()
    split = commands.index('trjconv', 1)
    first, second = commands[:split], commands[split:]
    assert first == ['trjconv', '-f', '../production_combined_1us.xtc', '-s', 'step7_production.tpr',
                     '-n', 'index.ndx', '-pbc', 'mol', '-ur', 'compact', '-center', '-o', 'processed.xtc']
    assert second == ['trjconv', '-f', 'processed.xtc', '-s', 'step7_production.tpr',
                      '-n', 'index.ndx', '-skip', '100', '-o', 'processed_every100.xtc']
    assert (analysis / 'selections.log').read_text() == 'SOLU\nSYSTEM\nSYSTEM\n'
    assert (analysis.parent / 'production_combined_1us.xtc').read_bytes() == raw


@pytest.mark.parametrize('key,value,message', [
    ('STRIDE', '0', 'positive integer'), ('STRIDE', '-2', 'positive integer'),
    ('TPR', 'view.gro', 'production .tpr'), ('PREVIEW', 'processed.xtc', 'must differ'),
    ('PROCESSED', '../bad.xtc', 'inside this folder'), ('INPUT', '../missing.xtc', 'unreadable'),
    ('GMX', 'missing_gmx_executable', 'Cannot find'),
])
def test_invalid_settings(analysis, key, value, message):
    if value == 'view.gro':
        (analysis / value).write_text('gro fixture')
    change(analysis, key, value)
    result = run(analysis, '--dry-run')
    assert result.returncode != 0 and message in result.stderr
    assert not (analysis / 'commands.log').exists()


@pytest.mark.parametrize('extra,message', [('UNKNOWN=x', 'Unknown setting'),
                                           ('STRIDE=2', 'Duplicate setting')])
def test_bad_keys(analysis, extra, message):
    with (analysis / 'instrcution.txt').open('a') as handle:
        handle.write(extra + '\n')
    result = run(analysis, '--dry-run')
    assert result.returncode != 0 and message in result.stderr


@pytest.mark.parametrize('filename', ['processed.xtc', 'processed_every100.xtc'])
def test_output_protection(analysis, filename):
    (analysis / filename).write_bytes(b'keep this')
    result = run(analysis)
    assert result.returncode != 0 and 'Output already exists' in result.stderr
    assert (analysis / filename).read_bytes() == b'keep this'
    assert not (analysis / 'commands.log').exists()


def test_dangling_output_symlink_protected(analysis):
    (analysis / 'processed.xtc').symlink_to(analysis / 'nonexistent')
    assert run(analysis, '--dry-run').returncode != 0
    assert (analysis / 'processed.xtc').is_symlink()


def test_gromacs_failure_stops_before_second_command(analysis):
    result = run(analysis, FAIL_GMX='true')
    assert result.returncode != 0
    assert (analysis / 'commands.log').read_text().splitlines().count('trjconv') == 1
    assert not (analysis / 'processed_every100.xtc').exists()


def test_settings_are_never_executed(analysis):
    change(analysis, 'INPUT', '$(touch UNEXPECTED).xtc')
    assert run(analysis, '--dry-run').returncode != 0
    assert not (analysis / 'UNEXPECTED').exists()


def test_contact_configuration_preserves_output_location():
    from iphasimulator.analysis_contacts import load_config
    config = load_config(CONTACTS / 'GK13_P3HO_4.yaml')
    assert config.output_root == ROOT / 'examples/output/enzyme_contacts'
    assert config.preview_frames == 31


def test_contact_runner_help_from_arbitrary_directory(tmp_path):
    environment = dict(os.environ)
    environment.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable, str(CONTACTS / 'run_enzyme_contacts.py'), '--help'],
                            cwd=tmp_path, env=environment, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '--full' in result.stdout and '--sample-interval-ns' in result.stdout


def test_contact_notebook_sources_and_relocated_paths():
    notebook = json.loads((CONTACTS / 'enzyme_contacts.ipynb').read_text())
    assert notebook['nbformat'] == 4
    for cell in notebook['cells']:
        if cell['cell_type'] == 'code':
            compile(''.join(cell['source']), 'enzyme_contacts.ipynb', 'exec')
    source = '\n'.join(''.join(c['source']) for c in notebook['cells'])
    assert 'src/md_simulation_scripts/enzyme_contacts/GK13_P3HO_4.yaml' in source
    assert 'examples/enzyme_contacts_GK13_P3HO_4.yaml' not in source
    assert (CONTACTS / '../../../docs/enzyme_contacts.md').resolve().is_file()
    assert not (ROOT / 'src/md_simulation_scripts/01_GK13_PHO4_trajactory_process.ipynb').exists()
    assert not (ROOT / 'examples/enzyme_trajectory_GK13_P3HO_4.yaml').exists()
