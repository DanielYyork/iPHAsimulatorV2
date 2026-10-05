"""Check scientific integration and window-limited reuse without private MD data."""
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


NOTEBOOK = (Path(__file__).resolve().parents[1] / 'src/md_simulation_scripts/'
            'enzyme_pha_analysis/enzyme_p3ho4_binding_modes.ipynb')


@pytest.fixture(scope='module')
def analysis():
    nb = json.loads(NOTEBOOK.read_text())
    scope = {'np': np, 'pd': pd, 'Path': Path, 'hashlib': hashlib,
             'CONTACT_CUTOFF_A': 3.5, 'RANDOM_SEED': 42,
             'BOOTSTRAP_REPEATS': 2000, 'DISTANCE_TOLERANCE_A': 2e-5}
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            source = ''.join(cell['source'])
            tree = ast.parse(source)
            if 'def contact_duration(' in source:
                functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
                exec(compile(ast.Module(body=functions, type_ignores=[]), str(NOTEBOOK), 'exec'), scope)
    return scope


def test_notebook_code_syntax_and_numbering():
    nb = json.loads(NOTEBOOK.read_text())
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            ast.parse(''.join(cell['source']))
    assert nb['metadata']['binding_modes']['figures'] == [7, 8]


def test_time_weighting_and_multiple_residues(analysis):
    t = np.array([0., 3., 10.])
    d = np.array([[2., 5.], [2., 5.], [5., 2.]])
    durations = analysis['contact_duration'](t, d, 3.5)
    np.testing.assert_allclose(durations, [6.5, 3.5])
    # Per-frame fractions would give 2/3 and 1/3 rather than these time fractions.
    assert not np.isclose(durations[0] / 10, 2/3)


def test_clip_boundaries_and_additivity(analysis):
    t = np.array([0., 3., 10.]); d = np.array([2., 2., 5.])
    whole = analysis['contact_duration'](t, d, 3.5)
    pieces = []
    for start, end in [(0, 4), (4, 10)]:
        ct, cd = analysis['clip_series'](t, d, start, end)
        pieces.append(analysis['contact_duration'](ct, cd, 3.5))
    assert np.isclose(sum(pieces), whole)
    with pytest.raises(ValueError, match='outside'):
        analysis['clip_series'](t, d, -1, 10)


def test_equal_cutoff_and_invalid_samples(analysis):
    assert analysis['contact_duration']([0, 10], np.array([3.5, 3.5]), 3.5) == 10
    for times, values in [([0, 0], [2, 3]), ([0, 1], [2, np.nan]), ([0, 1], [-1, 2])]:
        with pytest.raises(ValueError):
            analysis['contact_duration'](times, np.array(values), 3.5)


def test_fallback_reads_only_selected_windows_and_brackets(analysis):
    g = pd.DataFrame({'frame': range(11), 'time_ns': np.arange(11.)})
    chosen = analysis['required_frames'](g, {'A': (1.5, 3.5), 'B': (7, 9)})
    assert chosen.frame.tolist() == [1, 2, 3, 4, 7, 8, 9]
    assert not {0, 5, 6, 10}.intersection(chosen.frame)


def test_blocks_integrate_whole_intervals(analysis):
    t = np.arange(0., 31.); v = 2*t+3
    centers, values = analysis['block_values'](t, v, 10)
    np.testing.assert_allclose(centers, [5, 15, 25])
    np.testing.assert_allclose(values, [13, 33, 53])
    assert np.isclose(values.mean(), analysis['time_mean'](t, v))
    with pytest.raises(ValueError, match='divide'):
        analysis['block_values'](t, v, 12)
    np.testing.assert_allclose(analysis['block_interval'](np.ones((10, 2))*7), [[7, 7], [7, 7]])


def test_residue_schema_and_minima_are_required(analysis):
    r = pd.DataFrame({'distance_column': ['resindex_0_SER1_A', 'resindex_1_GLY2_A'],
                      'resindex': [0, 1]})
    d = pd.DataFrame({'frame': [0, 1], 'time_ns': [0., 1.],
                      'resindex_0_SER1_A': [2., 8.], 'resindex_1_GLY2_A': [5., 3.]})
    g = pd.DataFrame({'frame': [0, 1], 'time_ns': [0., 1.],
                      'enzyme_pha_min_distance_A': [2., 3.]})
    analysis['check_distances'](d, r, g)
    with pytest.raises(ValueError, match='schema'):
        analysis['check_distances'](g, r, g)
    wrong = d.copy(); wrong.loc[1, 'resindex_1_GLY2_A'] = 4
    with pytest.raises(ValueError, match='reproduce'):
        analysis['check_distances'](wrong, r, g)
