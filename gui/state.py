#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Streamlit session-state management for the iPHAsimulatorV2 GUI.
"""

from copy import deepcopy

import streamlit as st


SESSION_DEFAULTS = {
    "sequence": [],
    "preview_PHA": None,
    "sequence_history": [],
    "sequence_revision": 0,

    # Current OpenMM workflow
    "openmm_steps": [],
    "openmm_workflow_name": "Test",
    "openmm_run_name": "Test",
    "openmm_loaded_workflow_path": None,

    # Generated OpenMM script
    "generated_openmm_script": None,
    "generated_openmm_script_path": None,
    "generated_openmm_system_name": None,
    "generated_openmm_system_type": None,
    "generated_openmm_workflow_name": None,
}


def initialise_session_state():
    """
    Initialise any missing Streamlit session-state values.
    """

    for key, default_value in SESSION_DEFAULTS.items():
        if key not in st.session_state:
            st.session_state[key] = deepcopy(
                default_value
            )

def set_polymer_sequence(sequence):
    """Apply one user edit and preserve the previous state for session undo."""
    updated = list(sequence)
    previous = list(st.session_state.sequence)
    if updated == previous:
        return

    history = list(st.session_state.get("sequence_history", []))
    history.append((previous, st.session_state.preview_PHA))
    st.session_state.sequence_history = history[-50:]
    st.session_state.sequence = updated
    st.session_state.preview_PHA = updated[-1] if updated else None
    st.session_state.sequence_revision = st.session_state.get("sequence_revision", 0) + 1


def undo_polymer_sequence():
    """Restore the sequence and preview from before the most recent edit."""
    history = list(st.session_state.get("sequence_history", []))
    if not history:
        return
    sequence, preview = history.pop()
    st.session_state.sequence_history = history
    st.session_state.sequence = list(sequence)
    st.session_state.preview_PHA = preview
    st.session_state.sequence_revision = st.session_state.get("sequence_revision", 0) + 1


def clear_polymer_sequence():
    """
    Clear the selected polymer sequence and preview.
    """

    set_polymer_sequence([])


def clear_generated_openmm_script():
    """
    Clear the generated OpenMM script information.
    """

    st.session_state.generated_openmm_script = None
    st.session_state.generated_openmm_script_path = None
    st.session_state.generated_openmm_system_name = None
    st.session_state.generated_openmm_system_type = None
    st.session_state.generated_openmm_workflow_name = None


def clear_openmm_workflow():
    """
    Clear the current OpenMM workflow and generated-script information.
    """

    st.session_state.openmm_steps = []

    st.session_state.openmm_workflow_name = "Test"
    st.session_state.openmm_run_name = "Test"

    st.session_state.openmm_loaded_workflow_path = None

    clear_generated_openmm_script()
