#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Main Streamlit entry point for iPHAsimulatorV2.

Run this application from the project root with:

    streamlit run pha_gui.py
"""

import streamlit as st

from gui.config import (
    IPHASIMULATOR_PYTHON,
    PROJECT_ROOT,
    STRUCTURE_DATABASE,
)
from gui.data_loading import load_gui_data
from gui.state import initialise_session_state
from gui.styles import (
    apply_styles,
    configure_page,
    render_header,
)
from gui.tabs.build_console_tab import (
    render_build_console_tab,
)
from gui.tabs.molecular_preview_tab import (
    render_molecular_preview_tab,
)
from gui.tabs.openmm_builder_tab import (
    render_openmm_builder_tab,
)
from gui.tabs.polymer_builder_tab import (
    render_polymer_builder_tab,
)
from gui.tabs.system_viewer_tab import (
    render_system_viewer_tab,
)
from gui.tabs.md_system_builder_tab import (
    render_md_system_builder_tab,
)

from gui.tabs.analysis_tab import (
    render_analysis_tab,
)


# ==========================================================
# Application initialisation
# ==========================================================

configure_page()
apply_styles()
initialise_session_state()
render_header()


# ==========================================================
# Load shared GUI data
# ==========================================================

try:
    gui_data = load_gui_data()

except Exception as error:
    st.error(
        "Could not load the iPHAsimulator databases."
    )

    st.code(
        str(error)
    )

    st.stop()


# ==========================================================
# Sidebar
# ==========================================================

with st.sidebar:
    st.markdown(
        "## ⚙️ Control Centre"
    )

    st.markdown(
        '<div class="card">',
        unsafe_allow_html=True,
    )

    st.write(
        "Project root"
    )

    st.code(
        str(PROJECT_ROOT)
    )

    st.write(
        "Structure database"
    )

    st.code(
        str(STRUCTURE_DATABASE)
    )

    st.write(
        "iPHAsimulator Python"
    )

    st.code(
        str(IPHASIMULATOR_PYTHON)
    )

    st.markdown(
        "</div>",
        unsafe_allow_html=True,
    )

    # ======================================================
    # Database summary
    # ======================================================

    st.markdown(
        "## 📊 Database"
    )

    st.metric(
        "Available monomers",
        len(
            gui_data.available_phas
        ),
    )

    st.metric(
        "Registered MD systems",
        len(
            gui_data.md_systems_df
        ),
    )

    st.metric(
        "Current polymer length",
        len(
            st.session_state.sequence
        ),
    )

    st.metric(
        "Unique selected monomers",
        len(
            set(
                st.session_state.sequence
            )
        ),
    )

    st.metric(
        "OpenMM workflow steps",
        len(
            st.session_state.openmm_steps
        ),
    )

    # ======================================================
    # Environment checks
    # ======================================================

    with st.expander(
        "Environment status"
    ):
        if PROJECT_ROOT.exists():
            st.success(
                "Project directory found."
            )

        else:
            st.error(
                "Project directory not found."
            )

        if STRUCTURE_DATABASE.exists():
            st.success(
                "Structure database found."
            )

        else:
            st.error(
                "Structure database not found."
            )

        if IPHASIMULATOR_PYTHON.exists():
            st.success(
                "iPHAsimulator Python found."
            )

        else:
            st.error(
                "iPHAsimulator Python not found."
            )

    # ======================================================
    # Database tables
    # ======================================================

    with st.expander(
        "Show registered monomers"
    ):
        if gui_data.mainchain_df.empty:
            st.info(
                "No registered monomers were found."
            )

        else:
            st.dataframe(
                gui_data.mainchain_df,
                use_container_width=True,
                hide_index=True,
            )

    with st.expander(
        "Show MD-system registry"
    ):
        if gui_data.md_systems_df.empty:
            st.info(
                "No MD systems have been registered."
            )

        else:
            st.dataframe(
                gui_data.md_systems_df,
                use_container_width=True,
                hide_index=True,
            )

    # ======================================================
    # Manual refresh
    # ======================================================

    st.divider()

    if st.button(
        "🔄 Refresh application data",
        use_container_width=True,
        key="refresh_application_data",
    ):
        st.cache_data.clear()
        st.rerun()


# ==========================================================
# Main application tabs
# ==========================================================

(
    polymer_builder_tab,
    molecular_preview_tab,
    build_console_tab,
    md_system_builder_tab,
    system_viewer_tab,
    openmm_builder_tab,
    analysis_tab,
) = st.tabs(
    [
        "🧱 Polymer Builder",
        "🔬 Molecular Preview",
        "🖥 Build Console",
        "🤖 MD System Builder",
        "🧬 MD System Viewer",
        "⚛️ OpenMM Script Builder",
        "Analysis"
    ]
)


# ==========================================================
# Polymer Builder
# ==========================================================

with polymer_builder_tab:
    render_polymer_builder_tab(
        gui_data
    )


# ==========================================================
# Molecular Preview
# ==========================================================

with molecular_preview_tab:
    render_molecular_preview_tab(
        gui_data
    )


# ==========================================================
# Polymer Build Console
# ==========================================================

with build_console_tab:
    render_build_console_tab(
        gui_data
    )
    
# ==========================================================
# Build MD system 
# ==========================================================

with md_system_builder_tab:
    render_md_system_builder_tab()

# ==========================================================
# MD System Viewer
# ==========================================================

with system_viewer_tab:
    render_system_viewer_tab(
        gui_data
    )


# ==========================================================
# OpenMM Workflow Builder
# ==========================================================

with openmm_builder_tab:
    render_openmm_builder_tab(
        gui_data
    )
    
# ==========================================================
# Analysis
# ==========================================================

with analysis_tab:
    render_analysis_tab(
        gui_data
    )
