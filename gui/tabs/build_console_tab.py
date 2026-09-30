#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Polymer Build Console tab for the iPHAsimulatorV2 Streamlit GUI.

This tab launches the backend PHAPolymerBuilder in the configured
iPHAsimulator environment and displays the captured subprocess output.
"""

import streamlit as st

from gui.config import (
    IPHASIMULATOR_PYTHON,
    PROJECT_ROOT,
)
from gui.models import GUIData
from gui.styles import (
    render_info_box,
    render_success_box,
)
from gui.subprocess_helpers import (
    build_polymer_subprocess,
)


def render_build_console_tab(
    gui_data: GUIData,
) -> None:
    """
    Render the iPHAsimulator polymer-build console.

    Parameters
    ----------
    gui_data : GUIData
        Shared GUI data. It is currently retained to keep the same public
        function signature as the other tab-rendering functions.
    """

    del gui_data

    st.markdown(
        "## 🖥 Polymer Build Console"
    )

    render_info_box(
        "This tab launches PHAPolymerBuilder in the configured "
        "iphasimulator environment. The Streamlit process remains "
        "responsible only for the graphical interface."
    )

    st.divider()

    # ======================================================
    # Current build target
    # ======================================================

    sequence = list(
        st.session_state.sequence
    )

    if sequence:
        unique_units = []

        for unit in sequence:
            if unit not in unique_units:
                unique_units.append(unit)

        if len(unique_units) == 1:
            build_target = (
                f"P{unique_units[0]}_"
                f"{len(sequence)}"
            )

            build_mode = "Homopolymer"

        else:
            build_target = (
                "Custom copolymer"
            )

            build_mode = "Copolymer"

        st.markdown(
            "### Current Build Target"
        )

        summary_columns = st.columns(3)

        summary_columns[0].metric(
            "Chain length",
            len(sequence),
        )

        summary_columns[1].metric(
            "Unique monomers",
            len(unique_units),
        )

        summary_columns[2].metric(
            "Build mode",
            build_mode,
        )

        st.markdown(
            '<div class="card">',
            unsafe_allow_html=True,
        )

        st.write(
            f"**Target:** `{build_target}`"
        )

        st.write(
            "**Sequence:** "
            f"`{' - '.join(sequence)}`"
        )

        st.markdown(
            "</div>",
            unsafe_allow_html=True,
        )

    else:
        st.warning(
            "No polymer sequence has been selected. "
            "Create one in the Polymer Builder tab first."
        )

    st.divider()

    # ======================================================
    # Build action
    # ======================================================

    build_clicked = st.button(
        "🚀 Build polymer now",
        use_container_width=True,
        key="build_console_build_polymer",
        disabled=not sequence,
    )

    if build_clicked:
        st.info(
            "Launching polymer build in iPHAsimulator..."
        )

        progress = st.progress(0)
        status = st.empty()
        result = None

        with st.spinner(
            "Building polymer with PHAPolymerBuilder..."
        ):
            try:
                status.write(
                    "Validating build request..."
                )

                progress.progress(10)

                if not sequence:
                    raise ValueError(
                        "Cannot build an empty polymer sequence."
                    )

                status.write(
                    "Preparing iPHAsimulator subprocess..."
                )

                progress.progress(25)

                result = build_polymer_subprocess(
                    sequence
                )

                progress.progress(85)

                status.write(
                    "Build subprocess finished."
                )

            except Exception as error:
                progress.progress(100)

                st.error(
                    "Could not launch the polymer-build subprocess."
                )

                st.code(
                    str(error)
                )

        if result is not None:
            progress.progress(100)

            if result.returncode == 0:
                st.success(
                    "Polymer build completed successfully."
                )

                st.balloons()

                render_success_box(
                    "PHAPolymerBuilder returned successfully."
                )

            else:
                st.error(
                    "Polymer build failed."
                )

            # ==================================================
            # Subprocess summary
            # ==================================================

            result_columns = st.columns(3)

            result_columns[0].metric(
                "Return code",
                result.returncode,
            )

            result_columns[1].metric(
                "STDOUT characters",
                len(
                    result.stdout or ""
                ),
            )

            result_columns[2].metric(
                "STDERR characters",
                len(
                    result.stderr or ""
                ),
            )

            # ==================================================
            # Captured output
            # ==================================================

            if result.stdout:
                st.subheader(
                    "STDOUT"
                )

                st.code(
                    result.stdout
                )

            else:
                st.caption(
                    "The subprocess did not produce STDOUT."
                )

            if result.stderr:
                st.subheader(
                    (
                        "STDERR / warnings"
                        if result.returncode == 0
                        else "STDERR"
                    )
                )

                st.code(
                    result.stderr
                )

            else:
                st.caption(
                    "The subprocess did not produce STDERR."
                )

    st.divider()

    # ======================================================
    # Command context
    # ======================================================

    st.markdown(
        "### Current Command Context"
    )

    context_columns = st.columns(2)

    with context_columns[0]:
        st.write(
            "iPHAsimulator Python:"
        )

        st.code(
            str(IPHASIMULATOR_PYTHON)
        )

        if IPHASIMULATOR_PYTHON.exists():
            st.success(
                "iPHAsimulator Python was found."
            )

        else:
            st.error(
                "iPHAsimulator Python was not found."
            )

    with context_columns[1]:
        st.write(
            "Project working directory:"
        )

        st.code(
            str(PROJECT_ROOT)
        )

        if PROJECT_ROOT.exists():
            st.success(
                "Project directory was found."
            )

        else:
            st.error(
                "Project directory was not found."
            )

    st.caption(
        "The build runs in a separate Python process so that "
        "iPHAsimulator dependencies do not need to be installed in "
        "the Streamlit GUI environment."
    )
