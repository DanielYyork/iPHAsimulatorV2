#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis tab for the iPHAsimulatorV2 Streamlit GUI.

This tab allows the user to:

- select a registered molecular-dynamics system
- discover available analysis workflows
- select an analysis workflow
- select a completed simulation / replica
- configure analysis options
- execute the selected analysis workflow
- inspect completed analysis outputs

Analysis workflows are discovered automatically from:

    src/iphasimulator/analysis/<analysis_name>/workflow.py

Scientific analysis workflows are executed inside the dedicated
iphasimulator Python environment rather than the pha_gui environment.
"""

from __future__ import annotations

from pathlib import Path
import json
from numbers import Integral

import streamlit as st

from gui.models import GUIData
from gui.tabs.analysis_comparison import render_analysis_comparison
from gui.styles import (
    render_info_box,
    render_warning_box,
)

from gui.subprocess_helpers import (
    run_python_script_with_iphasimulator,
)


# =============================================================================
# Analysis workflow discovery
# =============================================================================

def get_analysis_root() -> Path:
    """
    Return the src/iphasimulator/analysis directory.
    """

    project_root = (
        Path(__file__)
        .resolve()
        .parents[2]
    )

    analysis_root = (
        project_root
        / "src"
        / "iphasimulator"
        / "analysis"
    )

    return analysis_root


def discover_analysis_workflows() -> dict[str, Path]:
    """
    Discover available analysis workflow.py files.

    Any workflow.py found beneath:

        src/iphasimulator/analysis/

    is treated as an available analysis workflow.

    Returns
    -------
    dict[str, pathlib.Path]
        Mapping of workflow identifier to workflow.py path.
    """

    analysis_root = (
        get_analysis_root()
    )

    workflows = {}


    if not analysis_root.is_dir():

        return workflows


    for workflow_path in sorted(
        analysis_root.rglob(
            "workflow.py"
        )
    ):

        relative_parent = (
            workflow_path
            .parent
            .relative_to(
                analysis_root
            )
        )


        workflow_name = (
            str(
                relative_parent
            )
            .replace(
                "/",
                "_",
            )
            .replace(
                "\\",
                "_",
            )
        )


        workflows[
            workflow_name
        ] = workflow_path


    return workflows


# =============================================================================
# System selection
# =============================================================================

def _render_system_selection(
    gui_data: GUIData,
):
    """
    Render MD-system selection controls.

    Parameters
    ----------
    gui_data : GUIData
        Shared GUI data.

    Returns
    -------
    tuple
        selected_system_name, selected_system_type
    """

    md_systems_df = (
        gui_data.md_systems_df
    )


    st.markdown(
        "### 1. Registered MD System"
    )


    if md_systems_df.empty:

        render_warning_box(
            "No molecular-dynamics systems are registered yet."
        )

        return (
            None,
            None,
        )


    # -------------------------------------------------------------------------
    # System-type filter
    # -------------------------------------------------------------------------

    available_system_types = sorted(
        str(system_type)
        for system_type in (
            md_systems_df[
                "system_type"
            ]
            .dropna()
            .unique()
        )
        if str(system_type).strip()
    )


    type_filter = st.selectbox(
        "System type filter",
        [
            "All",
            *available_system_types,
        ],
        key="analysis_system_type_filter",
    )


    if type_filter == "All":

        filtered_systems_df = (
            md_systems_df.copy()
        )

    else:

        filtered_systems_df = (
            md_systems_df[
                md_systems_df[
                    "system_type"
                ]
                == type_filter
            ]
            .copy()
        )


    if filtered_systems_df.empty:

        render_warning_box(
            "No systems match the selected system type."
        )

        return (
            None,
            None,
        )


    # -------------------------------------------------------------------------
    # System selection
    # -------------------------------------------------------------------------

    selected_system_name = st.selectbox(
        "Prepared MD system",
        filtered_systems_df[
            "system_name"
        ].tolist(),
        key="analysis_selected_system",
    )


    selected_system_row = (
        filtered_systems_df[
            filtered_systems_df[
                "system_name"
            ]
            == selected_system_name
        ]
        .iloc[0]
    )


    selected_system_type = str(
        selected_system_row[
            "system_type"
        ]
    )


    # -------------------------------------------------------------------------
    # Display selected system
    # -------------------------------------------------------------------------

    st.write(
        f"**System:** `{selected_system_name}`"
    )

    st.write(
        f"**System type:** `{selected_system_type}`"
    )


    try:

        system_files = (
            gui_data.paths.get_md_system_files(
                system_name=(
                    selected_system_name
                ),
                system_type=(
                    selected_system_type
                ),
            )
        )


        system_directory = (
            system_files[
                "system_dir"
            ]
        )


        with st.expander(
            "Show system path"
        ):

            st.code(
                str(
                    system_directory
                )
            )


    except Exception as error:

        render_warning_box(
            "The selected system is registered, but its "
            "system directory could not be resolved."
        )

        st.code(
            str(
                error
            )
        )


    return (
        selected_system_name,
        selected_system_type,
    )


# =============================================================================
# Analysis workflow selection
# =============================================================================

def _render_analysis_workflow_selection():
    """
    Render controls for selecting an available analysis workflow.

    Returns
    -------
    tuple
        selected_workflow_name, selected_workflow_path
    """

    st.markdown(
        "### 2. Analysis Workflow"
    )


    workflows = (
        discover_analysis_workflows()
    )


    if not workflows:

        render_warning_box(
            "No analysis workflows were found beneath "
            "src/iphasimulator/analysis/."
        )

        return (
            None,
            None,
        )


    workflow_names = sorted(
        workflows
    )


    selected_workflow_name = (
        st.selectbox(
            "Available analysis workflow",
            workflow_names,
            key="analysis_selected_workflow",
            format_func=lambda name: (
                name
                .replace(
                    "_",
                    " ",
                )
                .title()
            ),
        )
    )


    selected_workflow_path = (
        workflows[
            selected_workflow_name
        ]
    )


    st.write(
        "**Selected workflow:** "
        f"`{selected_workflow_name.replace('_', ' ').title()}`"
    )


    with st.expander(
        "Show workflow file"
    ):

        st.code(
            str(
                selected_workflow_path
            )
        )


    return (
        selected_workflow_name,
        selected_workflow_path,
    )


# =============================================================================
# Simulation selection
# =============================================================================

def _render_simulation_selection(
    gui_data: GUIData,
    selected_system_name,
    selected_system_type,
):
    """
    Render available simulation replicas for the selected MD system.

    Returns
    -------
    str or None
        Selected simulation directory name.
    """

    st.markdown(
        "### 3. Simulation / Replica"
    )


    if (
        selected_system_name is None
        or selected_system_type is None
    ):

        render_info_box(
            "Select a molecular-dynamics system first."
        )

        return None


    try:

        system_files = (
            gui_data.paths.get_md_system_files(
                system_name=(
                    selected_system_name
                ),
                system_type=(
                    selected_system_type
                ),
            )
        )


        simulations_directory = (
            Path(
                system_files[
                    "simulations_dir"
                ]
            )
        )


    except Exception as error:

        render_warning_box(
            "Could not resolve the simulations directory "
            "for the selected system."
        )

        st.code(
            str(
                error
            )
        )

        return None


    if not simulations_directory.is_dir():

        render_warning_box(
            "The selected system does not currently have "
            "a simulations directory."
        )

        return None


    simulation_directories = sorted(
        path
        for path in (
            simulations_directory.iterdir()
        )
        if path.is_dir()
    )


    if not simulation_directories:

        render_warning_box(
            "No simulation replicas were found for "
            "the selected system."
        )

        return None


    selected_simulation = (
        st.selectbox(
            "Simulation replica",
            simulation_directories,
            format_func=lambda path: (
                path.name
            ),
            key="analysis_selected_simulation",
        )
    )


    with st.expander(
        "Show simulation directory"
    ):

        st.code(
            str(
                selected_simulation
            )
        )


    return (
        selected_simulation.name
    )


# =============================================================================
# Workflow module resolution
# =============================================================================

def _workflow_path_to_module_name(
    workflow_path: Path,
) -> str:
    """
    Convert an analysis workflow filepath into its Python module name.

    Example
    -------

    src/iphasimulator/analysis/tg_analysis/workflow.py

    becomes

    iphasimulator.analysis.tg_analysis.workflow
    """

    workflow_path = (
        Path(workflow_path)
        .expanduser()
        .resolve()
    )


    project_root = (
        Path(__file__)
        .resolve()
        .parents[2]
    )


    src_directory = (
        project_root
        / "src"
    )


    relative_path = (
        workflow_path
        .relative_to(
            src_directory
        )
        .with_suffix("")
    )


    module_name = (
        ".".join(
            relative_path.parts
        )
    )


    return module_name


# =============================================================================
# Selected simulation path
# =============================================================================

def _get_simulation_directory(
    gui_data: GUIData,
    selected_system_name,
    selected_system_type,
    selected_simulation_name,
):
    """
    Resolve the selected simulation directory.
    """

    if (
        selected_system_name is None
        or selected_system_type is None
        or selected_simulation_name is None
    ):

        return None


    try:

        system_files = (
            gui_data.paths.get_md_system_files(
                system_name=(
                    selected_system_name
                ),
                system_type=(
                    selected_system_type
                ),
            )
        )


        simulations_directory = (
            Path(
                system_files[
                    "simulations_dir"
                ]
            )
        )


        simulation_directory = (
            simulations_directory
            / selected_simulation_name
        )


        return simulation_directory


    except Exception:

        return None


# =============================================================================
# Analysis output directory
# =============================================================================

def _get_analysis_output_directory(
    gui_data: GUIData,
    selected_system_name,
    selected_system_type,
    selected_simulation_name,
    selected_workflow_name,
):
    """
    Determine the expected replica-level output directory.

    Analysis workflows currently follow:

        <simulation>/analysis/<workflow_name>/
    """

    simulation_directory = (
        _get_simulation_directory(
            gui_data=(
                gui_data
            ),
            selected_system_name=(
                selected_system_name
            ),
            selected_system_type=(
                selected_system_type
            ),
            selected_simulation_name=(
                selected_simulation_name
            ),
        )
    )


    if (
        simulation_directory is None
        or selected_workflow_name is None
    ):

        return None


    return (
        simulation_directory
        / "analysis"
        / selected_workflow_name
    )


# =============================================================================
# Existing analysis detection
# =============================================================================

def _analysis_summary_path(
    analysis_output_directory,
):
    """
    Return the expected analysis-summary path.
    """

    if analysis_output_directory is None:

        return None


    return (
        Path(
            analysis_output_directory
        )
        / "analysis_summary.json"
    )


def _analysis_already_exists(
    analysis_output_directory,
) -> bool:
    """
    Determine whether a completed analysis summary exists.
    """

    summary_path = (
        _analysis_summary_path(
            analysis_output_directory
        )
    )


    return (
        summary_path is not None
        and summary_path.is_file()
    )


# =============================================================================
# Load analysis summary
# =============================================================================

def _load_analysis_summary(
    analysis_output_directory,
):
    """
    Load analysis_summary.json if available.
    """

    summary_path = (
        _analysis_summary_path(
            analysis_output_directory
        )
    )


    if (
        summary_path is None
        or not summary_path.is_file()
    ):

        return None


    try:

        with open(
            summary_path,
            "r",
            encoding="utf-8",
        ) as handle:

            return json.load(
                handle
            )


    except Exception:

        return None


# =============================================================================
# Workflow execution
# =============================================================================

def run_selected_analysis_workflow(
    workflow_path,
    selected_system_name,
    selected_simulation_name,
    generate_figures=True,
    analysis_stride=50,
):
    """
    Run the selected analysis workflow inside the iphasimulator environment.

    The Streamlit GUI itself runs inside the pha_gui environment, so the
    scientific analysis is launched as a separate Python process using the
    iphasimulator environment.
    """

    if isinstance(analysis_stride, bool) or not isinstance(analysis_stride, Integral) or analysis_stride < 1:
        raise ValueError("Analysis stride must be a positive whole number.")
    analysis_stride = int(analysis_stride)

    workflow_path = (
        Path(workflow_path)
        .expanduser()
        .resolve()
    )


    # =========================================================================
    # Resolve project paths
    # =========================================================================

    project_root = (
        Path(__file__)
        .resolve()
        .parents[2]
    )


    temp_directory = (
        project_root
        / "temp"
        / "gui_analysis"
    )


    temp_directory.mkdir(
        parents=True,
        exist_ok=True,
    )


    runner_script_path = (
        temp_directory
        / "analysis_runner.py"
    )


    # =========================================================================
    # Resolve selected workflow module
    # =========================================================================

    module_name = (
        _workflow_path_to_module_name(
            workflow_path
        )
    )


    # =========================================================================
    # Build runner script
    # =========================================================================

    runner_script = f'''#!/usr/bin/env python3

"""
Automatically generated iPHAsimulator analysis runner.
"""

import importlib


MODULE_NAME = {module_name!r}
SYSTEM_NAME = {selected_system_name!r}
SIMULATION_NAME = {selected_simulation_name!r}
GENERATE_FIGURES = {generate_figures!r}
ANALYSIS_STRIDE = {analysis_stride!r}


module = importlib.import_module(
    MODULE_NAME
)


if hasattr(
    module,
    "run_analysis",
):

    result = module.run_analysis(
        system_name=SYSTEM_NAME,
        simulation_name=SIMULATION_NAME,
        generate_figures=GENERATE_FIGURES,
        analysis_stride=ANALYSIS_STRIDE,
    )


elif hasattr(
    module,
    "run_tg_analysis",
):

    result = module.run_tg_analysis(
        system_name=SYSTEM_NAME,
        simulation_name=SIMULATION_NAME,
        generate_figures=GENERATE_FIGURES,
        analysis_stride=ANALYSIS_STRIDE,
    )


else:

    raise AttributeError(
        "Selected analysis workflow does not expose a supported "
        "entry point. Expected run_analysis(...) or run_tg_analysis(...)."
    )


print()
print("=" * 80)
print("GUI ANALYSIS RUNNER COMPLETE")
print("=" * 80)


if isinstance(
    result,
    dict,
):

    for key, value in result.items():

        print(
            f"{{key}}: {{value}}"
        )

else:

    print(
        result
    )
'''


    # =========================================================================
    # Save runner
    # =========================================================================

    runner_script_path.write_text(
        runner_script,
        encoding="utf-8",
    )


    # =========================================================================
    # Execute using iphasimulator environment
    # =========================================================================

    result = (
        run_python_script_with_iphasimulator(
            runner_script_path
        )
    )


    return result


# =============================================================================
# Analysis summary display
# =============================================================================

def _render_completed_analysis_summary(
    analysis_summary,
):
    """
    Display useful values from a completed analysis summary.

    The function is intentionally defensive so that other workflow types can
    still use the Analysis tab even if they do not contain Tg-specific data.
    """

    if not isinstance(
        analysis_summary,
        dict,
    ):

        return


    st.markdown(
        "### Analysis Results"
    )


    # -------------------------------------------------------------------------
    # Tg result
    # -------------------------------------------------------------------------

    tg_data = (
        analysis_summary.get(
            "tg",
            {},
        )
    )


    pca_data = (
        analysis_summary.get(
            "pca",
            {},
        )
    )


    dbscan_data = (
        analysis_summary.get(
            "dbscan",
            {},
        )
    )


    tg_value = (
        tg_data.get(
            "tg_K"
        )
        if isinstance(
            tg_data,
            dict,
        )
        else None
    )


    selected_components = (
        pca_data.get(
            "selected_components"
        )
        if isinstance(
            pca_data,
            dict,
        )
        else None
    )


    selected_min_samples = (
        dbscan_data.get(
            "selected_min_samples"
        )
        if isinstance(
            dbscan_data,
            dict,
        )
        else None
    )


    median_noise_fraction = (
        dbscan_data.get(
            "median_noise_fraction"
        )
        if isinstance(
            dbscan_data,
            dict,
        )
        else None
    )


    metrics = []


    if tg_value is not None:

        metrics.append(
            (
                "Estimated Tg",
                f"{float(tg_value):.2f} K",
            )
        )


    if selected_components is not None:

        metrics.append(
            (
                "Selected PCs",
                str(
                    selected_components
                ),
            )
        )


    if selected_min_samples is not None:

        metrics.append(
            (
                "DBSCAN min_samples",
                str(
                    selected_min_samples
                ),
            )
        )


    if median_noise_fraction is not None:

        metrics.append(
            (
                "Median noise fraction",
                f"{float(median_noise_fraction):.3f}",
            )
        )


    if metrics:

        metric_columns = (
            st.columns(
                len(
                    metrics
                )
            )
        )


        for column, (
            label,
            value,
        ) in zip(
            metric_columns,
            metrics,
        ):

            column.metric(
                label,
                value,
            )


    sampling = analysis_summary.get("sampling", {})
    saved_stride = sampling.get("stride") if isinstance(sampling, dict) else None
    if saved_stride is not None:
        st.caption(f"Saved results used an analysis stride of {saved_stride} saved frame(s).")

    with st.expander(
        "Show full analysis summary"
    ):

        st.json(
            analysis_summary
        )


# =============================================================================
# Main analysis tab
# =============================================================================

def render_analysis_tab(
    gui_data: GUIData,
) -> None:
    """
    Render the iPHAsimulator analysis tab.

    The tab allows the user to:

    - select a registered MD system
    - select an available analysis workflow
    - select a simulation replica
    - configure analysis options
    - execute the selected analysis workflow
    - inspect completed analysis outputs
    """

    # =========================================================================
    # Header
    # =========================================================================

    st.markdown(
        "## 📊 Analysis"
    )


    render_info_box(
        "Select a prepared molecular-dynamics system, "
        "choose an analysis workflow, then run or inspect one replica "
        "or compare results already saved for several replicas."
    )


    st.divider()


    # =========================================================================
    # System and workflow
    # =========================================================================

    system_column, workflow_column = (
        st.columns(
            [
                1,
                1,
            ]
        )
    )


    with system_column:

        (
            selected_system_name,
            selected_system_type,
        ) = (
            _render_system_selection(
                gui_data
            )
        )


    with workflow_column:

        (
            selected_workflow_name,
            selected_workflow_path,
        ) = (
            _render_analysis_workflow_selection()
        )


    st.divider()


    analysis_view = st.radio(
        "Analysis view",
        ["Run / inspect one replica", "Compare saved replicas"],
        horizontal=True,
        key="analysis_view",
    )

    if analysis_view == "Compare saved replicas":
        render_analysis_comparison(
            gui_data, selected_system_name, selected_system_type, selected_workflow_name,
        )
        return

    # =========================================================================
    # Simulation
    # =========================================================================

    selected_simulation_name = (
        _render_simulation_selection(
            gui_data=(
                gui_data
            ),
            selected_system_name=(
                selected_system_name
            ),
            selected_system_type=(
                selected_system_type
            ),
        )
    )


    st.divider()


    # =========================================================================
    # Analysis options
    # =========================================================================

    st.markdown(
        "### 4. Analysis Options"
    )

    analysis_stride = st.number_input(
        "Analysis stride (saved frames)",
        min_value=1, value=50, step=1,
        key="analysis_stride",
        help=(
            "Use every Nth saved trajectory frame. 1 uses every frame; "
            "50 uses every 50th frame. Larger strides analyse fewer frames, "
            "usually reducing runtime and memory use, but can change the results. "
            "This does not change the simulation timestep or saved trajectory."
        ),
    )
    st.caption(
        "Next run: use every saved frame (stride 1)."
        if analysis_stride == 1 else
        f"Next run: analyse one in every {analysis_stride} saved frames."
    )


    generate_figures = (
        st.toggle(
            "Generate analysis figures",
            value=True,
            help=(
                "Generate and save the figures associated with "
                "the selected analysis workflow."
            ),
            key="analysis_generate_figures",
        )
    )


    if generate_figures:

        st.caption(
            "Figures will be generated and saved with the analysis outputs."
        )

    else:

        st.caption(
            "The analysis will run without generating new figures."
        )


    # =========================================================================
    # Resolve output directory
    # =========================================================================

    analysis_output_directory = (
        _get_analysis_output_directory(
            gui_data=(
                gui_data
            ),
            selected_system_name=(
                selected_system_name
            ),
            selected_system_type=(
                selected_system_type
            ),
            selected_simulation_name=(
                selected_simulation_name
            ),
            selected_workflow_name=(
                selected_workflow_name
            ),
        )
    )


    existing_analysis = (
        _analysis_already_exists(
            analysis_output_directory
        )
    )


    st.divider()


    # =========================================================================
    # Current analysis selection
    # =========================================================================

    st.markdown(
        "### 5. Current Analysis Selection"
    )


    selection_columns = (
        st.columns(
            [
                1,
                1,
                1,
                1,
            ]
        )
    )


    selection_columns[0].write(
        "**System**"
    )

    selection_columns[0].code(
        (
            selected_system_name
            if selected_system_name
            is not None
            else "Not selected"
        )
    )


    selection_columns[1].write(
        "**Simulation / Replica**"
    )

    selection_columns[1].code(
        (
            selected_simulation_name
            if selected_simulation_name
            is not None
            else "Not selected"
        )
    )


    selection_columns[2].write(
        "**Analysis workflow**"
    )

    selection_columns[2].code(
        (
            selected_workflow_name
            if selected_workflow_name
            is not None
            else "Not selected"
        )
    )


    selection_columns[3].write(
        "**Generate figures**"
    )

    selection_columns[3].code(
        str(
            generate_figures
        )
    )


    # =========================================================================
    # Output location / existing analysis
    # =========================================================================

    if analysis_output_directory is not None:

        st.write(
            "**Analysis output directory:**"
        )

        st.code(
            str(
                analysis_output_directory
            )
        )


        if existing_analysis:

            st.success(
                "A completed analysis already exists for this "
                "workflow and simulation."
            )

            st.caption(
                "Changing the stride applies to the next run. Re-run analysis "
                "to update the results in this replica's existing analysis folder."
            )


            if not generate_figures:

                st.warning(
                    "Existing figure files are not automatically deleted. "
                    "Running with figure generation disabled will prevent new "
                    "figures from being created, but figures from an earlier "
                    "analysis may remain in the output directory."
                )

        else:

            st.info(
                "No existing completed analysis was found for "
                "this workflow and simulation."
            )


    st.divider()


    # =========================================================================
    # Run analysis
    # =========================================================================

    st.markdown(
        "### 6. Run Analysis"
    )


    ready_to_run = (
        selected_system_name
        is not None
        and selected_system_type
        is not None
        and selected_simulation_name
        is not None
        and selected_workflow_name
        is not None
        and selected_workflow_path
        is not None
    )


    if not ready_to_run:

        render_info_box(
            "Select a system, simulation replica, "
            "and analysis workflow before running the analysis."
        )


    if existing_analysis:

        run_button_label = (
            "🔄 Re-run analysis"
        )

    else:

        run_button_label = (
            "▶️ Run analysis"
        )


    run_clicked = (
        st.button(
            run_button_label,
            use_container_width=True,
            disabled=(
                not ready_to_run
            ),
            key="run_analysis_workflow",
        )
    )


    # =========================================================================
    # Execute analysis
    # =========================================================================

    if run_clicked:

        st.info(
            "Running the selected analysis workflow "
            "inside the iphasimulator environment."
        )


        progress = (
            st.progress(
                0
            )
        )


        status = (
            st.empty()
        )


        result = None


        try:

            status.write(
                "Launching analysis workflow..."
            )


            progress.progress(
                20
            )


            with st.spinner(
                "Running analysis workflow..."
            ):

                result = (
                    run_selected_analysis_workflow(
                        workflow_path=(
                            selected_workflow_path
                        ),
                        selected_system_name=(
                            selected_system_name
                        ),
                        selected_simulation_name=(
                            selected_simulation_name
                        ),
                        generate_figures=(
                            generate_figures
                        ),
                        analysis_stride=analysis_stride,
                    )
                )


            progress.progress(
                90
            )


            status.write(
                "Analysis process finished."
            )


        except Exception as error:

            progress.progress(
                100
            )


            st.error(
                "Could not launch the analysis workflow."
            )


            st.code(
                str(
                    error
                )
            )


            return


        progress.progress(
            100
        )


        if result is None:

            st.error(
                "The analysis process did not return a result."
            )

            return


        # =====================================================================
        # Successful / failed process
        # =====================================================================

        if result.returncode == 0:

            status.write(
                "Analysis complete."
            )


            st.success(
                "Analysis completed successfully."
            )

        else:

            status.write(
                "Analysis failed."
            )


            st.error(
                "Analysis workflow failed."
            )


        # =====================================================================
        # Load structured analysis output
        # =====================================================================

        if result.returncode == 0:

            analysis_summary = (
                _load_analysis_summary(
                    analysis_output_directory
                )
            )


            if analysis_summary is not None:

                st.divider()


                _render_completed_analysis_summary(
                    analysis_summary
                )

            else:

                st.info(
                    "The workflow completed successfully, but no "
                    "analysis_summary.json file was found."
                )


        # =====================================================================
        # Process / debugging information
        # =====================================================================

        st.divider()


        st.markdown(
            "### Process Information"
        )


        result_columns = (
            st.columns(
                [
                    1,
                    1,
                    1,
                ]
            )
        )


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


        # ---------------------------------------------------------------------
        # STDOUT
        # ---------------------------------------------------------------------

        if result.stdout:

            with st.expander(
                "Console output",
                expanded=(
                    result.returncode
                    != 0
                ),
            ):

                st.code(
                    result.stdout,
                    language="text",
                )


        # ---------------------------------------------------------------------
        # STDERR
        # ---------------------------------------------------------------------

        if result.stderr:

            with st.expander(
                "Warnings / errors",
                expanded=(
                    result.returncode
                    != 0
                ),
            ):

                st.code(
                    result.stderr,
                    language="text",
                )


    # =========================================================================
    # Existing analysis preview
    # =========================================================================

    elif existing_analysis:

        existing_summary = (
            _load_analysis_summary(
                analysis_output_directory
            )
        )


        if existing_summary is not None:

            st.divider()


            st.markdown(
                "### Existing Analysis"
            )


            render_info_box(
                "The results below were loaded from the existing "
                "analysis for the currently selected simulation."
            )


            _render_completed_analysis_summary(
                existing_summary
            )
