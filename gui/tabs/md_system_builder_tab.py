#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MD System Builder tab for the iPHAsimulatorV2 Streamlit GUI.

Scientific build operations are executed as subprocesses using the dedicated
iPHAsimulator environment. The Streamlit GUI therefore does not directly
import AmberTools, OpenMM, ACPYPE, Polyply, or the scientific builder modules.
"""

from __future__ import annotations

import re

import streamlit as st

from gui.config import (
    IPHASIMULATOR_PYTHON,
    STRUCTURE_DATABASE,
)
from gui.md_system_helpers import (
    build_dry_system_subprocess,
    build_melt_system_subprocess,
    build_solvated_ions_system_subprocess,
    build_solvated_system_subprocess,
    get_built_polymer_choices,
    get_incomplete_polymer_builds,
    get_salt_name,
    get_water_model_settings,
)


_RESULT_MARKER = "MD_SYSTEM_BUILD_RESULT"


def _extract_result_json(stdout: str) -> str | None:
    """
    Extract the JSON result printed by the scientific backend.

    The subprocess prints:

        MD_SYSTEM_BUILD_RESULT
        {
            ...
        }

    Parameters
    ----------
    stdout : str
        Captured subprocess standard output.

    Returns
    -------
    str or None
        JSON text following the result marker.
    """

    if _RESULT_MARKER not in stdout:
        return None

    _, result_text = stdout.split(
        _RESULT_MARKER,
        maxsplit=1,
    )

    result_text = result_text.strip()

    return result_text or None


def _display_process_output(
    process,
) -> None:
    """
    Display captured backend stdout and stderr.
    """

    st.markdown(
        "### Backend Process"
    )

    metric_columns = st.columns(3)

    metric_columns[0].metric(
        "Return code",
        process.returncode,
    )

    metric_columns[1].metric(
        "Standard output",
        f"{len(process.stdout or '')} characters",
    )

    metric_columns[2].metric(
        "Standard error",
        f"{len(process.stderr or '')} characters",
    )

    if process.stdout:
        with st.expander(
            "Backend standard output",
            expanded=process.returncode != 0,
        ):
            st.code(
                process.stdout,
                language="text",
            )

    if process.stderr:
        with st.expander(
            "Backend standard error",
            expanded=True,
        ):
            st.code(
                process.stderr,
                language="text",
            )


def _display_success_summary(
    *,
    proposed_system_name: str,
    process,
) -> None:
    """
    Display the successful build summary.
    """

    st.success(
        f"Successfully built `{proposed_system_name}`."
    )

    st.write(
        f"**Structure database:** `{STRUCTURE_DATABASE}`"
    )

    result_json = _extract_result_json(
        process.stdout or ""
    )

    if result_json:
        with st.expander(
            "Backend build result",
            expanded=False,
        ):
            st.code(
                result_json,
                language="json",
            )

    st.info(
        "The completed system has been registered in "
        "`structure_database/md_systems.csv`."
    )


def _render_environment_status() -> None:
    """
    Show the scientific subprocess environment configuration.
    """

    environment_exists = (
        IPHASIMULATOR_PYTHON.exists()
        and IPHASIMULATOR_PYTHON.is_file()
    )

    status_columns = st.columns(2)

    status_columns[0].write(
        "**Scientific Python**"
    )

    status_columns[0].code(
        str(IPHASIMULATOR_PYTHON)
    )

    status_columns[1].write(
        "**Environment status**"
    )

    if environment_exists:
        status_columns[1].success(
            "✓ iPHAsimulator environment found"
        )

    else:
        status_columns[1].error(
            "✗ iPHAsimulator environment not found"
        )


def render_md_system_builder_tab() -> None:
    """
    Render the MD System Builder tab.
    """

    st.markdown(
        "## 🧪 MD System Builder"
    )

    st.markdown(
        """
Build simulation-ready molecular dynamics systems from previously generated
PHA polymers.

All scientific operations run as subprocesses in the dedicated
`iphasimulator` environment.
"""
    )

    with st.expander(
        "📘 How this page works",
        expanded=False,
    ):
        st.markdown(
            """
### Workflow

1. Select a polymer type and enter its chain length. Only completed builds can be used.
2. Select the required MD system type.
3. Configure the system-building parameters.
4. Press **Build MD System**.
5. The GUI launches the scientific backend in the separate
   `iphasimulator` environment.
6. The completed system is stored in the structure database and added to
   `md_systems.csv`.

### Available system types

- **Dry PHA** — a single polymer chain in a periodic box.
- **Solvated PHA** — a single polymer chain solvated with water.
- **Solvated PHA with ions** — a water-solvated chain containing a target
  salt concentration.
- **Polymer melt** — one or more chains packed into an amorphous bulk system.
"""
        )

    _render_environment_status()

    st.divider()

    polymers = get_built_polymer_choices(
        structure_database=STRUCTURE_DATABASE,
    )

    incomplete = get_incomplete_polymer_builds(STRUCTURE_DATABASE)
    if incomplete:
        with st.expander("Polymers with incomplete build files"):
            st.caption("A SMILES entry or PDB alone does not mean a polymer build completed. "
                       "Check its build log and rebuild before using it here.")
            for name, reason in incomplete.items():
                st.write(f"**{name}**: {reason}.")

    if not polymers:
        st.warning(
            "No complete built polymers were found."
        )

        st.info(
            "Build a polymer using the Polymer Builder and Build Console "
            "before creating an MD system."
        )

        return

    summary_columns = st.columns(3)

    summary_columns[0].metric(
        "Available polymer types",
        len(polymers),
    )

    summary_columns[1].metric(
        "Build environment",
        "iphasimulator",
    )

    summary_columns[2].metric(
        "Database",
        STRUCTURE_DATABASE.name,
    )

    st.divider()

    system_type = st.selectbox(
        "MD system type",
        [
            "Dry PHA",
            "Solvated PHA",
            "Solvated PHA with ions",
            "Polymer melt",
        ],
        key="md_builder_system_type",
    )

    polymer_type = st.selectbox(
        "Polymer type", list(polymers), key="md_builder_polymer_type",
    )
    lengths = polymers[polymer_type]
    st.caption("Available lengths: " + ", ".join(str(length) for length in lengths))
    length = int(st.number_input(
        "Polymer length (monomer units)", min_value=1, value=min(lengths), step=1,
        key=f"md_builder_length_{polymer_type}",
    ))
    polymer_name = lengths.get(length)
    if polymer_name is None:
        st.error(f"{polymer_type} with {length} monomer units has not been built successfully. "
                 "Build this length in the Polymer Builder first, or choose an available length.")
        st.button("🔨 Build MD System", disabled=True, key="md_builder_build_button")
        return

    build_function = None
    build_arguments = {}
    proposed_system_name = None

    # ======================================================
    # Dry PHA
    # ======================================================

    if system_type == "Dry PHA":
        st.markdown(
            "### Dry single-chain system"
        )

        option_columns = st.columns(2)

        with option_columns[0]:
            forcefield = st.selectbox(
                "Force field",
                ["gaff2"],
                key="md_builder_dry_forcefield",
            )

        with option_columns[1]:
            box_radius = st.number_input(
                "Box radius (Å)",
                min_value=5.0,
                value=20.0,
                step=1.0,
                key="md_builder_dry_box_radius",
            )

        proposed_system_name = (
            f"{polymer_name}_dry"
        )

        build_function = (
            build_dry_system_subprocess
        )

        build_arguments = {
            "polymer_name": polymer_name,
            "forcefield": forcefield,
            "box_radius": float(box_radius),
        }

    # ======================================================
    # Solvated PHA
    # ======================================================

    elif system_type == "Solvated PHA":
        st.markdown(
            "### Solvated single-chain system"
        )

        option_columns = st.columns(3)

        with option_columns[0]:
            forcefield = st.selectbox(
                "Force field",
                ["gaff2"],
                key="md_builder_solvated_forcefield",
            )

        with option_columns[1]:
            water_model = st.selectbox(
                "Water model",
                ["TIP3P"],
                key="md_builder_solvated_water",
            )

        with option_columns[2]:
            box_radius = st.number_input(
                "Solvent radius (Å)",
                min_value=5.0,
                value=20.0,
                step=1.0,
                key="md_builder_solvated_radius",
            )

        water_leaprc, water_box = (
            get_water_model_settings(
                water_model
            )
        )

        proposed_system_name = (
            f"{polymer_name}_solvated"
        )

        build_function = (
            build_solvated_system_subprocess
        )

        build_arguments = {
            "polymer_name": polymer_name,
            "forcefield": forcefield,
            "water_leaprc": water_leaprc,
            "water_box": water_box,
            "box_radius": float(box_radius),
        }

    # ======================================================
    # Solvated PHA with ions
    # ======================================================

    elif system_type == "Solvated PHA with ions":
        st.markdown(
            "### Solvated and ionised single-chain system"
        )

        first_row = st.columns(3)

        with first_row[0]:
            forcefield = st.selectbox(
                "Force field",
                ["gaff2"],
                key="md_builder_ions_forcefield",
            )

        with first_row[1]:
            water_model = st.selectbox(
                "Water model",
                ["TIP3P"],
                key="md_builder_ions_water",
            )

        with first_row[2]:
            box_radius = st.number_input(
                "Solvent radius (Å)",
                min_value=5.0,
                value=20.0,
                step=1.0,
                key="md_builder_ions_radius",
            )

        second_row = st.columns(3)

        with second_row[0]:
            positive_ion = st.selectbox(
                "Positive ion",
                [
                    "K+",
                    "Na+",
                ],
                key="md_builder_positive_ion",
            )

        with second_row[1]:
            negative_ion = st.selectbox(
                "Negative ion",
                [
                    "Cl-",
                ],
                key="md_builder_negative_ion",
            )

        with second_row[2]:
            ion_concentration = st.number_input(
                "Ion concentration (mol L⁻¹)",
                min_value=0.0,
                value=0.15,
                step=0.05,
                format="%.3f",
                key="md_builder_ion_concentration",
            )

        water_leaprc, water_box = (
            get_water_model_settings(
                water_model
            )
        )

        salt = get_salt_name(
            positive_ion,
            negative_ion,
        )

        concentration_label = (
            str(float(ion_concentration))
            .rstrip("0")
            .rstrip(".")
            .replace(".", "_")
        )

        proposed_system_name = (
            f"{polymer_name}_solvated_"
            f"{salt}_{concentration_label}"
        )

        build_function = (
            build_solvated_ions_system_subprocess
        )

        build_arguments = {
            "polymer_name": polymer_name,
            "forcefield": forcefield,
            "water_leaprc": water_leaprc,
            "water_box": water_box,
            "box_radius": float(box_radius),
            "salt": salt,
            "positive_ion": positive_ion,
            "negative_ion": negative_ion,
            "ion_concentration": float(
                ion_concentration
            ),
        }

    # ======================================================
    # Polymer melt
    # ======================================================

    elif system_type == "Polymer melt":
        st.markdown(
            "### Amorphous polymer melt"
        )

        st.caption(
            "The current GUI version builds a melt containing one polymer "
            "type. Multi-component melts can be added later."
        )

        option_columns = st.columns(2)

        with option_columns[0]:
            number_of_polymers = st.number_input(
                "Number of chains",
                min_value=1,
                value=25,
                step=1,
                key="md_builder_melt_chains",
            )

        with option_columns[1]:
            density = st.number_input(
                "Target density (kg m⁻³)",
                min_value=100.0,
                value=750.0,
                step=10.0,
                key="md_builder_melt_density",
            )

        proposed_system_name = (
            f"{int(number_of_polymers)}_"
            f"{polymer_name}_melt"
        )

        build_function = (
            build_melt_system_subprocess
        )

        build_arguments = {
            "polymer_names": [
                polymer_name
            ],
            "number_of_polymers": [
                int(number_of_polymers)
            ],
            "density": float(density),
        }

    # ======================================================
    # Build target
    # ======================================================

    st.divider()

    st.markdown(
        "### Build Target"
    )

    if proposed_system_name:
        st.code(
            proposed_system_name
        )

    st.write(
        f"**Output database:** `{STRUCTURE_DATABASE}`"
    )

    environment_ready = (
        IPHASIMULATOR_PYTHON.exists()
        and IPHASIMULATOR_PYTHON.is_file()
    )

    build_clicked = st.button(
        "🔨 Build MD System",
        type="primary",
        use_container_width=True,
        disabled=not environment_ready,
        key="md_builder_build_button",
    )

    # ======================================================
    # Run backend subprocess
    # ======================================================

    if build_clicked:
        if build_function is None:
            st.error(
                "No MD system build function was selected."
            )

            return

        try:
            with st.spinner(
                f"Building {proposed_system_name} in the "
                "iPHAsimulator environment..."
            ):
                process = build_function(
                    **build_arguments
                )

            st.session_state[
                "md_builder_last_process"
            ] = process

            st.session_state[
                "md_builder_last_system_name"
            ] = proposed_system_name

        except Exception as error:
            st.error(
                "The GUI could not launch the scientific backend."
            )

            st.exception(
                error
            )

            return

    # ======================================================
    # Last build result
    # ======================================================

    process = st.session_state.get(
        "md_builder_last_process"
    )

    last_system_name = st.session_state.get(
        "md_builder_last_system_name"
    )

    if process is not None:
        st.divider()

        if process.returncode == 0:
            _display_success_summary(
                proposed_system_name=(
                    last_system_name
                    or "MD system"
                ),
                process=process,
            )

        else:
            st.error(
                f"The backend failed with return code "
                f"`{process.returncode}`."
            )

        _display_process_output(
            process
        )