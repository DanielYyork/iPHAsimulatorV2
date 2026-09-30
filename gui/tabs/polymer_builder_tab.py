#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Polymer Builder tab for the iPHAsimulatorV2 Streamlit GUI.

This tab allows the user to:

- browse registered PHA monomers
- filter available monomers
- construct an ordered polymer sequence
- remove or clear sequence entries
- repeat the current sequence
- inspect sequence details
- preview the generated polymer name
- preview the generated polymer SMILES
"""

import streamlit as st

from gui.models import GUIData
from gui.polymer_helpers import (
    draw_monomer,
    generate_polymer_smiles_from_sequence,
    get_polymer_name,
    repeat_sequence,
)
from gui.state import (
    clear_polymer_sequence,
    set_polymer_sequence,
    undo_polymer_sequence,
)
from gui.styles import render_warning_box
from pathlib import Path

import streamlit as st

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def _get_monomer_information(
    pha_type,
    monomer_smiles,
    gui_data,
):
    """
    Collect calculated and database-backed information for one monomer.

    Parameters
    ----------
    pha_type : str
        Registered PHA type, for example ``3HB``.

    monomer_smiles : dict
        Mapping between PHA type and monomer SMILES.

    gui_data : GUIData
        Shared GUI data containing the filepath manager and monomer table.

    Returns
    -------
    dict
        Monomer properties and parameterisation status.
    """

    smiles = monomer_smiles.get(
        pha_type,
        "",
    )

    molecule = Chem.MolFromSmiles(
        smiles
    )

    if molecule is None:
        molecular_formula = "Unavailable"
        molecular_weight = None
        stereocentre_count = None

    else:
        molecular_formula = CalcMolFormula(
            molecule
        )

        molecular_weight = Descriptors.MolWt(
            molecule
        )

        stereocentres = Chem.FindMolChiralCenters(
            molecule,
            includeUnassigned=True,
        )

        stereocentre_count = len(
            stereocentres
        )

    residue_code = "Unavailable"

    matching_rows = gui_data.mainchain_df[
        gui_data.mainchain_df["PHA_type"] == pha_type
    ]

    if not matching_rows.empty:
        first_row = matching_rows.iloc[0]

        if "residue_code" in first_row.index:
            value = first_row["residue_code"]

            if value is not None and str(value).strip():
                residue_code = str(value).strip()

    parameter_files = (
        gui_data.paths.get_PHA_monomer_unit_files(
            pha_type
        )
    )

    head_prepin_exists = Path(
        parameter_files["head_prepin"]
    ).is_file()

    mainchain_prepin_exists = Path(
        parameter_files["mainchain_prepin"]
    ).is_file()

    tail_prepin_exists = Path(
        parameter_files["tail_prepin"]
    ).is_file()

    frcmod_exists = Path(
        parameter_files["frcmod"]
    ).is_file()

    monomer_units_ready = all(
        [
            head_prepin_exists,
            mainchain_prepin_exists,
            tail_prepin_exists,
        ]
    )

    polymer_builder_ready = (
        monomer_units_ready
        and frcmod_exists
    )

    return {
        "smiles": smiles,
        "formula": molecular_formula,
        "molecular_weight": molecular_weight,
        "stereocentre_count": stereocentre_count,
        "residue_code": residue_code,
        "head_prepin_exists": head_prepin_exists,
        "mainchain_prepin_exists": mainchain_prepin_exists,
        "tail_prepin_exists": tail_prepin_exists,
        "frcmod_exists": frcmod_exists,
        "monomer_units_ready": monomer_units_ready,
        "polymer_builder_ready": polymer_builder_ready,
    }

def _render_monomer_card(
    pha_type,
    monomer_smiles,
    monomer_information,
    add_count=1,
):
    """
    Render a monomer information popover with a quick-add button below it.

    The popover contains calculated molecular properties
    and parameterisation status. Adding a monomer does not require opening it.
    """

    with st.popover(
        pha_type,
        use_container_width=True,
    ):
        st.markdown(
            f"## {pha_type}"
        )

        # ==================================================
        # Molecular structure
        # ==================================================

        try:
            monomer_image = draw_monomer(
                pha_type=pha_type,
                monomer_smiles=monomer_smiles,
            )

        except Exception as error:
            monomer_image = None

            st.error(
                "Could not generate the monomer preview."
            )

            st.code(
                str(error)
            )

        if monomer_image is not None:
            st.image(
                monomer_image,
                use_container_width=True,
            )

        else:
            st.warning(
                "Structure preview unavailable."
            )

        # ==================================================
        # Calculated molecular properties
        # ==================================================

        st.markdown(
            "### Molecular Information"
        )

        property_columns = st.columns(2)

        with property_columns[0]:
            st.write(
                "**Formula**"
            )

            st.code(
                monomer_information["formula"]
            )

            st.write(
                "**Residue code**"
            )

            st.code(
                monomer_information["residue_code"]
            )

        with property_columns[1]:
            st.write(
                "**Molecular weight**"
            )

            molecular_weight = (
                monomer_information[
                    "molecular_weight"
                ]
            )

            if molecular_weight is None:
                st.code(
                    "Unavailable"
                )

            else:
                st.code(
                    f"{molecular_weight:.3f} g mol⁻¹"
                )

            st.write(
                "**Stereocentres**"
            )

            stereocentre_count = (
                monomer_information[
                    "stereocentre_count"
                ]
            )

            if stereocentre_count is None:
                st.code(
                    "Unavailable"
                )

            else:
                st.code(
                    str(stereocentre_count)
                )

        # ==================================================
        # SMILES
        # ==================================================

        st.markdown(
            "### Monomer SMILES"
        )

        if monomer_information["smiles"]:
            st.code(
                monomer_information["smiles"]
            )

        else:
            st.warning(
                "No monomer SMILES is available."
            )

        # ==================================================
        # Parameterisation status
        # ==================================================

        st.markdown(
            "### Parameterisation Status"
        )

        status_rows = [
            (
                "Head prepin",
                monomer_information[
                    "head_prepin_exists"
                ],
            ),
            (
                "Mainchain prepin",
                monomer_information[
                    "mainchain_prepin_exists"
                ],
            ),
            (
                "Tail prepin",
                monomer_information[
                    "tail_prepin_exists"
                ],
            ),
            (
                "Force-field parameters",
                monomer_information[
                    "frcmod_exists"
                ],
            ),
        ]

        for status_name, status_value in status_rows:
            if status_value:
                st.success(
                    f"✓ {status_name}"
                )

            else:
                st.warning(
                    f"✗ {status_name}"
                )

        if monomer_information["polymer_builder_ready"]:
            st.success(
                "✓ Ready for polymer building"
            )

        else:
            st.error(
                "This monomer is not yet fully ready "
                "for polymer building."
            )

    # Status and adding remain visible while the detailed card is closed.
    ready = monomer_information["polymer_builder_ready"]
    status_class = "monomer-ready" if ready else "monomer-incomplete"
    status_text = "✓ Ready" if ready else "! Needs files"
    st.markdown(
        f'<span class="monomer-status {status_class}">{status_text}</span>',
        unsafe_allow_html=True,
    )
    st.button(
        "＋ Add" if add_count == 1 else f"＋ Add ×{add_count}",
        key=f"add_monomer_{pha_type}",
        help=f"Append {add_count} × {pha_type} to the polymer sequence.",
        on_click=set_polymer_sequence,
        args=(list(st.session_state.sequence) + [pha_type] * add_count,),
    )


def _render_sequence_editor(sequence, available_phas):
    """Edit a selected position without rendering a form for every unit."""
    if not sequence:
        return
    with st.expander("Edit individual monomers"):
        # New widget identities after an edit prevent stale positions or choices.
        revision = st.session_state.sequence_revision
        position = st.selectbox(
            "Position to edit",
            range(len(sequence)),
            format_func=lambda i: f"{i + 1}. {sequence[i]}",
            key=f"sequence_edit_position_{revision}",
        )
        key = f"sequence_edit_{revision}_{position}"
        replace_column, move_column, delete_column = st.columns(3)
        with replace_column:
            choices = list(dict.fromkeys([sequence[position], *available_phas]))
            replacement = st.selectbox("Replace with", choices, key=f"{key}_replacement")
            replaced = list(sequence)
            replaced[position] = replacement
            st.button(
                "Replace unit", key="sequence_replace",
                disabled=replacement == sequence[position],
                on_click=set_polymer_sequence, args=(replaced,),
            )
        with move_column:
            destination = st.number_input(
                "Move to position", min_value=1, max_value=len(sequence),
                value=position + 1, step=1, key=f"{key}_destination",
                help="The selected monomer will occupy this position in the final sequence.",
            )
            moved = list(sequence)
            moved.insert(destination - 1, moved.pop(position))
            st.button(
                "Move unit", key="sequence_move", disabled=destination == position + 1,
                on_click=set_polymer_sequence, args=(moved,),
            )
        with delete_column:
            st.write(f"Selected unit: **{sequence[position]}**")
            st.caption(f"Delete position {position + 1}. Later units move left.")
            st.button(
                "Delete unit", key="sequence_delete",
                on_click=set_polymer_sequence,
                args=(list(sequence[:position]) + list(sequence[position + 1:]),),
            )
        st.caption("Undo restores the sequence before an edit, including a deletion or clear.")


def render_polymer_builder_tab(
    gui_data: GUIData,
) -> None:
    """
    Render the interactive PHA polymer-builder tab.

    Parameters
    ----------
    gui_data : GUIData
        Shared GUI data containing the available PHA monomers and the
        monomer-SMILES lookup.
    """
    # ======================================================
    # Page guidance
    # ======================================================

    st.markdown(
        "## 🧱 Polymer Builder"
    )

    st.markdown(
        """
<div class="info-box">
Construct homopolymer and copolymer sequences interactively before
generating parameterised PHA structures.
</div>
""",
        unsafe_allow_html=True,
    )

    with st.expander(
        "📘 How to use this page",
        expanded=False,
    ):
        st.markdown(
            """
### Workflow

**1. Browse the registered monomers**

Use the search box to filter the available PHA monomers.

Click a monomer name to open its information panel. This shows:

- the monomer structure
- molecular formula and molecular weight
- stereocentre count
- residue code
- monomer SMILES
- parameterisation status

**2. Build a polymer sequence**

Press **＋ Add** below a monomer name to append that
monomer to the current sequence. Set **Units per add** to append several
copies in one click. **Ready** means the required parameter files were found.

Continue adding monomers in the required order. This can be used to
construct:

- homopolymers
- alternating copolymers
- block copolymers
- custom copolymer sequences

**3. Edit the sequence**

Use **Remove last** to remove the most recently added monomer.

Use **Clear** to delete the complete sequence.

Open **Edit individual monomers** to replace, move or delete a chosen position.
Use **Undo** to reverse the last change, including adding, repeating or clearing.

Use **Repeat sequence** to duplicate the current sequence one or more
times. This is useful for quickly generating longer chains or repeating
copolymer patterns.

**4. Review the generated polymer**

The polymer name and polymer SMILES are updated automatically whenever
the sequence changes.

The ordered sequence can also be inspected in the sequence-details
table.

**5. Preview and build the polymer**

Open the **Molecular Preview** tab to inspect the generated polymer
structure.

Open the **Build Console** tab when the sequence is ready to generate
the parameterised polymer using iPHAsimulator.
"""
        )

        st.markdown(
            "### Tips"
        )

        st.markdown(
            """
- A green parameterisation status means that the required monomer-unit
  and force-field files were found.
- Monomers that are not fully parameterised may not build successfully.
- The order in which monomers are added determines the final copolymer
  sequence.
- Build short test polymers before generating very long chains.
"""
        )

    st.markdown(
        "### iPHAsimulator workflow"
    )

    workflow_columns = st.columns(6)

    workflow_steps = [
        ("1", "Choose monomers"),
        ("2", "Build sequence"),
        ("3", "Preview polymer"),
        ("4", "Build polymer"),
        ("5", "Prepare MD system"),
        ("6", "Run simulation"),
    ]

    for column, (step_number, step_name) in zip(
        workflow_columns,
        workflow_steps,
    ):
        with column:
            st.markdown(
                f"""
<div class="workflow-stage">
    <div class="workflow-stage-number">
        {step_number}
    </div>
    <div class="workflow-stage-label">
        {step_name}
    </div>
</div>
""",
                unsafe_allow_html=True,
            )

    st.divider()

    available_phas = gui_data.available_phas
    monomer_smiles = gui_data.monomer_smiles

    sequence = st.session_state.sequence

    polymer_name = get_polymer_name(
        sequence
    )

    polymer_smiles = generate_polymer_smiles_from_sequence(
        sequence=sequence,
        smiles_lookup=monomer_smiles,
    )

    # ======================================================
    # Polymer summary metrics
    # ======================================================

    top_metrics = st.columns([1, 1, 1.6, 1.6])

    top_metrics[0].metric(
        "Sequence length",
        len(sequence),
    )

    top_metrics[1].metric(
        "Unique monomers",
        len(set(sequence)),
    )

    if not sequence:
        polymer_mode = "None"

    elif len(set(sequence)) == 1:
        polymer_mode = "Homopolymer"

    else:
        polymer_mode = "Copolymer"

    top_metrics[2].metric(
        "Mode",
        polymer_mode,
    )

    top_metrics[3].metric(
        "Build environment",
        "iphasimulator",
    )

    # ======================================================
    # Available PHA monomers
    # ======================================================

    st.markdown(
        "### Available PHA Monomer Units"
    )

    search_column, count_column = st.columns([3, 1])
    search = search_column.text_input(
        "Filter monomers",
        placeholder=(
            "Try 3HB, 4HB, phenyl, fluor..."
        ),
        key="polymer_builder_monomer_search",
    )

    add_count = count_column.number_input(
        "Units per add", min_value=1, max_value=1000, value=1, step=1,
        key="polymer_builder_add_count",
        help="Applies to every Add button. A whole batch can be undone in one step.",
    )
    only_ready = st.checkbox("Show only ready monomers", key="polymer_builder_only_ready")
    st.caption("Readiness checks the required files. Open a monomer card to see which files are missing.")

    search_text = (
        search
        .strip()
        .lower()
    )

    filtered_phas = [
        pha_type
        for pha_type in available_phas
        if search_text in pha_type.lower()
    ]

    monomer_information = {
        pha_type: _get_monomer_information(pha_type, monomer_smiles, gui_data)
        for pha_type in filtered_phas
    }
    if only_ready:
        filtered_phas = [
            pha_type for pha_type in filtered_phas
            if monomer_information[pha_type]["polymer_builder_ready"]
        ]

    if not filtered_phas:
        st.info(
            "No registered PHA monomers match the current filter."
        )

    else:
        number_of_columns = 6

        monomer_columns = st.columns(
            number_of_columns
        )

        for index, pha_type in enumerate(
            filtered_phas
        ):
            with monomer_columns[
                index % number_of_columns
            ]:
                _render_monomer_card(
                    pha_type=pha_type,
                    monomer_smiles=monomer_smiles,
                    monomer_information=monomer_information[pha_type],
                    add_count=add_count,
                )

    # ======================================================
    # Current sequence
    # ======================================================

    st.markdown(
        "### Current Sequence"
    )

    if sequence:
        sequence_chip_html = "".join(
            (
                '<span class="sequence-chip">'
                f"{index + 1}. {unit}"
                "</span>"
            )
            for index, unit
            in enumerate(sequence)
        )

        st.markdown(
            sequence_chip_html,
            unsafe_allow_html=True,
        )

    else:
        render_warning_box(
            "Click ＋ Add below a monomer name to begin "
            "building a polymer sequence."
        )

    # ======================================================
    # Sequence controls
    # ======================================================

    control_columns = st.columns([1, 1, 1, 1, 1.5])

    with control_columns[0]:
        st.button(
            "↩️ Remove last",
            use_container_width=True,
            key="polymer_builder_remove_last",
            disabled=not sequence,
            on_click=set_polymer_sequence, args=(list(sequence[:-1]),),
        )

    with control_columns[1]:
        st.button(
            "🧹 Clear",
            use_container_width=True,
            key="polymer_builder_clear",
            disabled=not sequence,
            on_click=clear_polymer_sequence,
        )

    with control_columns[2]:
        st.button(
            "Undo", key="polymer_builder_undo", use_container_width=True,
            disabled=not st.session_state.sequence_history,
            help="Undo the last edit. The last 50 changes are kept for this session.",
            on_click=undo_polymer_sequence,
        )

    with control_columns[3]:
        repeat_count = st.number_input(
            "Repeat sequence",
            min_value=1,
            max_value=20,
            value=1,
            step=1,
            key="polymer_builder_repeat_count",
        )

    with control_columns[4]:
        st.button(
            "🔁 Apply repeat",
            use_container_width=True,
            key="polymer_builder_apply_repeat",
            disabled=not sequence or repeat_count == 1,
            on_click=set_polymer_sequence,
            args=(repeat_sequence(sequence, repeat_count),),
        )

    _render_sequence_editor(sequence, available_phas)

    # ======================================================
    # Sequence details
    # ======================================================

    if sequence:
        with st.expander(
            "Show sequence details"
        ):
            sequence_rows = []

            for index, pha_type in enumerate(
                sequence,
                start=1,
            ):
                sequence_rows.append(
                    {
                        "Position": index,
                        "PHA type": pha_type,
                        "SMILES": monomer_smiles.get(
                            pha_type,
                            "",
                        ),
                    }
                )

            st.dataframe(
                sequence_rows,
                use_container_width=True,
                hide_index=True,
            )

    # ======================================================
    # Generated build target
    # ======================================================

    st.divider()

    if polymer_name is None:
        st.info(
            "A polymer build target will appear once "
            "a sequence has been created."
        )

    else:
        st.markdown(
            "### Generated Build Target"
        )

        st.code(
            polymer_name
        )

    # ======================================================
    # Generated polymer SMILES
    # ======================================================

    if polymer_smiles:
        st.markdown(
            "### Generated Polymer SMILES"
        )

        st.code(
            polymer_smiles
        )

        st.caption(
            "This SMILES is generated from the ordered "
            "monomer sequence selected above."
        )

    # ======================================================
    # Build summary
    # ======================================================

    if sequence:
        st.markdown(
            "### Build Summary"
        )

        st.markdown(
            '<div class="card">',
            unsafe_allow_html=True,
        )

        st.write(
            f"**Polymer name:** `{polymer_name}`"
        )

        st.write(
            f"**Chain length:** `{len(sequence)}`"
        )

        st.write(
            "**Unique PHA units:** "
            f"`{len(set(sequence))}`"
        )

        st.write(
            "**Sequence type:** "
            f"`{polymer_mode}`"
        )

        st.write(
            "**Ordered sequence:** "
            f"`{' - '.join(sequence)}`"
        )

        st.markdown(
            "</div>",
            unsafe_allow_html=True,
        )
