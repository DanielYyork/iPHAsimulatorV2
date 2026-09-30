#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import json
import os
import subprocess

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

from src.iphasimulator.openmmscript_builder import OpenMMScriptBuilder


STRUCTURE_DATABASE = Path("structure_database")
RESIDUE_CODES_CSV = STRUCTURE_DATABASE / "residue_codes.csv"

MD_SCRIPT_DIR = Path("src/md_simulation_scripts")
MD_SCRIPT_DIR.mkdir(parents=True, exist_ok=True)

AMBERTOOLS_PYTHON = Path.home() / "miniconda3/envs/AmberTools23/bin/python"


st.set_page_config(
    page_title="iPHAsimulatorV2",
    page_icon="🧬",
    layout="wide",
)


st.markdown(
    """
<style>
.stApp {
    background: linear-gradient(135deg, #0f172a 0%, #111827 45%, #172554 100%);
    color: #e5e7eb;
}
.main-title {
    font-size: 3.2rem;
    font-weight: 900;
    background: linear-gradient(90deg, #67e8f9, #a78bfa, #f472b6);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    margin-bottom: -0.4rem;
}
.subtitle {
    color: #cbd5e1;
    font-size: 1.1rem;
    margin-bottom: 1.5rem;
}
.card {
    background: rgba(15, 23, 42, 0.78);
    border: 1px solid rgba(148, 163, 184, 0.28);
    border-radius: 22px;
    padding: 1.2rem;
    box-shadow: 0 12px 35px rgba(0,0,0,0.25);
}
.sequence-chip {
    display: inline-block;
    padding: 0.35rem 0.7rem;
    margin: 0.2rem;
    border-radius: 999px;
    background: linear-gradient(90deg, #0891b2, #7c3aed);
    color: white;
    font-weight: 700;
    font-size: 0.9rem;
}
.workflow-chip {
    display: inline-block;
    padding: 0.4rem 0.8rem;
    margin: 0.25rem;
    border-radius: 999px;
    background: linear-gradient(90deg, #0f766e, #2563eb);
    color: white;
    font-weight: 700;
    font-size: 0.9rem;
}
.small-muted {
    color: #94a3b8;
    font-size: 0.9rem;
}
.good-box {
    border-left: 5px solid #22c55e;
    background: rgba(34,197,94,0.12);
    padding: 0.8rem 1rem;
    border-radius: 12px;
}
.warn-box {
    border-left: 5px solid #f59e0b;
    background: rgba(245,158,11,0.12);
    padding: 0.8rem 1rem;
    border-radius: 12px;
}
.info-box {
    border-left: 5px solid #38bdf8;
    background: rgba(56,189,248,0.12);
    padding: 0.8rem 1rem;
    border-radius: 12px;
}
div.stButton > button {
    border-radius: 999px;
    font-weight: 700;
    border: 1px solid rgba(103,232,249,0.45);
    background: rgba(15,23,42,0.75);
    color: #e0f2fe;
}
div.stButton > button:hover {
    border-color: #f472b6;
    color: white;
    transform: scale(1.02);
}
</style>
""",
    unsafe_allow_html=True,
)


st.markdown(
    '<div class="main-title">🧬 iPHAsimulatorV2</div>',
    unsafe_allow_html=True,
)
st.markdown(
    '<div class="subtitle">Interactive PHA polymer construction, OpenMM workflow design, and script generation.</div>',
    unsafe_allow_html=True,
)


@st.cache_data
def load_available_PHA_monomers(residue_codes_csv):
    df = pd.read_csv(residue_codes_csv)
    mainchain_df = df[df["component"] == "mainchain"].copy()
    mainchain_df = mainchain_df.dropna(subset=["PHA_type", "smiles"])
    mainchain_df = mainchain_df[
        mainchain_df["smiles"].astype(str).str.strip() != ""
    ]

    available_phas = sorted(mainchain_df["PHA_type"].unique())
    monomer_smiles = {
        row["PHA_type"]: row["smiles"]
        for _, row in mainchain_df.iterrows()
    }

    return available_phas, monomer_smiles, mainchain_df


try:
    available_phas, monomer_smiles, mainchain_df = load_available_PHA_monomers(
        RESIDUE_CODES_CSV
    )
except Exception as error:
    st.error("Could not load available PHA monomers.")
    st.code(str(error))
    st.stop()


# ==========================================================
# Session state
# ==========================================================

if "sequence" not in st.session_state:
    st.session_state.sequence = []

if "preview_PHA" not in st.session_state:
    st.session_state.preview_PHA = None

if "openmm_steps" not in st.session_state:
    st.session_state.openmm_steps = []

if "generated_openmm_script" not in st.session_state:
    st.session_state.generated_openmm_script = None

if "generated_openmm_script_path" not in st.session_state:
    st.session_state.generated_openmm_script_path = None


# ==========================================================
# Polymer helper functions
# ==========================================================

def get_polymer_name(sequence):
    if not sequence:
        return None

    unique_units = []
    for unit in sequence:
        if unit not in unique_units:
            unique_units.append(unit)

    if len(unique_units) == 1:
        return f"P{unique_units[0]}_{len(sequence)}"

    return (
        "co_"
        + "_".join([f"P{unit}" for unit in unique_units])
        + f"_custom_{len(sequence)}"
    )


def generate_polymer_smiles_from_sequence(sequence, smiles_lookup):
    if not sequence:
        return None

    polymer_smiles = ""
    for i, pha_type in enumerate(sequence):
        smiles = smiles_lookup[pha_type]
        polymer_smiles += smiles if i == len(sequence) - 1 else smiles[:-1]

    return polymer_smiles


def draw_monomer(pha_type):
    mol = Chem.MolFromSmiles(monomer_smiles[pha_type])
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=(420, 300), legend=pha_type)


def draw_polymer(polymer_smiles, polymer_name):
    mol = Chem.MolFromSmiles(polymer_smiles)
    if mol is None:
        return None
    return Draw.MolToImage(mol, size=(1200, 360), legend=polymer_name)


def build_polymer_subprocess(sequence):
    build_code = r"""
import sys
import json
from src.iphasimulator.build_pha import PHAPolymerBuilder

root_dir = sys.argv[1]
sequence = json.loads(sys.argv[2])

builder = PHAPolymerBuilder(root_dir)

unique_units = []
for unit in sequence:
    if unit not in unique_units:
        unique_units.append(unit)

if len(unique_units) == 1:
    PHA_type = unique_units[0]
    length = len(sequence)
    print(f"Building homopolymer: P{PHA_type}_{length}")
    output = builder.build_PHA_polymer(PHA_type=PHA_type, length=length)

else:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if len(unique_units) > len(letters):
        raise ValueError("Too many unique monomer types.")

    unit_to_letter = {unit: letters[i] for i, unit in enumerate(unique_units)}
    pattern = "".join(unit_to_letter[unit] for unit in sequence)
    length = len(sequence)

    print("Building copolymer.")
    print("PHA types:", unique_units)
    print("Pattern:", pattern)
    print("Length:", length)

    output = builder.build_PHA_copolymer(
        PHA_types=unique_units,
        pattern=pattern,
        length=length,
    )

print("Build complete.")
print(output)
"""

    if not AMBERTOOLS_PYTHON.exists():
        raise FileNotFoundError(
            f"Could not find AmberTools23 Python:\n{AMBERTOOLS_PYTHON}"
        )

    ambertools_bin = AMBERTOOLS_PYTHON.parent
    env = dict(**os.environ)
    env["PATH"] = f"{ambertools_bin}:{env.get('PATH', '')}"

    command = [
        str(AMBERTOOLS_PYTHON),
        "-c",
        build_code,
        str(STRUCTURE_DATABASE),
        json.dumps(sequence),
    ]

    return subprocess.run(
        command,
        cwd=Path.cwd(),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


# ==========================================================
# OpenMM workflow helper functions
# ==========================================================

def add_workflow_step(step):
    st.session_state.openmm_steps.append(step)


def move_workflow_step_up(index):
    if index <= 0:
        return
    steps = st.session_state.openmm_steps
    steps[index - 1], steps[index] = steps[index], steps[index - 1]


def move_workflow_step_down(index):
    steps = st.session_state.openmm_steps
    if index >= len(steps) - 1:
        return
    steps[index + 1], steps[index] = steps[index], steps[index + 1]


def remove_workflow_step(index):
    st.session_state.openmm_steps.pop(index)


def build_openmm_script_builder(
    polymer_names,
    number_of_polymers,
    run_name,
):
    builder = OpenMMScriptBuilder(
        polymer_names=polymer_names,
        number_of_polymers=number_of_polymers,
        run_name=run_name,
    )

    for step in st.session_state.openmm_steps:
        method = step["method"]

        if method == "minimize_energy":
            builder.add_minimization()

        elif method == "basic_NVT":
            builder.add_basic_NVT(
                total_steps=step["total_steps"],
                temp=step["temp"],
                filename=step["filename"],
                save_restart=step["save_restart"],
                restart_name=step["restart_name"],
            )

        elif method == "basic_NPT":
            builder.add_basic_NPT(
                total_steps=step["total_steps"],
                temp=step["temp"],
                pressure=step["pressure"],
                filename=step["filename"],
                save_restart=step["save_restart"],
                restart_name=step["restart_name"],
            )

        elif method == "anneal_NVT":
            builder.add_anneal_NVT(
                start_temp=step["start_temp"],
                max_temp=step["max_temp"],
                cycles=step["cycles"],
                quench_rate=step["quench_rate"],
                steps_per_cycle=step["steps_per_cycle"],
                filename=step["filename"],
                save_restart=step["save_restart"],
                restart_name=step["restart_name"],
            )

        elif method == "thermal_ramp":
            builder.add_thermal_ramp(
                heating=step["heating"],
                ensemble=step["ensemble"],
                start_temp=step["start_temp"],
                max_temp=step["max_temp"],
                quench_rate=step["quench_rate"],
                total_steps=step["total_steps"],
                pressure=step["pressure"],
                filename=step["filename"],
                save_restart=step["save_restart"],
                restart_name=step["restart_name"],
            )

    return builder

def get_next_md_script_path(run_name):

    safe_run_name = run_name.strip().replace(" ", "_")

    if safe_run_name == "":

        safe_run_name = "OpenMM_Run"

    counter = 1

    while True:

        script_path = MD_SCRIPT_DIR / f"{safe_run_name}_{counter:02d}.py"

        if not script_path.exists():

            return script_path

        counter += 1

def run_python_script_with_ambertools(script_path):
    if not AMBERTOOLS_PYTHON.exists():
        raise FileNotFoundError(
            f"Could not find AmberTools23 Python:\n{AMBERTOOLS_PYTHON}"
        )

    ambertools_bin = AMBERTOOLS_PYTHON.parent
    env = dict(**os.environ)
    env["PATH"] = f"{ambertools_bin}:{env.get('PATH', '')}"

    command = [
        str(AMBERTOOLS_PYTHON),
        str(script_path),
    ]

    return subprocess.run(
        command,
        cwd=Path.cwd(),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def workflow_step_label(step):
    method = step["method"]

    if method == "minimize_energy":
        return "Minimization"

    if method == "basic_NVT":
        return f"Basic NVT · {step['temp']} K · {step['total_steps']} steps"

    if method == "basic_NPT":
        return (
            f"Basic NPT · {step['temp']} K · "
            f"{step['pressure']} atm · {step['total_steps']} steps"
        )

    if method == "anneal_NVT":
        return (
            f"Anneal NVT · {step['start_temp']}→{step['max_temp']} K · "
            f"{step['cycles']} cycles"
        )

    if method == "thermal_ramp":
        direction = "Heat" if step["heating"] else "Cool"
        return (
            f"{direction} ramp · {step['ensemble']} · "
            f"{step['start_temp']}→{step['max_temp']} K · "
            f"{step['quench_rate']} K"
        )

    return method


polymer_name = get_polymer_name(st.session_state.sequence)
polymer_smiles = generate_polymer_smiles_from_sequence(
    st.session_state.sequence,
    monomer_smiles,
)


# ==========================================================
# Sidebar
# ==========================================================

with st.sidebar:
    st.markdown("## ⚙️ Control Centre")
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.write("Structure database")
    st.code(str(STRUCTURE_DATABASE))
    st.write("AmberTools Python")
    st.code(str(AMBERTOOLS_PYTHON))
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("## 📊 Database")
    st.metric("Available monomers", len(available_phas))
    st.metric("Current length", len(st.session_state.sequence))
    st.metric("Unique units", len(set(st.session_state.sequence)))
    st.metric("OpenMM steps", len(st.session_state.openmm_steps))

    with st.expander("Show monomer table"):
        monomer_table_html = mainchain_df.to_html(index=False, escape=True)
        st.markdown(
            f'<div style="overflow-x: auto;">{monomer_table_html}</div>',
            unsafe_allow_html=True,
        )


tab_builder, tab_preview, tab_logs, tab_openmm = st.tabs(
    [
        "🧱 Polymer Builder",
        "🔬 Molecular Preview",
        "🖥 Build Console",
        "⚛️ OpenMM Script Builder",
    ]
)


# ==========================================================
# Polymer Builder tab
# ==========================================================

with tab_builder:
    top_metrics = st.columns(4)
    top_metrics[0].metric("Sequence length", len(st.session_state.sequence))
    top_metrics[1].metric("Unique monomers", len(set(st.session_state.sequence)))
    top_metrics[2].metric(
        "Mode",
        "Homo" if len(set(st.session_state.sequence)) <= 1 else "Co",
    )
    top_metrics[3].metric("Build env", "AmberTools23")

    st.markdown("### Available PHA Monomer Units")

    search = st.text_input(
        "Filter monomers",
        placeholder="Try 3HB, 4HB, phenyl, fluor...",
    )

    filtered_phas = [
        pha for pha in available_phas
        if search.lower() in pha.lower()
    ]

    cols = st.columns(6)
    for i, pha in enumerate(filtered_phas):
        with cols[i % 6]:
            if st.button(f"➕ {pha}", use_container_width=True):
                st.session_state.sequence.append(pha)
                st.session_state.preview_PHA = pha
                st.rerun()

    st.markdown("### Current Sequence")

    if st.session_state.sequence:
        chip_html = "".join(
            [
                f'<span class="sequence-chip">{unit}</span>'
                for unit in st.session_state.sequence
            ]
        )
        st.markdown(chip_html, unsafe_allow_html=True)
    else:
        st.markdown(
            '<div class="warn-box">Click monomer buttons above to start building a polymer sequence.</div>',
            unsafe_allow_html=True,
        )

    control_cols = st.columns([1, 1, 1, 2])

    with control_cols[0]:
        if st.button("↩️ Remove last", use_container_width=True):
            if st.session_state.sequence:
                st.session_state.sequence.pop()
                st.rerun()

    with control_cols[1]:
        if st.button("🧹 Clear", use_container_width=True):
            st.session_state.sequence = []
            st.session_state.preview_PHA = None
            st.rerun()

    with control_cols[2]:
        duplicate_n = st.number_input(
            "Repeat current sequence",
            min_value=1,
            max_value=20,
            value=1,
        )

    with control_cols[3]:
        if st.button("🔁 Apply repeat", use_container_width=True):
            if st.session_state.sequence:
                original = list(st.session_state.sequence)
                st.session_state.sequence = original * int(duplicate_n)
                st.rerun()

    st.divider()

    if polymer_name:
        st.markdown("### Generated Build Target")
        st.code(polymer_name)

    if polymer_smiles:
        st.markdown("### Generated Polymer SMILES")
        st.code(polymer_smiles)


# ==========================================================
# Preview tab
# ==========================================================

with tab_preview:
    col_a, col_b = st.columns([1, 2])

    with col_a:
        st.markdown("### Selected Monomer")

        if st.session_state.preview_PHA is None:
            st.info("Click a monomer to preview it.")
        else:
            pha = st.session_state.preview_PHA
            img = draw_monomer(pha)

            if img is None:
                st.error(f"Could not draw {pha}")
            else:
                st.image(img)
                st.markdown(f"### `{pha}`")
                st.code(monomer_smiles[pha])

    with col_b:
        st.markdown("### Polymer Preview")

        if polymer_smiles is None:
            st.info("No polymer sequence yet.")
        else:
            img = draw_polymer(polymer_smiles, polymer_name)

            if img is None:
                st.warning("RDKit could not render this polymer SMILES.")
            else:
                st.image(img, use_container_width=True)

            st.markdown("### Polymer Summary")
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.write(f"**Name:** `{polymer_name}`")
            st.write(f"**Length:** `{len(st.session_state.sequence)}`")
            st.write(
                f"**Unique monomers:** `{', '.join(sorted(set(st.session_state.sequence)))}`"
            )
            st.write(f"**SMILES characters:** `{len(polymer_smiles)}`")
            st.markdown("</div>", unsafe_allow_html=True)


# ==========================================================
# Build Console tab
# ==========================================================

with tab_logs:
    st.markdown("### Build Polymer")

    st.markdown(
        """
<div class="small-muted">
This launches the backend polymer builder in the AmberTools23 environment. 
The GUI environment only runs Streamlit and RDKit rendering.
</div>
""",
        unsafe_allow_html=True,
    )

    build_clicked = st.button("🚀 Build polymer now", use_container_width=True)

    if build_clicked:
        if not st.session_state.sequence:
            st.error("No monomers selected.")
        else:
            st.info("Launching polymer build in AmberTools23 environment...")

            progress = st.progress(0)
            status = st.empty()

            with st.spinner("Building polymer with PHAPolymerBuilder..."):
                try:
                    status.write("Preparing subprocess...")
                    progress.progress(15)

                    result = build_polymer_subprocess(st.session_state.sequence)

                    progress.progress(85)
                    status.write("Build subprocess finished.")

                except Exception as error:
                    progress.progress(100)
                    st.error("Could not launch build subprocess.")
                    st.code(str(error))
                    st.stop()

            progress.progress(100)

            if result.returncode == 0:
                st.success("Polymer build complete.")
                st.balloons()

                st.markdown(
                    '<div class="good-box">Backend returned success.</div>',
                    unsafe_allow_html=True,
                )

                if result.stdout:
                    st.subheader("STDOUT")
                    st.code(result.stdout)

            else:
                st.error("Polymer build failed.")

                if result.stdout:
                    st.subheader("STDOUT")
                    st.code(result.stdout)

                if result.stderr:
                    st.subheader("STDERR")
                    st.code(result.stderr)

    st.divider()

    st.markdown("### Current command context")
    st.write("AmberTools Python:")
    st.code(str(AMBERTOOLS_PYTHON))
    st.write("Working directory:")
    st.code(str(Path.cwd()))


# ==========================================================
# OpenMM Script Builder tab
# ==========================================================

with tab_openmm:
    st.markdown("## ⚛️ OpenMM Simulation Script Builder")

    st.markdown(
        """
<div class="info-box">
Build an ordered OpenMM workflow here. The GUI writes a normal Python script using 
<code>OpenMMScriptBuilder</code>, then you can inspect it or run it with AmberTools23.
</div>
""",
        unsafe_allow_html=True,
    )

    st.divider()

    settings_col, workflow_col = st.columns([1, 1.4])

    with settings_col:
        st.markdown("### System Selection")

        polymer_names_input = st.text_input(
            "Polymer names",
            value="P3HB_10",
            help="Comma-separated, e.g. P3HB_10 or P3HB_10,P4HB_10",
        )

        polymer_counts_input = st.text_input(
            "Number of polymers",
            value="25",
            help="Comma-separated counts matching polymer names, e.g. 25 or 25,25",
        )

        run_name_input = st.text_input(
            "Run name",
            value="Test",
            help="Used to create directories such as Test_01, Test_02, Tg_01.",
        )

        st.write("Generated scripts will be saved in:")
        st.code(str(MD_SCRIPT_DIR))

        st.markdown("### Add Workflow Step")

        step_type = st.selectbox(
            "Step type",
            [
                "minimize_energy",
                "basic_NVT",
                "basic_NPT",
                "anneal_NVT",
                "thermal_ramp",
            ],
        )

        if step_type == "minimize_energy":
            st.markdown(
                '<div class="small-muted">No parameters required.</div>',
                unsafe_allow_html=True,
            )

            if st.button("➕ Add minimization", use_container_width=True):
                add_workflow_step({"method": "minimize_energy"})
                st.rerun()

        elif step_type == "basic_NVT":
            total_steps = st.number_input(
                "Total steps",
                min_value=1,
                value=3000,
                step=1000,
                key="add_nvt_steps",
            )
            temp = st.number_input(
                "Temperature / K",
                min_value=0,
                value=300,
                step=10,
                key="add_nvt_temp",
            )
            filename = st.text_input(
                "Filename label",
                value="NVT",
                key="add_nvt_filename",
            )
            save_restart = st.checkbox(
                "Save restart",
                value=False,
                key="add_nvt_restart",
            )
            restart_name = st.text_input(
                "Restart name",
                value="",
                key="add_nvt_restart_name",
            )

            if st.button("➕ Add Basic NVT", use_container_width=True):
                add_workflow_step(
                    {
                        "method": "basic_NVT",
                        "total_steps": int(total_steps),
                        "temp": float(temp),
                        "filename": filename,
                        "save_restart": bool(save_restart),
                        "restart_name": restart_name or None,
                    }
                )
                st.rerun()

        elif step_type == "basic_NPT":
            total_steps = st.number_input(
                "Total steps",
                min_value=1,
                value=3000,
                step=1000,
                key="add_npt_steps",
            )
            temp = st.number_input(
                "Temperature / K",
                min_value=0,
                value=300,
                step=10,
                key="add_npt_temp",
            )
            pressure = st.number_input(
                "Pressure / atm",
                min_value=0.0,
                value=1.0,
                step=0.1,
                key="add_npt_pressure",
            )
            filename = st.text_input(
                "Filename label",
                value="NPT",
                key="add_npt_filename",
            )
            save_restart = st.checkbox(
                "Save restart",
                value=False,
                key="add_npt_restart",
            )
            restart_name = st.text_input(
                "Restart name",
                value="",
                key="add_npt_restart_name",
            )

            if st.button("➕ Add Basic NPT", use_container_width=True):
                add_workflow_step(
                    {
                        "method": "basic_NPT",
                        "total_steps": int(total_steps),
                        "temp": float(temp),
                        "pressure": float(pressure),
                        "filename": filename,
                        "save_restart": bool(save_restart),
                        "restart_name": restart_name or None,
                    }
                )
                st.rerun()

        elif step_type == "anneal_NVT":
            start_temp = st.number_input(
                "Start temperature / K",
                min_value=0,
                value=300,
                step=10,
                key="add_anneal_start",
            )
            max_temp = st.number_input(
                "Max temperature / K",
                min_value=0,
                value=700,
                step=10,
                key="add_anneal_max",
            )
            cycles = st.number_input(
                "Cycles",
                min_value=1,
                value=5,
                step=1,
                key="add_anneal_cycles",
            )
            quench_rate = st.number_input(
                "Temperature increment / K",
                min_value=1,
                value=10,
                step=1,
                key="add_anneal_quench",
            )
            steps_per_cycle = st.number_input(
                "Steps per cycle",
                min_value=1,
                value=500000,
                step=10000,
                key="add_anneal_steps",
            )
            filename = st.text_input(
                "Filename label",
                value="anneal_NVT",
                key="add_anneal_filename",
            )
            save_restart = st.checkbox(
                "Save restart",
                value=False,
                key="add_anneal_restart",
            )
            restart_name = st.text_input(
                "Restart name",
                value="",
                key="add_anneal_restart_name",
            )

            if st.button("➕ Add Anneal NVT", use_container_width=True):
                add_workflow_step(
                    {
                        "method": "anneal_NVT",
                        "start_temp": float(start_temp),
                        "max_temp": float(max_temp),
                        "cycles": int(cycles),
                        "quench_rate": float(quench_rate),
                        "steps_per_cycle": int(steps_per_cycle),
                        "filename": filename,
                        "save_restart": bool(save_restart),
                        "restart_name": restart_name or None,
                    }
                )
                st.rerun()

        elif step_type == "thermal_ramp":
            heating_choice = st.radio(
                "Ramp direction",
                ["Heating", "Cooling"],
                horizontal=True,
                key="add_ramp_direction",
            )
            heating = heating_choice == "Heating"

            ensemble = st.selectbox(
                "Ensemble",
                ["NVT", "NPT"],
                key="add_ramp_ensemble",
            )
            start_temp = st.number_input(
                "Start temperature / K",
                min_value=0,
                value=300 if heating else 700,
                step=10,
                key="add_ramp_start",
            )
            max_temp = st.number_input(
                "Target temperature / K",
                min_value=0,
                value=700 if heating else 140,
                step=10,
                key="add_ramp_max",
            )
            quench_rate = st.number_input(
                "Temperature increment / K",
                min_value=1,
                value=10,
                step=1,
                key="add_ramp_quench",
            )
            total_steps = st.number_input(
                "Total steps",
                min_value=1,
                value=100000,
                step=10000,
                key="add_ramp_steps",
            )
            pressure = st.number_input(
                "Pressure / atm",
                min_value=0.0,
                value=1.0,
                step=0.1,
                key="add_ramp_pressure",
            )
            filename = st.text_input(
                "Filename label",
                value="thermal_ramp",
                key="add_ramp_filename",
            )
            save_restart = st.checkbox(
                "Save restart",
                value=False,
                key="add_ramp_restart",
            )
            restart_name = st.text_input(
                "Restart name",
                value="",
                key="add_ramp_restart_name",
            )

            if st.button("➕ Add Thermal Ramp", use_container_width=True):
                add_workflow_step(
                    {
                        "method": "thermal_ramp",
                        "heating": bool(heating),
                        "ensemble": ensemble,
                        "start_temp": float(start_temp),
                        "max_temp": float(max_temp),
                        "quench_rate": float(quench_rate),
                        "total_steps": int(total_steps),
                        "pressure": float(pressure),
                        "filename": filename,
                        "save_restart": bool(save_restart),
                        "restart_name": restart_name or None,
                    }
                )
                st.rerun()

    with workflow_col:
        st.markdown("### Current OpenMM Workflow")

        if st.session_state.openmm_steps:
            chip_html = "".join(
                [
                    f'<span class="workflow-chip">{i + 1}. {workflow_step_label(step)}</span>'
                    for i, step in enumerate(st.session_state.openmm_steps)
                ]
            )
            st.markdown(chip_html, unsafe_allow_html=True)
        else:
            st.markdown(
                '<div class="warn-box">No OpenMM workflow steps added yet. Start with minimization.</div>',
                unsafe_allow_html=True,
            )

        st.divider()

        for index, step in enumerate(st.session_state.openmm_steps):
            with st.expander(
                f"Step {index + 1}: {workflow_step_label(step)}",
                expanded=True,
            ):
                st.json(step)

                button_cols = st.columns(4)

                with button_cols[0]:
                    if st.button(
                        "⬆️ Up",
                        key=f"openmm_up_{index}",
                        use_container_width=True,
                    ):
                        move_workflow_step_up(index)
                        st.rerun()

                with button_cols[1]:
                    if st.button(
                        "⬇️ Down",
                        key=f"openmm_down_{index}",
                        use_container_width=True,
                    ):
                        move_workflow_step_down(index)
                        st.rerun()

                with button_cols[2]:
                    if st.button(
                        "🗑 Remove",
                        key=f"openmm_remove_{index}",
                        use_container_width=True,
                    ):
                        remove_workflow_step(index)
                        st.rerun()

                with button_cols[3]:
                    if st.button(
                        "📋 Duplicate",
                        key=f"openmm_duplicate_{index}",
                        use_container_width=True,
                    ):
                        st.session_state.openmm_steps.insert(
                            index + 1,
                            dict(step),
                        )
                        st.rerun()

        st.divider()

        action_cols = st.columns(3)

        with action_cols[0]:
            if st.button("🧹 Clear workflow", use_container_width=True):
                st.session_state.openmm_steps = []
                st.session_state.generated_openmm_script = None
                st.session_state.generated_openmm_script_path = None
                st.rerun()

        with action_cols[1]:
            generate_clicked = st.button(
                "📝 Generate script",
                use_container_width=True,
            )

        with action_cols[2]:
            run_clicked = st.button(
                "🚀 Run generated script",
                use_container_width=True,
            )

        if generate_clicked:
            try:
                polymer_names = [
                    item.strip()
                    for item in polymer_names_input.split(",")
                    if item.strip()
                ]

                number_of_polymers = [
                    int(item.strip())
                    for item in polymer_counts_input.split(",")
                    if item.strip()
                ]

                if len(polymer_names) != len(number_of_polymers):
                    raise ValueError(
                        "Polymer names and polymer counts must have the same length."
                    )

                if not st.session_state.openmm_steps:
                    raise ValueError("No OpenMM workflow steps have been added.")

                if st.session_state.openmm_steps[0]["method"] != "minimize_energy":
                    raise ValueError(
                        "For now, the first OpenMM workflow step must be minimization."
                    )

                script_builder = build_openmm_script_builder(
                    polymer_names=polymer_names,
                    number_of_polymers=number_of_polymers,
                    run_name=run_name_input,
                )

                output_script = get_next_md_script_path(run_name_input)

                script_text = script_builder.to_script()
                output_script = script_builder.write_script(output_script)

                st.session_state.generated_openmm_script = script_text
                st.session_state.generated_openmm_script_path = str(output_script)

                st.success("Generated OpenMM simulation script.")
                st.code(str(output_script))

            except Exception as error:
                st.error("Could not generate OpenMM script.")
                st.code(str(error))

        if st.session_state.generated_openmm_script is not None:
            st.markdown("### Generated Script Preview")
            st.code(
                st.session_state.generated_openmm_script,
                language="python",
            )

        if run_clicked:
            if st.session_state.generated_openmm_script_path is None:
                st.error("Generate a script before running it.")
            else:
                script_path = Path(st.session_state.generated_openmm_script_path)

                if not script_path.exists():
                    st.error(f"Generated script does not exist:\n{script_path}")
                else:
                    st.info("Running generated OpenMM script with AmberTools23...")

                    progress = st.progress(0)
                    status = st.empty()

                    with st.spinner("Running OpenMM script..."):
                        try:
                            status.write("Launching generated script...")
                            progress.progress(20)

                            result = run_python_script_with_ambertools(script_path)

                            progress.progress(90)
                            status.write("Script finished.")

                        except Exception as error:
                            progress.progress(100)
                            st.error("Could not run generated script.")
                            st.code(str(error))
                            st.stop()

                    progress.progress(100)

                    if result.returncode == 0:
                        st.success("Generated OpenMM script finished successfully.")
                        st.balloons()

                        if result.stdout:
                            st.subheader("STDOUT")
                            st.code(result.stdout)

                    else:
                        st.error("Generated OpenMM script failed.")

                        if result.stdout:
                            st.subheader("STDOUT")
                            st.code(result.stdout)

                        if result.stderr:
                            st.subheader("STDERR")
                            st.code(result.stderr)
