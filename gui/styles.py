#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visual styling for the iPHAsimulatorV2 Streamlit GUI.
"""

import streamlit as st
from pathlib import Path


APP_CSS = """
<style>
.stApp {
    background: linear-gradient(
        135deg,
        #0f172a 0%,
        #111827 45%,
        #172554 100%
    );
    color: #e5e7eb;
}

.main-title {
    font-size: 3.2rem;
    font-weight: 900;
    background: linear-gradient(
        90deg,
        #67e8f9,
        #a78bfa,
        #f472b6
    );
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
    box-shadow: 0 12px 35px rgba(0, 0, 0, 0.25);
}

.sequence-chip {
    display: inline-block;
    padding: 0.35rem 0.7rem;
    margin: 0.2rem;
    border-radius: 999px;
    background: linear-gradient(
        90deg,
        #0891b2,
        #7c3aed
    );
    color: white;
    font-weight: 700;
    font-size: 0.9rem;
}

.workflow-chip {
    display: inline-block;
    padding: 0.4rem 0.8rem;
    margin: 0.25rem;
    border-radius: 999px;
    background: linear-gradient(
        90deg,
        #0f766e,
        #2563eb
    );
    color: white;
    font-weight: 700;
    font-size: 0.9rem;
}

.system-chip {
    display: inline-block;
    padding: 0.4rem 0.8rem;
    margin: 0.25rem;
    border-radius: 999px;
    background: linear-gradient(
        90deg,
        #7c3aed,
        #db2777
    );
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
    background: rgba(34, 197, 94, 0.12);
    padding: 0.8rem 1rem;
    border-radius: 12px;
}

.warn-box {
    border-left: 5px solid #f59e0b;
    background: rgba(245, 158, 11, 0.12);
    padding: 0.8rem 1rem;
    border-radius: 12px;
}

.info-box {
    border-left: 5px solid #38bdf8;
    background: rgba(56, 189, 248, 0.12);
    padding: 0.8rem 1rem;
    border-radius: 12px;
}

.error-box {
    border-left: 5px solid #ef4444;
    background: rgba(239, 68, 68, 0.12);
    padding: 0.8rem 1rem;
    border-radius: 12px;
}

div.stButton > button {
    border-radius: 999px;
    font-weight: 700;
    border: 1px solid rgba(103, 232, 249, 0.45);
    background: rgba(15, 23, 42, 0.75);
    color: #e0f2fe;
}

div.stButton > button:hover {
    border-color: #f472b6;
    color: white;
    transform: scale(1.02);
}

/* Compact monomer actions; other application buttons keep their usual size. */
[class*="st-key-add_monomer_"] button {
    min-height: 1.75rem;
    padding: 0.125rem 0.55rem;
}

[class*="st-key-add_monomer_"] button p {
    font-size: 0.8rem;
}

.monomer-status {
    display: inline-block;
    padding: 0.1rem 0.45rem;
    border-radius: 999px;
    font-size: 0.72rem;
    font-weight: 600;
}
.monomer-ready {
    background: rgba(34, 197, 94, 0.15);
    color: #86efac;
}
.monomer-incomplete {
    background: rgba(245, 158, 11, 0.15);
    color: #fcd34d;
}

.monomer-hover-container {
    position: relative;
    display: block;
    width: 100%;
    margin-bottom: 0.25rem;
}

.monomer-hover-name {
    width: 100%;
    box-sizing: border-box;
    padding: 0.55rem 0.75rem;
    border-radius: 14px;
    border: 1px solid rgba(103, 232, 249, 0.38);
    background: rgba(15, 23, 42, 0.72);
    color: #e0f2fe;
    text-align: center;
    font-weight: 750;
    cursor: help;
    transition:
        border-color 0.15s ease,
        background 0.15s ease,
        transform 0.15s ease;
}

.monomer-hover-name:hover {
    border-color: #f472b6;
    background: rgba(30, 41, 59, 0.95);
    transform: translateY(-1px);
}

.monomer-hover-tooltip {
    visibility: hidden;
    opacity: 0;
    position: absolute;
    z-index: 9999;
    left: 50%;
    bottom: calc(100% + 10px);
    transform: translateX(-50%) translateY(5px);

    width: 340px;
    max-width: 80vw;
    padding: 0.8rem;

    background: rgba(15, 23, 42, 0.98);
    border: 1px solid rgba(167, 139, 250, 0.7);
    border-radius: 18px;
    box-shadow: 0 16px 45px rgba(0, 0, 0, 0.5);

    pointer-events: none;

    transition:
        opacity 0.15s ease,
        transform 0.15s ease,
        visibility 0.15s ease;
}

.monomer-hover-container:hover .monomer-hover-tooltip {
    visibility: visible;
    opacity: 1;
    transform: translateX(-50%) translateY(0);
}

.monomer-tooltip-title {
    color: #f8fafc;
    font-size: 1rem;
    font-weight: 800;
    text-align: center;
    margin-bottom: 0.4rem;
}

.monomer-tooltip-image {
    display: block;
    width: 100%;
    height: auto;
    border-radius: 12px;
    background: white;
}

.monomer-tooltip-smiles {
    margin-top: 0.55rem;
    padding: 0.4rem;
    border-radius: 8px;
    background: rgba(30, 41, 59, 0.9);
    color: #cbd5e1;
    font-family: monospace;
    font-size: 0.7rem;
    line-height: 1.2;
    overflow-wrap: anywhere;
}

.monomer-tooltip-text {
    color: #cbd5e1;
    text-align: center;
    padding: 1rem;
}
.workflow-stage {
    min-height: 105px;
    padding: 0.8rem 0.45rem;
    border: 1px solid rgba(103, 232, 249, 0.28);
    border-radius: 16px;
    background: rgba(15, 23, 42, 0.72);
    text-align: center;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    gap: 0.45rem;
}

.workflow-stage-number {
    width: 34px;
    height: 34px;
    border-radius: 50%;
    background: linear-gradient(
        135deg,
        #0891b2,
        #7c3aed
    );
    color: white;
    font-weight: 800;
    display: flex;
    align-items: center;
    justify-content: center;
}

.workflow-stage-label {
    color: #e2e8f0;
    font-size: 0.82rem;
    font-weight: 700;
    line-height: 1.2;
}

</style>
"""


def configure_page():
    """
    Configure the Streamlit page.
    """

    st.set_page_config(
        page_title="iPHAsimulatorV2",
        page_icon="🧬",
        layout="wide",
    )


def apply_styles():
    """
    Apply the custom application CSS.
    """

    st.markdown(
        APP_CSS,
        unsafe_allow_html=True,
    )


from pathlib import Path

import streamlit as st


def render_header():
    """
    Render the application title and subtitle.
    """

    st.markdown(
        '<div class="main-title">🧬 iPHAsimulatorV2</div>',
        unsafe_allow_html=True,
    )

    st.markdown(
        """
<div class="subtitle">
Interactive PHA polymer construction, registered MD-system selection,
molecular visualisation, OpenMM workflow design, and script generation.
</div>
""",
        unsafe_allow_html=True,
    )


def render_info_box(message):
    st.markdown(
        f'<div class="info-box">{message}</div>',
        unsafe_allow_html=True,
    )


def render_success_box(message):
    st.markdown(
        f'<div class="good-box">{message}</div>',
        unsafe_allow_html=True,
    )


def render_warning_box(message):
    st.markdown(
        f'<div class="warn-box">{message}</div>',
        unsafe_allow_html=True,
    )


def render_error_box(message):
    st.markdown(
        f'<div class="error-box">{message}</div>',
        unsafe_allow_html=True,
    )
