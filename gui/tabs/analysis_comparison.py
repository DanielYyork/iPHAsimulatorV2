"""Compare existing replica outputs using lightweight GUI charts and saved images."""

import json
import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from gui.analysis_results import (
    comparison_table,
    discover_saved_analyses,
    finite_number,
    load_noise_heatmap,
    load_pca_spectrum,
    load_temperature_response,
    saved_tanh_fit,
)


@st.cache_data(show_spinner=False, max_entries=48)
def _load_table(path, modified_ns, size, kind):
    """Invalidate cached results when a saved file changes or is replaced."""
    del modified_ns, size
    readers = {
        "response": load_temperature_response,
        "spectrum": load_pca_spectrum,
        "noise": load_noise_heatmap,
    }
    return readers[kind](path)


def _read_table(result, kind):
    specifications = {
        "response": ("temperature_response_file", "tg/temperature_response.csv"),
        "spectrum": ("pca_spectrum_file", "pca/pca_spectrum_by_chain.csv"),
        "noise": ("conformational_state_file", "conformational_state_by_chain.csv"),
    }
    try:
        path = result.output_path(*specifications[kind])
        stat = path.stat()
        return _load_table(str(path), stat.st_mtime_ns, stat.st_size, kind)
    except (OSError, ValueError, KeyError, pd.errors.ParserError) as error:
        st.info(f"{result.name}: {kind} data unavailable. {error}")
        return None


def _replica_colour(names):
    # Distinct colours remain legible against the GUI's dark background.
    palette = ["#38bdf8", "#fb923c", "#c084fc", "#4ade80"]
    return alt.Color(
        "Replica:N", scale=alt.Scale(domain=names, range=palette[:len(names)]),
        legend=alt.Legend(title="Replica"),
    )


def _render_tg(results):
    st.markdown("#### Temperature response and saved Tg fits")
    observations, fitted = [], []
    observables = {result.section("tg").get("observable", "unspecified") for result in results}
    if len(observables) > 1:
        st.warning("These replicas use different response definitions. Their responses are shown separately.")
    for result in results:
        frame = _read_table(result, "response")
        if frame is None:
            continue
        frame = frame[["Temperature (K)", "Response"]].assign(Replica=result.name)
        observations.append(frame)
        temperatures = np.linspace(frame["Temperature (K)"].min(), frame["Temperature (K)"].max(), 250)
        fit = saved_tanh_fit(result, temperatures)
        if fit is not None:
            fitted.append(pd.DataFrame({"Temperature (K)": temperatures, "Response": fit, "Replica": result.name}))
        else:
            st.caption(f"{result.name}: showing measured responses only; a supported saved fit is unavailable.")
    if not observations:
        return
    names = [result.name for result in results]
    observed = pd.concat(observations, ignore_index=True)
    fits = pd.concat(fitted, ignore_index=True) if fitted else None

    def plot(data, fit_data):
        points = alt.Chart(data).mark_point(size=38, filled=True).encode(
            x=alt.X("Temperature (K):Q", scale=alt.Scale(zero=False)),
            y=alt.Y("Response:Q", title="Saved temperature response", scale=alt.Scale(zero=False)),
            color=_replica_colour(names),
            tooltip=["Replica:N", alt.Tooltip("Temperature (K):Q", format=".1f"), alt.Tooltip("Response:Q", format=".3f")],
        )
        chart = points
        if fit_data is not None and not fit_data.empty:
            lines = alt.Chart(fit_data).mark_line().encode(
                x="Temperature (K):Q", y="Response:Q", color=_replica_colour(names),
            )
            chart = lines + points
        return chart.properties(height=350).interactive()

    if len(observables) == 1:
        st.altair_chart(plot(observed, fits), use_container_width=True)
    else:
        for result in results:
            data = observed[observed.Replica == result.name]
            if not data.empty:
                st.write(f"**{result.name}** — {result.section('tg').get('observable', 'unspecified response')}")
                fit_data = fits[fits.Replica == result.name] if fits is not None else None
                st.altair_chart(plot(data, fit_data), use_container_width=True)
    st.caption("Points are saved responses. Lines reconstruct saved fit parameters; no analysis is rerun.")
    if "historical_mean_dbscan_cluster_id" in observables:
        st.caption("This historical response averages numerical cluster IDs. Its value can change when clusters are renamed, so the fitted Tg is an operational estimate.")


def _render_pca(results, scope):
    view = st.selectbox(
        "PCA view", ["Variance per component", "Cumulative variance", "Chain projections"],
        key=f"{scope}_pca_view",
    )
    if view == "Chain projections":
        _render_chain_projections(results, scope)
        return
    frames = []
    for result in results:
        frame = _read_table(result, "spectrum")
        if frame is not None:
            frames.append(frame.assign(Replica=result.name))
    if not frames:
        return
    frame = pd.concat(frames, ignore_index=True)
    column = "Explained Variance" if view == "Variance per component" else "Cumulative Variance"
    frame["Variance (%)"] = 100 * frame[column]
    chart = alt.Chart(frame).mark_line(point=True).encode(
        x=alt.X("PC:Q", title="Principal component", axis=alt.Axis(tickMinStep=1)),
        y=alt.Y("Variance (%):Q", title=f"Median {view.lower()} (%)"),
        color=_replica_colour([result.name for result in results]),
        tooltip=["Replica:N", "PC:Q", alt.Tooltip("Variance (%):Q", format=".2f")],
    ).properties(height=350).interactive()
    st.altair_chart(chart, use_container_width=True)
    st.caption("Each curve is the median across chains within one replica. The selected PC counts are shown in the table above.")


def _render_chain_projections(results, scope):
    suffixes = {"Temperature": "_pc1_pc2_temperature.png", "DBSCAN clusters": "_pc1_pc2_dbscan.png"}
    colour = st.radio("Colour points by", list(suffixes), horizontal=True, key=f"{scope}_pca_colour")
    suffix = suffixes[colour]
    files = {}
    for result in results:
        # A run that disabled figures can leave older images behind. Do not show those.
        if result.summary.get("figures_generated") is False:
            files[result.name] = {}
        else:
            files[result.name] = {
                path.name[:-len(suffix)]: path
                for path in (result.directory / "figures" / "pca_by_chain").glob(f"*{suffix}")
                if path.is_file() and path.resolve().is_relative_to(result.directory.resolve())
            }
    chains = sorted({chain for mapping in files.values() for chain in mapping})
    if not chains:
        st.info("No saved chain projections are available for these replicas. The variance views can still use the saved PCA tables.")
        return
    chain = st.selectbox("Chain to display", chains, key=f"{scope}_chain")
    st.caption("Each chain and replica has its own PCA basis. Compare patterns, rather than matching coordinate positions or cluster colours across plots.")
    for start in range(0, len(results), 2):
        for column, result in zip(st.columns(2), results[start:start + 2]):
            with column:
                st.write(f"**{result.name}**")
                path = files[result.name].get(chain)
                if path:
                    st.image(str(path), use_container_width=True)
                else:
                    st.info("This chain projection was not saved for this replica.")


def _render_heatmaps(results):
    st.markdown("#### Noise fraction by chain and temperature")
    frames = {}
    with st.spinner("Reading saved chain labels..."):
        for result in results:
            frame = _read_table(result, "noise")
            if frame is not None:
                frames[result.name] = frame
    if not frames:
        return
    temperatures = sorted({value for frame in frames.values() for value in frame["Temperature (K)"]})
    chains = sorted({str(value) for frame in frames.values() for value in frame["Chain"]})
    ticks = temperatures[::max(1, len(temperatures) // 8)]
    if temperatures[-1] not in ticks:
        ticks.append(temperatures[-1])
    st.caption("All panels use the same temperature, chain and 0–1 colour scales. Blank cells have no saved observations. Noise is a DBSCAN classification, not a physical phase assignment.")
    for start in range(0, len(results), 2):
        for column, result in zip(st.columns(2), results[start:start + 2]):
            with column:
                st.write(f"**{result.name}**")
                frame = frames.get(result.name)
                if frame is None:
                    st.info("No heatmap data available.")
                    continue
                frame = frame.assign(Chain=frame.Chain.astype(str))
                # Pass the table directly so Streamlit includes its Arrow dataset.
                chart = alt.Chart(frame).mark_rect().encode(
                    x=alt.X("Temperature (K):O", sort=temperatures, scale=alt.Scale(domain=temperatures), axis=alt.Axis(values=ticks, labelAngle=-45)),
                    y=alt.Y("Chain:N", sort=chains, scale=alt.Scale(domain=chains)),
                    color=alt.Color("Noise fraction:Q", scale=alt.Scale(domain=[0, 1], scheme="viridis")),
                    tooltip=["Chain:N", "Temperature (K):Q", alt.Tooltip("Noise fraction:Q", format=".1%"), alt.Tooltip("Frames:Q", format=".0f")],
                ).properties(height=max(250, min(700, len(chains) * 16)))
                st.altair_chart(chart, use_container_width=True)


def render_analysis_comparison(gui_data, system_name, system_type, workflow):
    st.markdown("### Compare saved replicas")
    st.caption("Choose existing analyses to compare their results. This view does not launch simulations or rerun analysis.")
    if not system_name or not system_type or not workflow:
        st.info("Select a system and an analysis workflow first.")
        return
    try:
        directory = gui_data.paths.get_md_system_files(
            system_name=system_name, system_type=system_type,
        )["simulations_dir"]
        results, problems = discover_saved_analyses(directory, workflow)
    except (OSError, ValueError, KeyError) as error:
        st.warning(f"Could not find saved analyses: {error}")
        return
    for problem in problems:
        st.warning(problem)
    if not results:
        st.info("No saved analysis summaries were found for this system and workflow. Run an analysis in the single-replica view first.")
        return
    scope = f"comparison_{system_type}_{system_name}_{workflow}"
    by_name = {result.name: result for result in results}
    names = list(by_name)
    selection_key = f"{scope}_replicas"
    if selection_key in st.session_state:
        previous = st.session_state[selection_key]
        valid = [name for name in previous if name in by_name]
        if valid != previous:
            st.session_state[selection_key] = valid
    selected = st.multiselect(
        "Replicas to compare", names, default=names[:2], max_selections=4,
        key=selection_key, help="Select up to four replicas. Results are read from their individual analysis folders.",
    )
    if not selected:
        st.info("Choose at least one saved replica. Select two or more for a comparison.")
        return
    chosen = [by_name[name] for name in selected]
    table = comparison_table(chosen)
    st.dataframe(table, hide_index=True, use_container_width=True, column_config={
        "Tg (K)": st.column_config.NumberColumn(format="%.2f"),
        "Median noise (%)": st.column_config.NumberColumn(format="%.1f"),
    })
    st.download_button("Download comparison table", table.to_csv(index=False),
                       file_name=f"{system_name}_replica_comparison.csv", mime="text/csv", key=f"{scope}_download")
    with st.expander("Compare analysis settings"):
        for result in chosen:
            st.write(f"**{result.name}**")
            st.json({name: result.summary.get(name) for name in ["stage", "sampling", "temperature_assignment", "pca", "dbscan", "tg"]})
    protocols = {
        json.dumps({name: result.summary.get(name) for name in ["stage", "temperature_assignment"]}, sort_keys=True)
        for result in chosen
    }
    strides = {finite_number(result.section("sampling").get("stride")) for result in chosen}
    if len(protocols) > 1 or len(strides) > 1:
        st.info("These replicas have different recorded temperature protocols, stages or sampling strides. Review the settings when interpreting their differences.")
    view = st.radio("Comparison view", ["Tg fits", "PCA", "Heatmaps"], horizontal=True, key=f"{scope}_view")
    if view == "Tg fits":
        _render_tg(chosen)
    elif view == "PCA":
        _render_pca(chosen, scope)
    else:
        _render_heatmaps(chosen)
