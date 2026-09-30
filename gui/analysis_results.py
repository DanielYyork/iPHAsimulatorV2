"""Read saved replica results without importing or running scientific workflows."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SavedAnalysis:
    name: str
    directory: Path
    summary: dict

    def section(self, name):
        value = self.summary.get(name, {})
        return value if isinstance(value, dict) else {}

    def output_path(self, key, default):
        """Resolve recorded relative outputs inside this replica's analysis folder."""
        value = self.section("outputs").get(key, default)
        if not isinstance(value, str):
            raise ValueError(f"Invalid output filename for {key}.")
        path = (self.directory / value).resolve()
        if not path.is_relative_to(self.directory.resolve()):
            raise ValueError("Saved output points outside this analysis folder.")
        return path


def discover_saved_analyses(simulations_directory, workflow):
    """Return valid summaries and readable diagnostics for incomplete results."""
    root = Path(simulations_directory)
    results, problems = [], []
    if not root.is_dir():
        return results, problems
    for replica in sorted(root.iterdir()):
        if not replica.is_dir():
            continue
        directory = replica / "analysis" / workflow
        summary_path = directory / "analysis_summary.json"
        if not summary_path.is_file():
            continue
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            if not isinstance(summary, dict) or not summary:
                raise ValueError("Expected a non-empty summary object.")
            results.append(SavedAnalysis(replica.name, directory, summary))
        except (OSError, ValueError) as error:
            problems.append(f"{replica.name}: could not read its saved summary ({error}).")
    return results, problems


def finite_number(value):
    """Keep missing or invalid summary values missing, rather than treating them as zero."""
    if isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def comparison_table(results):
    return pd.DataFrame([
        {
            "Replica": result.name,
            "Tg (K)": finite_number(result.section("tg").get("tg_K")),
            "PCs": finite_number(result.section("pca").get("selected_components")),
            "min_samples": finite_number(result.section("dbscan").get("selected_min_samples")),
            "Median noise (%)": _noise_percent(result),
            "Stride": finite_number(result.section("sampling").get("stride")),
            "Chains": finite_number(result.section("sampling").get("polymer_chains")),
            "Frames / chain": finite_number(result.section("sampling").get("frames_per_chain")),
        }
        for result in results
    ])


def _noise_percent(result):
    value = finite_number(result.section("dbscan").get("median_noise_fraction"))
    return None if value is None else value * 100


def _numeric_columns(frame, columns):
    for column in columns:
        if column not in frame:
            raise ValueError(f"Missing column: {column}.")
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if frame.empty or not np.isfinite(frame[columns].to_numpy(dtype=float)).all():
        raise ValueError("The saved table contains empty or invalid numerical data.")
    return frame


def load_temperature_response(path):
    frame = pd.read_csv(path)
    return _numeric_columns(frame, ["Temperature (K)", "Response"]).sort_values("Temperature (K)")


def saved_tanh_fit(result, temperatures):
    """Reconstruct only the model explicitly identified by the saved summary."""
    section = result.section("tg")
    if section.get("observable") != "historical_mean_dbscan_cluster_id":
        return None
    values = [finite_number(section.get(name)) for name in ("C", "s", "d")]
    if any(value is None for value in values):
        return None
    C, s, d = values
    if C <= 0 or s <= 0:
        return None
    return C / 2 * (1 - np.tanh(s * np.asarray(temperatures) - d)) - 1


def load_pca_spectrum(path):
    frame = pd.read_csv(path)
    _numeric_columns(frame, ["PC", "Explained Variance", "Cumulative Variance"])
    if ((frame["PC"] < 1) | (frame["PC"] % 1 != 0)).any():
        raise ValueError("PC indices must be positive integers.")
    variance = frame[["Explained Variance", "Cumulative Variance"]]
    if ((variance < 0) | (variance > 1 + 1e-8)).any().any():
        raise ValueError("Explained-variance ratios must lie between zero and one.")
    return frame.groupby("PC", as_index=False)[
        ["Explained Variance", "Cumulative Variance"]
    ].median().sort_values("PC")


def load_noise_heatmap(path):
    """Aggregate counts in chunks so large frame tables do not fill GUI memory."""
    columns = ["Chain", "Temperature (K)", "DBSCAN Label"]
    totals = None
    for chunk in pd.read_csv(path, usecols=columns, chunksize=100_000):
        if chunk.empty:
            continue
        _numeric_columns(chunk, ["Temperature (K)", "DBSCAN Label"])
        if chunk["Chain"].isna().any():
            raise ValueError("Some observations have no chain identifier.")
        labels = chunk["DBSCAN Label"]
        if ((labels < -1) | (labels % 1 != 0)).any():
            raise ValueError("DBSCAN labels must be -1 or non-negative integers.")
        chunk["Noise count"] = labels.eq(-1).astype(int)
        counts = chunk.groupby(["Chain", "Temperature (K)"]).agg(
            **{"Noise count": ("Noise count", "sum"), "Frames": ("DBSCAN Label", "size")}
        )
        totals = counts if totals is None else totals.add(counts, fill_value=0)
    if totals is None or totals.empty:
        raise ValueError("No saved conformational observations were found.")
    totals["Noise fraction"] = totals["Noise count"] / totals["Frames"]
    return totals.reset_index().sort_values(["Chain", "Temperature (K)"])
