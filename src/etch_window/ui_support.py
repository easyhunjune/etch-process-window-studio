from __future__ import annotations

from collections.abc import Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

from .data_access import RESPONSE_COLUMNS, load_doe_and_models  # noqa: F401
from .design import FACTOR_NAMES, FACTOR_SPECS
from .rsm import FullQuadraticRSM
from .window import ProcessSpecs

RESPONSE_LABELS = {
    "etch_rate": "Si 식각 속도 (µm/min)",
    "selectivity": "Si/SiO₂ 선택비",
    "nonuniformity": "식각 속도 비균일도 (%)",
    "anisotropy": "이방성",
    "bias_voltage": "DC bias voltage (V)",
}
PLOT_RESPONSE_LABELS = {
    "etch_rate": "Si etch rate (µm/min)",
    "selectivity": "Si/SiO₂ selectivity",
    "nonuniformity": "Etch-rate non-uniformity (%)",
}


def process_window_frame(result: Mapping[str, object]) -> pd.DataFrame:
    actual = np.asarray(result["actual"])
    data: dict[str, np.ndarray] = {
        f"{name}_{FACTOR_SPECS[name].unit}": actual[:, index]
        for index, name in enumerate(FACTOR_NAMES)
    }
    data.update(
        {
            "predicted_etch_rate_um_min": np.asarray(result["etch_rate"]).ravel(),
            "predicted_selectivity": np.asarray(result["selectivity"]).ravel(),
            "predicted_nonuniformity_pct": np.asarray(
                result["nonuniformity"]
            ).ravel(),
            "nominally_meets_specs": np.asarray(
                result["nominal_feasible"]
            ).ravel(),
            "confidence_state": np.asarray(
                result["confidence_state"]
            ).ravel(),
            # app.py downloads only meets_all_specs rows, which are supported,
            # confirmed points, so this column is always "confirmed" in that
            # subset. Other values become relevant if that filter changes.
            "supported_confidence_state": np.asarray(
                result["supported_confidence_state"]
            ).ravel(),
            "inside_design_hull": np.asarray(
                result["inside_design_hull"]
            ).ravel(),
            "inside_leverage_envelope": np.asarray(
                result["inside_leverage_envelope"]
            ).ravel(),
            "inside_design_region": np.asarray(
                result["inside_design_region"]
            ).ravel(),
            "prediction_leverage": np.asarray(
                result["prediction_leverage"]
            ).ravel(),
            "meets_all_specs": np.asarray(result["feasible"]).ravel(),
        }
    )
    return pd.DataFrame(data)


def process_window_figure(
    result: Mapping[str, object],
    specs: ProcessSpecs,
) -> Figure:
    first_name, second_name = result["vary"]
    first = FACTOR_SPECS[first_name]
    second = FACTOR_SPECS[second_name]
    x = np.asarray(result["first_grid"])
    y = np.asarray(result["second_grid"])
    feasible = np.asarray(result["feasible"], dtype=bool)
    inside = np.asarray(result["inside_design_region"], dtype=bool)
    state = np.asarray(result["confidence_state"])
    panels = (
        ("etch_rate", "Si etch rate", "µm/min", specs.min_etch_rate_um_min, "≥"),
        ("selectivity", "Si/SiO₂ selectivity", "ratio", specs.min_selectivity, "≥"),
        (
            "nonuniformity",
            "Non-uniformity",
            "%",
            specs.max_nonuniformity_pct,
            "≤",
        ),
    )
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for axis, (key, title, unit, threshold, operator) in zip(
        axes.flat[:3], panels, strict=True
    ):
        values = np.asarray(result[key])
        contour = axis.contourf(x, y, values, levels=18, cmap="viridis")
        # Each panel draws the boundary of its own spec. Drawing nominal_feasible here
        # would put the same three-condition line on all three panels, so the white line
        # would not mean what the panel's colour scale says it means.
        if values.min() < threshold < values.max():
            axis.contour(
                x,
                y,
                values,
                levels=[threshold],
                colors="white",
                linewidths=2,
            )
        if (~inside).any():
            axis.contourf(
                x,
                y,
                (~inside).astype(float),
                levels=[0.5, 1.5],
                colors="none",
                hatches=["////"],
                alpha=0,
            )
        figure.colorbar(contour, ax=axis, label=unit)
        axis.set_title(f"{title} · spec {operator} {threshold:g}")

    mask_axis = axes.flat[3]
    state_code = np.zeros(state.shape, dtype=int)
    state_code[state == "uncertain"] = 1
    state_code[state == "confirmed"] = 2
    state_code[~inside] = -1
    mask_axis.pcolormesh(
        x,
        y,
        state_code,
        shading="auto",
        cmap=ListedColormap(["#f1f3f5", "#d9dee5", "#f5b942", "#00a67d"]),
        vmin=-1,
        vmax=2,
    )
    if (~inside).any():
        mask_axis.contourf(
            x,
            y,
            (~inside).astype(float),
            levels=[0.5, 1.5],
            colors="none",
            hatches=["////"],
            alpha=0,
        )
    mask_axis.set_title(
        "Joint prediction state "
        f"(confirmed in design region: {100 * feasible.mean():.1f}%)"
    )
    for axis in axes.flat:
        axis.set_xlabel(f"{first.display_name} ({first.unit})")
        axis.set_ylabel(f"{second.display_name} ({second.unit})")
    return figure


def factor_contribution_figure(factors: pd.DataFrame) -> Figure:
    responses = ("etch_rate", "selectivity", "nonuniformity")
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    for axis, response in zip(axes, responses, strict=True):
        subset = factors[factors["response"] == response]
        axis.bar(
            subset["factor"],
            subset["normalized_partial_ss_pct"],
            color="#3478c0",
        )
        axis.set_title(PLOT_RESPONSE_LABELS[response])
        axis.set_ylabel("Group-deletion partial SS share (%)")
        axis.tick_params(axis="x", rotation=45)
        axis.grid(axis="y", alpha=0.2)
    return figure


def main_effect_figure(
    models: Mapping[str, FullQuadraticRSM],
) -> Figure:
    responses = ("etch_rate", "selectivity", "nonuniformity")
    coded_axis = np.linspace(-2.0, 2.0, 201)
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for axis, response in zip(axes, responses, strict=True):
        for factor_index, factor_name in enumerate(FACTOR_NAMES):
            center_line = np.zeros((coded_axis.size, len(FACTOR_NAMES)))
            center_line[:, factor_index] = coded_axis
            axis.plot(
                coded_axis,
                models[response].predict(center_line),
                label=factor_name,
                linewidth=2,
            )
        axis.axvline(0, color="#777777", linestyle="--", linewidth=0.8)
        axis.set_title(PLOT_RESPONSE_LABELS[response])
        axis.set_xlabel("Coded level (other factors = 0)")
        axis.set_ylabel("Predicted response")
        axis.grid(alpha=0.2)
    axes[-1].legend(title="Varied factor", loc="best")
    return figure
