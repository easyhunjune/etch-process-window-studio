from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.contour import ContourSet
from scipy.interpolate import RegularGridInterpolator

from etch_window.ui_support import (
    factor_contribution_figure,
    load_doe_and_models,
    main_effect_figure,
    process_window_figure,
    process_window_frame,
)
from etch_window.window import ProcessSpecs, evaluate_process_window

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_ui_support_builds_dynamic_window_and_figures() -> None:
    frame, models = load_doe_and_models(PROJECT_ROOT)
    assert len(frame) == 32
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.4, 8.0, 4.5),
        vary=("pressure", "rf_power"),
        fixed_actual={"sf6": 30.0, "o2": 10.0, "chf3": 12.0},
        resolution=41,
    )
    output = process_window_frame(result)
    assert len(output) == 41**2
    assert output["meets_all_specs"].dtype == bool
    assert 0 < output["nominally_meets_specs"].mean() < 1
    # This slice currently exercises both interval outcomes it can produce;
    # equality prevents either state from disappearing unnoticed.
    assert set(output["confidence_state"]) == {"uncertain", "rejected"}
    assert "supported_confidence_state" in output.columns
    unsupported = ~output["inside_design_region"]
    assert unsupported.any()
    assert (
        output.loc[unsupported, "supported_confidence_state"] == "unsupported"
    ).all()

    window_figure = process_window_figure(result, ProcessSpecs(0.4, 8.0, 4.5))
    effects_figure = main_effect_figure(models)
    factors = pd.read_csv(
        PROJECT_ROOT / "data" / "processed" / "factor_sensitivity.csv"
    )
    contribution_figure = factor_contribution_figure(factors)
    assert len(window_figure.axes) >= 4
    assert len(effects_figure.axes) == 3
    assert len(contribution_figure.axes) == 3
    plt.close("all")


def test_each_response_panel_draws_its_own_spec_boundary() -> None:
    _, models = load_doe_and_models(PROJECT_ROOT)
    specs = ProcessSpecs(0.4, 8.0, 4.5)
    result = evaluate_process_window(
        models,
        specs,
        vary=("pressure", "rf_power"),
        fixed_actual={"sf6": 30.0, "o2": 10.0, "chf3": 12.0},
        resolution=41,
    )

    figure = process_window_figure(result, specs)
    panels = (
        ("Si etch rate", "etch_rate", specs.min_etch_rate_um_min),
        ("Si/SiO", "selectivity", specs.min_selectivity),
        ("Non-uniformity", "nonuniformity", specs.max_nonuniformity_pct),
    )
    first = np.asarray(result["first_grid"])
    second = np.asarray(result["second_grid"])
    for title_prefix, key, threshold in panels:
        axis = next(a for a in figure.axes if a.get_title().startswith(title_prefix))
        contours = [
            child
            for child in axis.get_children()
            if isinstance(child, ContourSet) and len(child.levels) == 1
        ]
        # One white line per panel, at that panel's own spec — not the boundary of the
        # three specs combined, which would be identical across the three panels.
        assert [float(c.levels[0]) for c in contours] == [threshold]

        # The level alone would still pass if the line were drawn from another response,
        # so check the line actually sits where this response equals its spec.
        surface = RegularGridInterpolator(
            (second[:, 0], first[0, :]), np.asarray(result[key])
        )
        vertices = np.vstack([path.vertices for path in contours[0].get_paths()])
        # An empty vertex array would make the comparison below vacuously true, which is
        # exactly what happens when the level is drawn from a response that never reaches
        # it — the case this check exists to catch.
        assert vertices.size > 0
        on_line = surface(np.column_stack([vertices[:, 1], vertices[:, 0]]))
        assert on_line == pytest.approx(threshold, abs=1e-6)
    plt.close("all")


def test_ui_window_contains_no_extrapolated_points() -> None:
    _, models = load_doe_and_models(PROJECT_ROOT)
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.4, 8.0, 4.5),
        resolution=21,
    )
    coded = np.asarray(result["coded"])
    assert coded.min() >= -2.0
    assert coded.max() <= 2.0
    assert not np.asarray(result["inside_design_region"]).all()
