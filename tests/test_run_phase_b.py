from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from etch_window.data_access import load_doe_and_models
from etch_window.design import FACTOR_NAMES
from etch_window.rsm import FullQuadraticRSM
from etch_window.window import ProcessSpecs, evaluate_process_window
from scripts.run_phase_b import (
    COEFFICIENT_TOLERANCE,
    PAPER_RESPONSE_COLUMNS,
    coefficient_tables,
    fit_models,
    metric_table,
    window_table,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_window_table_marks_points_outside_design_support() -> None:
    _, models = load_doe_and_models(PROJECT_ROOT)
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.4, 8.0, 4.5),
        vary=("sf6", "o2"),
        resolution=21,
    )
    output = window_table(result)

    assert "supported_prediction_confidence_state" in output.columns
    unsupported = ~output["inside_design_region"]
    assert unsupported.any()
    assert (
        output.loc[unsupported, "supported_prediction_confidence_state"]
        == "unsupported"
    ).all()


def test_coefficient_comparison_matches_headlines_and_processed_csv() -> None:
    frame = pd.read_csv(PROJECT_ROOT / "data/raw/legtenberg_1995_ccd.csv")
    paper = pd.read_csv(
        PROJECT_ROOT / "data/raw/legtenberg_1995_paper_coefficients.csv"
    )
    _, models = fit_models(frame)
    _, comparison = coefficient_tables(models, paper)

    assert len(comparison) == 84
    mismatches = comparison.loc[~comparison["within_rounding_tolerance"]]
    assert len(mismatches) == 18
    assert mismatches.groupby("response").size().to_dict() == {
        "anisotropy": 14,
        "etch_rate": 2,
        "selectivity": 1,
        "bias_voltage": 1,
    }
    for response, symbol, fitted, printed in (
        ("etch_rate", "b3", -0.0044, 0.004),
        ("etch_rate", "b4", -0.0340, -0.043),
        ("selectivity", "b3", 0.1417, -0.14),
        ("bias_voltage", "b34", -6.5625, 7.0),
    ):
        row = comparison.loc[
            (comparison["response"] == response)
            & (comparison["paper_symbol"] == symbol)
        ]
        assert len(row) == 1
        assert row.iloc[0]["fitted_coefficient"] == pytest.approx(fitted, abs=0.0001)
        assert row.iloc[0]["paper_coefficient"] == pytest.approx(printed, abs=0.0001)

    saved = pd.read_csv(
        PROJECT_ROOT / "data/processed/paper_coefficient_comparison.csv",
        keep_default_na=False,
    )
    keys = ["response", "paper_symbol"]
    comparison = comparison.sort_values(keys).reset_index(drop=True)
    saved = saved.sort_values(keys).reset_index(drop=True)
    assert comparison.columns.tolist() == saved.columns.tolist()
    for column in (*keys, "within_rounding_tolerance", "note"):
        assert comparison[column].tolist() == saved[column].tolist()
    for column in (
        "fitted_coefficient",
        "paper_coefficient",
        "delta",
        "absolute_delta",
        "rounding_tolerance",
    ):
        assert comparison[column].to_numpy() == pytest.approx(
            saved[column].to_numpy(), rel=1e-8, abs=1e-12
        )


def test_metric_table_matches_processed_csv() -> None:
    frame = pd.read_csv(PROJECT_ROOT / "data/raw/legtenberg_1995_ccd.csv")
    paper = pd.read_csv(
        PROJECT_ROOT / "data/raw/legtenberg_1995_paper_coefficients.csv"
    )
    _, models = fit_models(frame)
    metrics = metric_table(models, paper).sort_values("response").reset_index(drop=True)
    saved = pd.read_csv(PROJECT_ROOT / "data/processed/model_metrics.csv")
    saved = saved.sort_values("response").reset_index(drop=True)

    assert len(metrics) == len(saved) == 5
    assert metrics.columns.tolist() == saved.columns.tolist()
    for column in ("response", "model_status"):
        assert metrics[column].tolist() == saved[column].tolist()
    for column in metrics.columns.difference(["response", "model_status"]):
        assert metrics[column].to_numpy() == pytest.approx(
            saved[column].to_numpy(), rel=1e-8, nan_ok=True
        )


@pytest.mark.parametrize("sign", [1, -1])
def test_coefficient_rounding_tolerance_includes_exact_boundary(
    doe_frame, sign
) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    # A zero response gives exactly zero coefficients, avoiding fit roundoff
    # when constructing deltas exactly on the comparison boundary.
    model = FullQuadraticRSM.fit(coded, np.zeros(len(coded)))
    models = {response: model for response in PAPER_RESPONSE_COLUMNS}
    paper = pd.DataFrame({"paper_symbol": ["b0", "b1"]})
    for response, column in PAPER_RESPONSE_COLUMNS.items():
        tolerance = COEFFICIENT_TOLERANCE[response]
        paper[column] = [-sign * tolerance, -sign * (tolerance + 1e-9)]

    _, comparison = coefficient_tables(models, paper)
    for response in PAPER_RESPONSE_COLUMNS:
        rows = comparison.loc[comparison["response"] == response].set_index(
            "paper_symbol"
        )
        tolerance = COEFFICIENT_TOLERANCE[response]
        assert rows.loc["b0", "delta"] == sign * tolerance
        assert rows.loc["b1", "delta"] == sign * (tolerance + 1e-9)
        assert bool(rows.loc["b0", "within_rounding_tolerance"]) is True
        assert bool(rows.loc["b1", "within_rounding_tolerance"]) is False


def test_etch_rate_coefficient_tolerance_is_fixed() -> None:
    assert COEFFICIENT_TOLERANCE["etch_rate"] == 0.0015
