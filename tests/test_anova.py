from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from etch_window.anova import (
    factor_sensitivity,
    lack_of_fit_anova,
    overall_anova,
    partial_term_anova,
)
from etch_window.data_access import RESPONSE_COLUMNS
from etch_window.design import COEFFICIENT_NAMES, FACTOR_NAMES
from etch_window.rsm import FullQuadraticRSM

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_anova_partition_and_partial_contributions(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    response = doe_frame["si_etch_rate_um_min"].to_numpy()
    model = FullQuadraticRSM.fit(coded, response)
    overall = overall_anova(model)
    assert overall["regression_ss"] + overall["error_ss"] == pytest.approx(
        overall["total_ss"]
    )
    assert overall["df_model"] == 20
    assert overall["df_residual"] == 11
    assert 0 <= overall["p_value"] <= 1

    terms = partial_term_anova(coded, response, model)
    assert len(terms) == 20
    assert sum(row["normalized_partial_ss_pct"] for row in terms) == pytest.approx(
        100.0
    )
    assert all(row["partial_ss"] >= 0 for row in terms)


def test_factor_group_sensitivity_has_all_factors(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    response = doe_frame["etch_rate_nonuniformity_pct"].to_numpy()
    model = FullQuadraticRSM.fit(coded, response)
    rows = factor_sensitivity(coded, response, model)
    assert [row["factor"] for row in rows] == list(FACTOR_NAMES)
    assert sum(row["normalized_partial_ss_pct"] for row in rows) == pytest.approx(
        100.0
    )


def test_lack_of_fit_partitions_residual_error(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    response = doe_frame["si_etch_rate_um_min"].to_numpy()
    model = FullQuadraticRSM.fit(coded, response)
    result = lack_of_fit_anova(coded, response, model)
    assert result["lack_of_fit_ss"] + result["pure_error_ss"] == pytest.approx(
        model.sse
    )
    assert result["lack_of_fit_df"] == 6
    assert result["pure_error_df"] == 5
    assert result["f_value"] == pytest.approx(11.98, rel=0.01)
    assert result["p_value"] == pytest.approx(0.0077, rel=0.03)


@pytest.mark.parametrize(
    ("response_name", "expected_f", "expected_p"),
    [
        ("selectivity", 71.04, 0.000112),
        ("nonuniformity", 1.60, 0.312),
    ],
)
def test_other_responses_lack_of_fit_partitions_residual_error(
    doe_frame, response_name, expected_f, expected_p
) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    response = doe_frame[RESPONSE_COLUMNS[response_name]].to_numpy()
    model = FullQuadraticRSM.fit(coded, response)
    result = lack_of_fit_anova(coded, response, model)
    assert result["lack_of_fit_ss"] + result["pure_error_ss"] == pytest.approx(
        model.sse
    )
    assert result["lack_of_fit_df"] == 6
    assert result["pure_error_df"] == 5
    assert result["f_value"] == pytest.approx(expected_f, rel=0.01)
    assert result["p_value"] == pytest.approx(expected_p, rel=0.03)


def test_factor_sensitivity_matches_headlines_and_processed_csv(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    saved = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "factor_sensitivity.csv")
    expected_leaders = {
        "etch_rate": ("rf_power", 71.087252, 8.075319e-10),
        "selectivity": ("pressure", 39.437741, 0.0020782837),
        "nonuniformity": ("chf3", 36.904829, 0.110977582),
    }
    assert (saved["df"] == 6).all()
    for response_name, response_column in RESPONSE_COLUMNS.items():
        response = doe_frame[response_column].to_numpy()
        model = FullQuadraticRSM.fit(coded, response)
        rows = factor_sensitivity(coded, response, model)
        assert len(rows) == 5
        for row in rows:
            assert row["df"] == 6
            saved_row = saved.loc[
                (saved["response"] == response_name)
                & (saved["factor"] == row["factor"])
            ]
            assert len(saved_row) == 1
            for column in ("normalized_partial_ss_pct", "p_value"):
                assert row[column] == pytest.approx(
                    saved_row.iloc[0][column], rel=1e-8, abs=0.0
                )
        if response_name in expected_leaders:
            factor, contribution, p_value = expected_leaders[response_name]
            leader = max(rows, key=lambda row: row["normalized_partial_ss_pct"])
            assert leader["factor"] == factor
            assert leader["normalized_partial_ss_pct"] == pytest.approx(
                contribution, abs=0.01
            )
            assert leader["p_value"] == pytest.approx(p_value, rel=0.05, abs=0.0)


@pytest.mark.parametrize("response_name", RESPONSE_COLUMNS)
def test_partial_term_f_matches_coefficient_covariance_and_csv(
    doe_frame, response_name
) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    response = doe_frame[RESPONSE_COLUMNS[response_name]].to_numpy()
    model = FullQuadraticRSM.fit(coded, response)
    rows = partial_term_anova(coded, response, model)
    saved = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "anova_terms.csv")
    saved = saved.loc[saved["response"] == response_name].set_index("term")
    mse = model.sse / model.df_residual

    assert len(rows) == len(saved) == 20
    assert [row["term"] for row in rows] == list(COEFFICIENT_NAMES[1:])
    for row in rows:
        assert row["df"] == 1
        j = COEFFICIENT_NAMES.index(row["term"])
        expected_f = model.coefficients[j] ** 2 / (mse * model.gram_inverse[j, j])
        assert row["f_value"] == pytest.approx(expected_f, rel=1e-8, abs=0.0)
        assert row["f_value"] == pytest.approx(
            saved.loc[row["term"], "f_value"], rel=1e-8, abs=0.0
        )
