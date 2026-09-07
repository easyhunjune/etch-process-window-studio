from __future__ import annotations

import pytest

from etch_window.anova import (
    factor_sensitivity,
    lack_of_fit_anova,
    overall_anova,
    partial_term_anova,
)
from etch_window.design import FACTOR_NAMES
from etch_window.rsm import FullQuadraticRSM


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
