from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import t as t_distribution

from etch_window.design import COEFFICIENT_NAMES, FACTOR_NAMES, build_design_matrix
from etch_window.rsm import FullQuadraticRSM


def test_design_matrix_matches_paper_order(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    matrix = build_design_matrix(coded)
    assert matrix.shape == (32, len(COEFFICIENT_NAMES))
    assert np.linalg.matrix_rank(matrix) == 21
    first = coded[0]
    assert matrix[0, 1:6] == pytest.approx(first)
    assert matrix[0, 6:11] == pytest.approx(first**2)
    assert matrix[0, 11:] == pytest.approx(
        [first[i] * first[j] for i in range(4) for j in range(i + 1, 5)]
    )


def test_etch_rate_model_reproduces_paper_metrics(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    model = FullQuadraticRSM.fit(
        coded, doe_frame["si_etch_rate_um_min"].to_numpy()
    )
    assert model.coefficient_dict()["b0"] == pytest.approx(0.434, abs=0.001)
    assert model.coefficient_dict()["b5"] == pytest.approx(0.147, abs=0.001)
    assert model.r_squared == pytest.approx(0.99, abs=0.005)
    # The transcribed table reproduces R², but not the paper's printed s=0.01.
    assert model.residual_standard_error == pytest.approx(0.025994, abs=1e-6)


def test_predictions_equal_training_projection(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    response = doe_frame["si_sio2_selectivity"].to_numpy()
    model = FullQuadraticRSM.fit(coded, response)
    residual = response - model.predict(coded)
    assert float(residual @ residual) == pytest.approx(model.sse)


def test_leverage_and_prediction_intervals_expand_uncertainty(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    response = doe_frame["si_sio2_selectivity"].to_numpy()
    model = FullQuadraticRSM.fit(coded, response)
    center = np.zeros((1, len(FACTOR_NAMES)))
    corner = np.full((1, len(FACTOR_NAMES)), 2.0)
    assert model.leverage(corner)[0] > model.leverage(center)[0]
    lower, upper = model.prediction_interval(center)
    assert lower.shape == upper.shape == (1,)
    assert lower[0] < model.predict(center)[0] < upper[0]


def test_prediction_interval_half_width_matches_the_closed_form(doe_frame) -> None:
    """Pin the width, not just the ordering.

    ``lower < prediction < upper`` holds for a confidence interval too, and the
    difference between the two is the whole point: a future observation carries
    the residual variance on top of the estimation variance, so the ``1 +`` in
    ``MSE * (1 + leverage)`` must be there.
    """
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    model = FullQuadraticRSM.fit(
        coded, doe_frame["si_etch_rate_um_min"].to_numpy()
    )
    probes = np.vstack(
        [
            np.zeros(len(FACTOR_NAMES)),
            np.full(len(FACTOR_NAMES), 1.0),
            np.eye(len(FACTOR_NAMES))[0] * 2.0,
        ]
    )
    for confidence in (0.90, 0.95, 0.99):
        lower, upper = model.prediction_interval(probes, confidence=confidence)
        mean_squared_error = model.sse / model.df_residual
        expected = t_distribution.ppf(
            0.5 + confidence / 2.0, model.df_residual
        ) * np.sqrt(mean_squared_error * (1.0 + model.leverage(probes)))
        assert (upper - lower) / 2.0 == pytest.approx(expected)
        assert (lower + upper) / 2.0 == pytest.approx(model.predict(probes))

    # A confidence interval on the mean response would drop the leading 1 and be
    # strictly narrower everywhere; make sure that is not what we computed.
    confidence_half_width = t_distribution.ppf(0.975, model.df_residual) * np.sqrt(
        (model.sse / model.df_residual) * model.leverage(probes)
    )
    lower, upper = model.prediction_interval(probes)
    assert ((upper - lower) / 2.0 > confidence_half_width).all()


def test_prediction_interval_rejects_impossible_confidence(doe_frame) -> None:
    coded = doe_frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy()
    model = FullQuadraticRSM.fit(
        coded, doe_frame["si_etch_rate_um_min"].to_numpy()
    )
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="between zero and one"):
            model.prediction_interval(coded, confidence=bad)
