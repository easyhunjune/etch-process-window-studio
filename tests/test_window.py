from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pytest

import etch_window.window as window_module
from etch_window.data_access import load_doe_and_models
from etch_window.design import FACTOR_NAMES, actual_to_coded
from etch_window.rsm import FullQuadraticRSM
from etch_window.window import (
    ProcessSpecs,
    evaluate_process_window,
    supported_confidence_state,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def real_doe_models() -> dict[str, FullQuadraticRSM]:
    _, models = load_doe_and_models(PROJECT_ROOT)
    return models


def _synthetic_model(intercept: float, sf6_linear: float = 0.0) -> FullQuadraticRSM:
    coded = np.vstack(
        [
            np.array(values, dtype=float)
            for values in itertools.product((-1.0, 1.0), repeat=5)
        ]
        + [
            np.eye(5)[index] * sign * 2.0
            for index in range(5)
            for sign in (-1.0, 1.0)
        ]
        + [np.zeros(5) for _ in range(6)]
    )
    response = intercept + sf6_linear * coded[:, 0]
    return FullQuadraticRSM.fit(coded, response)


def test_process_window_applies_all_specs() -> None:
    models = {
        "etch_rate": _synthetic_model(0.5, 0.1),
        "selectivity": _synthetic_model(10.0),
        "nonuniformity": _synthetic_model(3.0),
    }
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.5, 8.0, 4.0),
        resolution=9,
    )
    nominal = np.asarray(result["nominal_feasible"])
    assert nominal.shape == (9, 9)
    assert not nominal[:, :4].any()
    assert nominal[:, 4:].all()
    assert np.asarray(result["feasible"]).shape == (9, 9)


def test_process_window_masks_factorwise_valid_but_unsupported_corner() -> None:
    models = {
        "etch_rate": _synthetic_model(0.5),
        "selectivity": _synthetic_model(10.0),
        "nonuniformity": _synthetic_model(3.0),
    }
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.4, 8.0, 4.0),
        resolution=9,
    )
    inside_hull = np.asarray(result["inside_design_hull"])
    inside_leverage = np.asarray(result["inside_leverage_envelope"])
    assert not inside_hull[0, 0]
    assert not inside_hull[-1, -1]
    assert not inside_leverage[-1, -1]
    assert inside_hull[4, 4]


def test_feasible_is_false_wherever_the_design_cannot_support_the_point() -> None:
    """Check the values of ``feasible``, not only its shape.

    These models are exactly constant and clear every spec, so the prediction
    intervals collapse and ``definitely_pass`` is true across the whole grid.
    Anything ``feasible`` rules out here it rules out for the support reason
    alone, which is what makes the two claims below separable.
    """
    models = {
        "etch_rate": _synthetic_model(0.5),
        "selectivity": _synthetic_model(10.0),
        "nonuniformity": _synthetic_model(3.0),
    }
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.4, 8.0, 4.0),
        resolution=9,
    )
    feasible = np.asarray(result["feasible"], dtype=bool)
    inside = np.asarray(result["inside_design_region"], dtype=bool)
    assert (np.asarray(result["confidence_state"]) == "confirmed").all()

    assert not feasible[0, 0]
    assert not feasible[-1, -1]
    assert feasible[4, 4]
    assert (feasible == inside).all()
    # Not a vacuous comparison: the grid has to contain both kinds of point.
    assert 0 < feasible.mean() < 1


def test_feasible_is_false_where_the_interval_straddles_a_spec() -> None:
    """The other half of ``feasible``: support alone is not enough."""
    models = {
        "etch_rate": _synthetic_model(0.5, 0.1),
        "selectivity": _synthetic_model(10.0),
        "nonuniformity": _synthetic_model(3.0),
    }
    result = evaluate_process_window(
        models,
        # Set the etch-rate floor above every value the sloped model produces.
        ProcessSpecs(5.0, 8.0, 4.0),
        resolution=9,
    )
    assert not np.asarray(result["nominal_feasible"], dtype=bool).any()
    assert not np.asarray(result["feasible"], dtype=bool).any()
    assert np.asarray(result["inside_design_region"], dtype=bool).any()


def test_unsupported_points_are_labelled_in_the_supported_state() -> None:
    models = {
        "etch_rate": _synthetic_model(0.5),
        "selectivity": _synthetic_model(10.0),
        "nonuniformity": _synthetic_model(3.0),
    }
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.4, 8.0, 4.0),
        resolution=9,
    )
    state = np.asarray(result["confidence_state"])
    supported = np.asarray(result["supported_confidence_state"])
    inside = np.asarray(result["inside_design_region"], dtype=bool)

    # These models are exactly constant, so the prediction intervals collapse and
    # every point clears the specs -- including the corners the design cannot
    # support. That is the case a raw-CSV reader would misread.
    assert state[0, 0] == "confirmed"
    assert not inside[0, 0]
    assert supported[0, 0] == "unsupported"

    assert (supported[~inside] == "unsupported").all()
    assert (supported[inside] == state[inside]).all()


def test_conversion_rejects_extrapolation() -> None:
    center = np.array([[30.0, 10.0, 12.0, 100.0, 100.0]])
    assert actual_to_coded(center) == pytest.approx(np.zeros((1, 5)))
    outside = center.copy()
    outside[0, FACTOR_NAMES.index("sf6")] = 51.0
    with pytest.raises(ValueError, match="outside the DOE range"):
        actual_to_coded(outside)


def test_real_doe_distinguishes_nominal_from_interval_feasibility(
    real_doe_models: dict[str, FullQuadraticRSM],
) -> None:
    """Pin the 95% prediction-interval decision used by feasible."""
    result = evaluate_process_window(
        real_doe_models,
        ProcessSpecs(0.4, 8.0, 4.5),
        vary=("sf6", "o2"),
        resolution=81,
    )
    nominal = np.asarray(result["nominal_feasible"], dtype=bool)
    feasible = np.asarray(result["feasible"], dtype=bool)
    inside = np.asarray(result["inside_design_region"], dtype=bool)
    state = np.asarray(result["confidence_state"])

    assert nominal.size == 81**2
    assert 0 < nominal.mean() < 1
    # Real residual variance makes many supported mean predictions pass while
    # their 95% prediction intervals do not clear all three specifications.
    assert (inside & nominal & ~feasible).any()
    assert feasible.sum() == 0
    assert (state == "rejected").sum() > 0

    summary = json.loads(
        (PROJECT_ROOT / "data" / "processed" / "phase_b_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert nominal.mean() == pytest.approx(
        summary["window"]["nominal_grid_fraction"]
    )
    assert feasible.mean() == pytest.approx(
        summary["window"]["feasible_grid_fraction"]
    )


def test_etch_rate_must_clear_the_prediction_lower_bound(
    real_doe_models: dict[str, FullQuadraticRSM],
) -> None:
    # Make the other two specifications non-binding so an etch-rate interval
    # that straddles 0.4 is an observable uncertain point, not a hidden detail.
    result = evaluate_process_window(
        real_doe_models,
        ProcessSpecs(0.4, -1e6, 1e6),
        vary=("sf6", "o2"),
        resolution=81,
    )
    lower = np.asarray(result["etch_rate_prediction_lower"])
    upper = np.asarray(result["etch_rate_prediction_upper"])
    state = np.asarray(result["confidence_state"])
    straddles = (lower < 0.4) & (upper >= 0.4)

    assert straddles.any()
    assert (state[straddles] == "uncertain").all()


def test_anisotropy_adjustment_precedes_supported_state_derivation(
    real_doe_models: dict[str, FullQuadraticRSM],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # At the default 95% level the primary models have no confirmed points.
    # A 50% interval keeps the real DOE residuals while making the anisotropy
    # downgrade observable, so moving the derived-state call earlier is caught.
    base = evaluate_process_window(
        real_doe_models,
        ProcessSpecs(0.4, 8.0, 4.5),
        vary=("sf6", "o2"),
        resolution=81,
        confidence=0.5,
    )
    observed_inputs: list[np.ndarray] = []
    original = window_module.supported_confidence_state

    def record_supported_input(
        state: np.ndarray,
        inside: np.ndarray,
    ) -> np.ndarray:
        observed_inputs.append(np.asarray(state).copy())
        return original(state, inside)

    monkeypatch.setattr(
        window_module, "supported_confidence_state", record_supported_input
    )
    constrained = evaluate_process_window(
        real_doe_models,
        ProcessSpecs(0.4, 8.0, 4.5, min_anisotropy=0.5),
        vary=("sf6", "o2"),
        resolution=81,
        confidence=0.5,
    )
    base_state = np.asarray(base["confidence_state"])
    constrained_state = np.asarray(constrained["confidence_state"])
    inside = np.asarray(constrained["inside_design_region"], dtype=bool)
    supported = np.asarray(constrained["supported_confidence_state"])
    downgraded = (base_state == "confirmed") & np.isin(
        constrained_state, ["uncertain", "rejected"]
    )

    assert downgraded.any()
    assert (base_state != constrained_state).any()
    assert inside.any()
    assert (supported[inside] == constrained_state[inside]).all()
    # This records call order directly; output masking outside DOE support
    # would otherwise hide that the derivation happened before anisotropy.
    assert len(observed_inputs) == 1
    assert (observed_inputs[0] == constrained_state).all()


def test_supported_confidence_state_preserves_the_full_unsupported_label() -> None:
    result = supported_confidence_state(
        np.array(["confirmed", "uncertain"]),
        np.array([True, False]),
    )
    assert result.tolist() == ["confirmed", "unsupported"]
    assert result[1] == "unsupported"
    assert len(result[1]) == len("unsupported")


def test_supported_confidence_state_rejects_mismatched_shapes() -> None:
    with pytest.raises(ValueError, match="share a shape"):
        supported_confidence_state(
            np.array(["confirmed", "uncertain"]),
            np.array([[True, False]]),
        )
