from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
from scipy.spatial import Delaunay

from .design import FACTOR_NAMES, FACTOR_SPECS, actual_row, actual_to_coded
from .rsm import FullQuadraticRSM


@dataclass(frozen=True)
class ProcessSpecs:
    min_etch_rate_um_min: float
    min_selectivity: float
    max_nonuniformity_pct: float
    min_anisotropy: float | None = None

    def __post_init__(self) -> None:
        numeric = (
            self.min_etch_rate_um_min,
            self.min_selectivity,
            self.max_nonuniformity_pct,
        )
        if not np.isfinite(numeric).all():
            raise ValueError("Process specifications must be finite.")
        if self.max_nonuniformity_pct < 0:
            raise ValueError("Maximum non-uniformity cannot be negative.")
        if self.min_anisotropy is not None and not np.isfinite(self.min_anisotropy):
            raise ValueError("Minimum anisotropy must be finite.")


def supported_confidence_state(
    confidence_state: np.ndarray,
    inside_design_region: np.ndarray,
) -> np.ndarray:
    """Fold design-region support into the interval-only confidence state.

    ``confidence_state`` answers one question — does the pointwise prediction
    interval clear every spec? — so a grid point far outside the DOE support
    can read ``"confirmed"`` on its own. Every consumer inside this project
    pairs it with ``inside_design_region``, but the raw ``process_window_grid``
    export leaves the pairing to whoever reads the file. This view spells the
    unsupported points out instead of relying on that discipline.

    Support takes precedence over the interval result: every point outside the
    design region is overwritten with ``"unsupported"``, even when its original
    ``confidence_state`` was ``"rejected"``. Count interval-defined rejected
    points from the original state (or ``evaluate_process_window``'s
    ``confidence_state`` output), not from this supported view.
    """
    state = np.asarray(confidence_state)
    inside = np.asarray(inside_design_region, dtype=bool)
    if state.shape != inside.shape:
        raise ValueError("State and support masks must share a shape.")
    supported = np.array(state, dtype=object)
    supported[~inside] = "unsupported"
    return supported


def evaluate_process_window(
    models: Mapping[str, FullQuadraticRSM],
    specs: ProcessSpecs,
    *,
    vary: tuple[str, str] = ("sf6", "o2"),
    fixed_actual: Mapping[str, float] | None = None,
    resolution: int = 81,
    confidence: float = 0.95,
) -> dict[str, np.ndarray | tuple[int, int] | tuple[str, str]]:
    required = {"etch_rate", "selectivity", "nonuniformity"}
    missing = required - set(models)
    if missing:
        raise ValueError(f"Missing response model(s): {sorted(missing)}")
    if specs.min_anisotropy is not None and "anisotropy" not in models:
        raise ValueError("An anisotropy model is required by the specification.")
    if len(set(vary)) != 2 or any(name not in FACTOR_NAMES for name in vary):
        raise ValueError("vary must name two distinct DOE factors.")
    if resolution < 2:
        raise ValueError("resolution must be at least 2.")

    base = actual_row(fixed_actual)
    first_index = FACTOR_NAMES.index(vary[0])
    second_index = FACTOR_NAMES.index(vary[1])
    first_spec = FACTOR_SPECS[vary[0]]
    second_spec = FACTOR_SPECS[vary[1]]
    first_axis = np.linspace(first_spec.actual_min, first_spec.actual_max, resolution)
    second_axis = np.linspace(
        second_spec.actual_min, second_spec.actual_max, resolution
    )
    first_grid, second_grid = np.meshgrid(first_axis, second_axis)
    actual = np.repeat(base.reshape(1, -1), resolution**2, axis=0)
    actual[:, first_index] = first_grid.ravel()
    actual[:, second_index] = second_grid.ravel()
    coded = actual_to_coded(actual, check_bounds=True)

    etch_rate = models["etch_rate"].predict(coded)
    selectivity = models["selectivity"].predict(coded)
    nonuniformity = models["nonuniformity"].predict(coded)
    nominal_feasible = (
        (etch_rate >= specs.min_etch_rate_um_min)
        & (selectivity >= specs.min_selectivity)
        & (nonuniformity <= specs.max_nonuniformity_pct)
    )
    interval_bounds = {
        name: models[name].prediction_interval(coded, confidence=confidence)
        for name in required
    }
    definitely_pass = (
        (interval_bounds["etch_rate"][0] >= specs.min_etch_rate_um_min)
        & (interval_bounds["selectivity"][0] >= specs.min_selectivity)
        & (
            interval_bounds["nonuniformity"][1]
            <= specs.max_nonuniformity_pct
        )
    )
    definitely_fail = (
        (interval_bounds["etch_rate"][1] < specs.min_etch_rate_um_min)
        | (interval_bounds["selectivity"][1] < specs.min_selectivity)
        | (
            interval_bounds["nonuniformity"][0]
            > specs.max_nonuniformity_pct
        )
    )

    reference_design = models["etch_rate"].coded_design
    for name in required - {"etch_rate"}:
        if not np.array_equal(models[name].coded_design, reference_design):
            raise ValueError("All response models must use the same DOE design.")
    hull = Delaunay(np.unique(reference_design, axis=0))
    inside_design_hull = hull.find_simplex(coded) >= 0
    prediction_leverage = models["etch_rate"].leverage(coded)
    observed_max_leverage = float(
        models["etch_rate"].leverage(reference_design).max()
    )
    inside_leverage_envelope = (
        prediction_leverage <= observed_max_leverage + 1e-12
    )
    inside_design_region = inside_design_hull & inside_leverage_envelope
    confidence_state = np.full(coded.shape[0], "uncertain", dtype=object)
    confidence_state[definitely_fail] = "rejected"
    confidence_state[definitely_pass] = "confirmed"
    feasible = definitely_pass & inside_design_region

    output: dict[str, np.ndarray | tuple[int, int] | tuple[str, str]] = {
        "vary": vary,
        "shape": first_grid.shape,
        "first_axis": first_axis,
        "second_axis": second_axis,
        "first_grid": first_grid,
        "second_grid": second_grid,
        "actual": actual,
        "coded": coded,
        "etch_rate": etch_rate.reshape(first_grid.shape),
        "selectivity": selectivity.reshape(first_grid.shape),
        "nonuniformity": nonuniformity.reshape(first_grid.shape),
        "nominal_feasible": nominal_feasible.reshape(first_grid.shape),
        "confidence_state": confidence_state.reshape(first_grid.shape),
        "inside_design_hull": inside_design_hull.reshape(first_grid.shape),
        "inside_leverage_envelope": inside_leverage_envelope.reshape(
            first_grid.shape
        ),
        "inside_design_region": inside_design_region.reshape(first_grid.shape),
        "prediction_leverage": prediction_leverage.reshape(first_grid.shape),
        "observed_max_leverage": np.asarray(observed_max_leverage),
        "feasible": feasible.reshape(first_grid.shape),
    }
    for name, (lower, upper) in interval_bounds.items():
        output[f"{name}_prediction_lower"] = lower.reshape(first_grid.shape)
        output[f"{name}_prediction_upper"] = upper.reshape(first_grid.shape)
    if "anisotropy" in models:
        anisotropy = models["anisotropy"].predict(coded)
        output["anisotropy"] = anisotropy.reshape(first_grid.shape)
        if specs.min_anisotropy is not None:
            lower, upper = models["anisotropy"].prediction_interval(
                coded,
                confidence=confidence,
            )
            output["anisotropy_prediction_lower"] = lower.reshape(
                first_grid.shape
            )
            output["anisotropy_prediction_upper"] = upper.reshape(
                first_grid.shape
            )
            anisotropy_pass = lower >= specs.min_anisotropy
            anisotropy_fail = upper < specs.min_anisotropy
            state = np.asarray(output["confidence_state"]).ravel()
            state[anisotropy_fail] = "rejected"
            state[
                (~anisotropy_fail)
                & (~anisotropy_pass)
                & (state != "rejected")
            ] = "uncertain"
            output["confidence_state"] = state.reshape(first_grid.shape)
            output["nominal_feasible"] = np.asarray(
                output["nominal_feasible"]
            ) & (np.asarray(output["anisotropy"]) >= specs.min_anisotropy)
            output["feasible"] = np.asarray(output["feasible"]) & (
                anisotropy_pass.reshape(first_grid.shape)
            )
    # Derived last so it picks up the anisotropy adjustment above.
    output["supported_confidence_state"] = supported_confidence_state(
        np.asarray(output["confidence_state"]),
        inside_design_region.reshape(first_grid.shape),
    )
    return output
