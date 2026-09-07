from __future__ import annotations

from typing import Any

import numpy as np
from scipy.stats import f as f_distribution

from .design import COEFFICIENT_NAMES, FACTOR_NAMES, TERM_FACTORS, build_design_matrix
from .rsm import FullQuadraticRSM


def overall_anova(model: FullQuadraticRSM) -> dict[str, float | int]:
    regression_ss = max(model.sst - model.sse, 0.0)
    regression_ms = regression_ss / model.df_model
    error_ms = model.sse / model.df_residual
    f_value = regression_ms / error_ms if error_ms > 0 else float("inf")
    p_value = float(
        f_distribution.sf(f_value, model.df_model, model.df_residual)
    )
    return {
        "regression_ss": regression_ss,
        "error_ss": model.sse,
        "total_ss": model.sst,
        "df_model": model.df_model,
        "df_residual": model.df_residual,
        "regression_ms": regression_ms,
        "error_ms": error_ms,
        "f_value": f_value,
        "p_value": p_value,
    }


def lack_of_fit_anova(
    coded_values: np.ndarray,
    response: np.ndarray,
    model: FullQuadraticRSM,
) -> dict[str, float | int]:
    """Partition residual error into pure error and lack of fit.

    Pure error is estimated from replicated rows with identical coded factors.
    """
    coded = np.asarray(coded_values, dtype=float)
    y = np.asarray(response, dtype=float).reshape(-1)
    if coded.shape[0] != y.size:
        raise ValueError("Factor rows and response rows must match.")

    _, inverse = np.unique(coded, axis=0, return_inverse=True)
    group_count = int(inverse.max()) + 1
    pure_error_ss = 0.0
    for group_index in range(group_count):
        values = y[inverse == group_index]
        pure_error_ss += float(np.sum((values - values.mean()) ** 2))

    pure_error_df = y.size - group_count
    lack_of_fit_df = group_count - model.rank
    if pure_error_df <= 0:
        raise ValueError("Lack-of-fit ANOVA requires replicated design points.")
    if lack_of_fit_df <= 0:
        raise ValueError("The model has no lack-of-fit degrees of freedom.")

    lack_of_fit_ss = max(model.sse - pure_error_ss, 0.0)
    pure_error_ms = pure_error_ss / pure_error_df
    lack_of_fit_ms = lack_of_fit_ss / lack_of_fit_df
    f_value = lack_of_fit_ms / pure_error_ms if pure_error_ms > 0 else float("inf")
    p_value = float(
        f_distribution.sf(f_value, lack_of_fit_df, pure_error_df)
    )
    return {
        "residual_ss": model.sse,
        "residual_df": model.df_residual,
        "lack_of_fit_ss": lack_of_fit_ss,
        "lack_of_fit_df": lack_of_fit_df,
        "lack_of_fit_ms": lack_of_fit_ms,
        "pure_error_ss": pure_error_ss,
        "pure_error_df": pure_error_df,
        "pure_error_ms": pure_error_ms,
        "f_value": f_value,
        "p_value": p_value,
    }


def _reduced_sse(matrix: np.ndarray, y: np.ndarray, removed: list[int]) -> float:
    keep = [index for index in range(matrix.shape[1]) if index not in removed]
    reduced = matrix[:, keep]
    coefficients = np.linalg.lstsq(reduced, y, rcond=None)[0]
    residuals = y - reduced @ coefficients
    return float(residuals @ residuals)


def partial_term_anova(
    coded_values: np.ndarray,
    response: np.ndarray,
    model: FullQuadraticRSM,
) -> list[dict[str, Any]]:
    """Return extra sums of squares from deleting each term from the full model."""
    matrix = build_design_matrix(coded_values)
    y = np.asarray(response, dtype=float)
    error_ms = model.sse / model.df_residual
    rows: list[dict[str, Any]] = []
    for index, name in enumerate(COEFFICIENT_NAMES[1:], start=1):
        partial_ss = max(_reduced_sse(matrix, y, [index]) - model.sse, 0.0)
        f_value = partial_ss / error_ms if error_ms > 0 else float("inf")
        rows.append(
            {
                "term": name,
                "partial_ss": partial_ss,
                "df": 1,
                "f_value": f_value,
                "p_value": float(
                    f_distribution.sf(f_value, 1, model.df_residual)
                ),
            }
        )
    total_partial = sum(row["partial_ss"] for row in rows)
    for row in rows:
        row["normalized_partial_ss_pct"] = (
            100.0 * row["partial_ss"] / total_partial if total_partial else 0.0
        )
    return rows


def factor_sensitivity(
    coded_values: np.ndarray,
    response: np.ndarray,
    model: FullQuadraticRSM,
) -> list[dict[str, Any]]:
    """Group-deletion sensitivity; interactions are counted in each parent factor."""
    matrix = build_design_matrix(coded_values)
    y = np.asarray(response, dtype=float)
    error_ms = model.sse / model.df_residual
    rows: list[dict[str, Any]] = []
    for factor_index, factor_name in enumerate(FACTOR_NAMES):
        removed = [
            term_index
            for term_index, factors in enumerate(TERM_FACTORS)
            if factor_index in factors
        ]
        partial_ss = max(_reduced_sse(matrix, y, removed) - model.sse, 0.0)
        df = len(removed)
        f_value = (partial_ss / df) / error_ms if error_ms > 0 else float("inf")
        center_line = np.zeros((201, len(FACTOR_NAMES)))
        center_line[:, factor_index] = np.linspace(-2.0, 2.0, 201)
        predictions = model.predict(center_line)
        rows.append(
            {
                "factor": factor_name,
                "partial_ss": partial_ss,
                "df": df,
                "f_value": f_value,
                "p_value": float(
                    f_distribution.sf(f_value, df, model.df_residual)
                ),
                "centerline_prediction_range": float(np.ptp(predictions)),
            }
        )
    total_partial = sum(row["partial_ss"] for row in rows)
    for row in rows:
        row["normalized_partial_ss_pct"] = (
            100.0 * row["partial_ss"] / total_partial if total_partial else 0.0
        )
    return rows
