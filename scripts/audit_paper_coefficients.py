from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data" / "raw"
FACTOR_NAMES = ["sf6", "o2", "chf3", "pressure", "rf_power"]
COEFFICIENT_NAMES = [
    "b0",
    "b1",
    "b2",
    "b3",
    "b4",
    "b5",
    "b11",
    "b22",
    "b33",
    "b44",
    "b55",
    "b12",
    "b13",
    "b14",
    "b15",
    "b23",
    "b24",
    "b25",
    "b34",
    "b35",
    "b45",
]
TARGETS = [
    "si_etch_rate_um_min",
    "si_sio2_selectivity",
    "anisotropy",
    "bias_voltage_v",
]


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA_DIR / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def unidentifiable_terms(matrix: np.ndarray) -> list[str]:
    """Return the coefficients a rank-deficient design cannot separate.

    ``np.linalg.lstsq`` answers even when the design is rank deficient: it
    returns the minimum-norm solution out of the infinitely many that fit
    equally well. That number is well defined but arbitrary, so comparing it
    term by term against a paper says nothing. Terms with a non-zero loading in
    the null space are exactly the ones where that applies.
    """
    expected_columns = len(COEFFICIENT_NAMES)
    if matrix.shape[1] != expected_columns:
        raise ValueError(
            f"expected {expected_columns} columns (one per COEFFICIENT_NAMES entry), "
            f"got {matrix.shape[1]}"
        )
    _, singular_values, right_vectors = np.linalg.svd(matrix)
    tolerance = max(matrix.shape) * np.finfo(float).eps * singular_values[0]
    kept = int((singular_values > tolerance).sum())
    null_space = right_vectors[kept:]
    if null_space.size == 0:
        return []
    loading = np.abs(null_space).max(axis=0)
    return [
        name
        for name, value in zip(COEFFICIENT_NAMES, loading, strict=True)
        if value > 1e-8
    ]


def design_matrix(rows: list[dict[str, str]]) -> np.ndarray:
    coded = np.asarray(
        [
            [float(row[f"{factor}_code"]) for factor in FACTOR_NAMES]
            for row in rows
        ]
    )
    columns = [np.ones(len(rows))]
    columns.extend(coded[:, index] for index in range(5))
    columns.extend(coded[:, index] ** 2 for index in range(5))
    columns.extend(
        coded[:, left] * coded[:, right]
        for left in range(4)
        for right in range(left + 1, 5)
    )
    return np.column_stack(columns)


def main() -> None:
    rows = read_csv("legtenberg_1995_ccd.csv")
    paper_rows = read_csv("legtenberg_1995_paper_coefficients.csv")[:21]
    matrix = design_matrix(rows)
    print(f"design matrix: {matrix.shape}, rank={np.linalg.matrix_rank(matrix)}")
    unidentifiable = unidentifiable_terms(matrix)
    if unidentifiable:
        matrix_rank = np.linalg.matrix_rank(matrix)
        print(
            "  WARNING: the full design is rank deficient "
            f"({matrix_rank}/{matrix.shape[1]}). lstsq still returns the "
            "minimum-norm fit, so the affected terms are not separately "
            f"identifiable: {', '.join(unidentifiable)}."
        )
    for target in TARGETS:
        response = np.asarray([float(row[target]) for row in rows])
        fitted = np.linalg.lstsq(matrix, response, rcond=None)[0]
        paper = np.asarray([float(row[target]) for row in paper_rows])
        print(f"\n{target}")
        for name, fitted_value, paper_value in zip(
            COEFFICIENT_NAMES, fitted, paper, strict=True
        ):
            flag = " [not identifiable]" if name in unidentifiable else ""
            print(
                f"{name:>3}: transcribed_fit={fitted_value:>9.4f} "
                f"paper={paper_value:>9.4f} delta={fitted_value-paper_value:>9.4f}{flag}"
            )

    # Table 4.3 prints outward-sloped profiles with an asterisk but only gives
    # the magnitude below 1. The signed geometry implied by Eq. 4.2 is 2-A for
    # those rows. Check this explicitly instead of silently changing raw data.
    signed_anisotropy = np.asarray(
        [
            2.0 - float(row["anisotropy"])
            if row["anisotropy_outward_slope"] == "true"
            else float(row["anisotropy"])
            for row in rows
        ]
    )
    signed_unidentifiable = unidentifiable_terms(matrix)
    fitted = np.linalg.lstsq(matrix, signed_anisotropy, rcond=None)[0]
    paper = np.asarray([float(row["anisotropy"]) for row in paper_rows])
    print("\nanisotropy_with_outward_slope_sign_restored")
    if signed_unidentifiable:
        matrix_rank = np.linalg.matrix_rank(matrix)
        print(
            "  WARNING: the signed-anisotropy design is rank deficient "
            f"({matrix_rank}/{matrix.shape[1]}). lstsq still returns the "
            "minimum-norm fit, so the affected terms are not separately "
            f"identifiable: {', '.join(signed_unidentifiable)}."
        )
    for name, fitted_value, paper_value in zip(
        COEFFICIENT_NAMES, fitted, paper, strict=True
    ):
        flag = " [not identifiable]" if name in signed_unidentifiable else ""
        print(
            f"{name:>3}: signed_fit={fitted_value:>9.4f} "
            f"paper={paper_value:>9.4f} delta={fitted_value-paper_value:>9.4f}{flag}"
        )

    non_outward_mask = np.asarray(
        [row["anisotropy_outward_slope"] == "false" for row in rows]
    )
    filtered_matrix = matrix[non_outward_mask]
    unidentifiable = unidentifiable_terms(filtered_matrix)
    fitted = np.linalg.lstsq(
        filtered_matrix,
        np.asarray([float(row["anisotropy"]) for row in rows])[non_outward_mask],
        rcond=None,
    )[0]
    filtered_rank = np.linalg.matrix_rank(filtered_matrix)
    print(
        "\nanisotropy_excluding_outward_slope_rows "
        f"(rows={filtered_matrix.shape[0]}, rank={filtered_rank}"
        f"/{filtered_matrix.shape[1]})"
    )
    if unidentifiable:
        dropped = [
            row["run"]
            for row, keep in zip(rows, non_outward_mask, strict=True)
            if not keep
        ]
        print(
            f"  WARNING: dropping runs {', '.join(dropped)} leaves the design rank "
            f"deficient ({filtered_rank}/{filtered_matrix.shape[1]}). lstsq still "
            "returns a fit -- the minimum-norm one out of infinitely many that fit "
            "the retained rows equally well."
        )
        print(
            f"  {len(unidentifiable)} of {len(COEFFICIENT_NAMES)} terms are not "
            "separately identifiable here, so their `delta` column below is an "
            "artefact of that choice and is NOT evidence about the paper: "
            f"{', '.join(unidentifiable)}."
        )
        print(
            "  Only "
            f"{', '.join(n for n in COEFFICIENT_NAMES if n not in unidentifiable)} "
            "can be read as a genuine comparison."
        )
    for name, fitted_value, paper_value in zip(
        COEFFICIENT_NAMES, fitted, paper, strict=True
    ):
        flag = " [not identifiable]" if name in unidentifiable else ""
        print(
            f"{name:>3}: filtered_fit={fitted_value:>9.4f} "
            f"paper={paper_value:>9.4f} delta={fitted_value-paper_value:>9.4f}{flag}"
        )


if __name__ == "__main__":
    main()
