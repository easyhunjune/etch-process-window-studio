from __future__ import annotations

import numpy as np
import pytest

from scripts.audit_paper_coefficients import (
    COEFFICIENT_NAMES,
    design_matrix,
    read_csv,
    unidentifiable_terms,
)

CCD_ROWS = read_csv("legtenberg_1995_ccd.csv")


def test_full_ccd_design_identifies_every_term() -> None:
    matrix = design_matrix(CCD_ROWS)
    assert matrix.shape == (32, 21)
    assert np.linalg.matrix_rank(matrix) == 21
    assert unidentifiable_terms(matrix) == []


def test_dropping_outward_slope_rows_leaves_the_design_rank_deficient() -> None:
    """Pin the rank-deficiency the 2026-08-13 review could not verify.

    The anisotropy audit refits after excluding the eight outward-sloped runs.
    That subset is 24x21 but only rank 19, so a term-by-term comparison against
    the paper is not identifiable for most coefficients.
    """
    keep = [row["anisotropy_outward_slope"] == "false" for row in CCD_ROWS]
    dropped = [
        row["run"] for row, kept in zip(CCD_ROWS, keep, strict=True) if not kept
    ]
    assert dropped == ["3", "9", "11", "15", "17", "20", "21", "24"]

    matrix = design_matrix(CCD_ROWS)[np.asarray(keep)]
    assert matrix.shape == (24, 21)
    assert np.linalg.matrix_rank(matrix) == 19

    unidentifiable = unidentifiable_terms(matrix)
    assert len(unidentifiable) == 15
    identifiable = [n for n in COEFFICIENT_NAMES if n not in unidentifiable]
    assert identifiable == ["b0", "b1", "b5", "b11", "b55", "b15"]


def test_unidentifiable_terms_names_only_the_columns_in_the_null_space() -> None:
    # A design whose last column duplicates the first leaves exactly those two
    # inseparable -- and, importantly, leaves the other 19 alone.
    matrix = np.eye(len(COEFFICIENT_NAMES))
    matrix[:, -1] = matrix[:, 0]
    assert np.linalg.matrix_rank(matrix) == len(COEFFICIENT_NAMES) - 1
    assert set(unidentifiable_terms(matrix)) == {
        COEFFICIENT_NAMES[0],
        COEFFICIENT_NAMES[-1],
    }


def test_unidentifiable_terms_rejects_the_wrong_number_of_columns() -> None:
    matrix = np.zeros((32, len(COEFFICIENT_NAMES) - 1))
    with pytest.raises(ValueError, match=r"expected 21 columns.*got 20"):
        unidentifiable_terms(matrix)
