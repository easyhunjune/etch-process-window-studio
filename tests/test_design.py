from __future__ import annotations

import numpy as np

from etch_window.design import TERM_FACTORS, build_design_matrix


def test_design_matrix_columns_match_term_factors() -> None:
    coded = np.random.default_rng(0).normal(size=(50, 5))
    matrix = build_design_matrix(coded)

    assert matrix.shape == (50, 21)
    assert len(TERM_FACTORS) == 21
    for column, factors in enumerate(TERM_FACTORS):
        if not factors:
            expected = np.ones(50)
        elif len(factors) == 1:
            (factor,) = factors
            if 1 <= column <= 5:
                expected = coded[:, factor]
            else:
                assert 6 <= column <= 10
                expected = coded[:, factor] ** 2
        else:
            left, right = factors
            expected = coded[:, left] * coded[:, right]
        assert np.allclose(matrix[:, column], expected), (column, factors)
