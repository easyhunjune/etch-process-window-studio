from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import t as t_distribution

from .design import COEFFICIENT_NAMES, build_design_matrix


@dataclass(frozen=True)
class FullQuadraticRSM:
    coefficients: np.ndarray
    coded_design: np.ndarray
    gram_inverse: np.ndarray
    n_observations: int
    rank: int
    sse: float
    sst: float
    r_squared: float
    adjusted_r_squared: float
    residual_standard_error: float

    @classmethod
    def fit(
        cls,
        coded_values: np.ndarray,
        response: np.ndarray,
    ) -> FullQuadraticRSM:
        matrix = build_design_matrix(coded_values)
        y = np.asarray(response, dtype=float).reshape(-1)
        if matrix.shape[0] != y.size:
            raise ValueError("Factor rows and response rows must match.")
        if not np.isfinite(y).all():
            raise ValueError("Response values must be finite.")

        coefficients, _, rank, _ = np.linalg.lstsq(matrix, y, rcond=None)
        if rank != matrix.shape[1]:
            raise ValueError(
                f"Full quadratic design is rank deficient: {rank}/{matrix.shape[1]}."
            )
        fitted = matrix @ coefficients
        residuals = y - fitted
        sse = float(residuals @ residuals)
        centered = y - y.mean()
        sst = float(centered @ centered)
        r_squared = 1.0 - sse / sst if sst > 0 else 1.0
        df_residual = y.size - rank
        if df_residual <= 0:
            raise ValueError("The model needs residual degrees of freedom.")
        adjusted = 1.0 - (1.0 - r_squared) * (y.size - 1) / df_residual
        residual_standard_error = float(np.sqrt(sse / df_residual))
        return cls(
            coefficients=coefficients,
            coded_design=np.asarray(coded_values, dtype=float).copy(),
            gram_inverse=np.linalg.inv(matrix.T @ matrix),
            n_observations=y.size,
            rank=int(rank),
            sse=sse,
            sst=sst,
            r_squared=float(r_squared),
            adjusted_r_squared=float(adjusted),
            residual_standard_error=residual_standard_error,
        )

    @property
    def df_model(self) -> int:
        return self.rank - 1

    @property
    def df_residual(self) -> int:
        return self.n_observations - self.rank

    def predict(self, coded_values: np.ndarray) -> np.ndarray:
        return build_design_matrix(coded_values) @ self.coefficients

    def leverage(self, coded_values: np.ndarray) -> np.ndarray:
        """Return x'(X'X)^-1x for each requested design point."""
        matrix = build_design_matrix(coded_values)
        return np.einsum(
            "ij,jk,ik->i",
            matrix,
            self.gram_inverse,
            matrix,
        )

    def prediction_interval(
        self,
        coded_values: np.ndarray,
        *,
        confidence: float = 0.95,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return pointwise OLS prediction intervals for future observations."""
        if not 0.0 < confidence < 1.0:
            raise ValueError("confidence must lie between zero and one.")
        prediction = self.predict(coded_values)
        mean_squared_error = self.sse / self.df_residual
        standard_error = np.sqrt(
            mean_squared_error * (1.0 + self.leverage(coded_values))
        )
        critical = float(
            t_distribution.ppf(
                0.5 + confidence / 2.0,
                self.df_residual,
            )
        )
        margin = critical * standard_error
        return prediction - margin, prediction + margin

    def coefficient_dict(self) -> dict[str, float]:
        return {
            name: float(value)
            for name, value in zip(
                COEFFICIENT_NAMES, self.coefficients, strict=True
            )
        }
