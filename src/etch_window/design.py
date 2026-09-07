from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FactorSpec:
    name: str
    display_name: str
    unit: str
    center: float
    step: float
    coded_min: float = -2.0
    coded_max: float = 2.0

    @property
    def actual_min(self) -> float:
        return self.center + self.step * self.coded_min

    @property
    def actual_max(self) -> float:
        return self.center + self.step * self.coded_max


FACTOR_NAMES = ("sf6", "o2", "chf3", "pressure", "rf_power")
FACTOR_SPECS = {
    "sf6": FactorSpec("sf6", "SF₆ flow", "sccm", 30.0, 10.0),
    "o2": FactorSpec("o2", "O₂ flow", "sccm", 10.0, 4.0),
    "chf3": FactorSpec("chf3", "CHF₃ flow", "sccm", 12.0, 5.0),
    "pressure": FactorSpec("pressure", "Pressure", "mTorr", 100.0, 40.0),
    "rf_power": FactorSpec("rf_power", "RF power", "W", 100.0, 40.0),
}

COEFFICIENT_NAMES = (
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
)

TERM_FACTORS: tuple[tuple[int, ...], ...] = (
    (),
    (0,),
    (1,),
    (2,),
    (3,),
    (4,),
    (0,),
    (1,),
    (2,),
    (3,),
    (4,),
    (0, 1),
    (0, 2),
    (0, 3),
    (0, 4),
    (1, 2),
    (1, 3),
    (1, 4),
    (2, 3),
    (2, 4),
    (3, 4),
)


def _as_2d(values: np.ndarray | Sequence[Sequence[float]]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.ndim != 2 or array.shape[1] != len(FACTOR_NAMES):
        raise ValueError(
            f"Expected shape (n, {len(FACTOR_NAMES)}), got {array.shape}."
        )
    if not np.isfinite(array).all():
        raise ValueError("Factor values must be finite.")
    return array


def build_design_matrix(
    coded_values: np.ndarray | Sequence[Sequence[float]],
) -> np.ndarray:
    """Build the 21-column model matrix in the paper's coefficient order."""
    coded = _as_2d(coded_values)
    columns = [np.ones(coded.shape[0], dtype=float)]
    columns.extend(coded[:, index] for index in range(5))
    columns.extend(coded[:, index] ** 2 for index in range(5))
    columns.extend(
        coded[:, left] * coded[:, right]
        for left in range(4)
        for right in range(left + 1, 5)
    )
    return np.column_stack(columns)


def actual_to_coded(
    actual_values: np.ndarray | Sequence[Sequence[float]],
    *,
    check_bounds: bool = True,
) -> np.ndarray:
    actual = _as_2d(actual_values)
    centers = np.asarray([FACTOR_SPECS[name].center for name in FACTOR_NAMES])
    steps = np.asarray([FACTOR_SPECS[name].step for name in FACTOR_NAMES])
    coded = (actual - centers) / steps
    if check_bounds:
        lower = np.asarray([FACTOR_SPECS[name].coded_min for name in FACTOR_NAMES])
        upper = np.asarray([FACTOR_SPECS[name].coded_max for name in FACTOR_NAMES])
        outside = (coded < lower - 1e-12) | (coded > upper + 1e-12)
        if outside.any():
            row, column = np.argwhere(outside)[0]
            spec = FACTOR_SPECS[FACTOR_NAMES[column]]
            raise ValueError(
                f"{spec.name}={actual[row, column]:g} {spec.unit} is outside "
                f"the DOE range [{spec.actual_min:g}, {spec.actual_max:g}]."
            )
    return coded


def coded_to_actual(
    coded_values: np.ndarray | Sequence[Sequence[float]],
    *,
    check_bounds: bool = True,
) -> np.ndarray:
    coded = _as_2d(coded_values)
    if check_bounds and ((coded < -2.0 - 1e-12) | (coded > 2.0 + 1e-12)).any():
        raise ValueError("Coded values must remain inside the DOE range [-2, 2].")
    centers = np.asarray([FACTOR_SPECS[name].center for name in FACTOR_NAMES])
    steps = np.asarray([FACTOR_SPECS[name].step for name in FACTOR_NAMES])
    return centers + coded * steps


def actual_row(values: Mapping[str, float] | None = None) -> np.ndarray:
    supplied = values or {}
    unknown = set(supplied) - set(FACTOR_NAMES)
    if unknown:
        raise ValueError(f"Unknown factor(s): {sorted(unknown)}")
    return np.asarray(
        [supplied.get(name, FACTOR_SPECS[name].center) for name in FACTOR_NAMES],
        dtype=float,
    )
