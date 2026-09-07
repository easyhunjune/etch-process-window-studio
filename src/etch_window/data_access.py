from __future__ import annotations

from pathlib import Path

import pandas as pd

from .design import FACTOR_NAMES
from .rsm import FullQuadraticRSM

RESPONSE_COLUMNS = {
    "etch_rate": "si_etch_rate_um_min",
    "selectivity": "si_sio2_selectivity",
    "nonuniformity": "etch_rate_nonuniformity_pct",
    "anisotropy": "anisotropy",
    "bias_voltage": "bias_voltage_v",
}


def load_doe_and_models(
    project_root: Path,
) -> tuple[pd.DataFrame, dict[str, FullQuadraticRSM]]:
    frame = pd.read_csv(project_root / "data" / "raw" / "legtenberg_1995_ccd.csv")
    coded = frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy(float)
    models = {
        name: FullQuadraticRSM.fit(coded, frame[column].to_numpy(float))
        for name, column in RESPONSE_COLUMNS.items()
    }
    return frame, models
