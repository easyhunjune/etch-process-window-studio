from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture(scope="session")
def doe_frame() -> pd.DataFrame:
    path = Path(__file__).resolve().parents[1] / "data" / "raw" / (
        "legtenberg_1995_ccd.csv"
    )
    return pd.read_csv(path)
