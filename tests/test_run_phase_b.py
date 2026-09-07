from __future__ import annotations

from pathlib import Path

from etch_window.data_access import load_doe_and_models
from etch_window.window import ProcessSpecs, evaluate_process_window
from scripts.run_phase_b import window_table

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_window_table_marks_points_outside_design_support() -> None:
    _, models = load_doe_and_models(PROJECT_ROOT)
    result = evaluate_process_window(
        models,
        ProcessSpecs(0.4, 8.0, 4.5),
        vary=("sf6", "o2"),
        resolution=21,
    )
    output = window_table(result)

    assert "supported_prediction_confidence_state" in output.columns
    unsupported = ~output["inside_design_region"]
    assert unsupported.any()
    assert (
        output.loc[unsupported, "supported_prediction_confidence_state"]
        == "unsupported"
    ).all()
