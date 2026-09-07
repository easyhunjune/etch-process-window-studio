from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_streamlit_app_renders_required_phase_c_controls() -> None:
    app_path = Path(__file__).resolve().parents[1] / "app.py"
    app = AppTest.from_file(str(app_path), default_timeout=30).run()
    assert not app.exception
    assert app.title[0].value == "Etch Process Window Studio"
    assert "외삽은 유효하지 않음" in app.warning[0].value
    assert len(app.tabs) == 4
    assert len(app.number_input) == 3
    assert len(app.selectbox) >= 4
    assert len(app.metric) >= 5


def test_streamlit_spec_change_recalculates_window() -> None:
    app_path = Path(__file__).resolve().parents[1] / "app.py"
    app = AppTest.from_file(str(app_path), default_timeout=30).run()
    initial = next(
        metric.value
        for metric in app.metric
        if metric.label == "평균예측 충족률"
    )
    app.number_input(key="spec_min_rate").set_value(0.48)
    app.run()
    assert not app.exception
    updated = next(
        metric.value
        for metric in app.metric
        if metric.label == "평균예측 충족률"
    )
    assert updated != initial
