from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from etch_window import provenance
from etch_window.provenance import (
    DIGEST_FIELD,
    SUMMARY_RELATIVE_PATH,
    TRACKED_INPUT_FILES,
    current_input_digests,
    stale_inputs,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _build_project(root: Path, *, record_digests: bool = True) -> Path:
    """Copy the tracked inputs into a scratch tree with a matching summary file."""
    for name in TRACKED_INPUT_FILES:
        destination = root / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((PROJECT_ROOT / name).read_bytes())
    summary: dict[str, object] = {"dataset_rows": 32}
    if record_digests:
        summary[DIGEST_FIELD] = current_input_digests(root)
    summary_path = root / SUMMARY_RELATIVE_PATH
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    return root


def test_digests_cover_every_tracked_input_and_match_the_files(tmp_path) -> None:
    root = _build_project(tmp_path)
    digests = current_input_digests(root)
    assert set(digests) == set(TRACKED_INPUT_FILES)
    for name, digest in digests.items():
        expected = hashlib.sha256((root / name).read_bytes()).hexdigest()
        assert digest == expected


def test_tracked_inputs_include_raw_data_and_the_numeric_core() -> None:
    # A digest list covering only the raw CSVs would miss the other half of the
    # drift: editing rsm.py and skipping the rebuild.
    assert "data/raw/legtenberg_1995_ccd.csv" in TRACKED_INPUT_FILES
    assert "src/etch_window/rsm.py" in TRACKED_INPUT_FILES
    assert "src/etch_window/anova.py" in TRACKED_INPUT_FILES


def test_freshly_built_project_reports_nothing_stale(tmp_path) -> None:
    assert stale_inputs(_build_project(tmp_path)) == []


def test_repository_tables_match_the_current_inputs() -> None:
    # Guards the checked-in data/processed/ tables, not just the synthetic copy.
    assert stale_inputs(PROJECT_ROOT) == []


@pytest.mark.parametrize(
    "changed",
    [
        "data/raw/legtenberg_1995_ccd.csv",
        "src/etch_window/data_access.py",
        "src/etch_window/rsm.py",
        "scripts/run_phase_b.py",
    ],
)
def test_edited_input_is_reported_as_stale(tmp_path, changed: str) -> None:
    root = _build_project(tmp_path)
    target = root / changed
    target.write_bytes(target.read_bytes() + b"\n# edited after the rebuild\n")
    assert stale_inputs(root) == [changed]


def test_deleted_input_is_reported_rather_than_ignored(tmp_path) -> None:
    root = _build_project(tmp_path)
    (root / "src" / "etch_window" / "anova.py").unlink()
    with pytest.raises(FileNotFoundError):
        stale_inputs(root)


def test_missing_summary_raises(tmp_path) -> None:
    root = _build_project(tmp_path)
    (root / SUMMARY_RELATIVE_PATH).unlink()
    with pytest.raises(FileNotFoundError, match="phase_b_summary.json"):
        stale_inputs(root)


def test_summary_without_digest_block_raises(tmp_path) -> None:
    # Silently passing here would leave every pre-guard summary looking current.
    root = _build_project(tmp_path, record_digests=False)
    with pytest.raises(ValueError, match=DIGEST_FIELD):
        stale_inputs(root)


def test_app_stops_when_the_precomputed_tables_are_stale(monkeypatch) -> None:
    monkeypatch.setattr(
        provenance,
        "stale_inputs",
        lambda root: ["src/etch_window/rsm.py"],
    )
    app = AppTest.from_file(str(PROJECT_ROOT / "app.py"), default_timeout=30).run()
    assert not app.exception
    assert any("rsm.py" in message.value for message in app.error)
    # st.stop() must run before the tabs render, or the mismatched numbers are
    # on screen anyway and the message is decoration.
    assert len(app.tabs) == 0


def test_app_renders_normally_when_nothing_is_stale() -> None:
    app = AppTest.from_file(str(PROJECT_ROOT / "app.py"), default_timeout=30).run()
    assert not app.exception
    assert len(app.error) == 0
    assert len(app.tabs) == 4
