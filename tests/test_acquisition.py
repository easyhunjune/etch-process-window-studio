from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from scripts import validate_acquisition

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"


def test_require_raises_instead_of_asserting() -> None:
    validate_acquisition.require(True, "never raised")
    with pytest.raises(validate_acquisition.AcquisitionError, match="row count"):
        validate_acquisition.require(False, "bad row count")


def test_primary_validation_catches_a_broken_coded_to_actual_mapping(
    tmp_path,
    monkeypatch,
) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    rows = (RAW_DIR / "legtenberg_1995_ccd.csv").read_text(encoding="utf-8")
    # Run 1 sits at sf6 coded -1, which must read 20 sccm.
    broken = rows.replace(",20,", ",25,", 1)
    assert broken != rows
    (raw / "legtenberg_1995_ccd.csv").write_text(broken, encoding="utf-8")
    (raw / "legtenberg_1995_paper_coefficients.csv").write_text(
        (RAW_DIR / "legtenberg_1995_paper_coefficients.csv").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(validate_acquisition, "DATA_DIR", raw)
    with pytest.raises(
        validate_acquisition.AcquisitionError,
        match=r"Run 1 sf6: coded -1 maps to 20, but the table records 25",
    ):
        validate_acquisition.validate_primary()


def test_validation_still_fires_under_python_dash_o(tmp_path) -> None:
    """`python -O` strips `assert`, so the checks must not be assertions."""
    raw = tmp_path / "raw"
    raw.mkdir()
    lines = (RAW_DIR / "legtenberg_1995_ccd.csv").read_text(
        encoding="utf-8"
    ).splitlines(keepends=True)
    (raw / "legtenberg_1995_ccd.csv").write_text(
        "".join(lines[:-1]), encoding="utf-8"
    )
    program = textwrap.dedent(
        f"""
        import sys
        from pathlib import Path
        sys.path.insert(0, {str(PROJECT_ROOT)!r})
        from scripts import validate_acquisition as module
        module.DATA_DIR = Path({str(raw)!r})
        module.validate_primary()
        """
    )
    result = subprocess.run(
        [sys.executable, "-O", "-c", program],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "32 runs" in result.stderr


def test_manifest_allows_missing_private_pdfs_unless_required(
    tmp_path,
    monkeypatch,
) -> None:
    data = tmp_path / "data"
    data.mkdir()
    manifest = {
        "sources": [
            {
                "id": "paper",
                "local_pdf": "papers/source.pdf",
                "sha256": "abc123",
            }
        ]
    }
    (data / "source_manifest.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    monkeypatch.setattr(validate_acquisition, "PROJECT_ROOT", tmp_path)

    audit = validate_acquisition.validate_manifest()
    assert audit["paper"]["status"] == "not_present"
    assert audit["paper"]["expected_sha256"] == "abc123"

    with pytest.raises(FileNotFoundError, match="required but missing"):
        validate_acquisition.validate_manifest(require_pdfs=True)
