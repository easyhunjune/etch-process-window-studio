"""Digest bookkeeping that ties the precomputed tables to the inputs that made them.

The app reads two kinds of numbers side by side: the "Process window" tab refits
the raw DOE live, while the "논문 재현 대조" and "ANOVA" tabs read the CSVs
`scripts/run_phase_b.py` wrote earlier. Nothing forced those two to agree — edit
the raw CSV or the numeric core, skip the rebuild, and the tabs quietly disagree.
Recording the input digests at build time turns that into a visible mismatch.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

RAW_SOURCE_FILES = (
    "data/raw/legtenberg_1995_ccd.csv",
    "data/raw/legtenberg_1995_paper_coefficients.csv",
)
# The modules whose output lands in data/processed/. Data loading and response
# mapping live in data_access.py; ui_support.py and app.py only render tables and
# figures, so changing their display labels does not make the stored tables wrong.
NUMERIC_CORE_FILES = (
    "src/etch_window/data_access.py",
    "src/etch_window/design.py",
    "src/etch_window/rsm.py",
    "src/etch_window/anova.py",
    "src/etch_window/window.py",
    "scripts/run_phase_b.py",
)
TRACKED_INPUT_FILES = RAW_SOURCE_FILES + NUMERIC_CORE_FILES
SUMMARY_RELATIVE_PATH = "data/processed/phase_b_summary.json"
DIGEST_FIELD = "input_digests"
REBUILD_COMMAND = "python scripts/run_phase_b.py"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def current_input_digests(project_root: Path) -> dict[str, str]:
    """Hash every input the processed tables are derived from."""
    return {name: sha256_file(project_root / name) for name in TRACKED_INPUT_FILES}


def stale_inputs(project_root: Path) -> list[str]:
    """Return the tracked inputs that changed since the tables were last built.

    An empty list means the precomputed CSVs and a live refit read the same
    sources. Missing bookkeeping raises instead of returning an empty list: a
    guard that passes when it cannot check is not a guard.
    """
    summary_path = project_root / SUMMARY_RELATIVE_PATH
    if not summary_path.is_file():
        raise FileNotFoundError(
            f"{SUMMARY_RELATIVE_PATH} is missing. Run `{REBUILD_COMMAND}` first."
        )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    recorded = summary.get(DIGEST_FIELD)
    if not isinstance(recorded, dict) or not recorded:
        raise ValueError(
            f"{SUMMARY_RELATIVE_PATH} has no `{DIGEST_FIELD}` block, so the "
            "precomputed tables cannot be matched to their inputs. Run "
            f"`{REBUILD_COMMAND}` to record it."
        )
    current = current_input_digests(project_root)
    return sorted(
        name
        for name in set(recorded) | set(current)
        if recorded.get(name) != current.get(name)
    )
