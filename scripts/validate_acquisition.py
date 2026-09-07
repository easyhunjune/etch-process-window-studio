from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data" / "raw"


class AcquisitionError(ValueError):
    """Raised when the transcribed source data fails a structural check."""


def require(condition: bool, message: str) -> None:
    """Enforce a check that must survive `python -O`.

    These checks are the whole point of this script, and `assert` statements
    vanish under `-O`, which would leave it printing ``"status": "PASS"`` without
    having validated anything.
    """
    if not condition:
        raise AcquisitionError(message)


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA_DIR / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_primary() -> dict[str, object]:
    rows = read_csv("legtenberg_1995_ccd.csv")
    require(len(rows) == 32, f"CCD table must hold 32 runs, found {len(rows)}.")
    run_numbers = {int(row["run"]) for row in rows}
    require(
        run_numbers == set(range(1, 33)),
        f"Run numbers must be 1-32, found {sorted(run_numbers)}.",
    )

    mappings = {
        "sf6": {-2: 10, -1: 20, 0: 30, 1: 40, 2: 50},
        "o2": {-2: 2, -1: 6, 0: 10, 1: 14, 2: 18},
        "chf3": {-2: 2, -1: 7, 0: 12, 1: 17, 2: 22},
        "pressure": {-2: 20, -1: 60, 0: 100, 1: 140, 2: 180},
        "rf_power": {-2: 20, -1: 60, 0: 100, 1: 140, 2: 180},
    }
    actual_columns = {
        "sf6": "sf6_sccm",
        "o2": "o2_sccm",
        "chf3": "chf3_sccm",
        "pressure": "pressure_mtorr",
        "rf_power": "rf_power_w",
    }
    for row in rows:
        for factor, mapping in mappings.items():
            coded = int(row[f"{factor}_code"])
            actual = int(row[actual_columns[factor]])
            require(
                mapping[coded] == actual,
                f"Run {row['run']} {factor}: coded {coded:+d} maps to "
                f"{mapping[coded]}, but the table records {actual}.",
            )

    centers = [
        row
        for row in rows
        if all(
            int(row[f"{factor}_code"]) == 0
            for factor in ("sf6", "o2", "chf3", "pressure", "rf_power")
        )
    ]
    center_runs = [int(row["run"]) for row in centers]
    require(
        center_runs == [27, 28, 29, 30, 31, 32],
        f"Center replicates must be runs 27-32, found {center_runs}.",
    )
    surfaces = Counter(row["surface"] for row in rows)
    require(
        surfaces == {"smooth": 24, "rough": 8},
        f"Surface counts must be 24 smooth / 8 rough, found {dict(surfaces)}.",
    )
    outward = sum(row["anisotropy_outward_slope"] == "true" for row in rows)
    require(
        outward == 8,
        f"Table 4.3 marks 8 outward-sloped profiles, found {outward}.",
    )

    coefficients = read_csv("legtenberg_1995_paper_coefficients.csv")
    require(
        len(coefficients) == 23,
        f"Coefficient table must hold 21 terms plus s and R2, found "
        f"{len(coefficients)} rows.",
    )
    require(
        coefficients[-1]["paper_symbol"] == "R2",
        "The coefficient table must end with the R2 row, found "
        f"{coefficients[-1]['paper_symbol']!r}.",
    )
    return {
        "runs": len(rows),
        "center_replicates": len(centers),
        "responses": [
            "si_etch_rate_um_min",
            "si_sio2_selectivity",
            "anisotropy",
            "etch_rate_nonuniformity_pct",
        ],
        "paper_coefficient_rows": len(coefficients),
    }


def rounded_level_means(
    rows: list[dict[str, str]], factor: str, response: str
) -> dict[int, float]:
    values: dict[int, list[float]] = defaultdict(list)
    for row in rows:
        values[int(row[factor])].append(float(row[response]))
    return {level: round(sum(items) / len(items), 1) for level, items in values.items()}


def validate_rejected_taguchi_candidate() -> dict[str, object]:
    rows = read_csv("muttalib_2017_l9_audit.csv")
    require(len(rows) == 9, f"L9 array must hold 9 runs, found {len(rows)}.")

    printed_pairs = {
        (int(row["gas_ratio_pct"]), int(row["pressure_mtorr_as_printed"]))
        for row in rows
    }
    # The printed table collapses 9 runs onto 3 distinct (gas, pressure) pairs.
    # That collapse is the reason this source was rejected, so losing the check
    # would lose the documented rejection ground.
    require(
        len(printed_pairs) == 3,
        f"Printed table should collapse to 3 distinct pairs, found "
        f"{len(printed_pairs)}.",
    )

    inferred_pairs = {
        (int(row["gas_ratio_pct"]), int(row["pressure_mtorr_inferred_from_table4"]))
        for row in rows
    }
    require(
        len(inferred_pairs) == 9,
        f"Table 4 implies 9 distinct pairs, found {len(inferred_pairs)}.",
    )

    inferred_rate_means = rounded_level_means(
        rows, "pressure_mtorr_inferred_from_table4", "ta2o5_etch_rate_nm_min"
    )
    inferred_selectivity_means = rounded_level_means(
        rows, "pressure_mtorr_inferred_from_table4", "ta2o5_cr_selectivity"
    )
    require(
        inferred_rate_means == {15: 12.7, 30: 14.5, 45: 12.5},
        f"Inferred etch-rate level means changed: {inferred_rate_means}.",
    )
    require(
        inferred_selectivity_means == {15: 17.6, 30: 17.7, 45: 7.5},
        f"Inferred selectivity level means changed: {inferred_selectivity_means}.",
    )
    return {
        "runs": len(rows),
        "printed_unique_gas_pressure_pairs": len(printed_pairs),
        "inferred_unique_gas_pressure_pairs": len(inferred_pairs),
        "status": "rejected_due_to_internal_table_conflict_and_no_anova_or_regression",
    }


def validate_manifest(*, require_pdfs: bool = False) -> dict[str, dict[str, object]]:
    manifest_path = PROJECT_ROOT / "data" / "source_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    audit: dict[str, dict[str, object]] = {}
    for source in manifest["sources"]:
        pdf_path = PROJECT_ROOT / source["local_pdf"]
        if not pdf_path.is_file():
            if require_pdfs:
                raise FileNotFoundError(
                    f"Source PDF is required but missing: {source['local_pdf']}"
                )
            audit[source["id"]] = {
                "status": "not_present",
                "path": source["local_pdf"],
                "expected_sha256": source["sha256"],
            }
            continue
        actual = sha256(pdf_path)
        require(
            actual == source["sha256"],
            f"{source['local_pdf']} sha256 is {actual}, manifest expects "
            f"{source['sha256']}.",
        )
        audit[source["id"]] = {
            "status": "verified",
            "path": source["local_pdf"],
            "sha256": actual,
        }
    return audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--require-pdfs",
        action="store_true",
        help="Fail when local source PDFs are absent; public clones omit them.",
    )
    arguments = parser.parse_args()
    result = {
        "primary": validate_primary(),
        "rejected_taguchi_candidate": validate_rejected_taguchi_candidate(),
        "pdf_audit": validate_manifest(require_pdfs=arguments.require_pdfs),
        "status": "PASS",
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
