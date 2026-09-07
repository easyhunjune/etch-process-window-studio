from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from etch_window.anova import (  # noqa: E402
    factor_sensitivity,
    lack_of_fit_anova,
    overall_anova,
    partial_term_anova,
)
from etch_window.data_access import RESPONSE_COLUMNS  # noqa: E402
from etch_window.design import (  # noqa: E402
    COEFFICIENT_NAMES,
    FACTOR_NAMES,
    FACTOR_SPECS,
)
from etch_window.provenance import current_input_digests  # noqa: E402
from etch_window.rsm import FullQuadraticRSM  # noqa: E402
from etch_window.ui_support import process_window_figure  # noqa: E402
from etch_window.window import ProcessSpecs, evaluate_process_window  # noqa: E402

RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
REPORT_DIR = PROJECT_ROOT / "reports"
FIGURE_DIR = REPORT_DIR / "figures"

RESPONSE_LABELS = {
    "etch_rate": "Si etch rate (µm/min)",
    "selectivity": "Si/SiO₂ selectivity",
    "nonuniformity": "Etch-rate non-uniformity (%)",
    "anisotropy": "Anisotropy",
    "bias_voltage": "DC bias voltage (V)",
}
PAPER_RESPONSE_COLUMNS = {
    "etch_rate": "si_etch_rate_um_min",
    "selectivity": "si_sio2_selectivity",
    "anisotropy": "anisotropy",
    "bias_voltage": "bias_voltage_v",
}
COEFFICIENT_TOLERANCE = {
    "etch_rate": 0.0015,
    "selectivity": 0.015,
    "anisotropy": 0.0015,
    "bias_voltage": 0.75,
}
KNOWN_COEFFICIENT_ISSUES = {
    ("etch_rate", "b3"): "전사 적합과 논문 부호 불일치",
    ("selectivity", "b3"): "전사 적합과 논문 부호 불일치",
    ("bias_voltage", "b34"): "전사 적합과 논문 부호 불일치",
}


def fit_models(frame: pd.DataFrame) -> tuple[np.ndarray, dict[str, FullQuadraticRSM]]:
    coded = frame[[f"{name}_code" for name in FACTOR_NAMES]].to_numpy(float)
    models = {
        name: FullQuadraticRSM.fit(coded, frame[column].to_numpy(float))
        for name, column in RESPONSE_COLUMNS.items()
    }
    return coded, models


def coefficient_tables(
    models: dict[str, FullQuadraticRSM],
    paper: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    coefficient_rows = []
    comparison_rows = []
    paper_terms = paper[paper["paper_symbol"].isin(COEFFICIENT_NAMES)].copy()
    for response, model in models.items():
        for name, value in model.coefficient_dict().items():
            coefficient_rows.append(
                {
                    "response": response,
                    "paper_symbol": name,
                    "fitted_coefficient": value,
                    "source": (
                        "project_fit_only"
                        if response == "nonuniformity"
                        else "fit_from_transcribed_table"
                    ),
                }
            )

    for response, paper_column in PAPER_RESPONSE_COLUMNS.items():
        fitted = models[response].coefficient_dict()
        tolerance = COEFFICIENT_TOLERANCE[response]
        for record in paper_terms.to_dict("records"):
            symbol = record["paper_symbol"]
            paper_value = float(record[paper_column])
            delta = fitted[symbol] - paper_value
            comparison_rows.append(
                {
                    "response": response,
                    "paper_symbol": symbol,
                    "fitted_coefficient": fitted[symbol],
                    "paper_coefficient": paper_value,
                    "delta": delta,
                    "absolute_delta": abs(delta),
                    "rounding_tolerance": tolerance,
                    "within_rounding_tolerance": abs(delta) <= tolerance,
                    "note": KNOWN_COEFFICIENT_ISSUES.get((response, symbol), ""),
                }
            )
    return pd.DataFrame(coefficient_rows), pd.DataFrame(comparison_rows)


def metric_table(
    models: dict[str, FullQuadraticRSM],
    paper: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    paper_s = paper.loc[paper["paper_symbol"] == "s"].iloc[0]
    paper_r2 = paper.loc[paper["paper_symbol"] == "R2"].iloc[0]
    for response, model in models.items():
        paper_column = PAPER_RESPONSE_COLUMNS.get(response)
        printed_s = float(paper_s[paper_column]) if paper_column else np.nan
        printed_r2 = float(paper_r2[paper_column]) if paper_column else np.nan
        rows.append(
            {
                "response": response,
                "n": model.n_observations,
                "rank": model.rank,
                "df_model": model.df_model,
                "df_residual": model.df_residual,
                "r_squared": model.r_squared,
                "adjusted_r_squared": model.adjusted_r_squared,
                "residual_standard_error": model.residual_standard_error,
                "paper_r_squared": printed_r2,
                "r_squared_delta": (
                    model.r_squared - printed_r2 if paper_column else np.nan
                ),
                "paper_standard_deviation": printed_s,
                "standard_deviation_delta": (
                    model.residual_standard_error - printed_s
                    if paper_column
                    else np.nan
                ),
                "model_status": (
                    "exploratory_low_adjusted_r2"
                    if response == "nonuniformity"
                    and model.adjusted_r_squared < 0.5
                    else (
                        "auxiliary_transcription_mismatch"
                        if response == "anisotropy"
                        else (
                            "audit_only"
                            if response == "bias_voltage"
                            else "primary"
                        )
                    )
                ),
            }
        )
    return pd.DataFrame(rows)


def anova_tables(
    coded: np.ndarray,
    frame: pd.DataFrame,
    models: dict[str, FullQuadraticRSM],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    overall_rows = []
    lack_of_fit_rows = []
    term_rows = []
    factor_rows = []
    for response, model in models.items():
        y = frame[RESPONSE_COLUMNS[response]].to_numpy(float)
        overall_rows.append({"response": response, **overall_anova(model)})
        lack_of_fit_rows.append(
            {
                "response": response,
                **lack_of_fit_anova(coded, y, model),
            }
        )
        term_rows.extend(
            {"response": response, **row}
            for row in partial_term_anova(coded, y, model)
        )
        factor_rows.extend(
            {"response": response, **row}
            for row in factor_sensitivity(coded, y, model)
        )
    return (
        pd.DataFrame(overall_rows),
        pd.DataFrame(lack_of_fit_rows),
        pd.DataFrame(term_rows),
        pd.DataFrame(factor_rows),
    )


def window_table(
    result: dict[str, object],
) -> pd.DataFrame:
    actual = np.asarray(result["actual"])
    shape = tuple(result["shape"])
    rows = {
        FACTOR_SPECS[name].name: actual[:, index]
        for index, name in enumerate(FACTOR_NAMES)
    }
    rows.update(
        {
            "predicted_etch_rate_um_min": np.asarray(result["etch_rate"]).ravel(),
            "predicted_selectivity": np.asarray(result["selectivity"]).ravel(),
            "predicted_nonuniformity_pct": np.asarray(
                result["nonuniformity"]
            ).ravel(),
            "predicted_anisotropy_auxiliary": np.asarray(
                result["anisotropy"]
            ).ravel(),
            "nominally_meets_primary_specs": np.asarray(
                result["nominal_feasible"]
            ).ravel(),
            "prediction_confidence_state": np.asarray(
                result["confidence_state"]
            ).ravel(),
            # Prediction-interval state alone says nothing about DOE support, so
            # a reader who takes only that column would treat an unsupported
            # point as confirmed. This column carries the pairing itself.
            "supported_prediction_confidence_state": np.asarray(
                result["supported_confidence_state"]
            ).ravel(),
            "inside_design_hull": np.asarray(
                result["inside_design_hull"]
            ).ravel(),
            "inside_leverage_envelope": np.asarray(
                result["inside_leverage_envelope"]
            ).ravel(),
            "inside_design_region": np.asarray(
                result["inside_design_region"]
            ).ravel(),
            "prediction_leverage": np.asarray(
                result["prediction_leverage"]
            ).ravel(),
            "meets_all_primary_specs": np.asarray(result["feasible"]).ravel(),
            "grid_row": np.repeat(np.arange(shape[0]), shape[1]),
            "grid_column": np.tile(np.arange(shape[1]), shape[0]),
        }
    )
    return pd.DataFrame(rows)


def plot_process_window(
    result: dict[str, object],
    specs: ProcessSpecs,
    destination: Path,
) -> None:
    figure = process_window_figure(result, specs)
    figure.suptitle(
        "Center slice: CHF₃=12 sccm, pressure=100 mTorr, RF power=100 W",
        fontsize=13,
    )
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def plot_factor_sensitivity(factors: pd.DataFrame, destination: Path) -> None:
    responses = ["etch_rate", "selectivity", "nonuniformity"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for axis, response in zip(axes, responses, strict=True):
        subset = factors[factors["response"] == response]
        axis.bar(
            subset["factor"],
            subset["normalized_partial_ss_pct"],
            color="#3478c0",
        )
        axis.set_title(RESPONSE_LABELS[response])
        axis.set_ylim(0, max(50, subset["normalized_partial_ss_pct"].max() * 1.15))
        axis.tick_params(axis="x", rotation=45)
        axis.set_ylabel("Normalized partial SS (%)")
    fig.suptitle(
        "Factor-group deletion sensitivity (interactions overlap by definition)"
    )
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def plot_main_effects(
    models: dict[str, FullQuadraticRSM],
    destination: Path,
) -> None:
    responses = ["etch_rate", "selectivity", "nonuniformity"]
    coded_axis = np.linspace(-2.0, 2.0, 201)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for axis, response in zip(axes, responses, strict=True):
        for factor_index, factor_name in enumerate(FACTOR_NAMES):
            center_line = np.zeros((coded_axis.size, len(FACTOR_NAMES)))
            center_line[:, factor_index] = coded_axis
            axis.plot(
                coded_axis,
                models[response].predict(center_line),
                label=factor_name,
                linewidth=2,
            )
        axis.axvline(0.0, color="#888888", linewidth=0.8, linestyle="--")
        axis.set_title(RESPONSE_LABELS[response])
        axis.set_xlabel("Coded factor level (others fixed at 0)")
        axis.set_ylabel("Predicted response")
        axis.grid(alpha=0.2)
    axes[-1].legend(title="Varied factor", loc="best")
    fig.suptitle("Center-line main effects within the DOE range")
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def write_report(
    metrics: pd.DataFrame,
    comparison: pd.DataFrame,
    factors: pd.DataFrame,
    specs: ProcessSpecs,
    feasible_fraction: float,
    nominal_fraction: float,
    design_region_fraction: float,
) -> None:
    metric_lookup = metrics.set_index("response")
    discrepancy_count = int((~comparison["within_rounding_tolerance"]).sum())
    discrepancy_breakdown = (
        comparison.loc[~comparison["within_rounding_tolerance"]]
        .groupby("response")
        .size()
        .to_dict()
    )
    top_factors = {}
    for response in ("etch_rate", "selectivity", "nonuniformity"):
        subset = factors[factors["response"] == response]
        row = subset.loc[subset["normalized_partial_ss_pct"].idxmax()]
        top_factors[response] = (
            row["factor"],
            float(row["normalized_partial_ss_pct"]),
        )

    def metric_row(
        response: str,
        label: str,
        paper_value: str,
        status: str,
    ) -> str:
        values = metric_lookup.loc[response]
        return (
            f"| {label} | {values['r_squared']:.4f} | "
            f"{values['adjusted_r_squared']:.4f} | "
            f"{values['residual_standard_error']:.4f} | "
            f"{paper_value} | {status} |"
        )

    metric_rows = "\n".join(
        (
            metric_row(
                "etch_rate",
                "Etch rate",
                f"{metric_lookup.loc['etch_rate', 'paper_r_squared']:.2f}",
                "주 모델",
            ),
            metric_row(
                "selectivity",
                "Selectivity",
                f"{metric_lookup.loc['selectivity', 'paper_r_squared']:.2f}",
                "주 모델",
            ),
            metric_row("nonuniformity", "Non-uniformity", "—", "탐색 모델"),
            metric_row(
                "anisotropy",
                "Anisotropy",
                f"{metric_lookup.loc['anisotropy', 'paper_r_squared']:.2f}",
                "보조 모델",
            ),
            metric_row(
                "bias_voltage",
                "DC bias",
                f"{metric_lookup.loc['bias_voltage', 'paper_r_squared']:.2f}",
                "감사용",
            ),
        )
    )
    factor_lines = "\n".join(
        f"- {label} 최대 기여 그룹: `{top_factors[key][0]}` "
        f"({top_factors[key][1]:.1f}%)"
        for key, label in (
            ("etch_rate", "Etch rate"),
            ("selectivity", "Selectivity"),
            ("nonuniformity", "Non-uniformity"),
        )
    )

    report = f"""# Project 2 Phase B — RSM 수치 코어 검증

## 완료 범위

- 논문과 같은 coded factor 순서(SF6, O2, CHF3, pressure, RF power)로
  21항 full-quadratic OLS를 구현했다.
- 32행 CCD의 5개 응답을 적합하고 R², adjusted R², 잔차표준오차, 전체
  ANOVA, 항별 부분 제곱합, 요인 그룹 민감도를 산출했다.
- etch rate, selectivity, non-uniformity의 세 조건을 동시에 판정하는
  2차원 process-window 엔진을 구현했다.
- 인자별 범위에 더해 CCD 볼록껍질과 관측 레버리지 한계를 적용하고,
  점별 95% 예측구간으로 확정·불확실·탈락을 구분한다.
- 중심점 6회 반복으로 pure error와 lack-of-fit을 분리한다.

## 모델 재현 결과

| 응답 | R² | adjusted R² | 잔차표준오차 | 논문 R² | 판정 |
|---|---:|---:|---:|---:|---|
{metric_rows}

논문 계수와 반올림 허용오차를 벗어난 대조 항은 {discrepancy_count}개다.
응답별로 etch rate {discrepancy_breakdown.get("etch_rate", 0)}개,
selectivity {discrepancy_breakdown.get("selectivity", 0)}개,
anisotropy {discrepancy_breakdown.get("anisotropy", 0)}개,
DC bias {discrepancy_breakdown.get("bias_voltage", 0)}개다.
Phase A에서 확인한 부호 불일치 외에도 anisotropy 전사값은 논문의 계수와
R²를 충분히 재현하지 못한다. 특히 etch-rate R²는 재현되지만 전사값 기반
잔차표준오차는 {metric_lookup.loc["etch_rate", "residual_standard_error"]:.4f}로,
논문 표의 0.01과 다르다. 원자료를 바꾸지 않고 이 차이를 감사 결과로 보존했다.

non-uniformity 모델의 adjusted R²는
{metric_lookup.loc["nonuniformity", "adjusted_r_squared"]:.3f}로 낮다.
따라서 이 응답을 포함한 process window는 생산 레시피 추천이 아니라
포트폴리오용 탐색·설명 기능으로만 해석해야 한다.

## 민감도와 ANOVA 해석

{factor_lines}

기여율은 각 요인에 속한 선형항·제곱항·상호작용항을 full model에서 함께
제거했을 때 늘어나는 부분 제곱합을 정규화한 값이다. 상호작용항이 양쪽
요인 그룹에 포함되므로 고전적인 가산형 ANOVA 분해와 같지 않다.

원 논문에는 별도의 ANOVA 표가 없고 21개 회귀계수, 모델 표준편차, R²만
제시되어 있다. 따라서 ANOVA는 전사 데이터로 계산한 프로젝트 산출물이며,
논문과의 직접 대조는 계수·R² 표에서 수행했다.

## 데모 process window

- 고정값: CHF3 12 sccm, pressure 100 mTorr, RF power 100 W
- 탐색축: SF6 10–50 sccm, O2 2–18 sccm
- 예시 spec: etch rate ≥ {specs.min_etch_rate_um_min:g} µm/min,
  selectivity ≥ {specs.min_selectivity:g},
  non-uniformity ≤ {specs.max_nonuniformity_pct:g}%
- 평균예측 충족률: {100 * nominal_fraction:.1f}%
- 설계 지지영역 비율: {100 * design_region_fraction:.1f}%
- 95% 예측구간 확정률: {100 * feasible_fraction:.1f}%

예시 spec은 소프트웨어 동작을 재현하기 위한 데모 기준이며 논문이나 양산
승인 기준이 아니다. 사용자가 spec을 바꾸면 동일한 판정 엔진이 즉시 다시
계산할 수 있다.

## 산출물

- `data/processed/model_coefficients.csv`
- `data/processed/model_metrics.csv`
- `data/processed/paper_coefficient_comparison.csv`
- `data/processed/anova_overall.csv`
- `data/processed/anova_lack_of_fit.csv`
- `data/processed/anova_terms.csv`
- `data/processed/factor_sensitivity.csv`
- `data/processed/process_window_grid.csv`
  - `prediction_confidence_state`는 예측구간만으로 판정한 값이라 설계 지지영역
    밖에서도 `confirmed`가 될 수 있다. 이 컬럼만 읽는 소비자를 위해
    `supported_prediction_confidence_state`가 지지영역 밖을 `unsupported`로
    표시한다. `unsupported`가 `rejected`보다 우선하므로 지지영역 밖에서
    예측구간상 탈락한 점도 이 컬럼에서는 `unsupported`로 보인다. 확정 탈락
    개수는 `prediction_confidence_state`에서 세어야 한다.
- `data/processed/phase_b_summary.json`
- `reports/figures/process_window_sf6_o2.png`
- `reports/figures/factor_contribution.png`
- `reports/figures/main_effects.png`

## 재현

```powershell
$env:PYTHONPATH = "src"
python scripts/run_phase_b.py
python -m pytest -q
```
"""
    (REPORT_DIR / "PHASE_B_NUMERICAL_CORE.md").write_text(
        report, encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=81)
    arguments = parser.parse_args()
    for directory in (PROCESSED_DIR, REPORT_DIR, FIGURE_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(RAW_DIR / "legtenberg_1995_ccd.csv")
    paper = pd.read_csv(RAW_DIR / "legtenberg_1995_paper_coefficients.csv")
    coded, models = fit_models(frame)
    coefficients, comparison = coefficient_tables(models, paper)
    metrics = metric_table(models, paper)
    anova_overall, lack_of_fit, anova_terms, factors = anova_tables(
        coded,
        frame,
        models,
    )

    specs = ProcessSpecs(
        min_etch_rate_um_min=0.4,
        min_selectivity=8.0,
        max_nonuniformity_pct=4.5,
    )
    result = evaluate_process_window(
        models,
        specs,
        vary=("sf6", "o2"),
        resolution=arguments.resolution,
    )
    window = window_table(result)

    coefficients.to_csv(PROCESSED_DIR / "model_coefficients.csv", index=False)
    metrics.to_csv(PROCESSED_DIR / "model_metrics.csv", index=False)
    comparison.to_csv(
        PROCESSED_DIR / "paper_coefficient_comparison.csv", index=False
    )
    anova_overall.to_csv(PROCESSED_DIR / "anova_overall.csv", index=False)
    lack_of_fit.to_csv(PROCESSED_DIR / "anova_lack_of_fit.csv", index=False)
    anova_terms.to_csv(PROCESSED_DIR / "anova_terms.csv", index=False)
    factors.to_csv(PROCESSED_DIR / "factor_sensitivity.csv", index=False)
    window.to_csv(PROCESSED_DIR / "process_window_grid.csv", index=False)

    feasible_fraction = float(np.asarray(result["feasible"]).mean())
    nominal_fraction = float(np.asarray(result["nominal_feasible"]).mean())
    design_region_fraction = float(
        np.asarray(result["inside_design_region"]).mean()
    )
    summary = {
        "dataset_rows": len(frame),
        # Lets the app tell whether these tables still match the sources a live
        # refit would read. See etch_window/provenance.py.
        "input_digests": current_input_digests(PROJECT_ROOT),
        "design_rank": int(next(iter(models.values())).rank),
        "model_terms": len(COEFFICIENT_NAMES),
        "specifications": {
            "min_etch_rate_um_min": specs.min_etch_rate_um_min,
            "min_selectivity": specs.min_selectivity,
            "max_nonuniformity_pct": specs.max_nonuniformity_pct,
            "source": "demonstration_only_not_from_paper",
        },
        "window": {
            "vary": ["sf6", "o2"],
            "resolution": arguments.resolution,
            "feasible_grid_fraction": feasible_fraction,
            "nominal_grid_fraction": nominal_fraction,
            "design_region_grid_fraction": design_region_fraction,
            "definition": (
                "Pointwise 95% prediction intervals satisfy all three specs "
                "and the point lies inside both the DOE convex hull and "
                "observed leverage envelope."
            ),
            "fixed_actual": {
                "chf3_sccm": 12.0,
                "pressure_mtorr": 100.0,
                "rf_power_w": 100.0,
            },
        },
        "warnings": [
            "Non-uniformity adjusted R-squared is below 0.5.",
            "Anisotropy coefficients and R-squared do not reproduce the paper.",
            "Etch-rate R-squared reproduces the paper but residual s does not.",
            "No prediction outside the original DOE bounds is allowed.",
            "Factor-wise bounds alone are insufficient; the process window "
            "also applies convex-hull and leverage support checks.",
            "Prediction states use pointwise, not simultaneous, 95% intervals.",
        ],
    }
    (PROCESSED_DIR / "phase_b_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    plot_process_window(
        result, specs, FIGURE_DIR / "process_window_sf6_o2.png"
    )
    plot_factor_sensitivity(
        factors, FIGURE_DIR / "factor_contribution.png"
    )
    plot_main_effects(models, FIGURE_DIR / "main_effects.png")
    write_report(
        metrics,
        comparison,
        factors,
        specs,
        feasible_fraction,
        nominal_fraction,
        design_region_fraction,
    )

    print(f"rows={len(frame)}, terms={len(COEFFICIENT_NAMES)}, rank=21")
    print(f"feasible_grid_fraction={feasible_fraction:.6f}")
    print(
        "coefficient_discrepancies="
        f"{int((~comparison['within_rounding_tolerance']).sum())}"
    )
    print(f"report={REPORT_DIR / 'PHASE_B_NUMERICAL_CORE.md'}")


if __name__ == "__main__":
    main()
