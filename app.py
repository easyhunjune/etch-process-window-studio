from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from etch_window.design import FACTOR_NAMES, FACTOR_SPECS  # noqa: E402
from etch_window.provenance import (  # noqa: E402
    REBUILD_COMMAND,
    SUMMARY_RELATIVE_PATH,
    sha256_file,
    stale_inputs,
)
from etch_window.ui_support import (  # noqa: E402
    RESPONSE_LABELS,
    factor_contribution_figure,
    load_doe_and_models,
    main_effect_figure,
    process_window_figure,
    process_window_frame,
)
from etch_window.window import ProcessSpecs, evaluate_process_window  # noqa: E402

CITATION = (
    "R. Legtenberg et al., “Anisotropic Reactive Ion Etching of Silicon "
    "Using SF₆/O₂/CHF₃ Gas Mixtures,” J. Electrochem. Soc. 142(6), "
    "2020–2028 (1995)"
)


@st.cache_resource
def load_models(summary_digest: str) -> tuple[pd.DataFrame, dict]:
    del summary_digest  # Included only so a rebuilt summary invalidates the cache.
    return load_doe_and_models(PROJECT_ROOT)


@st.cache_data
def load_validation_tables(summary_digest: str) -> dict[str, pd.DataFrame]:
    del summary_digest  # Included only so a rebuilt summary invalidates the cache.
    processed = PROJECT_ROOT / "data" / "processed"
    names = (
        "model_metrics",
        "paper_coefficient_comparison",
        "anova_overall",
        "anova_lack_of_fit",
        "anova_terms",
        "factor_sensitivity",
    )
    return {name: pd.read_csv(processed / f"{name}.csv") for name in names}


def check_precomputed_tables_are_current() -> None:
    """Stop before showing live and precomputed numbers that came from different sources.

    "Process window" refits the raw DOE on every interaction; the other tabs read
    CSVs from an earlier `run_phase_b.py`. Rendering both without this check is how
    the two silently drift apart.

    Deliberately uncached: the seven tracked files hash in well under a
    millisecond, and a cached staleness check would keep reporting "current"
    after the very edit it exists to catch.
    """
    try:
        stale = stale_inputs(PROJECT_ROOT)
    except (FileNotFoundError, ValueError) as error:
        st.error(str(error))
        st.stop()
    if stale:
        st.error(
            "사전 계산 표가 현재 입력과 어긋납니다. 아래 파일이 마지막 "
            f"`{REBUILD_COMMAND}` 실행 이후 바뀌었으므로, 재실행 전까지 "
            "'논문 재현 대조'·'ANOVA' 탭은 'Process window' 탭과 다른 "
            "버전을 보여 줍니다.\n\n"
            + "\n".join(f"- `{name}`" for name in stale)
        )
        st.stop()


def format_factor(name: str) -> str:
    spec = FACTOR_SPECS[name]
    return f"{spec.display_name} ({spec.unit})"


def sidebar_controls() -> tuple[ProcessSpecs, tuple[str, str], dict[str, float], int]:
    st.sidebar.header("공정창 제어")
    st.sidebar.caption("세 조건을 바꾸면 모든 contour와 충족 영역이 갱신됩니다.")
    min_rate = st.sidebar.number_input(
        "최소 Si 식각 속도 (µm/min)",
        min_value=0.0,
        max_value=1.0,
        value=0.40,
        step=0.01,
        key="spec_min_rate",
    )
    min_selectivity = st.sidebar.number_input(
        "최소 Si/SiO₂ 선택비",
        min_value=0.0,
        max_value=25.0,
        value=8.0,
        step=0.5,
        key="spec_min_selectivity",
    )
    max_nonuniformity = st.sidebar.number_input(
        "최대 비균일도 (%)",
        min_value=0.0,
        max_value=15.0,
        value=4.5,
        step=0.1,
        key="spec_max_nonuniformity",
    )
    st.sidebar.divider()
    st.sidebar.subheader("탐색 평면")
    first = st.sidebar.selectbox(
        "가로축",
        FACTOR_NAMES,
        index=0,
        format_func=format_factor,
        key="first_factor",
    )
    second_options = tuple(name for name in FACTOR_NAMES if name != first)
    default_second = second_options.index("o2") if "o2" in second_options else 0
    second = st.sidebar.selectbox(
        "세로축",
        second_options,
        index=default_second,
        format_func=format_factor,
        key="second_factor",
    )

    fixed_actual: dict[str, float] = {}
    st.sidebar.subheader("고정 인자")
    for name in FACTOR_NAMES:
        if name in (first, second):
            continue
        factor = FACTOR_SPECS[name]
        fixed_actual[name] = st.sidebar.slider(
            format_factor(name),
            min_value=float(factor.actual_min),
            max_value=float(factor.actual_max),
            value=float(factor.center),
            step=float(factor.step / 10),
            key=f"fixed_{name}",
        )
    resolution = st.sidebar.select_slider(
        "격자 해상도",
        options=(41, 61, 81, 101, 121),
        value=81,
        help="높을수록 경계가 부드럽지만 계산량이 늘어납니다.",
    )
    specs = ProcessSpecs(
        min_etch_rate_um_min=min_rate,
        min_selectivity=min_selectivity,
        max_nonuniformity_pct=max_nonuniformity,
    )
    return specs, (first, second), fixed_actual, resolution


def process_window_tab(
    models: dict,
    specs: ProcessSpecs,
    vary: tuple[str, str],
    fixed_actual: dict[str, float],
    resolution: int,
) -> None:
    result = evaluate_process_window(
        models,
        specs,
        vary=vary,
        fixed_actual=fixed_actual,
        resolution=resolution,
    )
    feasible = np.asarray(result["feasible"], dtype=bool)
    nominal = np.asarray(result["nominal_feasible"], dtype=bool)
    inside = np.asarray(result["inside_design_region"], dtype=bool)
    uncertain = (
        np.asarray(result["confidence_state"]) == "uncertain"
    ) & inside
    columns = st.columns(5)
    columns[0].metric("95% PI 확정률", f"{100 * feasible.mean():.1f}%")
    columns[1].metric("평균예측 충족률", f"{100 * nominal.mean():.1f}%")
    columns[2].metric("불확실 격자", f"{uncertain.sum():,}")
    columns[3].metric("설계영역 내부", f"{100 * inside.mean():.1f}%")
    columns[4].metric("전체 격자", f"{feasible.size:,}")
    if feasible.any():
        st.success(
            "설계영역 내부에서 세 조건의 95% 예측구간을 모두 만족하는 "
            "확정 영역이 있습니다."
        )
    else:
        st.warning(
            "평균예측 충족점이 있더라도 95% 예측구간까지 보수적으로 "
            "통과하는 확정 격자는 없습니다."
        )

    figure = process_window_figure(result, specs)
    st.pyplot(figure, width="stretch")
    plt.close(figure)
    st.caption(
        "각 응답 패널의 흰 선은 그 응답의 평균예측이 자기 기준값과 만나는 경계입니다. "
        "세 기준을 동시에 만족하는 영역은 위의 평균예측 충족률로 확인하세요. "
        "종합 패널에서 초록은 95% 예측구간 확정, 노랑은 불확실, "
        "회색은 예측구간이 기준을 명백히 못 만족하는 확정 탈락(rejected), "
        "사선은 DOE 볼록껍질·레버리지 지지영역 밖입니다. "
        "구간은 점별(pointwise) 95%이며 동시신뢰영역은 아닙니다."
    )

    result_frame = process_window_frame(result)
    feasible_frame = result_frame[result_frame["meets_all_specs"]]
    st.download_button(
        "충족 격자 CSV 다운로드",
        feasible_frame.to_csv(index=False).encode("utf-8-sig"),
        file_name="feasible_process_window.csv",
        mime="text/csv",
        disabled=feasible_frame.empty,
    )


def reproduction_tab(tables: dict[str, pd.DataFrame]) -> None:
    metrics = tables["model_metrics"].copy()
    display_metrics = metrics[
        [
            "response",
            "r_squared",
            "adjusted_r_squared",
            "residual_standard_error",
            "paper_r_squared",
            "model_status",
        ]
    ].rename(
        columns={
            "response": "응답",
            "r_squared": "내 적합 R²",
            "adjusted_r_squared": "adjusted R²",
            "residual_standard_error": "잔차표준오차",
            "paper_r_squared": "논문 R²",
            "model_status": "사용 상태",
        }
    )
    st.subheader("모델 지표 대조")
    st.dataframe(
        display_metrics.style.format(
            {
                "내 적합 R²": "{:.4f}",
                "adjusted R²": "{:.4f}",
                "잔차표준오차": "{:.4f}",
                "논문 R²": "{:.2f}",
            },
            na_rep="—",
        ),
        width="stretch",
        hide_index=True,
    )

    comparison = tables["paper_coefficient_comparison"].copy()
    response = st.selectbox(
        "계수 대조 응답",
        ("etch_rate", "selectivity", "anisotropy", "bias_voltage"),
        format_func=lambda name: RESPONSE_LABELS[name],
        key="comparison_response",
    )
    subset = comparison[comparison["response"] == response].copy()
    subset["판정"] = np.where(
        subset["within_rounding_tolerance"], "일치", "확인 필요"
    )
    bad = int((~subset["within_rounding_tolerance"]).sum())
    st.metric("반올림 허용오차 밖 계수", f"{bad} / {len(subset)}")
    st.dataframe(
        subset[
            [
                "paper_symbol",
                "fitted_coefficient",
                "paper_coefficient",
                "delta",
                "판정",
                "note",
            ]
        ].rename(
            columns={
                "paper_symbol": "항",
                "fitted_coefficient": "내 적합",
                "paper_coefficient": "논문",
                "delta": "차이",
                "note": "감사 메모",
            }
        ),
        width="stretch",
        hide_index=True,
    )
    st.info(
        "원 논문에는 별도의 ANOVA 표가 없고 21개 회귀계수, 모델 "
        "표준편차와 R²가 제시되어 있습니다. 아래 ANOVA는 전사된 32행에서 "
        "계산한 프로젝트 결과이며, 논문과의 직접 대조는 위 표에서 수행합니다."
    )
    st.dataframe(
        tables["anova_overall"],
        width="stretch",
        hide_index=True,
    )
    st.subheader("Lack-of-fit 검정")
    st.dataframe(
        tables["anova_lack_of_fit"],
        width="stretch",
        hide_index=True,
    )
    st.caption(
        "중심점 6회 반복에서 pure error를 분리했습니다. 낮은 p값은 "
        "반올림·전사 오차를 포함해 2차식이 run 평균을 충분히 설명하지 "
        "못한다는 신호이지, 원인을 하나로 확정하는 증거는 아닙니다."
    )


def analysis_tab(models: dict, tables: dict[str, pd.DataFrame]) -> None:
    st.subheader("인자 그룹 기여도")
    st.caption(
        "선형·제곱·상호작용항을 요인별로 함께 제거한 부분 제곱합입니다. "
        "상호작용이 양쪽 그룹에 포함되므로 고전적인 가산형 분해는 아닙니다."
    )
    contribution = factor_contribution_figure(tables["factor_sensitivity"])
    st.pyplot(contribution, width="stretch")
    plt.close(contribution)
    st.dataframe(
        tables["factor_sensitivity"],
        width="stretch",
        hide_index=True,
    )

    st.subheader("중심선 main effect")
    st.caption("한 인자만 coded −2~+2로 변화시키고 나머지는 중심점에 고정했습니다.")
    effects = main_effect_figure(models)
    st.pyplot(effects, width="stretch")
    plt.close(effects)

    st.subheader("항별 유의성")
    term_response = st.selectbox(
        "항별 ANOVA 응답",
        ("etch_rate", "selectivity", "nonuniformity"),
        format_func=lambda name: RESPONSE_LABELS[name],
        key="anova_response",
    )
    terms = tables["anova_terms"]
    st.dataframe(
        terms[terms["response"] == term_response].sort_values("p_value"),
        width="stretch",
        hide_index=True,
    )


def data_tab(frame: pd.DataFrame) -> None:
    st.subheader("출처와 실험 설계")
    st.write(CITATION)
    factor_rows = [
        {
            "인자": spec.display_name,
            "단위": spec.unit,
            "최솟값": spec.actual_min,
            "중심값": spec.center,
            "최댓값": spec.actual_max,
        }
        for spec in FACTOR_SPECS.values()
    ]
    st.dataframe(pd.DataFrame(factor_rows), width="stretch", hide_index=True)
    st.write("32-run central composite design · 중심점 반복 6회 · 21항 2차 모델")
    with st.expander("전사 DOE 32행 보기"):
        st.dataframe(frame, width="stretch", hide_index=True)
    st.subheader("해석 제한")
    st.markdown(
        """
- 비균일도 모델 adjusted R²는 약 0.143이므로 생산 레시피 추천에 사용할 수 없습니다.
- anisotropy 전사값은 논문의 계수와 R²를 충분히 재현하지 못해 보조 응답으로만 둡니다.
- 예시 spec은 프로그램 동작 시연용이며 논문 또는 양산 승인 기준이 아닙니다.
- 여러 논문이나 synthetic response를 혼합하지 않았습니다.
"""
    )


def main() -> None:
    st.set_page_config(
        page_title="Etch Process Window Studio",
        page_icon="◫",
        layout="wide",
    )
    st.title("Etch Process Window Studio")
    st.caption("공개 plasma-etch DOE 재현 · RSM 검증 · 고객 spec 기반 공정창 탐색")
    st.warning(
        f"본 도구는 {CITATION}의 DOE 데이터 재현이며, "
        "해당 조건 범위 밖 외삽은 유효하지 않음"
    )

    check_precomputed_tables_are_current()
    summary_digest = sha256_file(PROJECT_ROOT / SUMMARY_RELATIVE_PATH)
    frame, models = load_models(summary_digest)
    tables = load_validation_tables(summary_digest)
    specs, vary, fixed_actual, resolution = sidebar_controls()
    process_tab, reproduction, analysis, data = st.tabs(
        ("Process window", "논문 재현 대조", "ANOVA · Main effect", "데이터 · 한계")
    )
    with process_tab:
        process_window_tab(
            models, specs, vary, fixed_actual, resolution
        )
    with reproduction:
        reproduction_tab(tables)
    with analysis:
        analysis_tab(models, tables)
    with data:
        data_tab(frame)


if __name__ == "__main__":
    main()
