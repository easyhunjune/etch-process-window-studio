# Project 2 Phase B — RSM 수치 코어 검증

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
| Etch rate | 0.9901 | 0.9720 | 0.0260 | 0.99 | 주 모델 |
| Selectivity | 0.8876 | 0.6832 | 2.7908 | 0.90 | 주 모델 |
| Non-uniformity | 0.6958 | 0.1428 | 1.6445 | — | 탐색 모델 |
| Anisotropy | 0.8371 | 0.5410 | 0.0995 | 0.95 | 보조 모델 |
| DC bias | 0.9931 | 0.9805 | 20.6847 | 0.99 | 감사용 |

논문 계수와 반올림 허용오차를 벗어난 대조 항은 18개다.
응답별로 etch rate 2개,
selectivity 1개,
anisotropy 14개,
DC bias 1개다.
Phase A에서 확인한 부호 불일치 외에도 anisotropy 전사값은 논문의 계수와
R²를 충분히 재현하지 못한다. 특히 etch-rate R²는 재현되지만 전사값 기반
잔차표준오차는 0.0260로,
논문 표의 0.01과 다르다. 원자료를 바꾸지 않고 이 차이를 감사 결과로 보존했다.

non-uniformity 모델의 adjusted R²는
0.143로 낮다.
따라서 이 응답을 포함한 process window는 생산 레시피 추천이 아니라
포트폴리오용 탐색·설명 기능으로만 해석해야 한다.

## 민감도와 ANOVA 해석

- Etch rate 최대 기여 그룹: `rf_power` (71.1%)
- Selectivity 최대 기여 그룹: `pressure` (39.4%)
- Non-uniformity 최대 기여 그룹: `chf3` (36.9%)

기여율은 각 요인에 속한 선형항·제곱항·상호작용항을 full model에서 함께
제거했을 때 늘어나는 부분 제곱합을 정규화한 값이다. 상호작용항이 양쪽
요인 그룹에 포함되므로 고전적인 가산형 ANOVA 분해와 같지 않다.

원 논문에는 별도의 ANOVA 표가 없고 21개 회귀계수, 모델 표준편차, R²만
제시되어 있다. 따라서 ANOVA는 전사 데이터로 계산한 프로젝트 산출물이며,
논문과의 직접 대조는 계수·R² 표에서 수행했다.

## 데모 process window

- 고정값: CHF3 12 sccm, pressure 100 mTorr, RF power 100 W
- 탐색축: SF6 10–50 sccm, O2 2–18 sccm
- 예시 spec: etch rate ≥ 0.4 µm/min,
  selectivity ≥ 8,
  non-uniformity ≤ 4.5%
- 평균예측 충족률: 46.4%
- 설계 지지영역 비율: 50.0%
- 95% 예측구간 확정률: 0.0%

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
