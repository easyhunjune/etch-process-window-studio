# Project 2 Phase A - 데이터 확보 결과

- 실행일: 2026-07-30
- 상태: 완료(3편 탐색 컷 규칙 적용)
- 주 논문: 선정
- 예비 논문: 적합 후보 없음
- 다음 단계: 주 논문 데이터만으로 Phase B 재현 진행

## 1. 선정 결과

주 논문은 Legtenberg et al. (1995)로 확정한다.

> R. Legtenberg, H. Jansen, M. de Boer, and M. Elwenspoek, "Anisotropic
> Reactive Ion Etching of Silicon Using SF6/O2/CHF3 Gas Mixtures," Journal of
> The Electrochemical Society, 142(6), 2020-2028, 1995.
> DOI: 10.1149/1.2044234

공개 학위논문의 Chapter 4가 위 저널 논문을 재수록한다. 해당 PDF의 85-89쪽에서
설계표, 실제 인자 수준, 32개 실험의 응답값, 논문 회귀계수를 직접 확인했다.

### 선정 기준 대조

| 기준 | 판정 | 근거 |
|---|---|---|
| 완전한 DOE 입력 조건과 응답값 | 통과 | 중앙합성회전 2차 설계 32행, Tables 4.1-4.3 |
| 응답 변수 2개 이상, etch rate 필수 | 통과 | etch rate, selectivity, anisotropy, non-uniformity |
| 논문 대조값 | 통과 | 2차 회귀계수 21개와 R2, Table 4.4 |
| 실측 데이터 | 통과 | Si etch depth, SiO2 thickness, 5-point etch depth 측정 절차 명시 |
| 외삽 방지 범위 | 통과 | 각 인자의 -2에서 +2 실제 수준이 Table 4.2에 명시 |

## 2. 확보한 주 데이터

- `data/raw/legtenberg_1995_ccd.csv`
  - 32개 실험행
  - coded/actual factor를 모두 보존
  - 응답: bias voltage, Si etch rate, Si/SiO2 selectivity, anisotropy,
    etch-rate non-uniformity, qualitative surface state
  - 논문의 별표는 `anisotropy_outward_slope`로 별도 보존
- `data/raw/legtenberg_1995_paper_coefficients.csv`
  - 논문의 full quadratic model 계수
  - 변수 순서: SF6, O2, CHF3, pressure, RF power
  - 논문 R2: etch rate 0.99, selectivity 0.90, anisotropy 0.95, DCB 0.99

## 3. 탈락 후보

### Muttalib et al. (2017)

9개 실험의 etch rate와 selectivity가 완전히 공개되어 있으나 예비 논문으로
채택하지 않았다.

- ANOVA 표와 회귀계수가 없다.
- Table 2의 pressure 열은 gas ratio 열과 완전히 결합되어 있어 L9 직교성이 없다.
- Table 4의 pressure-level 평균은 Table 2의 인쇄값으로 재현되지 않는다.
- Table 4 평균을 역으로 만족시키는 pressure 배치는 추론할 수 있지만, 이는 원문
  실측표의 명백한 정정 공지가 아니므로 재현의 기준 데이터로 사용할 수 없다.

감사 목적으로 인쇄값과 추론값을
`data/raw/muttalib_2017_l9_audit.csv`에 함께 보존했다. 추론값을 학습이나 모델
피팅의 정답으로 사용해서는 안 된다.

### Tipton et al. (1993)

20-run CCF RSM, 세 응답(et­ch rate, non-uniformity, selectivity)과 회귀계수는
제공하지만, 20개 run별 실측 응답표가 본문에 없다. 따라서 논문 계수의 독립 재현이
불가능해 예비 논문에서 제외했다.

## 4. 결정

가이드라인의 "3편 이상 탐색하지 않는다" 규칙에 따라 추가 검색을 중단한다.
예비 논문을 채우기 위해 기준을 낮추지 않는다. Phase B는 주 논문 하나로 진행하며,
논문 실측 32행 이외의 synthetic response를 만들지 않는다.

Process window의 초기 지원 응답은 다음과 같이 고정한다.

1. Si etch rate 최소값
2. Si/SiO2 selectivity 최소값
3. etch-rate non-uniformity 최대값

anisotropy는 추가 최소 spec으로 제공할 수 있으나, 별표가 붙은 outward-sloped
profile과 일반 undercut profile의 동일 수치가 같은 형상을 뜻하지 않는다는 제한을
UI에 표시해야 한다.

## 5. 재현성과 주의사항

- 모든 PDF SHA-256과 출처 URL은 `data/source_manifest.json`에 기록했다.
- `scripts/validate_acquisition.py`는 32행 설계, coded-to-actual 매핑, 중앙점
  6회와 후보 논문의 표 불일치를 검증한다. 원문 PDF는 공개 저장소에서
  제외하며, 별도로 받은 뒤 `--require-pdfs`를 사용하면 manifest의 SHA-256까지
  엄격하게 검증한다.
- `scripts/audit_paper_coefficients.py`는 전사한 32행으로 full quadratic OLS를
  적합해 Table 4.4와 항별로 비교한다.
- 논문 식 (4.3)은 PDF 추출에서 괄호가 손실되어 보일 수 있다. 실제 coding은
  `(phi_i - center_i) / step_i`이며 Table 4.1과 4.2의 대응으로 확인했다.
- 주 논문의 non-uniformity는 5개 위치 etch depth 측정에 기반한 값이다.
- 데이터 범위 밖 예측은 허용하지 않는다.

### Phase B로 넘길 계수 사전 감사 이슈

전사값에 full quadratic OLS를 적용한 결과, etch rate·selectivity·DCB 계수의
대부분은 논문 반올림 범위에서 일치했다. 다만 아래 항은 같은 크기의 반대 부호 또는
반올림으로 설명하기 어려운 차이를 보였다.

- etch rate `b3`(CHF3 linear): 전사값 적합 `-0.0044`, 논문 `+0.004`
- selectivity `b3`(CHF3 linear): 전사값 적합 `+0.1417`, 논문 `-0.14`
- DCB `b34`(CHF3 x pressure): 전사값 적합 `-6.5625`, 논문 `+7`
- anisotropy: 여러 계수에서 차이가 남음

앞의 세 항은 절댓값이 거의 같아 Table 4.4의 부호 인쇄 오류 가능성이 있다. 그러나
원저자 정정 공지가 없으므로 현 단계에서 논문 값을 고치지 않는다. anisotropy의
별표 행을 `2-A`로 바꾸거나 제외하는 두 해석도 Table 4.4를 복원하지 못했다.
Phase B에서는 etch rate·selectivity·non-uniformity를 먼저 구현하고, 이 차이를
재현 대조 탭에 공개한다. anisotropy는 계수 차이의 원인이 확인되기 전까지 보조
응답으로 취급한다.
