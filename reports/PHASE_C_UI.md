# Project 2 Phase C — Streamlit UI 검증

## 완료 범위

Phase B의 수치 코어를 변경하지 않고 Streamlit 인터페이스를 연결했다.
UI에서 raw DOE를 읽어 동일한 21항 모델을 적합하며, 논문 대조와 ANOVA는
Phase B가 생성한 검증 CSV를 표시한다.

## 화면 구성

### Process window

- etch rate 최소값, selectivity 최소값, non-uniformity 최대값 입력
- 5개 인자 중 contour의 가로축과 세로축 선택
- 나머지 세 인자는 논문 DOE 범위 안에서 slider로 고정
- 세 응답 contour와 동시 충족 영역 표시
- 충족률·격자 수 표시 및 충족 격자 CSV 다운로드

### 논문 재현 대조

- 응답별 내 적합 R², adjusted R², 잔차표준오차와 논문 R²
- 21개 항의 내 적합 계수·논문 계수·차이·허용오차 판정
- 원 논문에 ANOVA 표가 없다는 사실과 프로젝트 계산 ANOVA를 함께 표시

### ANOVA · Main effect

- 요인 그룹 제거 기반 정규화 부분 제곱합
- 중심점에서 한 인자만 coded −2~+2로 변화시킨 main-effect 곡선
- 응답별 항의 부분 F 검정과 p-value

### 데이터 · 한계

- 논문 서지사항과 5개 인자의 실제 단위 범위
- 32행 전사 DOE
- 낮은 비균일도 adjusted R², anisotropy 전사 불일치, 데모 spec의
  비생산 기준 성격을 명시

## 가이드라인 완료 기준

| 기준 | 결과 | 구현 위치 |
|---|---|---|
| 원 논문 계수/ANOVA 대조 | 통과 | `논문 재현 대조` 탭 |
| spec 3개 입력 후 window 표시 | 통과 | sidebar + `Process window` 탭 |
| 설계 지지영역·불확실성 표시 | 통과 | 사선 마스크 + 3상태 공정창 |
| core unit test 3개 이상 | 통과 | 전체 14개 테스트 |

원 논문은 회귀계수·모델 표준편차·R²만 제공하고 ANOVA 표는 제공하지
않는다. 따라서 논문 대조는 계수와 R²에서 수행하고, ANOVA는 전사 데이터로
계산한 프로젝트 결과임을 화면에서 구분한다.

## 검증 결과

- Python compile: PASS
- Phase A 출처·해시 검증: PASS
- Phase B 산출물 재생성: PASS
- pytest: 15 passed
- Streamlit AppTest: 예외 없음, 4개 탭과 필수 입력 확인
- AppTest spec 변경: 최소 식각 속도 변경 후 충족률 재계산 확인
- 실제 서버 health endpoint: HTTP 200, `ok`
- 실제 브라우저: 첫 화면·논문 대조 탭 렌더링 확인
- AppTest에서 최소 식각 속도 변경 후 평균예측 충족률 재계산 확인
- 브라우저 console error: 0개

## 실행

```powershell
cd project2
$env:PYTHONPATH = "src"
python -m streamlit run app.py
```

## 해석 원칙

이 도구는 논문 DOE의 재현·설명용 분석기다. 공정조건은 실험 범위 밖으로
입력할 수 없으며, 표시된 공정창은 양산 레시피 승인이나 장비 간 일반화를
의미하지 않는다. 비균일도 모델의 adjusted R²가 낮으므로 동시 충족 영역은
탐색적 시각화로만 사용한다.
