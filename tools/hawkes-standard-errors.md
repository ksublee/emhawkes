# Hawkes 모형의 최종 Hessian과 표준오차

2026-10-13에 기존 계산과 선택적인 안정화 계산을 비교했습니다.
Windows 11, R 4.6.0, maxLik 1.5-2.2 환경이며, 원래 예제의 seed,
관측 자료, 모수 공유, 초기 상태, 최적화 설정을 유지했습니다.

## 변경 내용

`hfit()`의 최적화는 그대로 실행합니다. 최종 단계에서만 수치 Hessian과
공분산을 다시 계산합니다. `coef()`, 로그우도, 저장된 optimizer gradient,
반환 코드, 반복 횟수는 변경하지 않습니다.

기본값 `final_hessian = "maxLik"`는 기존 최종 Hessian과 표준오차 계산을 유지합니다.
`final_hessian = "auto"`를 명시했을 때만 built-in Hawkes 우도에 안정화를 적용합니다.
사용자 Hessian은 보존하고 정보행렬의 역산 가능 여부를 확인합니다.
`finalHessian = FALSE`, BHHH, 사용자 정의 `mylogLik`, FRPP 등 residual 우도의
계산 선택은 유지합니다. `final_hessian = "maxLik"`는 이전 계산을 재현합니다.

수치 계산은 다음을 확인합니다.

- 자유 모수의 scale은 `pmax(abs(coef), .1)`입니다. 고정 모수는 미분에서 제외합니다.
- 중앙 차분의 상대 간격은 `.004`, `.002`, `.001`부터 시작합니다.
  필요하면 절반씩 줄여 최대 여섯 간격을 확인합니다.
- Richardson 외삽의 두 결과를 비교합니다. Hessian 상대 변화는 1%,
  공분산 변화는 2%, 각 표준오차 변화는 1% 이내여야 합니다.
  인접 중앙 차분의 표준오차 변화도 2% 이내여야 합니다.
- 인접한 세 중앙 차분 정보행렬은 모두 양의 정부호여야 합니다.
  Scale 적용 최소 고유값이 최대 절대 고유값의 `1e-10` 이하이면 역산하지 않습니다.
- 등식·지원하지 않는 제약, 또는 최대 차분 영역에 가까운 부등식 제약에서는
  일반 Wald 공분산을 제공하지 않습니다. 부등식의 판정은 제약식의 배율에 불변입니다.

통과하면 `maxLik` 객체의 최종 Hessian과 공분산을 갱신합니다.
`summary()`와 `vcov()`가 같은 공분산을 사용합니다.
통과하지 못하면 자유 모수의 공분산을 `NA`로 저장합니다.
고정 모수의 공분산은 0입니다. 음의 분산에 절댓값을 취하거나,
정보행렬의 작은 고유값을 임의로 올리지 않습니다.

계산 상태와 원인은 `fit_diagnostics()`로 확인합니다.
`status = "ok"`는 수치 곡률 검사를 통과했다는 뜻입니다.
전역 최적해, 식별성, 정상성, 신뢰구간의 포함률을 보장하는 표시가 아닙니다.

## 기존 추정값과 표준오차 비교

세 예제 모두 실제 `hfit()`을 다시 실행했습니다. 모든 계수, 로그우도,
optimizer gradient, 반환 코드, 반복 횟수의 차이가 **0**입니다.
기본 `maxLik` 계산과 명시적으로 요청한 안정화 계산을 각각 실행하여 비교했습니다.

| 예제 | 유지된 로그우도 | 반복 횟수 | 유지된 반환 코드 |
| --- | ---: | ---: | ---: |
| Multi-kernel | 6084.70174869 | 75 | 0 |
| Synchronized intensity | -260.62730006 | 73 | 0 |
| Bid-ask price | -1810.89162660 | 41 | 0 |

| 예제·모수 | 이전 표준오차 | 새 표준오차 |
| --- | ---: | ---: |
| Multi-kernel `alpha1.1` | `NaN` | 0.35773384 |
| Multi-kernel `beta1.1` | 0.03791102 | 1.63305102 |
| Synchronized `kappa` | `NaN` | 0.15486893 |
| Synchronized `eta` | `NaN` | 0.44037175 |
| Synchronized `nu` | `NaN` | 0.39985180 |
| Synchronized `beta` | `NaN` | 4.73851611 |
| Bid-ask `xi` | `NaN` | 0.83863463 |
| Bid-ask `beta1.1` | 0.35535 | 2.87451121 |

모든 새 표준오차는 유한합니다. 원래 fit에 후처리만 적용한 Hessian과
공분산은 안정화를 요청한 `hfit()` 결과와 `identical()`까지 일치했습니다.
별도로 계산한 중앙 차분 또는 이전 진단의 Richardson 결과와 비교한
최대 표준오차 상대 차이는 multi-kernel `2.96e-5`, synchronized `1.46e-6`,
bid-ask `4.63e-6`입니다. 이것은 수치 계산의 일치 검증이며,
실제 표본 분포의 표준편차나 신뢰구간 포함률 검증은 아닙니다.

수정 전 finite 표준오차도 부정확했습니다. 따라서 `NaN` 모수만 선택적으로
보정하지 않고 전체 자유 모수 정보행렬을 함께 검증했습니다.

## 남아 있는 추정 진단

Multi-kernel의 Newton decrement는 약 0.19824입니다.
그 제곱의 절반은 국소적인 이차 근사에서 기대되는 우도 개선량입니다.
현재 값은 약 0.01965여서 score 경고를 냅니다. 앞 진단의 엄격한 재추정에서도
약 0.01830의 개선을 확인했습니다. 새 코드가 자동으로 그 해를 선택하지는 않습니다.

원래 시작 순서의 엄격한 추정은 `6084.72004374`에 도달했고,
커널 그룹의 시작 순서를 바꾸면 `6084.98473958`도 발견했습니다.
커널별 고정 초기 상태를 바꾸지 않았으므로 단순 라벨 교환이 아닙니다.
Multi-kernel 추정에서는 커널 순서와 초기 상태를 포함한 여러 시작점 비교가 필요합니다.

Synchronized와 bid-ask의 Newton decrement는 약 `5.45e-5`, `0.00475`입니다.
세 예제의 scale 적용 정보행렬 조건수는 각각 약 3788, 1449, 180입니다.
조건수는 scale에 의존하므로 다른 scale로 계산한 이전 보고서 값과 구분해야 합니다.

정상 수렴 코드도 최적화 방법별로 다릅니다.
NR·BFGSR 계열은 1·2·8, BFGS·CG·Nelder–Mead는 0을 정상 종료로 봅니다.
SANN의 0은 설정된 반복 실행의 완료여서 수렴 여부를 단정하지 않습니다.
이 상태는 수치 곡률 상태와 별도로 기록합니다.
[maxLik NR 문서](https://search.r-project.org/CRAN/refmans/maxLik/html/maxNR.html),
[BFGS 문서](https://search.r-project.org/CRAN/refmans/maxLik/html/maxBFGS.html)

FRPP는 변경하지 않았습니다. 기존 fit을 후처리에 전달해도 전체 객체가
그대로 유지됩니다. 비매끄러운 residual 우도는 별도의 추론 검증이 필요합니다.

## 사용 방법과 재현

```r
fit <- hfit(h, inter_arrival, type, lambda_component0 = initial,
            final_hessian = "auto")
summary(fit)
fit_diagnostics(fit)[c("status", "scaled_condition", "newton_decrement", "warnings")]

# 기본 최종 Hessian과 표준오차 계산
legacy <- hfit(h, inter_arrival, type, lambda_component0 = initial,
               final_hessian = "maxLik")
```

세 예제를 두 방식으로 재현하고 독립 중앙 차분과 비교할 수 있습니다.
저장소 루트에서 실행합니다.

```sh
Rscript tools/reproduce-fit-diagnostics.R <output>/legacy multi_kernel,synchronized,bid_ask maxLik
Rscript tools/reproduce-fit-diagnostics.R <output>/auto multi_kernel,synchronized,bid_ask auto
Rscript tools/check-hawkes-standard-errors.R <output>
```

[수치 결과](hawkes-standard-error-results)에 모수별 표준오차,
차분 간격 검사, 추정값 보존 비교를 기록했습니다.
이전 `NaN`은 R의 CSV 출력에서 `NA`로 기록됩니다.
후처리 단독 시간은 multi-kernel 0.68초, synchronized 2.89초,
bid-ask 28.3초였습니다. 이 시간은 한 환경의 측정값입니다.
사용자 정의 `impact`를 R에서 평가하는 bid-ask에서는 추가 계산 비용이 큽니다.

## 검증

`R CMD build --no-manual`과 `R CMD check --no-manual`은 `Status: OK`입니다.
11개 테스트, 예제, vignette 검사 모두 통과했습니다.

새 회귀 테스트에는 해석적 Poisson 정보행렬, 고정 모수, 저장·복원,
사용자 Hessian, disabled/BHHH 선택, 방법별 수렴 코드, 실제 미식별 decay,
부등식 경계와 배율 불변성, FRPP·custom 우도 제외, 꺾이는 우도의 간격 민감성이 포함됩니다.

Portable 비교 도구는 저장된 세 예제에서 추정값 불변과 독립 중앙 차분의
표준오차 일치를 확인했습니다. 수정 전·후 재현 선택도 실제로 실행했습니다.
확장 예제의 실제 출력에는 Hawkes의 `NaN`이 없고 FRPP의 두 `NaN`만 남습니다.
문서의 링크·이미지를 검사했고, 같은 플롯의 기존 고해상도 그림을 유지했습니다.
## 표준오차가 계속 불안정할 때

수치 미분의 정확도와 통계적 식별성은 다른 문제입니다.
표준 Hawkes에서도 excitation이 약하면 decay 정보가 적고,
같은 사건 유형에 작용하는 커널의 decay가 가까우면 개별 성분을 분리하기 어렵습니다.
짧은 표본, 제약 경계, 근임계 설정, 모수에 threshold를 쓰는 callback도
일반 Wald 추론에 불리할 수 있습니다. 고정된 nonzero 초기 강도는 decay에 대한
정보를 남길 수 있으므로 excitation이 작다는 이유만으로 일괄 비식별 판정을 하지 않습니다.

권장 순서는 다음과 같습니다.

1. 여러 시작점, score, 간격 안정성을 먼저 확인합니다. 특히 문제 모수의
   profile likelihood에서 다른 모수를 다시 최적화하여 비대칭·평탄함·다중 최대를 봅니다.
2. 중요한 구간 추론에는 관측 설계와 초기 상태를 보존한 bootstrap을 비교합니다.
   Hawkes 연구에서도 Gaussian 구간의 유한 표본 포함률이 낮을 수 있고,
   fixed intensity bootstrap이 개선을 보인 설정이 있습니다.
   모든 state dependent 모형에 같은 성능을 보장하는 결과는 아닙니다.
   [Cavaliere 외, 2023](https://arxiv.org/abs/2104.03122)
3. 수치 Hawkes 경로에는 analytic 미분 또는 AD를 적용하고 수치 미분과 교차 검증합니다.
   지수 Hawkes의 analytic 미분은 알려져 있습니다. 사용자 callback과
   공유 모수 매핑의 미분은 별도 처리가 필요합니다.
   [Ozaki, 1979](https://www.ism.ac.jp/editsec/aism/pdf/031_1_0145.pdf)
4. 표본 길이, 약한 excitation, 근접 decay, 근임계 설정을 바꾸는 Monte Carlo에서
   bias, 실제 표준편차, 평균 추정 표준오차, 구간 포함률, 추정 실패율을 함께 비교합니다.

현재 우도는 마지막 사건에서 끝나고 `hsim(size=n)`은 사건 수를 기준으로 생성합니다.
고정 시간 관측창의 bootstrap에서는 사건 수가 달라질 수 있고 마지막 무사건 구간도
고려해야 합니다. 이 관측 설계를 지원하는 인터페이스와 marks·state dependence·
초기 상태의 재현 정책을 마련한 뒤 bootstrap을 구현해야 합니다.
경계나 약한 식별성은 bootstrap이나 profile의 일반적인 임계값으로도 자동 해결되지 않습니다.
[Self와 Liang, 1987](https://pages.stat.wisc.edu/~larget/Stat998/Fall2015/Self-Liang-1987.pdf)

현재 수정은 수치 Hessian 문제를 줄이고 불안정한 추론을 드러내는 단계입니다.
Profile 구간, bootstrap 구간, Monte Carlo 포함률은 아직 구현·검증하지 않았습니다.
