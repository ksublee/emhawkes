# Hawkes curvature stress checks

## Findings and change

Across 35 deterministic cases, stable curvature returns finite standard errors in
10 cases and withholds them in 25. The fitted coefficients, log likelihood,
returned optimizer gradient, convergence code, iteration count, final Hessian,
covariance, curvature status, reason and step checks are exactly unchanged.
No numerical thresholds, covariance formulas or optimizer settings were changed.

One warning defect was reproduced: when optimization stopped early and curvature
was also unavailable, only the unavailable-standard-error warning was emitted.
The existing nonconvergence warning was inside the successful-curvature branch.
It is now emitted for checked unsuccessful fits regardless of curvature status.
The additional `optimizer_not_converged` code is recorded in `warnings`.
Both deliberately stopped cases now emit both warnings. Notices in the other
33 cases are unchanged. The default `final_hessian = "maxLik"` is unchanged.

## Design and independent reference

All cases are univariate Hawkes models with one or two exponential kernels.
Amplitudes and decay parameters remain separate with `reduced = FALSE`.
Initial components are specified explicitly and treated as constants.
Numerical optimization and the `auto` engine are used in the study.

The fixed-decay comparison fits the same 400-event sample generated with
`mu = 0.5`, amplitudes `(0.18, 0.22)` and rates `(1, 1.01)`. The fitted rates
are then held at `(1, 1 + gap)` for gaps from zero to 0.5. This isolates numerical
conditioning; these fitted fixed rates are not all the data-generating rates.
These amplitude fits are unconstrained. Separate free-decay and weak-excitation
fits impose nonnegative amplitudes and positive baseline and decay bounds.
Inequality runs start strictly inside their bounds, as required by maxLik.
Redundant inequalities on a coefficient fixed at zero are omitted.

The independent R reference forms all earlier-event time differences directly.
For fixed rates, intensity is linear in the free baseline and amplitudes. If
`Z` contains their event-intensity derivatives and `lambda` the event intensities,
the observed information is `crossprod(Z / lambda)`. SVD of the scaled design
computes the reference covariance without first forming an ill-conditioned Gram
matrix. For free rates, the reference differentiates the pair sums and closed-form
exponential exposure integrals, including all cross derivatives and known initial
transients. It uses neither the package event recursion nor its score differences.

Among the 10 finite cases, the maximum relative standard-error difference from
this reference is **0.000444%**. This is a numerical comparison at the fitted
points, not an interval-coverage study. Elapsed times in the raw CSV are incidental
and should not be read as a performance benchmark.

## Results

| Cases | Count | Finite SE | Unavailable SE | Interpretation |
| --- | ---: | ---: | ---: | --- |
| Fixed decay gap: 0 to 0.5 | 7 | 5 | 2 | Coincident kernels or numerical rank deficiency are rejected. |
| Both decays free, three starting values | 3 | 0 | 3 | A collapsed pair or a zero-amplitude boundary is detected. |
| Only 1, 2 or 5 regular events | 3 | 0 | 3 | Unidentified curvature or an amplitude boundary is detected. |
| Weak excitation: 12/40/150 events and 12 seeds at 40 events | 15 | 3 | 12 | Some fits are interior; most reach a boundary in these samples. |
| Zero amplitude, initial excitation 0 or 30 | 2 | 1 | 1 | A known transient can identify decay despite zero jumps. |
| Active amplitude or baseline bound | 2 | 0 | 2 | Ordinary Wald standard errors are withheld. |
| Interior constrained fit | 1 | 1 | 0 | Constraints alone do not prevent curvature estimation. |
| Deliberate early optimizer stops | 2 | 0 | 2 | Both nonconvergence and unavailable curvature are now warned. |

For the two otherwise identical kernels driven by the same events:

- With equal fixed rates, transferring 0.1 between amplitudes leaves the likelihood
  exactly unchanged on the tested grid. Only their sum is identified.
- With rates `(1, 1.0001)`, amplitude standard errors are approximately 1,077,
  and the scaled information condition number is `9.61e9`. Stable returns the
  same large uncertainty as the independent SVD reference and warns about weak
  curvature. A transfer of 0.1 changes log likelihood by less than `5e-7`.
- At gap `1e-5`, the numerical rank criterion withholds standard errors.

Equal rates do not by themselves imply nonidentification for kernels with
different event drivers or other distinguishing structure. The flat-direction
result concerns otherwise identical kernels excited by the same univariate events.

For the 40-event weak-excitation sample with seed 2, the estimated amplitude is
`0.0590` with SE `0.1788`, and the estimated decay is `0.8228` with SE `1.1831`.
These finite errors describe substantial uncertainty. At amplitude zero, holding
the other fitted values fixed reduces log likelihood by only `0.1552`.
The supplied likelihood slices are conditional slices, not profile likelihoods;
other parameters are held fixed. They do not justify Wald coverage in short data.

With amplitude fixed at zero and no initial excitation, decay is unidentifiable.
With initial excitation fixed at 30, decay is estimated as `0.9150` with SE
`0.1323` in the selected 150-event realization. The independent analytic Hessian
agrees. No blanket rule rejects small or zero amplitudes or short samples.

## Reproduce and inspect

```sh
Rscript --vanilla tools/check-stable-stress.R BEFORE_LIB BEFORE_OUT
Rscript --vanilla tools/check-stable-stress.R AFTER_LIB AFTER_OUT BEFORE_OUT
```

The second command asserts exact numerical agreement with the baseline and checks
that optimizer failure is warned even with unavailable curvature.

See [case outcomes](stable-stress-results/cases.csv),
[baseline outcomes](stable-stress-results/before-cases.csv),
[coefficients and standard errors](stable-stress-results/coefficients.csv),
[unchanged-result assertions](stable-stress-results/unchanged-results.csv),
[conditional likelihood slices](stable-stress-results/conditional-slices.csv),
[emitted warnings](stable-stress-results/warnings.csv),
[manifest](stable-stress-results/manifest.txt), and
[session information](stable-stress-results/session-info.txt).

## Validation

Source build and `R CMD check --no-manual` pass with **Status: OK**, including
all **27 test scripts**, examples and vignette rebuilding. There are no errors,
warnings or notes. The checked source and tests match the final working files.

The new regression fails on baseline `ecedd5d` because the unavailable-curvature
fit omits `optimizer_not_converged`; it passes after the warning relocation.
Other regression cases cover SVD reference covariance, coincident and nearly
coincident rates, explicit R likelihood differences, active Hawkes constraints,
and the effect of known initial excitation. The weak-curvature unit case uses
rate gap `2e-4`, away from the numerical-rank cutoff; the 35-case study also
retains the more demanding `1e-4` case.

Parsed-source checks confirm that the runtime change only relocates the existing
convergence-warning block. Numerical calculations are unchanged. The baseline
source matches the prior isolated installation, and the updated Rd file passes
`tools::checkRd`. See [source verification](stable-stress-results/source-verification.txt).
