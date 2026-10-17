# Pre-release bug audit

Date: 2026-10-17. Reference: `0a88f47` on `main`.

This document records successive validation stages. Earlier references to
uncommitted changes or pending CI describe those stages. See the final section
for the current release-preparation results.

## Confirmed defects

1. `print.hreal()` and `summary.hreal()` overwrote the caller's `digits` option.
   With `digits = 12`, successful output left 7; interrupted output left 4 or 5.
   `n = 0` also displayed the initial row. Restore options on exit, validate row
   limits, and keep matrix formatting for zero or one displayed rows.
2. `hvol()` failed with omitted `eta` and returned `NaN` with omitted marks.
   It also silently used diagonal decay rates for unequal-rate rows and ignored
   unsupported impacts or mappings. Normalize unit marks, treat missing `eta`
   as zero, and reject models outside the bivariate moment formula's scope.
   Replace its internal deprecated inference call. Remove unused count-weighted
   quantities that no longer have count paths in the inference result.
3. `expected_tau()` integrated only finite-event density and missed probability
   at an infinite waiting time. With `mu = 0`, `beta = 1`, and total excitation
   1, it returned about 0.484829 instead of `Inf`. A Poisson rate of `1e-6`
   produced a divergent-integral error instead of `1e6`. Detect missing-event
   probability and integrate survival on log time, splitting at baseline,
   excitation and decay scales. Custom laws use this route when their CDF is
   available. Preserve finite-limit truncated moments and density-only custom
   integration; these retain the usual numerical quadrature limitations.

## Validation

- 180 independently randomized simulation designs: dimensions 1–5, kernels
  1–8, overlapping/reversed mappings, signed equal-rate components, numeric
  mark effects, sizes 1–101, and time units from `1e-5` to `1e5`. Every returned
  value, warning and final RNG state matched the R reference exactly.
- R and C++ likelihood replay reconstructed all four pre/post intensity fields
  in those designs with maximum absolute error zero. Full-history sums also
  agreed within `1e-12` relative comparison.
- 90 ordinary-versus-exponential-residual likelihood comparisons passed,
  including initial-only paths and missing event types.
- Independent compound Poisson variances check both `hvol()` branches. Fifty
  randomized unit-mark Hawkes models also matched the resolvent covariance
  formula in both branches; maximum absolute volatility error was `4.45e-16`.
  Five of these identities are retained as package regression tests.
  Forty existing marked-model volatility outputs matched the old implementation
  exactly; the internal deprecation warning is intentionally removed.
- Expected waiting times are checked against independent Poisson-mixture
  moments, including small baselines and subnormal rates, and against bounded
  and exponential residual laws. New tests also cover infinite waits, truncated
  moments, unsupported inputs and print cleanup on error.

The initial bug-fix source archive passed `R CMD check --as-cran --run-donttest` on
Windows with R 4.5.1: **Status: OK**, including all 37 test files, vignette
rebuilding, examples, and PDF/HTML manuals. A first-pass code-analysis NOTE
about two local functions sharing a name was removed before that check.
The archive SHA-256 is
`bcc89127070cae90f282274ee096ea1f1bf0e78d5e784195c80cb82cc969ad54`.
Cross-platform CI for these uncommitted fixes has not run yet.

## Diagnostic-message follow-up

- Correct the singular-initializer warning: it uses **zero excitation**, not
  `mu`. Keep the same fallback and explain how to supply an explicit state.
- Use one mapped-initialization error for simulation, fitting and both likelihood
  engines. Explain why initialization is required and show the matrix shape or
  the scalar alternative. Dimension errors report the conflicting sizes.
- Explain the specific unsupported feature in `hvol()` rather than listing all
  restrictions. Callback errors suggest optional defaults or captured values.
- Report returned/requested rows on early simulation stopping, including the
  initial row. Keep warning severity and R/C++ wording aligned, independent of
  display precision and scientific-notation options.
- Explain thinning NaNs as numerical representation failures and suggest time
  rescaling. Do not automatically rescale, retry or change the RNG path.
- Distinguish unsupported stability checking from a failed stability condition;
  the latter reports the spectral radius. Keep existing warning conditions.
- Suppress misleading automatic-initialization announcements for custom or
  invalid fits. Keep verbose output for valid Hawkes initialization.
- Clarify deprecation replacements and unavailable Wald-curvature warnings.
  Preserve machine-readable reason codes, fitting decisions and callback errors.

`tests/user-diagnostics.R` covers shared errors, fallback values, callback
evaluation counts, verbosity, deprecations and three native early-stop paths.
It also compares warnings and RNG state on thinning failure and checks that
large row counts print consistently with nondefault display settings.

The follow-up also reproduced an existing RNG defect in the pre-message-change
archive: `options(warn = 2)` let internal native warnings jump over C++ RNG
cleanup. Protected R warning evaluation and explicit simulation RNG scopes now
cover warning handlers that draw random values, call nested simulations or
abort. Regression checks compare complete results, original handler conditions,
subsequent paths and RNG states, including Box-Muller normal-cache behavior.
All nine native forced-GC fixtures and callback error/interrupt recovery checks
also passed locally. One existing error-message assertion now checks the shared
mapped-state reason and expected matrix shape; numerical tolerances are unchanged.

The follow-up source archive passed `R CMD check --as-cran --run-donttest` on
Windows/R 4.5.1: **Status: OK**, with all 38 test files, rebuilt vignettes,
examples and PDF/HTML manuals. Its SHA-256 is
`34afd5f7075e1da5346af9c74dff090c7ace2bd5eb76211450ddfd211e0ca7c7`.
Cross-platform and Valgrind CI for the follow-up remain pending.

Event-selection, decay and likelihood formulas are unchanged; warning
evaluation and RNG synchronization were corrected. These checks
do not establish finite-sample Wald coverage or validate every custom callback
and residual distribution. Statistical calibration experiments remain separate
from the release implementation.

## Constraint argument compatibility

`hfit()` accepts both `constraint` and `constraints` as equivalent spellings,
without deprecation warnings. The named-only singular alias preserves existing
positional arguments. Two non-NULL values are rejected before fitting; an
explicit NULL for the unused spelling is allowed.

`tests/constraint-alias.R` compares estimates, likelihoods, scores, Hessians,
covariance fields, convergence information, warnings, diagnostics and restored
initial states exactly, across R/C++ engines and default/stable curvature.
It also covers positional compatibility, reconstruction, explicit NULL values
and identical handling of unsupported analytic equality constraints.

The updated source archive passed `R CMD check --as-cran --run-donttest` on
Windows/R 4.5.1: **Status: OK**, with all 39 test files, rebuilt vignettes,
examples and PDF/HTML manuals. Its SHA-256 is
`423d768ade0d6ff200efcff2be352bdf8d31c2776acf578db6bad6fce198f9ab`.
Cross-platform CI for these uncommitted changes remains pending.

At that stage, [Base R optimizer fallback](optimizer-fallback-plan.md) was a
proposal with a small Poisson likelihood experiment. The subsequent work below
implements the independent provider.

## Optional stats optimizer

- Add `hfit(optimizer = "stats")`, using `optim()` or `constrOptim()` directly.
  The default maxLik path, its result classes, likelihood engine selection and
  callback binding are preserved. Make maxLik a suggested dependency rather
  than an imported namespace; requesting it without installation gives explicit
  guidance to select stats.
- Reuse likelihoods, scores, parameter sharing, initial-state handling and stable
  curvature checks. Support fixed coordinates, linear inequalities and box
  bounds; reject equality constraints, unsupported algorithms and maxLik controls.
  Return a separate stats fit class with standard result methods, reconstruction
  and recorded provider information. Default stats fits to stable curvature;
  raw curvature never supplies Wald covariance near a constraint boundary.
- Check independent Poisson and correlated-quadratic optima and information,
  observed-event baselines, derivative attributes and precedence, fixed-only
  constraints, all-fixed fits, boundary estimates and nonconvergence.
- In a fresh R process with only emhawkes, Rcpp and R's own packages, fitting,
  uncertainty, saving/reloading, summary and model reconstruction all passed.
  Poisoning maxLik entry points also left stats fitting and reconstruction working.
- `tools/stats-optimizer-validation.R` passed 48 multivariate fits with mappings,
  mark effects, reduced/full maps, numerical/analytic scores and both engines.
  Fitted likelihood replay errors were zero. Coefficient, likelihood and
  covariance comparisons with maxLik passed their stated criteria.
- Sixteen prior-versus-current maxLik fits matched exactly in estimates,
  likelihoods, gradients, Hessians, covariance fields, convergence information,
  classes and final RNG states. Existing simulation code is not changed here.

The first full check without maxLik exposed a missing dependency guard in an
existing test that called `hfit` through `do.call`. Guard only that maxLik-specific
section, retaining the file's independent likelihood tests. The same source
implementation had already passed the full check with maxLik. The final checks
below use the updated test archive in both environments.

The final archive passed `R CMD check --as-cran --run-donttest` on Windows/R
4.5.1 in **both environments: Status: OK**, with zero errors, warnings or notes.
With maxLik present, all 41 test files ran, together with rebuilt vignettes,
examples and PDF/HTML manuals. With maxLik absent,
`_R_CHECK_FORCE_SUGGESTS_=false` permitted the intentional missing suggestion;
maxLik-specific tests/examples were skipped, while independent likelihood,
simulation and stats-provider checks, vignette rebuilding and manuals passed.
Missing suggested-package and external-Rd-link checks were informational only.
The final archive SHA-256 is
`0438c0aa4ada2d2cf2550898a11f39fb777044a3c643bf1fab69c434fefb0ee6`.
Validation logs are in the task's `stats-optimizer/verified` and
`stats-optimizer/verified-without-maxLik` artifact directories. Cross-platform
GitHub checks have not run for these uncommitted changes.

## CRAN release preparation (2026-10-17)

The package implementation remains `4737f96`; this preparation changes release
metadata, documentation and CI artifact retention only.

- Confirmed CRAN's published version is 0.9.8. A complete CRAN source index query,
  without R-version/platform filtering, found no reverse Depends, Imports,
  LinkingTo, Suggests or Enhances dependencies.
- Removed unused `LazyData`, added the bug-report URL and checked all nine URLs.
  Updated README for optional maxLik and the independent stats provider;
  consolidated NEWS around final behavior and corrected the duplicate historical
  version heading using the published 0.9.8 source archive.
- Built `emhawkes_0.9.9.tar.gz` on Windows/R 4.5.1. Compared 93 included tracked
  files and all original DESCRIPTION fields with the working source; no
  differences. Build-added DESCRIPTION fields are expected. The 106 archive
  entries exclude developer tools, the website, CI files, submission notes and
  compiled/log/session artifacts.
- The same archive passed `R CMD check --as-cran --run-donttest --timings` with
  maxLik present and absent: **Status: OK**, zero errors, warnings and notes.
  Both checks ran all 42 test files, rebuilt vignettes and PDF/HTML manuals.
  In the absent-provider check, maxLik-specific sections were skipped and
  `_R_CHECK_FORCE_SUGGESTS_=false` permitted the intentionally missing suggestion.
- The exact README quick start converged normally in an isolated session with
  only emhawkes, Rcpp and R's own packages. maxLik was neither available nor
  loaded. The longest package example took 2.16 seconds in the full check.
- Existing implementation CI passed Windows release/R 4.0.0, macOS release and
  Ubuntu release/devel/oldrel-1, plus native GC/interruption/Valgrind checks.
  The observed current release was R 4.6.1 and devel was R 4.7.0 (2026-10-15).
- CI now runs donttest examples, retains successful check logs and exports the
  Ubuntu release source archive with its SHA-256, source commit and session.
  This supplies a candidate built with current R rather than local R 4.5.1.

Local archive SHA-256:
`852efad2598b5f9983e5ba9ee87983e0dd43053432a766cfe8bca79be5be3341`.
Logs, the CRAN package-index snapshot, archive comparison and isolated README
check are in the task's `cran-preflight` artifact directory.

Before actual submission, confirm the latest main CI and retained release
artifact, then upload that checked source archive and complete the maintainer's
CRAN confirmation email. No CRAN or win-builder upload was made in this audit.
Stable-curvature checks remain distinct from finite-sample Wald coverage;
weak-excitation interval calibration is not claimed by this release.
