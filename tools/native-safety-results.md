# Minimum R and native safety checks

The package declares R >= 4.0.0. The R CMD check workflow now includes that
exact version on Windows, in addition to the existing release, devel and
old-release jobs. The minimum-version job checks the running version explicitly
and runs the package tests and vignettes.

The separate Native safety workflow runs on Linux and checks:

- Forced garbage collection during nine native fixtures: both numeric
  simulation paths, callback simulation, observed baseline, baseline Jacobian,
  state-dependent likelihood, observed-baseline likelihood, numeric score and
  internal simulation warnings with nested RNG-using handlers or `warn = 2`.
  Collection occurs at every allocation for the small direct calls and every
  tenth allocation for the two public likelihood calls. Each result and final
  RNG state must exactly match execution without forced collection.
- Errors and interrupt conditions from both mark and impact callbacks. The
  condition, callback count, RNG state, subsequent draws and recovery simulation
  must match the R reference path.
- A real POSIX SIGINT during the native thinning rejection loop. A bounded,
  owned child process is interrupted; R replays the rejected exponential and
  uniform draws to verify RNG recovery at an interrupt checkpoint.
- Valgrind Memcheck on a debug build, with origin tracking and failures for
  invalid accesses or definitely lost memory. GC tests complement Memcheck,
  which does not identify every R-object protection error in a normal R build.

Local R 4.5.1 validation passed all eight fixtures with forced collection and
all four callback failure/interrupt cases. The full source-package check,
including 34 test files and vignette rebuilding, completed with Status: OK.
R 4.0.0, POSIX signals and Valgrind require the corresponding CI environments;
their logs and the `native-safety` artifact provide the run-specific evidence.
The [first Linux run](https://github.com/ksublee/emhawkes/actions/runs/38915442553)
passed forced GC, callback unwinding, actual SIGINT and Valgrind. RNG replay
matched 215,039 rejected exponential/uniform pairs before the interrupt.

The first R 4.0.0 run used Rcpp 1.0.8.3 and found an ambiguous constructor for
sensitivity parameter names. An explicit `Rcpp::as<Rcpp::CharacterVector>`
conversion resolves the overload without changing the parameter values.
The next run compiled successfully but exposed lost callback condition classes
and changed messages. Rcpp 1.0.8.3 makes unwind protection opt-in, so both
`Makevars` files explicitly enable `RCPP_USE_UNWIND_PROTECT`. This preserves
original R conditions and lets C++ destructors restore RNG state while unwinding.
Newer Rcpp versions already enable that mechanism by default; the R and Rcpp
version requirements are unchanged. See the
[Rcpp 1.0.8.3 headers](https://github.com/cran/Rcpp/blob/1.0.8.3/inst/include/Rcpp/r/headers.h)
and [evaluation implementation](https://github.com/cran/Rcpp/blob/1.0.8.3/inst/include/Rcpp/api/meat/Rcpp_eval.h).
The minimum-version compatibility fixes above leave the numerical simulation,
likelihood, callback and RNG algorithms unchanged.
They exercise representative paths and do not prove the absence of all native
memory defects.

The pre-release diagnostic follow-up found an additional existing failure:
direct `Rcpp::warning()` calls can jump past RNG cleanup when a warning becomes
an error. Simulation warnings now use protected R evaluation with explicit RNG
synchronization. Shared explicit simulation scopes also preserve nested handler
draws. Local R 4.5.1 passed all nine forced-GC fixtures, including warning
unwinding, and the existing callback failure/interrupt recovery cases.
`tests/user-diagnostics.R` checks all three simulation loops with warnings
promoted to errors, RNG-using/nested handlers and original handler errors,
including Box-Muller normal-cache recovery. These follow-up changes still need
their own cross-platform, POSIX-signal and Valgrind CI run.
The follow-up source package also passed the full local R CMD check with
38 test files, vignette rebuilding and PDF/HTML manuals: Status: OK.

Run from the repository root after installing the package:

```sh
Rscript --vanilla tools/check-native-safety.R gc native-safety-results/gc
Rscript --vanilla tools/check-native-interrupt.R native-safety-results/interrupt
```

The interruption check requires POSIX signals and Python 3. Neither is a new
package runtime dependency; these development tools are excluded from the
source package.
