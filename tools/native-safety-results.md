# Minimum R and native safety checks

The package declares R >= 4.0.0. The R CMD check workflow now includes that
exact version on Windows, in addition to the existing release, devel and
old-release jobs. The minimum-version job checks the running version explicitly
and runs the package tests and vignettes.

The separate Native safety workflow runs on Linux and checks:

- Forced garbage collection during eight native fixtures: both numeric
  simulation paths, callback simulation, observed baseline, baseline Jacobian,
  state-dependent likelihood, observed-baseline likelihood and numeric score.
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

These checks change no simulation, likelihood, callback or RNG implementation.
They exercise representative paths and do not prove the absence of all native
memory defects.

Run from the repository root after installing the package:

```sh
Rscript --vanilla tools/check-native-safety.R gc native-safety-results/gc
Rscript --vanilla tools/check-native-interrupt.R native-safety-results/interrupt
```

The interruption check requires POSIX signals and Python 3. Neither is a new
package runtime dependency; these development tools are excluded from the
source package.
