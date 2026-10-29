## emhawkes 1.0.2: no-long-double maintenance update

This update addresses the Additional issues/noLD failures reported by CRAN,
with deadline 2026-11-19. It is based on the released 1.0.1 source; unrelated
development features are not included.

The eight reported test failures were reproduced with the unchanged submission
tarball on R-devel 2026-10-26 r90766 configured with --disable-long-double.
The same R revision, compiler, package and dependency versions passed all
45 original tests when long-double support was enabled.

Native recursions now match the accumulator used by the running R build for
R-equivalent sums, row sums and cumulative sums. Precision is selected once
outside the event loops. Product rounding is preserved where R materializes
vectors before summing. Marked likelihood jumps also preserve a separate double multiplication and
addition on all builds, correcting a fused multiply-add discrepancy exposed by
the new regression on Linux compilers. Independent high-precision derivative
calculations remain unchanged.

The eight failing tests and an additional precision/initial-only RNG regression
test pass without relaxing existing numerical thresholds. Ordinary-build
before/after comparisons preserve 36 complete simulation paths and RNG states,
72 inference records, and two fitted coefficient/gradient/Hessian records exactly.
Both ordinary and noLD builds also pass all 46 tests with native fused
multiply-add optimization enabled. Public interfaces, native registration
signatures and saved-object layouts are unchanged.

Local R CMD check, including all 46 test files and vignette rebuilds:

- Windows / R 4.5.1: Status: OK.
- Ubuntu 24.04 / R-devel r90766, ordinary build: Status: OK.
- Ubuntu 24.04 / R-devel r90766, --disable-long-double: Status: OK.

Local checks use --no-manual --run-donttest and
_R_CHECK_SUGGESTS_ONLY_=false, with maxLik installed. The local PDF reference
manual was not compiled. A dedicated noLD workflow accompanies the existing
multi-platform and native-safety workflows.
