# CRAN noLD maintenance validation

## Reproduction

Baseline: the submitted emhawkes 1.0.1 tarball, SHA256
`01f5115204c3facd072120299c8fc4da137bd50af8d4b191e5921c486ca98453`,
source commit `092607523943f5db45424d8d7c3a48515826ccba`.
R-devel 2026-10-26 r90766 was built twice on Ubuntu 24.04 with GCC 13.3:
ordinary and --disable-long-double. Package/dependency versions match.
The noLD setting is verified through capabilities("long.double"), not inferred
from the compiler's long-double type. CRAN uses Fedora/GCC 16; the reproduced
eight failing test names match its log exactly.

| Source | R build | Tests passed | Failed |
|---|---|---:|---:|
| 1.0.1 | ordinary | 45 | 0 |
| 1.0.1 | noLD | 37 | 8 |
| 1.0.2 maintenance patch | ordinary | 46 | 0 |
| 1.0.2 maintenance patch | noLD | 46 | 0 |

## Repair

The native simulation and likelihood entry points select a typed accumulator
once, using the running R build's long-double capability. R-equivalent totals,
row totals, validation and residual cumulative sums follow that selection.
Double-precision products are materialized before accumulation to preserve R's
rounding and prevent fused multiply-add contraction. Derivative intermediates
that are not reproducing an R sum retain their existing extended precision.
The exported R/native function signatures and registration files are unchanged.

## Regression evidence

- Original eight noLD failures pass with their original thresholds.
- The added native-precision test covers cancellation, an initial-only call
  without creating .Random.seed, the first mapped simulation mismatch, marked
  callbacks and independent native likelihood reconstruction.
- Ordinary-build before/after fixtures are identical: 36 complete paths,
  final RNG and following draws, 72 complete inference records, and two fitted
  coefficient/gradient/Hessian records.
- Windows R 4.5.1 and both Linux builds pass R CMD check, including all 46 test
  files and vignette rebuilds, with --no-manual --run-donttest and
  _R_CHECK_SUGGESTS_ONLY_=false. PDF reference manuals are not part of this local run.

The dedicated no-long-double workflow builds the pinned R revision with the
relevant configure option and runs package tests and vignette checks. Standard
multi-platform and native-safety checks also run on maintenance branches.

This release excludes unrelated shared-kernel/FRPP development and preserves the
primary checkout's existing uncommitted work. No performance claim is based on
test-suite elapsed times.
