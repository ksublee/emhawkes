# Shared-kernel API publication

2026-10-18, Windows / R 4.5.1. GitHub development version: 1.0.1.9000.
The previously submitted 1.0.1 source archive is unchanged.

## Public boundary

Export `exp_kernel_spec()`, `as_exp_kernel()` and `exp_kernel_replay()`.
Numeric alpha/beta/optional eta, positive row-wise equal beta, multivariate and
multi-kernel mappings, observed marks and signed components are supported.
The initial row has zero integral and no event update. Initial state is explicit;
full component histories are optional. Baseline, counts, residual law,
initialization policy, likelihood, optimization and RNG belong to callers.
Advance/update helpers remain internal. See the website shared-kernel guide
and function help for the complete matrix/order/zero-interval contract.

## Exact publication source

Three mixed files (DESCRIPTION, NAMESPACE and NEWS) were staged selectively
from the previous committed source. Only the new API/version/registration
changes were included. The candidate was archived from that index into a
separate source directory before compilation and testing.

Pending simulate/observation-count changes, roxygen metadata edits, other
working documents and compiled artifacts were preserved outside this commit.
The 33 pre-existing R/src Git blobs exactly match the parent commit. The validation
archive has the same content after platform line-ending normalization. Existing event loops and emhawkes.engine semantics did not
change. No native ABI or third package was added.

## Local results

| Check | Result |
|---|---|
| emhawkes exact publication source | R CMD build and R CMD check --no-manual: Status OK; 47 test files and vignette rebuilding included |
| frpp with publication installation | Build and full check --no-manual: Status OK; eight test files and both vignettes/rebuild included |
| Existing versus publication calculations | 40 path/inference/final-RNG/following-draw records exactly equal |
| stats and maxLik fitting | Both records exactly equal: coefficients, likelihood attributes, raw Hessian, covariance, initial state, diagnostics and RNG |
| Kernel checks | 72 exact R/C++ state comparisons, all-past-event reconstruction and independent quadrature passed; integral tolerance remains 1e-12 |
| Default joint quick suite | 9/9 files; 14.34 s test loop, 15.22 s full command including parent R startup/report preparation |
| Web documentation | Full local build passed; new guide and both reference pages checked |

Checks used `_R_CHECK_SUGGESTS_ONLY_=false`. No numerical tolerance was widened.
No speed improvement or confidence-interval coverage is claimed. PDF manuals
were not compiled. The quick-suite elapsed observation is not a benchmark.

## Integration adjustments

The earlier quick runner also selected the separately pending simulate/nobs
work. Its default now selects all eight FRPP files plus the public kernel test;
`--extended` explicitly requests the two additional development tests. Missing
selected files produce an error listing their paths rather than being skipped.
FRPP's base emhawkes dependency remains compatible with its bridge; state
calculations detect the public development API. Guides now document installing
GitHub development emhawkes 1.0.1.9000 for state calculations.

The first full website build hit Windows PATH_MAX in the deeply nested isolated
source path. Building the same source from a shorter temporary path passed;
no package calculation or website-building implementation was changed.

GitHub checks and deployment are recorded separately after publication.
