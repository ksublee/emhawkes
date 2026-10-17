## Resubmission of emhawkes 1.0.0

The Windows R-devel pretest reported failures in fit-rds-session.R and
stats-isolated-session.R. Fresh sessions without methods attached could skip
hspec's S4 initializer, leaving its dimension unset; the isolated optimizer
test also used an unqualified new().

Import methods in the namespace and export the initialize and show methods.
Use methods::new() in isolated tests. Force methods to remain unattached in
the subprocess regression tests, covering construction, display, validation,
R/C++ likelihoods and saved-fit restoration with both optimizer providers.
Add Windows R-devel to the GitHub Actions check matrix.

The corrected source passed `R CMD check --as-cran --run-donttest` on
Windows/R 4.5.1 with 0 errors, 0 warnings and 0 notes, including all 43 test
files, rebuilt vignettes and PDF/HTML manuals. Saved fits from the original
1.0.0 candidate were restored and their likelihoods reproduced in fresh
sessions with methods unattached. The isolated stats test excludes maxLik.
Twelve simulation comparisons against the original candidate preserved all
outputs and final RNG states exactly.

## Update from emhawkes 0.9.8 to 1.0.0

This release adds registered Rcpp likelihood, inference and simulation engines,
an independent optional stats optimizer, fitted-model reconstruction and
curvature diagnostics. Unsupported models retain the R implementation. Rcpp
and a C++ compiler are now required for source installation; maxLik is optional
and remains the default optimizer provider.

Automatic initial intensity components now use unconditional means only for
supported stable constant linear models, otherwise zero. The selected rule is
reported. This corrects the previous univariate half-mean rule and can change
results; explicit initial states are preserved. NEWS documents this behavior
and the other numerical corrections. Stable Hessians do not imply nominal
finite-sample confidence-interval coverage.

## Checks

The implementation in commit 4737f96 passed:

- Windows, R 4.5.1: `R CMD check --as-cran --run-donttest`, including PDF/HTML
  manuals and rebuilt vignettes, with 0 errors, 0 warnings and 0 notes.
- GitHub Actions: Windows release and R 4.0.0; macOS release; Ubuntu release,
  oldrel-1 and devel. Release was R 4.6.1 and devel was R 4.7.0 (2026-10-15).
- Additional native garbage-collection, interruption and Valgrind checks.

The release-preparation source archive was checked again on Windows/R 4.5.1
with `--as-cran --run-donttest`, including rebuilt vignettes and PDF/HTML manuals:
0 errors, 0 warnings and 0 notes, both with and without maxLik installed.
For the latter check, `_R_CHECK_FORCE_SUGGESTS_=false` allowed the intentionally
missing suggestion; maxLik-specific sections were skipped. Independent stats,
likelihood, simulation and saved-result checks passed. The README example also
converged in a separate session with only emhawkes, Rcpp and R's own packages.

GitHub checks now run donttest examples and retain successful check results.
The Ubuntu release job additionally retains the checked source archive, SHA-256,
source commit and build session as the `cran-source-package` artifact.

## Reverse dependencies and URLs

CRAN's source package index was checked on 2026-10-17, including Depends,
Imports, LinkingTo, Suggests and Enhances. No reverse dependencies were found.
Package URLs were checked with urlchecker.

## Minimum-version CI note

The 0.9.9 development candidate (commit 94ea371) passed all six GitHub check
jobs. Windows/R 4.0.0 returned 0 errors, 0 warnings and 1 NOTE for generated
multi-architecture check outputs: emhawkes-Ex_i386.Rout, emhawkes-Ex_x64.Rout,
examples_i386, examples_x64, tests_i386 and tests_x64. Examples and tests passed
on both i386 and x64. These files are absent from the source archive. The other
five environments, including current release and devel, returned Status: OK.
The 1.0.0 version change preserves the checked R and C++ implementation.
