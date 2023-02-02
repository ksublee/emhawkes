## Resubmission

Fixed some bugs, improved efficiency in this version and several features have been added.

-   Slot `eta` is introduced which represents the constant part of `impact`.

-   The concepts of `rambda` and `rambda_component` are introduced. They are closely related to the right-continuous version of the intensity process.

-   For inference of intensity and goodness of fit, `infer_lambda` and `residual_process` functions are implemented.

-   The method `vol` to measure the volatility is introduced and this feature is currently experimental.

## Test environment

-   local Windows install, R 4.2.1
-   windows-x86_64-devel (r-devel)
-   ubuntu-gcc-release (r-release)
-   fedora-clang-devel (r-devel)

## R CMD check results

The following is the result of R CMD check performed via `devtools:check()` in local machine:

─ using R version 4.2.1 (2022-06-23 ucrt)\
─ using platform: x86_64-w64-mingw32 (64-bit)\
─ using session charset: UTF-8\
─ using options '--no-manual --as-cran'

── R CMD check results ──────── emhawkes 0.9.7 ──── Duration: 2m 10.2s

0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔

R CMD check succeeded

In addition, `rhub::check_for_cran()` is also performed with the following summary:

\## Test environments

-   R-hub windows-x86_64-devel (r-devel)
-   R-hub ubuntu-gcc-release (r-release)
-   R-hub fedora-clang-devel (r-devel)

\## R CMD check results

❯ On windows-x86_64-devel (r-devel) checking for detritus in the temp directory ... NOTE Found the following files/directories: 'lastMiKTeXException'

❯ On fedora-clang-devel (r-devel) checking HTML version of manual ... NOTE Skipping checking HTML validation: no command 'tidy' found

0 errors ✔ \| 0 warnings ✔ \| 2 notes ✖

First, according to R-hub [issue #503](https://github.com/r-hub/rhub/issues/503), this could be a bug in R-hub and can be disregarded.

Second note might also seem from the problem in R-hub test platform of Fedora [issue #548](https://github.com/r-hub/rhub/issues/548).
