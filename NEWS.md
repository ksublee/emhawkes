# emhawkes 0.9.6

Fixed some bugs, improved efficiency in this version and several features have been added.

## Breaking changes

* Slot `eta` is introduced which represents the constant part of `impact`.

* The concepts of `rambda` and `rambda_component` are introduced. They are closely related to the right-continuous version of the intensity process.

* For inference of intensity and goodness of fit, `infer_lambda` and `residual_process` functions are implemented.

* The method `vol` to measure the volatility is introduced and this feature is currently experimental.


# emhawkes 0.9.7

Fixed some bugs, improved efficiency in this version and several features have been added.

## Breaking changes

* The `lambda0` argument name used in previous versions has been changed to `lambda_component0` in this version. This is to clearly indicate the meaning of the argument and to avoid confusion.

* The name of method `vol` is changed to `hvol` and this feature is currently experimental.

* The Vignette file contains more examples and explanations.

