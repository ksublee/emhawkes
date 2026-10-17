# Realization of Hawkes Process

`hreal` is an S3 list returned by
[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md),
containing the following components:

- `hspec`: An S4 object of class
  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
  that specifies the parameter values.

- `inter_arrival`: The time intervals between consecutive events.

- `arrival`: The cumulative sum of `inter_arrival` times.

- `type`: An event-type vector. Values 1, ..., `dimens` identify
  intensity rows.

- `mark`: A vector of event marks.

- `N`: A cumulative event-count matrix, with one column per event type.

- `Nc`: A cumulative mark-total matrix, with one column per event type.

- `lambda`: The left-continuous intensity process.

- `lambda_component`: The component of the intensity process,
  \\\lambda\_{ij}\\, that excludes `mu`.

- `rambda`: The right-continuous intensity process.

- `rambda_component`: The right-continuous version of
  `lambda_component`.

Functions for printing `hreal` objects are provided.

## Usage

``` r
# S3 method for class 'hreal'
print(x, n = 20, ...)

# S3 method for class 'hreal'
summary(object, n = 20, ...)

# S3 method for class 'hreal'
as.matrix(x, ...)
```

## Arguments

- x:

  An S3 object of class `hreal`.

- n:

  The number of rows to display: a nonnegative integer or `Inf`.

- ...:

  Additional arguments passed to or from other methods.

- object:

  An S3 object of class `hreal`.

## Details

The first row is the initial state, with zero `inter_arrival`, `type`,
and `mark`. Initial counts and mark totals are stored in the first rows
of `N` and `Nc`. Subsequent counting-process rows contain values after
each event. Kernel-component columns are ordered by intensity row, then
kernel column.
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
with `infer = TRUE` returns an inference list, with likelihood and
residual fields, rather than an `hreal` simulation result.

## See also

[`hsim()`](https://ksublee.github.io/emhawkes/reference/hsim.md),
[`logLik()`](https://ksublee.github.io/emhawkes/reference/logLik.md)
