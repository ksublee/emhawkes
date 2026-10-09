# Realization of Hawkes Process

`hreal` is a list containing the following components:

- `hspec`: An S4 object of class
  [`hspec-class`](https://ksublee.github.io/emhawkes/reference/hspec-class.md)
  that specifies the parameter values.

- `inter_arrival`: The time intervals between consecutive events.

- `arrival`: The cumulative sum of `inter_arrival` times.

- `type`: An integer representing the type of event.

- `mark`: The size of the mark, providing additional information
  associated with the event.

- `N`: A counting process that tracks the number of events.

- `Nc`: A counting process that tracks the number of events, weighted by
  mark.

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

  The number of rows to display.

- ...:

  Additional arguments passed to or from other methods.

- object:

  An S3 object of class `hreal`.
