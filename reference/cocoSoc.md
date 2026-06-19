# Computes Scores for Various Models Maintaining a Common Sample

This function computes log, quadrtic and ranked probability scores for
Poisson and Generalized Poisson models.

## Usage

``` r
cocoSoc(
  data,
  models = "all",
  print.progress = TRUE,
  max_x_score = 50,
  julia = FALSE,
  ...
)
```

## Arguments

- data:

  A numeric vector containing the data to be used for modeling

- models:

  A character string specifying which models to use. Default is `"all"`,
  which uses both Poisson and GP models.

- print.progress:

  A logical value indicating whether to print progress messages
  (Default: `TRUE`).

- max_x_score:

  An integer which is used as the maximum count for the computation of
  the score (defaul: `50`)

- julia:

  if TRUE, `cocoSoc` is run with julia (default: FALSE)

- ...:

  Additional arguments to be passed to the `cocoReg` function.

## Value

A list of class `"cocoSoc"` containing:

- fits:

  A list of fitted model objects.

- scores_list:

  A list of score objects for each model.

- scores_df:

  A data frame containing the logarithmic, quadratic, and ranked
  probability scores for each model.

## Details

Supports model selection by computing score over a range of models while
maintaining a common sample and a common specification.

## Author

Manuel Huth
