# Residual Based Model Assessment Procedure

Calculates the (Pearson) residuals of a fitted model for model
evaluation purposes.

## Usage

``` r
cocoResid(coco, val.num = 1e-10)
```

## Arguments

- coco:

  An object of class "coco

- val.num:

  A non-negative real number that halts the calculation once the
  cumulative probability reaches 1-`val.num`

## Value

a list that includes the (Pearson) residuals, conditional expectations,
conditional variances, and information on the model specifications.

## Details

The Pearson residuals are computed as the scaled deviation of the
observed count from its conditional expectation given the relevant past
history, including covariates, if applicable. If a fitted model is
correctly specified, the Pearson residuals should exhibit mean zero,
variance one, and no significant serial correlation.

## Author

Manuel Huth
