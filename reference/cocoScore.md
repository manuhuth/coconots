# Scoring Rule Based Model Assessment Procedure

The function computes log, quadratic and ranked probability scores for
assessing relative performance of a fitted model.

## Usage

``` r
cocoScore(coco, max_x = 50, julia = FALSE)
```

## Arguments

- coco:

  An object of class coco

- max_x:

  An integer which is used as the maximum count for the computation of
  the score (default: `50`)

- julia:

  if TRUE, the scores are computed with julia (default: FALSE).

## Value

a list containing the log score, quadratic score and ranked probability
score.

## Details

Scoring rules assign a numerical score based on the predictive
distribution and the observed data to measure the quality of
probabilistic predictions. They are provided here as a model selection
tool and are computed as averages over the relevant set of (in-sample)
predictions. Scoring rules are, generally, negatively oriented penalties
that one seeks to minimize. The literature has developed a large number
of scoring rules and, unless there is a unique and clearly defined
underlying decision problem, there is no automatic choice of a (proper)
scoring rule to be used in any given situation. Therefore, the use of a
variety of scoring rules may be appropriate to take advantage of
specific emphases and strengths. Three proper scoring rules (for a
definition of the concept of propriety see Gneiting and Raftery, 2007),
which Jung, McCabe and Tremayne (2016) found to be particularly useful,
are implemented. For more information see the references listed below.

## References

Czado, C. and Gneitling, T. and Held, L. (2009) Predictive Model
Assessment for Count Data. *Biometrics*, **65**, 1254–1261.

Gneiting, T. and Raftery, A. E. (2007) Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association*, 102:359-378.

Jung, R. C., McCabe, B.P.M. and Tremayne, A.R. (2016). Model validation
and diagnostics. *In Handbook of Discrete Valued Time Series*. Edited by
Davis, R.A., Holan, S.H., Lund, R. and Ravishanker, N.. Boca Raton:
Chapman and Hall, pp. 189–218.

Jung, R. C. and Tremayne, A. R. (2011) Convolution-closed models for
count time series with applications. *Journal of Time Series Analysis*,
**32**, 268–280.

## Author

Manuel Huth

## Examples

``` r
lambda <- 1
alpha <- 0.4
set.seed(12345)
data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
fit <- cocoReg(order = 1, type = "Poisson", data = data)

# scoring rules - R implementation
score_r <- cocoScore(fit)
```
