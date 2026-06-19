# Probability Integral Transform Based Model Assessment Procedure

Computes the probability integral transform (PIT) and provides the
non-randomized PIT histogram for assessing absolute performance of a
fitted model as proposed by Czado et al. (2009).

## Usage

``` r
cocoPit(coco, J = 10, conf.alpha = 0.05, julia = FALSE)
```

## Arguments

- coco:

  An object of class coco

- J:

  Number of bins for the histogram (default: 10)

- conf.alpha:

  Significance level for the confidence intervals (default: 0.05)

- julia:

  if TRUE, the PIT is computed with julia (default: FALSE)

## Value

an object of class cocoPit. It contains the probability integral
transform values, p-value of the chi-square goodness of fit test and
information on the model specifications.

## Details

The adequacy of a distributional assumption for a model is assessed by
checking the cumulative non-randomized PIT distribution for uniformity.
A useful graphical device is the PIT histogram, which displays this
distribution to J equally spaced bins. We supplement the graph by
incorporating approximately \\100(1 - \alpha)\\\\ confidence intervals
obtained from a standard chi-square goodness-of-fit test of the null
hypothesis that the J bins of the histogram are drawn from a uniform
distribution. For details, see Jung, McCabe and Tremayne (2016).

## References

Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment
for count data. *Biometrics* **65**, 1254–61.

Jung, R. C., McCabe, B.P.M. and Tremayne, A.R. (2016). Model validation
and diagnostics. *In Handbook of Discrete Valued Time Series*. Edited by
Davis, R.A., Holan, S.H., Lund, R. and Ravishanker, N.. Boca Raton:
Chapman and Hall, pp. 189–218.

Jung, R. C. and Tremayne, A. R. (2011) Convolution-closed models for
count time series with applications. *Journal of Time Series Analysis*,
**32**, 3, 268–280.

## Author

Manuel Huth

## Examples

``` r
lambda <- 1
alpha <- 0.4
set.seed(12345)
data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
fit <- cocoReg(order = 1, type = "Poisson", data = data)

#PIT R implementation
pit_r <- cocoPit(fit)
plot(pit_r)
```
