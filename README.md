  <!-- badges: start -->
  [![R-CMD-check](https://github.com/manuhuth/coconots/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/manuhuth/coconots/actions/workflows/R-CMD-check.yaml)
  [![Codecov test coverage](https://codecov.io/gh/manuhuth/coconots/branch/main/graph/badge.svg)](https://app.codecov.io/gh/manuhuth/coconots?branch=main)
  <!-- badges: end -->

# coconots
Likelihood-based methods for model fitting, assessment and prediction analysis of some convolution-closed count time series model are provided. The marginal distribution can be Poisson or Generalized Poisson. Regression effects can be modelled via time varying innovation rates.

![alt text](https://github.com/manuhuth/coconots/blob/main/images/functionality.png?raw=true)

# Details
The package allows simulation of convolution-closed count time series models with the cocoSim
function. Model fitting is performed with the cocoReg routine. By passing a cocoReg-type object,
cocoForecast computes the probability mass of the one-step ahead forecast. cocoBoot, cocoPIT,
cocoScore, and cocoResid provide routines for model assessment.

## Model

![alt text](https://github.com/manuhuth/coconots/blob/main/images/dgp.png?raw=true)
