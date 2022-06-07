# coconots
Likelihood-based methods for model fitting, assessment and prediction analysis of some convo- lution-closed count time series model are provided. The marginal distribution can be Poisson or Generalized Poisson. Regression effects can be modelled via time varying innovation rates.

# Details
The package allows simulation of convolution-closed count time series models with the cocoSim
function. Model fitting is performed with the cocoReg routine. By passing a cocoReg-type object,
cocoForecast computes the probability mass of the one-step ahead forecast. cocoBoot, cocoPIT,
cocoScore, and cocoResid provide routines for model assessment.
