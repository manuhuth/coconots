\name{coconots-package}
\alias{Analysis of Convolution-closed Nonlinear Count Time Series}
\alias{coconots}
\docType{package}
\title{
  Analysis of Convolution-closed Nonlinear Count Time Series XXX
}
\description{
  We provide R functions to conduct Likelihood-based methods for model fitting, assessment and prediction analysis of some convolution-closed count time series model. The marginal distribution can be Poisson or Generalized Poisson. Regression effects can be modelled via time varying innovation rates.
}
\details{
The package allows simulation of convolution-closed count time series models with the \code{\link{cocoSim}} function. Model fitting is performed with the \code{\link{cocoReg}} routine. By passing a cocoReg-type object, \code{\link{cocoForecast}} computes the probability mass of the one-step ahead forecast. \code{\link{cocoBoot}}, \code{\link{cocoPit}}, \code{\link{cocoScore}}, and \code{\link{cocoResid}} provide routines for model assessment. 


The main usage of the package is illustrated within the \code{\link{cocoReg}} function chapter. For more details and examples of the functions see the respective sections within this vignette. 

}
\author{
Maintainer: Manuel Huth <manuel.huth@yahoo.com>
}
\references{
Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment for count data. \emph{Biometrics} \bold{65}, 1254--61.

Jung, R. C. and Tremayne, A. R. (2010) Convolution-closed models for count time series with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.

Jung, R. C. and McCabe, B. P. and Tremayne, A. R. (2016) Model validation and diagnostics. \emph{Handbook of discrete-valued time series}, 189--218.

Joe, H. (1996) Time series models with univariate margins in the convolution-closed infinitely divisible class. \emph{Journal of Applied Probability}, 664--677.

Tsay, R. S. (1992) Model checking via parametric bootstraps in time series analysis. \emph{Applied Statistics} \bold{41}, 1--15.


}
