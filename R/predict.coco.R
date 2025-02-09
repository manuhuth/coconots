#' @title K-Step Ahead Forecast Distributions
#' @description Computes the k-step ahead forecast (distributions) using the models in the coconots package. 
#' 
#' @param object An object that has been fitted previously, of class coco.
#' @param k The number of steps ahead for which the forecast should be computed (Default: 1).
#' @param number_simulations The number of simulation runs to compute (Default: 1000).
#' @param alpha Significance level used to construct the prediction intervals (Default: 0.05).
#' @param simulate_one_step_ahead If FALSE, the one-step ahead prediction is obtained using the analytical predictive distribution. If TRUE, bootstrapping is used.
#' @param max The maximum number of the forecast support for the plot. If NULL all values for which the cumulative distribution function is below 1- epsilon are used for the plot.
#' @param epsilon If max is NULL, epsilon determines the range of the support that is used by subsequent automatic plotting using R's plot() function.
#' @param xcast An optional matrix of covariate values for the forecasting. If `NULL`, the function assumes no covariates.
#' @param decimals Number of decimal places for the forecast probabilities
#' @param julia  if TRUE, the estimate is predicted with \proglang{julia} (Default: FALSE).
#' @param ... Optional arguments.
#'
#' @return A list of frequency tables. Each table represents a k-step ahead forecast frequency distribution based on the simulation runs.
#' 
#' @details Returns forecasts for each mass point of the k-step ahead
#' distribution for the fitted model. The exact predictive distributions for
#' one-step ahead predictions for
#' the models included here are provided in Jung and Tremayne (2011), maximum
#' likelihood estimates replace the true model parameters. 
#' For k>1 forecast distributions are estimated using a parametric bootstrap. See Jung and Tremanye (2006).
#' Out-of-sample values
#' for covariates can be provided, if necessary.
#' 
#' 
#' for k > 1 
#' 
#' @examples
#' length <- 500
#' pars <- c(1, 0.4)
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = pars, length = length)
#' fit <- cocoReg(order = 1, type = "Poisson", data = data, julia_installed=TRUE)
#' forecast <- predict(fit, k=1, julia = TRUE, simulate_one_step_ahead = F)
#' plot(forecast[1]]) #plot one-step ahead forecast distribution
#' 
#' @references 
#' Jung, R.C. and Tremayne, A. R. (2011) Convolution-closed models for count time series with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 3, 268--280.
#' 
#' Jung, R.C. and Tremayne, A.R. (2006) Coherent forecasting in integer time series models. 
#'  \emph{International Journal of Forecasting} \bold{22}, 223--238
#' @rdname predict.coco
#' @export
predict.coco <- function(object, k=1, number_simulations=1000,
                         alpha=0.05,
                         simulate_one_step_ahead = FALSE,
                         max=NULL, epsilon=1e-8, xcast=NULL,
                         decimals = 4,
                         julia=FALSE,...) {
  
  if ((k == 1) & (!is.matrix(xcast)) & (!is.null(xcast))){
    xcast <- matrix(xcast, nrow=1)
  }
  
  if ((k==1) & (!simulate_one_step_ahead)){
    return(cocoForecastOneStep(object, max=max, epsilon=epsilon, xcast=xcast,
                               alpha=alpha,
                    decimals=decimals, julia=julia))
  } else if ((k > 1) & (!simulate_one_step_ahead))  {
    
    k_steps <- cocoForecastKSteps(object, k=k, number_simulations=number_simulations,
                       covariates=xcast, alpha=alpha,
                       decimals=decimals, 
                       julia=julia)
    k_steps[[1]] <- cocoForecastOneStep(object, max=max, epsilon=epsilon, xcast=xcast,
                                        alpha=alpha,
                                        decimals=decimals, julia=julia)[[1]]
    return(k_steps)
  } else {
    return(cocoForecastKSteps(object, k=k, number_simulations=number_simulations,
                              covariates=xcast, alpha=alpha,
                              decimals=decimals, 
                              julia=julia))
  }
}