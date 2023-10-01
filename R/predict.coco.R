#' @title K-Step Ahead Forecast Bootstrapping
#' @description Computes the k-step ahead forecast using the models in the coconots package. 
#' 
#' @param object An object that has been fitted previously, of class coco.
#' @param k The number of steps ahead for which the forecast should be computed. Defaults to 3.
#' @param number_simulations The number of simulation runs to compute. Defaults to 500.
#' @param alpha Level of confidence that is used to construct the prediction intervals.
#' @param simulate_one_step_ahead If FALSE, the one-step ahead prediciton is obtained using the analytical distribution. If TRUE, bootstrapping is used.
#' @param max The maximum number of the forecast support for the plot. If NULL all values for which the cumulative distribution function is below 1- epsilon are used for the plot.
#' @param epsilon If max is NULL, epsilon determines the range of the support that is used by subsequent automatic plotting using R's plot() function.
#' @param xcast An optional matrix of covariate values for the forecasting. If `NULL`, the function assumes no covariates.
#' @param decimals Number of decimal places for the forecast probabilities
#' @param julia  if TRUE, the estimate is predicted with Julia.
#'
#' @return A list of frequency tables. Each table represents a k-step ahead forecast frequency distribution based on the simulation runs.
#' 
#' @details Returns forecasts for each mass point of the k-step ahead
#' distribution for the fitted model. The exact predictive distributions for
#' one-step ahead predicitons for
#' the models included here are provided in Jung and Tremayne (2011), maximum
#' likelihood estimates replace the true model parameters. Out-of-sample values
#' for covariates can be provided, if necessary.
#' @rdname predict.coco
#' @export
predict.coco <- function(object, k=1, number_simulations=1000,
                         alpha=0.05,
                         simulate_one_step_ahead = FALSE,
                         max=NULL, epsilon=1e-8, xcast=NULL,
                         decimals = 4,
                         julia=FALSE) {
  
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