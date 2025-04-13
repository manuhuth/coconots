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
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#' forecast <- predict(fit, k=1, simulate_one_step_ahead = FALSE)
#' plot(forecast[[1]]) #plot one-step ahead forecast distribution
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


#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoForecast <- function(object, breaks=NULL, width=0.1, ...){
  
  pl <- ggplot2::ggplot(mapping = ggplot2::aes(x = object$x, y = object$densities_plot)) +
    ggplot2::geom_bar(stat="identity", position="dodge", width=width) + 
    ggplot2::labs(title = "Probability Forecast", x = "Support", y = "Probability") +
    ggplot2::xlim(c(0, max(object$x))) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20)) 
  
  pl
}

#' @export
plot.cocoForecast <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  suppressWarnings({print(p)})
}

#' @export
plot.cocoForecast <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  suppressWarnings({print(p)})
}

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' @export
autoplot.cocoForecastCollection <- function(object, forecast_type="mode",
                                            n_lags_plot=30, ...){
  return_forecast <- function(i){
    L <- object[[i]]
    c(df_data[nrow(df_data), "t"]+i, NA, L[[forecast_type]], L[["lower"]], L[["upper"]])
  }
  
  relevant_data <- object[[1]]$data[max(1,(length(object[[1]]$data)-n_lags_plot+1)):length(object[[1]]$data)]
  df_data <- cbind(max(1,(length(object[[1]]$data)-n_lags_plot+1)):length(object[[1]]$data), relevant_data, NA, NA, NA)
  colnames(df_data) <- c("t", "x", "forecast", "lower", "upper")
  
  df <- as.data.frame(t(sapply(1:length(object), return_forecast)))
  colnames(df) <- c("t", "x", "forecast", "lower", "upper")
  
  df_plot <- rbind(df_data, df)
  
  df_transition <- df_plot[(nrow(df_plot)-length(object)):nrow(df_plot) ,]
  df_transition[1, "forecast"] <- df_transition[1, "x"] 
  
  pl <- ggplot2::ggplot(mapping = ggplot2::aes(x = df_plot[, "t"], y = df_plot[, "x"])) +
    ggplot2::geom_point(ggplot2::aes(color="Data")) + ggplot2::geom_line() + 
    ggplot2::geom_point(mapping = ggplot2::aes(x = df[, "t"], y = df[, "forecast"],
                                               color = "Prediction")) +
    ggplot2::geom_line(mapping = ggplot2::aes(x = df_transition[, "t"],
                                              y = df_transition[, "forecast"]),
                       color = "steelblue", linetype = "dotted") +
    ggplot2::geom_errorbar(mapping = ggplot2::aes(
      ymin = df_plot[, "lower"],
      ymax = df_plot[, "upper"]), color = "steelblue") +
    ggplot2::labs(title = paste0(forecast_type, " forecast"), x = "Time", y = "Count (x)") +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::scale_color_manual(name='',
                                breaks=c('Data', 'Prediction'),
                                values=c('Data'='black', 'Prediction'='steelblue'))
  
  pl
}

#' @export
plot.cocoForecastCollection <- function(x, ...) {
  p <- autoplot(
    x,
    ...
  )
  suppressWarnings({print(p)})
}

#'@export
print.cocoForecastCollection <- function(x, ...) {
  suppressWarnings({print(autoplot(x, ...,))})
  invisible(x)
}