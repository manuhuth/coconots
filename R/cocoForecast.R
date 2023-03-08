#' @title One-Step Ahead Forecast Distribution
#' @description Computes the one-step ahead forecast distribution for the models included in the coconots package. 
#' @param coco An object of class coco
#' @param max The maximum number of the forecast support for the plot. If NULL all values for which the cumulative distribution function is below 1- epsilon are used for the plot.
#' @param epsilon If max is NULL, epsilon determines the range of the support that is used by subsequent automatic plotting using R's plot() function.
#' @param xcast A vector of covariate values for forecasting 
#' @param decimals Number of decimal places for the forecast probabilities
#' @param julia  if TRUE, the estimate is predicted with Julia.
#' @return A list containing the probability mass, mode, and median of the forecast
#' @details Returns forecasts for each mass point of the one-step ahead
#' distribution for the fitted model. The exact predictive distributions for
#' the models included here are provided in Jung and Tremayne (2011), maximum
#' likelihood estimates replace the true model parameters. Out-of-sample values
#' for covariates can be provided, if necessary.
#' 
#' Point forecasts are provided by returning the median and the mode of the predictive distribution. 
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
#' #julia_installed = TRUE ensures that the fit object
#' #is compatible with the julia cocoForecast implementation 
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#'
#' #median, mode, and density forecasts - R implementation
#' forecast_r <- cocoForecast(fit)
#' @export


cocoForecast <- function(coco, max=NULL, epsilon=1e-5, xcast=NULL,
                         decimals = 4, julia=FALSE) {
  
  seasonality <- c(1, 2) #will be used as argument in future versions
  
  if (is.null(max)){
    max_use <- 60
  } else{
    max_use <- max
  }
  
  if (!is.null(coco$julia_reg) & julia){
    addJuliaFunctions()
    coco_forecast <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("cocoPredict", coco$julia_reg, 0:max_use, xcast))
    densities <- coco_forecast$values[[4]]
    if (is.null(max)){
      cumulative <- cumsum(densities)
      index_use <- min(which(cumulative >= 1-epsilon))
      densities <- densities[1:index_use]
    }
    
    mode <- coco_forecast$values[[1]]
    median <- coco_forecast$values[[2]]
  } else {
    if (length(seasonality == 1)) {
      seasonality = c(seasonality, seasonality + 1)
    }
    
    data <- coco$ts
    y <- data[length(data) - seasonality[1] + 1]
    z <- data[length(data) - seasonality[2] + 1]
    parameter <-coco$par
    
    if ( !is.null(coco$cov) ){ 
      number_covariates <- ncol(coco$cov) 
      betas <- parameter[(length(parameter)-number_covariates+1):length(parameter)]
      parameter <- utils::head(parameter, -number_covariates)
      dot_product <- betas %*% c(xcast)
      lambda <- exp(dot_product)
      parameter <- c(lambda, parameter)
    }
    
    type <- coco$type
    order <- coco$order
    
    densities <- c()
    
    x <- 0:max_use
    
    if (order == 1){
      if (type == "Poisson"){parameter <- c(parameter,0)}
      for (number in x) { 
        densities[number+1] <- dGP1(x=number, y=y, par=parameter)
      }
    } 
    
    if (order == 2){
      if (type == "Poisson"){parameter <- c(parameter,0)}
      for (number in x) { 
        densities[number+1] <- dGP2(x=number, y=y, z=z, par=parameter)
      }
    }
    
    mode <- match(max(densities), densities) - 1
    distribution_function <- cumsum(densities)
    median <- min(which(distribution_function >= 0.5)) - 1
    if (is.null(max)){
      cumulative <- cumsum(densities)
      index_use <- min(which(cumulative >= 1-epsilon))
      densities <- densities[1:index_use]
    }
  }
  
  densities_plot <- round(densities, decimals)
  
  x <- 0:(length(densities)-1)
  
  out <- list("plot" = plot, "density" = densities, "mode" = mode,
              "median" = median, "densities_plot" = densities_plot, "x" = x)
  class(out) <- "cocoForecast"
  return(out) 
}
