#' @title Forecast for COCO models
#' @description Computes the one-step ahead forecast distribution for the models included in the coconots package. 
#' @param coco An object of class coco.fit or coco.fit.c
#' @param max The maximum number of the forecast support for the plot. If NULL all values for which the cumulative distribution function is below 1- epsilon are used for the plot.
#' @param epsilon If max is NULL, epsilon determines how big the support of the forecast is for the plot.
#' @param xcast A vector of covariate values for forecasting (only required for coco.fit.c class)
#' @param plot A logical value indicating whether to plot the forecast
#' @param title Plot title
#' @param xlab X-axis label for the plot
#' @param ylab Y-axis label for the plot
#' @param width_bars Width of bars in the plot
#' @param seasonality A vector of two integers indicating the seasonalities of the time series
#' @param decimals Number of decimal places for the forecast probabilities
#' @param julia  if TRUE, the estimate is predicted with Julia.
#' @return A list containing the plot, density, mode, and median of the forecast
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)$data
#' #julia_installed = TRUE ensures that the fit object
#' #is compatible with the julia cocoForecast implementation 
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#'
#' #median, mode, and density forecasts - R implementation
#' forecast_r <- cocoForecast(fit)
#' @export


cocoForecast <- function(coco, max=NULL, epsilon=1e-5, xcast=NULL, plot = FALSE, title = "Probability mass",
             xlab="x", ylab="Probabilities", width_bars = 0.04, seasonality=c(1,2), decimals = 4, julia=FALSE) {
  
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
    x <- 0:(length(densities)-1)
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
    
    if (methods::is(coco, "coco.fit.c")) { 
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
  
  if (isTRUE(plot)) {
  plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = densities_plot)) + ggplot2::geom_bar(stat="identity", position="dodge", width=width_bars) + 
    ggplot2::labs(title = title, x = xlab, y = ylab) + ggplot2::scale_x_continuous(breaks=x) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20)) 
  }
  
  out <- list("plot" = plot, "density" = densities, "mode" = mode, "median" = median)
  
  return(out) 
}
