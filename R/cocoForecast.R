#' @title Forecast for COCO models
#' @description Computes forecast for a COCO model. The function can handle both Poisson and GP models, with order 1 or 2. The function also has options for the maximum number of forecast steps, covariates for forecasting, the plot, plot title, axis labels, width of bars, seasonality, and decimal precision.
#' @param coco An object of class coco.fit or coco.fit.c
#' @param max The maximum number of forecast steps
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


cocoForecast <- function(coco, max=15, xcast=NULL, plot = FALSE, title = "Probability mass",
             xlab="x", ylab="Probabilities", width_bars = 0.04, seasonality=c(1,2), decimals = 4, julia=FALSE) {
  
  if (!is.null(coco$julia_reg) & julia){
    addJuliaFunctions()
    coco_forecast <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("cocoPredict", coco$julia_reg, 0:max, xcast))
    densities <- coco_forecast$values[[4]]
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
    
    x <- 0:max
    
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
