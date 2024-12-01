cocoForecastKSteps <- function(fit, k=3, number_simulations=1000, alpha=0.05, covariates=NULL,
                               decimals=4, 
                               julia=FALSE){
  if (julia){

    forecasts <- cocoForecastKStepsJulia(fit, k=k, number_simulations = number_simulations,
                                   covariates=covariates)
    
  } else {
    forecasts <- cocoForecastKStepsRCPP(fit, k=k, number_simulations = number_simulations,
                                  covariates=covariates)
  }
  
  make_class <- function(i){
    densities <- forecasts[[i]][,"frequency"]
    x <- as.numeric(forecasts[[i]][,"value"])
    
    if (min(x) > 0){
      x_fill <- 0:(min(x)-1)
      dens <- rep(0, (min(x)))
      x <- c(x_fill, x)
      densities <- c(dens, densities)
    }
    
    mode <- match(max(densities), densities) - 1
    distribution_function <- cumsum(densities)
    median <- min(which(distribution_function >= 0.5)) - 1
    densities_plot <- round(densities, decimals)
    
    out <- list("density" = densities, "mode" = mode, mean=sum(densities*x),
                "median" = median, "densities_plot" = densities_plot, "x" = x,
                "k" = i,
                "lower"=min(which(distribution_function >= alpha/2)) - 1,
                "upper"=min(which(distribution_function >= 1-alpha/2)) - 1,
                "data"=fit$ts)
    class(out) <- "cocoForecast"
    return(out)
  }
  
  list_out <- lapply(1:length(forecasts), make_class)
  class(list_out) <- "cocoForecastCollection"
  return(list_out)
}