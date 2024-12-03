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
    df_temp <- forecasts[[i]]
    df_temp$value <- as.numeric(df_temp$value)
    
    # Generate the full sequence of values
    all_values <- data.frame(value = 0:max(df_temp$value))
    
    # Merge with the original data, filling missing frequencies with 0
    df_complete <- merge(all_values, df_temp, by = "value", all.x = TRUE)
    df_complete$frequency[is.na(df_complete$frequency)] <- 0
    
    
    
    densities <- df_complete[,"frequency"]
    x <- df_complete[,"value"]
    
    
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