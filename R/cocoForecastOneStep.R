cocoForecastOneStep <- function(coco, max=NULL, epsilon=1e-12, xcast=NULL,
                                alpha=0.05,
                         decimals = 4, julia=FALSE) {
  
  seasonality <- c(1, 2) #will be used as argument in future versions
  
  
  if (!is.vector(xcast) & !is.null(xcast)){
    if (nrow(xcast) > 1){
      xcast <- xcast[1, ]
    }
  }
  
  if (is.null(max)){
    max_use <- 60
  } else{
    max_use <- max
  }
  
  if (!is.null(coco$julia_reg) & julia){
    addJuliaFunctions()
    coco_forecast <- JuliaConnectoR::juliaGet( JuliaConnectoR::juliaCall("cocoPredictOneStep", coco$julia_reg, 0:max_use, xcast))
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
    
    if ( !is.null(xcast) ){ 
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
  
  
  distribution_function <- cumsum(densities)
  densities_plot <- round(densities, decimals)
  
  x <- 0:(length(densities)-1)
  
  out <- list("density" = densities, "mode" = mode, mean=sum(densities*as.numeric(x)),
              "median" = median, "densities_plot" = densities_plot, "x" = x,
              "lower"=min(which(distribution_function >= alpha/2)) - 1,
              "upper"=min(which(distribution_function >= 1-alpha/2)) - 1,
              "data"=coco$ts)
  class(out) <- "cocoForecast"
  
  list_out <- list(out)
  class(list_out) <- "cocoForecastCollection"
  return(list_out) 
}
