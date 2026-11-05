#' @importFrom JuliaConnectoR juliaGet
#' @importFrom utils head
cocoForecastOneStep <- function(coco, max=NULL, epsilon=1e-12, xcast=NULL,
                                alpha=0.05,
                         decimals = 4, julia=FALSE) {
  
  seasonality <- c(1, 2) #will be used as argument in future versions
  
  
  if (!is.vector(xcast) & !is.null(xcast)){
    if (nrow(xcast) > 1){
      xcast <- xcast[1, ]
    }
  }
  
  # the predictive pmf on 0..max needs factorials up to max, and the
  # forecast conditions on the last `order` observations
  coco_check_support_limit(max, "max")
  coco_check_count_limit(coco$ts)
  coco_check_length(length(coco$ts), coco$order, "predict",
                    sprintf("an order-%d model conditions on the last %d observation(s)",
                            coco$order, coco$order))

  if (is.null(max)){
    max_use <- 60
  } else{
    max_use <- max
  }
  
  use_julia <- !is.null(coco$julia_reg) & julia
  if (use_julia){
    addJuliaFunctions()
    if ((is.matrix(xcast))) {
      xcast <- c(xcast)
    }
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
      lambda <- applyLinkFunction(dot_product, coco$link_function)
      parameter <- c(lambda, parameter)
    }
    
    type <- coco$type
    order <- coco$order
    if (type == "Poisson"){parameter <- c(parameter,0)}
    if (order == 1){
      validate_gp1_par(parameter[1], parameter[2], parameter[3])
    } else {
      validate_gp2_par(parameter[1], parameter[2], parameter[3], parameter[4], parameter[5])
    }
  }
  
  # Predictive pmf on the support 0:m.
  pmf_on_support <- function(m){
    if (use_julia){
      return(JuliaConnectoR::juliaGet(JuliaConnectoR::juliaCall("cocoPredictOneStep",
                                                                coco$julia_reg, 0:m, xcast)))
    }
    if (order == 1){
      return(as.numeric(dGP1h_vec(m, y, parameter[1], parameter[2], parameter[3])))
    }
    as.numeric(dGP2h_vec(m, y, z, parameter[1], parameter[2],
                         parameter[3], parameter[4], parameter[5]))
  }
  get_densities <- function(res){
    if (use_julia) getJuliaValue(res, "probabilities") else res
  }
  
  # Without an explicit max, the support is doubled until it carries at least
  # 1 - epsilon of the mass (heavy GP tails need more than the initial 60
  # points), until a non-finite density appears (the recursions can overflow
  # for counts close to the limit) or until max_support is reached. The
  # support ends at the largest supported count COCO_MAX_COUNT on both
  # backends (the C++ factorials overflow above it).
  max_support <- COCO_MAX_COUNT
  res <- pmf_on_support(max_use)
  if (is.null(max)){
    repeat {
      densities <- get_densities(res)
      finite <- is.finite(densities)
      if (!all(finite) || max_use >= max_support ||
          sum(densities) >= 1 - epsilon) break
      max_use <- min(2 * max_use, max_support)
      res <- pmf_on_support(max_use)
    }
  }
  densities <- get_densities(res)
  if (!is.null(max)){
    coco_check_finite_pmf(densities, "predictive pmf")
  }
  
  if (is.null(max)){
    # keep the leading run of finite densities, then truncate at 1 - epsilon
    nonfinite <- which(!is.finite(densities))
    if (length(nonfinite) > 0){
      densities <- densities[seq_len(max(1, min(nonfinite) - 1))]
    }
  }
  
  if (use_julia){
    mode <- getJuliaValue(res, "prediction_mode")
    median <- getJuliaValue(res, "prediction_median")
  } else {
    mode <- match(max(densities), densities) - 1
    distribution_function <- cumsum(densities)
    median <- min(which(distribution_function >= 0.5)) - 1
  }
  
  if (is.null(max)){
    cumulative <- cumsum(densities)
    reached <- which(cumulative >= 1-epsilon)
    index_use <- if (length(reached) > 0) min(reached) else length(densities)
    densities <- densities[1:index_use]
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
