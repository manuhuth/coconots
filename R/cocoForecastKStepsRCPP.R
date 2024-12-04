cocoForecastKStepsRCPP <- function(fit, k=1, number_simulations=500, covariates=NULL){
  
  init <- fit$ts[length(fit$ts)]
  if (fit$order == 2){
    init <- fit$ts[(length(fit$ts)-1):length(fit$ts)]
  }
  
  run_sim <- function(i){ 
    #cannot be cocosim as the init argument is not available in cocoSim
    if (!is.null(covariates)){
      output <- cocoSim_cov(fit$type, order = fit$order, par = fit$par, 
                          size= k+length(init), xreg = covariates, 
                          seasonality = c(1, 2), init = init, link_function=fit$link_function)$data
    } else {
      output <- cocoSim_base(fit$type, order = fit$order, par = fit$par, 
                            size= k+length(init),
                            seasonality = c(1, 2), init = init)$data
    }
  }
  
  if (k > 1) {
    out_matrix <- t(sapply(1:number_simulations, run_sim))
  } else if (k==1){
    out_matrix <- matrix(sapply(1:number_simulations, run_sim), ncol=1)
  }
  
  freq_table <- function(i){
    df <- as.data.frame(table(out_matrix[,i]) / number_simulations)
    colnames <- c("value", "frequency")
    colnames(df) <- colnames
    df["value"] <- as.numeric(levels(df[, "value"]))

    return(df)
  }
  
  return(lapply(1:ncol(out_matrix), freq_table))
  
}
