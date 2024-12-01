cocoForecastKStepsRCPP <- function(fit, k=1, number_simulations=500, covariates=NULL){
  
  init <- fit$ts[length(fit$ts)]
  if (fit$order == 2){
    init <- fit$ts[(length(fit$ts)-1):length(fit$ts)]
  }
  
  run_sim <- function(i){ 
    output <- cocoSim(type = fit$type, order = fit$order, par = fit$par, length =  k+length(init),
            xreg = covariates,
            julia=FALSE, link_function=fit$link_function)
  }
  
  if (k > 1) {
    out_matrix <- t(sapply(1:number_simulations, run_sim))
  } else if (k==1){
    out_matrix <- as.matrix(sapply(1:number_simulations, run_sim)[1, ])
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
