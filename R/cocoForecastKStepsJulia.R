cocoForecastKStepsJulia <- function(fit, k=3, number_simulations=500, covariates=NULL){
  dict_out <- JuliaConnectoR::juliaCall("cocoPredictKsteps", fit$julia_reg,k, number_simulations, covariates)
  return_frequencies <- function(i){as.data.frame(dict_out[[paste0("prediction_", i)]])}
  return(lapply(1:dict_out[["length"]], return_frequencies))
}