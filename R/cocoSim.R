#' @title Simulation of time series data 
#' @description 
#' The cocoSim_base function generates a time series of a specified innovation distribution, size, order, and parameters. It also allows for the option of specifying the seasonality range and an initial set of data. 
#' The function checks for valid input of the type, order, parameters, seasonality, and initial data before generating the time series.
#' @param type character, either "Poisson" or "GP" indicating the type of the innovation distribution
#' @param order integer, either 1 or 2 indicating the order of the model
#' @param par numeric vector, the parameters of the model, the number of elements in the vector depends on the type and order specified
#' @param length integer, the number of observations in the generated time series
#' @param xreg data.frame, data frame of control variables
#' @param seasonality integer vector, the range of seasonality, default is c(1,2)
#' @param init numeric vector, initial data to use, default is NULL
#' @return list containing 'time' and 'data' which are the computation time and generated time series respectively
#' @author Manuel Huth
#' @export

cocoSim <- function(type, order, par, length, xreg = NULL, seasonality = c(1, 2), init = NULL) {
  if (is.null(xreg)) {
    if (order == 2){
      
      if (par[2] + par[4] > 1) {
        stop("The condition alpha_1 + alpha_3 < 1 ist not satisfied")
      }
      
      if ((par[2] < 0) | (par[3] < 0) | (par[4] < 0) | (par[2] >= 1) | (par[3] >= 1) | (par[4] >= 1)) {
        stop("The alpha parameters must be within the unit interval")
      }
      

      if ((type == "GP") & ((par[5] < 0) | (par[5] >= 1))) {
          stop("eta must be within the unit interval")
      }


    } else {
      if ((par[2] < 0) | (par[2] >= 1))  {
        stop("alpha must be within the unit interval")
      }
      if ((type == "GP") & ((par[3] < 0) | (par[3] >= 1))){
        stop("eta must be within the unit interval")
      }
    }
    
    if (is.null(init)) {
      length_burn_in <- 200
    } else {length_burn_in <- 0}
    
    size <- length + length_burn_in
    output <- cocoSim_base(type = type, order = order, par = par, size = size, seasonality = seasonality, init = init)
    output$data <- output$data[(length_burn_in+1):(length+length_burn_in)]
  } else {
    if (order == 2){
      
      if (par[1] + par[3] > 1) {
        stop("The condition alpha_1 + alpha_3 < 1 ist not satisfied")
      }
      
      if ((par[1] < 0) | (par[2] < 0) | (par[3] < 0) | (par[1] >= 1) | (par[2] >= 1) | (par[3] >= 1)) {
        stop("The alpha parameters must be within the unit interval")
      }

      if ((type == "GP") & ((par[4] < 0) | (par[4] >= 1))){
        stop("eta must be within the unit interval")
      }
      
    } else {
      if ((par[1] < 0) | (par[1] >= 1))  {
        stop("alpha must be within the unit interval")
      }
      if ((type == "GP") & ((par[2] < 0) | (par[2] >= 1))){
        stop("eta must be within the unit interval")
      }
    }
    output <- cocoSim_cov(
      type = type, order = order, par = par, size = length, xreg = xreg, seasonality = seasonality, init = init
    )
  }

  return(output)
}
