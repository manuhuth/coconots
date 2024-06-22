#' @title Simulation of Count Time Series
#' @description The function generates a time series of low counts from the (G)PAR model class for a specified
#'  innovation distribution, sample size, lag order,
#' and parameter values. 
#' @param type character, either "Poisson" or "GP" indicating the type of the innovation distribution
#' @param order integer, either 1 or 2 indicating the order of the model
#' @param par numeric vector, the parameters of the model, the number of elements
#'  in the vector depends on the type and order specified. 
#' @param length integer, the number of observations in the generated time series
#' @param xreg data.frame, data frame of control variables
#' @param init numeric vector, initial data to use, default is NULL. See details
#' for more information on the usage.
#' @param julia If TRUE, the Julia implementation is used. In this case, init is ignored but it might be faster.
#' @param link_function Specifies the link function for the conditional mean of the innovation (\eqn{\lambda}). The default is `log`, but other available options include `identity` and `relu`. This parameter is applicable only when covariates are used. Note that using the `identity` link function may result in \eqn{\lambda} becoming negative. To prevent this, ensure all covariates are positive and restrict the parameter \eqn{\beta} to positive values.
#' @param julia_seed Seed for the Julia implementation. Only used if Julia equals TRUE.
#' @return a vector of the simulated time series. 
#' @details The function checks for valid input of the type, order, parameters, and initial data
#' before generating the time series.
#'
#'The init parameter allows users to set a custom burn-in period
#'for the simulation. By default, when simulating with covariates, no burn-in
#'period is specified since there is no clear choice on the covariates.
#'However, the init argument gives users the flexibility to select an
#'appropriate burn-in period for the covariate case. One way to do this is to
#'simulate a time series using \code{\link{cocoSim}} with appropriate covariates and pass the
#'resulting time series to the
#'init argument of a new \code{\link{cocoSim}} run so that the first time series is used as
#'the burn-in period. 
#'If init is not specified for the covariate case, a warning will be returned
#'to prompt the user to specify a custom burn-in period. This helps ensure that
#'the simulation accurately captures the dynamics of the system being modelled.
#' @author Manuel Huth
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' 
#' # Simulate using the RCPP implementation
#' data_rcpp <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
#' # Simulate using the Julia implementation
#' data_julia <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
#' @export

cocoSim <- function(type, order, par, length, xreg = NULL, init = NULL, 
                    julia=FALSE, julia_seed = NULL, link_function="log") {
  seasonality <- c(1, 2) #will be used as argument in future versions
  
  
  
  if (is.null(init)) {
    length_burn_in <- 200
    init_add <- 0
  } else {
    length_burn_in <- 0
    init_add <- order
  }
  
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
    
    if (julia){
      if (!is.null(julia_seed)){
        setJuliaSeed(julia_seed)
      }
      return(cocoSimJulia(type, order, par, length, xreg, link_function))
    }
    
    size <- length + length_burn_in + init_add
    output <- cocoSim_base(type = type, order = order, par = par, size = size,
                           seasonality = seasonality, init = init)
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
    
    if (julia){
      setJuliaSeed(julia_seed)
      return(cocoSimJulia(type, order, par, length, xreg, link_function))
    }
    size <- length + init_add
    output <- cocoSim_cov(
      type = type, order = order, par = par, size = size, xreg = xreg,
      seasonality = seasonality, init = init, link_function=link_function
    )
    
  }

  return(output$data)
}
