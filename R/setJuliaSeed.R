#' @title Set Seed for Julia's Random Number Generator
#' @description Sets the seed for Julia's random number generator to ensure reproducibility.
#' @param julia_seed An integer seed value to be passed to Julia's random number generator.
#' @details This function initializes the necessary Julia functions and sets the random seed for Julia. 
#' If the provided seed is NULL, the function does nothing.
#' @author Manuel Huth
#' @export
setJuliaSeed <- function(julia_seed){
  addJuliaFunctions()
  if (!is.null(julia_seed)){
    JuliaConnectoR::juliaEval(paste0('Random.seed!(Int64(', julia_seed, '))', collapse=""))
  }
}