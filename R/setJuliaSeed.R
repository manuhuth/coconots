#' @title Set Seed for \proglang{julia}'s Random Number Generator
#' @description Sets the seed for \proglang{julia}'s random number generator to ensure reproducibility.
#' @param julia_seed An integer seed value to be passed to \proglang{julia}'s random number generator.
#' @details This function initializes the necessary \proglang{julia} functions and sets the random seed for \proglang{julia}. 
#' If the provided seed is NULL, the function does nothing.
#' @importFrom JuliaConnectoR juliaEval
#' @author Manuel Huth
#' @export
setJuliaSeed <- function(julia_seed){
  addJuliaFunctions()
  if (!is.null(julia_seed)){
    JuliaConnectoR::juliaEval(paste0('Random.seed!(Int64(', julia_seed, '))', collapse=""))
  }
}