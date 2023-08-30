#' @export
setJuliaSeed <- function(julia_seed){
  addJuliaFunctions()
  if (!is.null(julia_seed)){
    JuliaConnectoR::juliaEval(paste0('Random.seed!(Int64(', julia_seed, '))', collapse=""))
  }
}