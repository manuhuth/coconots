setJuliaSeed <- function(julia_seed){
  addJuliaFunctions()
  if (!is.null(julia_seed)){
    JuliaConnectoR::juliaEval(paste0('Random.seed!(', julia_seed, ')', collapse=""))
  }
}