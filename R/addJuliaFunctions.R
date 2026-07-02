#' @importFrom JuliaConnectoR juliaEval
addJuliaFunctions <- function(){
  set_julia_threads()

  if (JuliaConnectoR::juliaEval("isdefined(Main, :cocoReg)")) {
  } else {
    #----------------Add Julia functions-----------------------------------
    JuliaConnectoR::juliaEval('
      using Coconots
      using Random
    ')
    #-------------------------------------
  }
}