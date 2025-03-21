#' @importFrom JuliaConnectoR juliaEval
addJuliaFunctions <- function(){
  
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