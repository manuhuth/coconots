#' @importFrom JuliaConnectoR juliaCall
cocoSimJulia <- function(type, order, par, length, xreg, link_function){
  addJuliaFunctions()
  
  if (is.null(xreg)){
    par_use <- c(par[2:length(par)], par[1])
  } else {
    par_use <- par
  }
  
  return(JuliaConnectoR::juliaCall("cocoSim", type, order, par_use, length, xreg, link_function))
}