cocoSimJulia <- function(type, order, par, length, xreg){
  addJuliaFunctions()
  
  if (is.null(xreg)){
    par_use <- c(par[2:length(par)], par[1])
  } else {
    par_use <- par
  }
  
  return(JuliaConnectoR::juliaCall("cocoSim", type, order, par_use, length, xreg))
}