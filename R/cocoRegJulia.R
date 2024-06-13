cocoRegJulia <- function(type, order, data, xreg, start, link_function, lower_bound_covariates){
  addJuliaFunctions()
  return(JuliaConnectoR::juliaCall("cocoReg", type, order, data, xreg, start, link_function, lower_bound_covariates))
}