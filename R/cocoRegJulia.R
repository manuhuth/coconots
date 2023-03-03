cocoRegJulia <- function(type, order, data, xreg, start){
  addJuliaFunctions()
  return(JuliaConnectoR::juliaCall("cocoReg", type, order, data, xreg,start))
}