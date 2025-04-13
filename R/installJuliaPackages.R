#' @title installJuliaPackages
#' @description checks for needed \proglang{julia} packages and installs them if not installed.
#' @importFrom JuliaConnectoR juliaEval
#' @return no return value, called to install \proglang{julia} packages in \proglang{julia}.
#' @export
installJuliaPackages <- function(){
  
  strings1 <- c('"Coconots"')
  strings2 <- c('Pkg.add("Coconots")')
  for (i in 1:length(strings1)){
    if (!JuliaConnectoR::juliaEval(paste0(strings1[i], ' in keys(Pkg.project().dependencies)'))){
      #JuliaConnectoR.utils::install_julia_packages(strings2[i])
      JuliaConnectoR::juliaEval("using Pkg")
      JuliaConnectoR::juliaEval(paste0(strings2[i]) )
    }
  }
  JuliaConnectoR::juliaEval('
      using Coconots
      using Random
    ')
}