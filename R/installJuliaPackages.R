#' @title installJuliaPackages
#' @description checks for needed Julia packages and installs them if not installed.
#' @return no return value, called to install Julia packages in Julia.
#' @export
installJuliaPackages <- function(){
  
  strings1 <- c('"Random"', '"Distributions"', 
    '"ForwardDiff"', '"Optim"', '"StatsBase"', '"LineSearches"', '"LinearAlgebra"')
  strings2 <- c("Random", "Distributions", 
                "ForwardDiff", "Optim", "StatsBase", "LineSearches", "LinearAlgebra")
  for (i in 1:length(strings1)){
    if (!JuliaConnectoR::juliaEval(paste0(strings1[i], ' in keys(Pkg.project().dependencies)'))){
      #JuliaConnectoR.utils::install_julia_packages(strings2[i])
      JuliaConnectoR::juliaEval("using Pkg")
      JuliaConnectoR::juliaEval(paste0('Pkg.add("', strings2[i], '")') )
    }
  }
}