#' @title installJuliaPackages
#' @description checks for needed \proglang{julia} packages and installs them if not installed.
#' @importFrom JuliaConnectoR juliaEval
#' @return no return value, called to install \proglang{julia} packages in \proglang{julia}.
#' @export
installJuliaPackages <- function(){
  set_julia_threads()

  min_version <- "0.2.0" # minimum Coconots.jl version the R bridge supports
  JuliaConnectoR::juliaEval("using Pkg")
  installed <- JuliaConnectoR::juliaEval(
    'v = [string(d.version) for d in values(Pkg.dependencies()) if d.name == "Coconots" && d.is_direct_dep]; isempty(v) ? "" : first(v)')
  if (installed == "" || numeric_version(installed) < min_version) {
    if (installed != "") message("Updating Coconots.jl ", installed, " to >= ", min_version, ".")
    JuliaConnectoR::juliaEval(sprintf('Pkg.add(name = "Coconots", version = "%s")',
                                      sub("\\.0$", "", min_version)))
  }
  JuliaConnectoR::juliaEval('
      using Coconots
      using Random
    ')
}

