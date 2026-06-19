# Skip a test unless a working Julia installation is available.
#
# The Julia backend (via JuliaConnectoR and the Coconots.jl package) is an
# optional system dependency. Guarding the Julia-dependent tests with this
# helper lets `R CMD check` pass on machines without Julia -- including CRAN
# and the main R-CMD-check matrix -- while the dedicated `julia-tests`
# workflow installs Julia and exercises these tests in full.
skip_if_no_julia <- function() {
  testthat::skip_if_not_installed("JuliaConnectoR")
  ok <- tryCatch(
    isTRUE(JuliaConnectoR::juliaSetupOk()),
    error = function(e) FALSE
  )
  if (!ok) {
    testthat::skip("Julia is not available")
  }
  invisible(TRUE)
}
