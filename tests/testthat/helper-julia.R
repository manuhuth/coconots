# Skip a test unless the Julia backend is genuinely usable.
#
# The Julia backend needs BOTH a working Julia install (via JuliaConnectoR)
# AND the Coconots.jl package. Checking only that Julia exists is not enough:
# some CI images (e.g. GitHub's ubuntu-latest) ship with Julia but without
# Coconots.jl, which would make the Julia tests fail rather than skip.
#
# Guarding on the real requirement lets `R CMD check` pass anywhere Julia or
# Coconots.jl is missing (CRAN, the main R-CMD-check matrix), while the
# dedicated julia-tests workflow installs Coconots.jl and runs these in full.
#
# Setting COCONOTS_SKIP_JULIA_TESTS=true forces a skip without touching Julia
# at all (used by the pure R-CMD-check matrix for a deterministic, fast run).
skip_if_no_julia <- function() {
  if (identical(tolower(Sys.getenv("COCONOTS_SKIP_JULIA_TESTS")), "true")) {
    testthat::skip("Julia integration tests disabled via COCONOTS_SKIP_JULIA_TESTS")
  }
  testthat::skip_if_not_installed("JuliaConnectoR")
  have_backend <- tryCatch(
    isTRUE(JuliaConnectoR::juliaSetupOk()) &&
      isTRUE(JuliaConnectoR::juliaEval('Base.find_package("Coconots") !== nothing')),
    error = function(e) FALSE
  )
  if (!have_backend) {
    testthat::skip("Julia with the Coconots.jl package is not available")
  }
  invisible(TRUE)
}
