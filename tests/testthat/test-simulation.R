#-----------------Simulation (RCPP backend): regression tests---------------------

# Former R implementation of rgenpois(): inversion of the cumulative
# HMMpa::dgenpois from x = 0 upwards.
ref_rgenpois_inv <- function(u, lambda, eta) {
  x <- 0
  cdf <- HMMpa::dgenpois(0, lambda, eta)
  while (cdf < u) {
    x <- x + 1
    cdf <- cdf + HMMpa::dgenpois(x, lambda, eta)
  }
  x
}

# Covariates with intercept 0.1 and lambda_t = 200 at row t0 (log link).
spike_xreg <- function(n, t0) {
  spike <- numeric(n)
  spike[t0] <- 1
  cbind(1, spike)
}
spike_beta <- c(log(0.1), log(2000))

test_that("cocoSim with covariates drives observation t by covariate row t", {
  n <- 20
  cases <- list(
    list(type = "Poisson", order = 1, dyn = c(0.01)),
    list(type = "GP", order = 1, dyn = c(0.01, 0.1)),
    list(type = "Poisson", order = 2, dyn = c(0.01, 0.01, 0.01)),
    list(type = "GP", order = 2, dyn = c(0.01, 0.01, 0.01, 0.1))
  )
  for (cs in cases) {
    for (t0 in c(1, cs$order, 10, n)) {
      xreg <- spike_xreg(n, t0)
      set.seed(100 + t0)
      x <- suppressWarnings(cocoSim(type = cs$type, order = cs$order,
                                    par = c(cs$dyn, spike_beta), length = n,
                                    xreg = xreg))
      expect_length(x, n)
      expect_equal(which.max(x), t0, info = paste(cs$type, cs$order, t0))
      expect_gt(x[t0], 100)

      # with init the output observation t is driven by row t as well
      set.seed(200 + t0)
      xi <- cocoSim(type = cs$type, order = cs$order, par = c(cs$dyn, spike_beta),
                    length = n, xreg = xreg, init = rep(0, cs$order))
      expect_length(xi, n)
      expect_equal(which.max(xi), t0, info = paste("init", cs$type, cs$order, t0))
    }
  }
})

test_that("order-2 GP simulation with covariates draws the start values from GP(lambda_t, eta)", {
  n <- 10
  xreg <- cbind(1, seq(0, 1, length.out = n))
  beta <- c(log(1.5), 0.5)
  eta <- 0.6
  lambdas <- exp(as.vector(xreg %*% beta))
  set.seed(31)
  x <- suppressWarnings(cocoSim(type = "GP", order = 2,
                                par = c(0.2, 0.1, 0.2, eta, beta), length = n,
                                xreg = xreg))
  # each start value consumes one uniform (rgenpois(n = 1, ...)), in time order
  set.seed(31)
  u <- stats::runif(2)
  expect_equal(x[1:2], c(ref_rgenpois_inv(u[1], lambdas[1], eta),
                         ref_rgenpois_inv(u[2], lambdas[2], eta)))
})

test_that("the C++ GP inversion reproduces the former R inversion", {
  set.seed(5)
  u <- stats::runif(30)
  lam <- rep(c(0.7, 3, 3, 6), length.out = 30)
  for (eta in c(0, 0.3, 0.7)) {
    expect_equal(coconots:::rgenpoisInv(u, lam, eta),
                 mapply(ref_rgenpois_inv, u, lam, MoreArgs = list(eta = eta)))
  }
})

test_that("the GP inversion terminates when the cumulative probability saturates below u", {
  u_hi <- c(1, 1 - .Machine$double.eps / 2, 1)
  for (eta in c(0, 0.4)) {
    x <- coconots:::rgenpoisInv(u_hi, c(2, 2, 5), eta)
    expect_true(all(is.finite(x)))
    expect_true(all(x >= 0 & x < 1000))
    # the returned value lies in the far tail
    expect_gt(x[1], 10)
  }
  # the cumulative probability of the HMMpa recursion stays near 0 at a large
  # lambda (denormal pmf values): informative error instead of an endless loop
  expect_error(coconots:::rgenpoisInv(0.5, 1000, 0), "did not converge")
})

test_that("the thinning inversion terminates at the support maximum", {
  # u = 1: the accumulated thinning probability never exceeds u, the former
  # loops never ended; now the largest feasible value is drawn
  expect_equal(coconots:::simGP1(1, 0.5, 0.2, 3L, 1L, 1L, c(5L, 0L, 0L),
                                 c(1, 1, 1), c(0L, 0L, 0L)),
               c(5L, 5L, 5L))
  expect_equal(coconots:::simGP2(1, 0.3, 0.1, 0.2, 0.2, 4L, 2L, 1L, 2L,
                                 c(2L, 3L, 0L, 0L), rep(1, 4), rep(0L, 4)),
               c(2L, 3L, 5L, 8L))
  xreg <- matrix(1, 2, 1)
  expect_equal(coconots:::simGP1cov(0.5, 0.2, 0, 3L, 1L, 1L, c(5L, 0L, 0L), xreg,
                                    c(1, 1, 1), c(0L, 0L), "log"),
               c(5L, 5L, 5L))
  expect_equal(coconots:::simGP2cov(0.3, 0.1, 0.2, 0.2, 0, 4L, 2L, 1L, 2L,
                                    c(2L, 3L, 0L, 0L), xreg, c(1, 1), c(0L, 0L),
                                    "log"),
               c(2L, 3L, 5L, 8L))
})

test_that("cocoSim stops for non-positive lambda_t under the identity link", {
  n <- 30
  xreg <- cbind(1, c(rep(1, n - 1), -10))
  beta <- c(1, 0.5)  # lambda_n = 1 - 5 < 0
  for (type in c("Poisson", "GP")) {
    for (order in 1:2) {
      dyn <- if (order == 1) 0.2 else c(0.2, 0.1, 0.2)
      if (type == "GP") dyn <- c(dyn, 0.2)
      expect_error(
        suppressWarnings(cocoSim(type = type, order = order, par = c(dyn, beta),
                                 length = n, xreg = xreg,
                                 link_function = "identity")),
        "finite and positive")
    }
  }
})

test_that("bootAcfCov replicates are the cocoSim series", {
  n <- 40
  xreg <- cbind(1, sin(2 * pi * (1:n) / 10))
  cases <- list(
    list(type = "Poisson", order = 1, dyn = c(0.3), alpha = c(0.3, 0, 0), eta = 0),
    list(type = "GP", order = 2, dyn = c(0.2, 0.1, 0.2, 0.3), alpha = c(0.2, 0.1, 0.2),
         eta = 0.3)
  )
  beta <- c(0.5, 0.4)
  nlags <- 5
  for (cs in cases) {
    lambdas <- exp(as.vector(xreg %*% beta))
    set.seed(77)
    ac <- coconots:::bootAcfCov(as.integer(cs$order), cs$type == "GP", cs$alpha[1],
                                cs$alpha[2], cs$alpha[3], cs$eta, beta, xreg, lambdas,
                                lambdas[1:cs$order], as.integer(n), as.integer(n),
                                as.integer(nlags), 2L, "log")
    set.seed(77)
    for (b in 1:2) {
      x <- suppressWarnings(cocoSim(type = cs$type, order = cs$order,
                                    par = c(cs$dyn, beta), length = n, xreg = xreg))
      expect_equal(ac[, b], stats::acf(x, lag.max = nlags, plot = FALSE)$acf[-1],
                   info = paste(cs$type, cs$order, b))
    }
  }
})

test_that("k-step forecast paths with covariates use xcast row j at step j", {
  k <- 4
  xcast <- cbind(1, c(0, 0, 1, 0))  # lambda = 200 at step 3
  for (type in c("Poisson", "GP")) {
    for (order in 1:2) {
      dyn <- if (order == 1) 0.01 else c(0.01, 0.01, 0.01)
      if (type == "GP") dyn <- c(dyn, 0.1)
      fit <- list(type = type, order = order, par = c(dyn, spike_beta),
                  link_function = "log")
      init <- rep(1L, order)
      set.seed(9)
      paths <- coconots:::forecastPathsRCPP(fit, k, 2, xcast, init)
      expect_true(all(apply(paths, 1, which.max) == 3))
      # the fast path equals the per-path cocoSim_cov(init = ...) fallback
      set.seed(9)
      for (s in 1:2) {
        ref <- coconots:::cocoSim_cov(type, order = order, par = fit$par,
                                      size = k + order, xreg = xcast,
                                      init = init, link_function = "log")$data
        expect_equal(as.numeric(paths[s, ]), as.numeric(ref),
                     info = paste(type, order, s))
      }
    }
  }
})
