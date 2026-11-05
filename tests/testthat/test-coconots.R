#-----------------cgeck install julia packages----
test_that("insatll_Julia_packages", {
  skip_on_cran()
  skip_if_no_julia()
  expect_no_error(installJuliaPackages())
})
#-----------------No covariates, constrained--------------------------------
test_that("GP2Works", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3) #lambda, alpha1, alpha2, alpha3, eta
  set.seed(123499)
  data <- cocoSim(order = 2, type = "GP", par = par, length = length)
  cocoSim(order = 2, type = "GP", par = par, length = length, julia = TRUE,
          julia_seed = 1234992)
  
  fit <- cocoReg(order = 2, type = "GP", data = data, julia_installed = TRUE)
  for (i in c(TRUE, FALSE)){
    cocoPit(fit, julia=i)
    cocoResid(fit)
    summary(fit, julia=i)
    cocoScore(fit, julia=i)
    predict(fit, k=2, julia=i)
    cocoBoot(fit, rep.Bootstrap = 6, julia=i)
  }
  expect_no_error(fit)
})

test_that("Poisson2Works", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(12347)
  data <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  cocoSim(order = 2, type = "Poisson", par = par, length = length, 
                      julia = TRUE,
                      julia_seed = 12342399)
  
  fit <- cocoReg(order = 2, type = "Poisson", data = data, julia_installed = TRUE)
  for (i in c(TRUE, FALSE)){
    cocoPit(fit, julia=i)
    cocoResid(fit)
    summary(fit, julia=i)
    cocoScore(fit, julia=i)
    predict(fit, k=3, julia=i)
    cocoBoot(fit, rep.Bootstrap = 6, julia=i)
  }
  expect_no_error(fit)
})

test_that("GP1Works", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5, 0.2, 0.2)
  set.seed(12341)
  data <- cocoSim(order = 1, type = "GP", par = par, length = length)
  cocoSim(order = 1, type = "GP", par = par, length = length, 
          julia = TRUE,
          julia_seed = 123423299)
  
  fit <- cocoReg(order = 1, type = "GP", data = data, julia_installed = TRUE)
  for (i in c(TRUE, FALSE)){
    cocoPit(fit, julia=i)
    cocoResid(fit)
    summary(fit, julia=i)
    cocoScore(fit, julia=i)
    predict(fit, k=4, julia=i)
    cocoBoot(fit, rep.Bootstrap = 6, julia=i)
  }
  expect_no_error(fit)
})

test_that("Poisson1Works", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5, 0.2)
  set.seed(12345)
  data <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  
  cocoSim(order = 1, type = "Poisson", par = par, length = length, 
          julia = TRUE,
          julia_seed = 1234232199)
  fit <- cocoReg(order = 1, type = "Poisson", data = data, julia_installed = TRUE)
  for (i in c(TRUE, FALSE)){
    cocoPit(fit, julia=i)
    cocoResid(fit)
    summary(fit, julia=i)
    cocoScore(fit, julia=i)
    predict(fit, k=4, julia=i)
    cocoBoot(fit, rep.Bootstrap = 6, julia=i)
  }
  expect_no_error(fit)
})

#-----------------Covariates, constrained--------------------------------
test_that("Poisson1Works_cov", {
  skip_on_cran()
  skip_if_no_julia()
 ##Poisson1 model with covariates
 length <- 300
 period <- 50
 sin <- sin(2*pi/period*(1:length))
 cos <- cos(2*pi/period*(1:length))
 cov <- cbind(sin, cos)
 par <- c(0.8, 0.2, -0.2)
 set.seed(1234)
 data <- cocoSim(order = 1, type = "Poisson", par = par,
                     xreg = cov, length = length)
 
 cocoSim(order = 1, type = "Poisson", par = par, xreg = cov, length = length, 
         julia = TRUE,
         julia_seed = 123499)
 fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov, 
                julia_installed = TRUE)
 for (i in c(TRUE, FALSE)){
   cocoPit(fit, julia=i)
   cocoResid(fit)
   summary(fit, julia=i)
   cocoScore(fit, julia=i)
   predict(fit, k=2, xcast = cov[1:2,], julia=i)
   cocoBoot(fit, rep.Bootstrap = 6, julia=i)
 }
 expect_no_error(fit)

})

test_that("GP1Works_cov", {
  skip_on_cran()
  skip_if_no_julia()
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.3, 0.2, 0.2, -0.2)
  set.seed(1234567)
  data <- cocoSim(order = 1, type = "GP", par = par,
                      xreg = cov, length = length)
  
  cocoSim(order = 1, type = "GP", par = par, xreg = cov, length = length, 
          julia = TRUE,
          julia_seed = 123499)
  fit <- cocoReg(order = 1, type = "GP", data = data, xreg = cov, julia_installed = TRUE)
  for (i in c(TRUE, FALSE)){
    cocoPit(fit, julia=i)
    cocoResid(fit)
    summary(fit, julia=i)
    cocoScore(fit, julia=i)
    predict(fit, k=2, xcast = cov[1:2,], julia=i)
    cocoBoot(fit, rep.Bootstrap = 6, julia=i)
  }
  expect_no_error(fit)
})

test_that("Poisson2Works_cov", {
  skip_on_cran()
  skip_if_no_julia()
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, -0.2)
  set.seed(12345678)
  data <- cocoSim(order = 2, type = "Poisson", par = par,
                      xreg = cov, length = length)
  
  cocoSim(order = 2, type = "Poisson", par = par, xreg = cov, length = length, 
          julia = TRUE,
          julia_seed = 123499)
  fit <- cocoReg(order = 2, type = "Poisson", data = data, xreg = cov, julia_installed = TRUE)
  for (i in c(TRUE, FALSE)){
    cocoPit(fit, julia=i)
    cocoResid(fit)
    summary(fit, julia=i)
    cocoScore(fit, julia=i)
    predict(fit, k=2, xcast = cov[1:2,], julia=i)
    cocoBoot(fit, rep.Bootstrap = 6, julia=i)
  }
  expect_no_error(fit)
})

test_that("GP2Works_cov", {
  skip_on_cran()
  skip_if_no_julia()
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, 0.2, -0.2)
  set.seed(123456798)
  data <- cocoSim(order = 2, type = "GP", par = par,
                      xreg = cov, length = length)
  
  cocoSim(order = 2, type = "GP", par = par, xreg = cov, length = length, 
          julia = TRUE,
          julia_seed = 123499)
  fit <- cocoReg(order = 2, type = "GP", data = data, xreg = cov, julia_installed = TRUE)
  for (i in c(TRUE, FALSE)){
    cocoPit(fit, julia=i)
    cocoResid(fit)
    summary(fit, julia=i)
    cocoScore(fit, julia=i)
    predict(fit, k=2, xcast = cov[1:2,], julia=i)
    cocoBoot(fit, rep.Bootstrap = 6, julia=i)
  }
  expect_no_error(fit)
})

#-----------------No covariates, unconstrained--------------------------------
test_that("GP2Works_constr", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(123499)
  data <- cocoSim(order = 2, type = "GP", par = par, length = length)
  
  fit <- cocoReg(order = 2, type = "GP", data = data, constrained.optim = FALSE)
  expect_no_error(fit)
})

test_that("Poisson2Works_constr", {
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(12347)
  data <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  
  fit <- cocoReg(order = 2, type = "Poisson", data = data, constrained.optim = FALSE)
  expect_no_error(fit)
})

test_that("GP1Works_constr", {
  length <- 300
  par <- c(0.5, 0.2, 0.2)
  set.seed(12341)
  data <- cocoSim(order = 1, type = "GP", par = par, length = length)
  
  fit <- cocoReg(order = 1, type = "GP", data = data, constrained.optim = FALSE)
  expect_no_error(fit)
})

test_that("Poisson1Works_constr", {
  length <- 300
  par <- c(0.5, 0.2)
  set.seed(12345)
  data <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  
  fit <- cocoReg(order = 1, type = "Poisson", data = data, constrained.optim = FALSE)
  expect_no_error(fit)
})

#-----------------Covariates, unconstrained--------------------------------
test_that("Poisson1Works_cov_constr", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.8, 0.2, -0.2)
  set.seed(1234)
  data <- cocoSim(order = 1, type = "Poisson", par = par,
                      xreg = cov, length = length)
  
  fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov, constrained.optim = FALSE)
  expect_no_error(fit)
})

test_that("GP1Works_cov_constr", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.8, 0.2, 0.2, -0.2)
  set.seed(1234567)
  data <- cocoSim(order = 1, type = "GP", par = par,
                      xreg = cov, length = length)
  
  fit <- cocoReg(order = 1, type = "GP", data = data, xreg = cov, constrained.optim = FALSE)
  expect_no_error(fit)
})

test_that("GP2Works_cov_uncosntrained", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, 0.2, -0.2)
  set.seed(123456798)
  data <- cocoSim(order = 2, type = "GP", par = par,
                      xreg = cov, length = length)
  
  fit <- cocoReg(order = 2, type = "GP", data = data, xreg = cov,
                 constrained.optim = FALSE)
  expect_no_error(fit)
})

test_that("Poisson2Works_cov_constr", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, -0.2)
  set.seed(12345678)
  data <- cocoSim(order = 2, type = "Poisson", par = par,
                      xreg = cov, length = length)
  
  fit <- cocoReg(order = 2, type = "Poisson", data = data, xreg = cov, constrained.optim = FALSE)
  expect_no_error(fit)
})



#-----------------Simulate initial time series-------------
test_that("initialize_time_series", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(123499)
  data <- cocoSim(order = 2, type = "GP", par = par, length = length)
  
  data2 <- cocoSim(order = 2, type = "GP", par = par, length = length,
                      init = data)
  expect_no_error(data)
})
#----------------Custom starting values covariates------------------------------
test_that("custom_start_values_GP2", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(122133499)
  data <- cocoSim(order = 2, type = "GP", par = par, length = length)
  
  fit <- cocoReg(order = 2, type = "GP", data = data, start = par)
  expect_no_error(fit)
})

test_that("custom_start_values_Poisson2", {
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(13423499)
  data <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  
  fit <- cocoReg(order = 2, type = "Poisson", data = data, start = par)
  expect_no_error(fit)
})

test_that("custom_start_values_GP1", {
  length <- 300
  par <- c(0.5,0.2,0.3)
  set.seed(12323499)
  data <- cocoSim(order = 1, type = "GP", par = par, length = length)
  
  fit <- cocoReg(order = 1, type = "GP", data = data, start = par)
  expect_no_error(fit)
})

test_that("custom_start_values_Poisson1", {
  length <- 300
  par <- c(0.5,0.2)
  set.seed(12349923)
  data <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  
  fit <- cocoReg(order = 1, type = "Poisson", data = data, start = par)
  expect_no_error(fit)
})
#----------------Custom starting values covariates------------------------------


test_that("GP1Works_cov_2", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.8, 0.2, 0.2, -0.2)
  set.seed(1234567)
  data <- cocoSim(order = 1, type = "GP", par = par,
                      xreg = cov, length = length)
  
  fit <- cocoReg(order = 1, type = "GP", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  summary(fit)
  cocoScore(fit)
  predict(fit, k=2, xcast = cov[1:2,])
  cocoBoot(fit, rep.Bootstrap = 6)
  expect_no_error(fit)
})

test_that("Poisson2Works_cov_2", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, -0.2)
  set.seed(12345678)
  data <- cocoSim(order = 2, type = "Poisson", par = par,
                      xreg = cov, length = length)
  
  fit <- cocoReg(order = 2, type = "Poisson", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  summary(fit)
  cocoScore(fit)
  predict(fit, k=2, xcast = cov[1:2,])
  cocoBoot(fit, rep.Bootstrap = 6)
  expect_no_error(fit)
})

test_that("GP2Works_cov_2", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, 0.2, -0.2)
  set.seed(123456798)
  data <- cocoSim(order = 2, type = "GP", par = par,
                      xreg = cov, length = length)
  
  fit <- cocoReg(order = 2, type = "GP", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  summary(fit)
  cocoScore(fit)
  predict(fit, k=2, xcast = cov[1:2,])
  cocoBoot(fit, rep.Bootstrap = 6)
  expect_no_error(fit)
})
#----------------Wrong models no covariates-----------------------------------
test_that("wrong_modelsGP2", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(123499)
  data <- cocoSim(order = 2, type = "GP", par = par, length = length)
  
  for (julia in c(TRUE, FALSE)){
    fit <- cocoReg(order = 1, type = "Poisson", data = data, julia = julia)
    fit <- cocoReg(order = 2, type = "Poisson", data = data, julia = julia)
    fit <- cocoReg(order = 1, type = "GP", data = data, julia = julia)
  }
  expect_no_error(fit)
})

test_that("wrong_modelsPoisson2", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(12347)
  data <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  
  for (julia in c(TRUE, FALSE)){
    fit <- cocoReg(order = 1, type = "Poisson", data = data,  julia = julia)
    fit <- cocoReg(order = 1, type = "GP", data = data,  julia = julia)
    fit <- cocoReg(order = 2, type = "GP", data = data,  julia = julia)
  }
  expect_no_error(fit)
  
})

test_that("wrong_modelsGP1", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5, 0.2, 0.2)
  set.seed(12341)
  data <- cocoSim(order = 1, type = "GP", par = par, length = length)
  
  for (julia in c(TRUE, FALSE)){
    fit <- cocoReg(order = 1, type = "Poisson", data = data,  julia = julia)
    fit <- cocoReg(order = 2, type = "Poisson", data = data,  julia = julia)
    fit <- cocoReg(order = 2, type = "GP", data = data,  julia = julia)
  }
  expect_no_error(fit)
})

test_that("wrong_modelsPoisson1", {
  skip_on_cran()
  skip_if_no_julia()
  length <- 300
  par <- c(0.5, 0.2)
  set.seed(12345)
  data <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  
  for (julia in c(TRUE, FALSE)){
    fit <- cocoReg(order = 2, type = "Poisson", data = data,  julia = julia)
    fit <- cocoReg(order = 1, type = "GP", data = data,  julia = julia)
    fit <- cocoReg(order = 2, type = "GP", data = data,  julia = julia)
  }
  expect_no_error(fit)

})
test_that("julia thread default respects an existing JULIA_NUM_THREADS", {
  old <- Sys.getenv("JULIA_NUM_THREADS", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("JULIA_NUM_THREADS") else
    Sys.setenv(JULIA_NUM_THREADS = old))

  Sys.setenv(JULIA_NUM_THREADS = "3")
  coconots:::set_julia_threads()
  expect_identical(Sys.getenv("JULIA_NUM_THREADS"), "3")

  Sys.unsetenv("JULIA_NUM_THREADS")
  coconots:::set_julia_threads()
  set_value <- Sys.getenv("JULIA_NUM_THREADS")
  if (!identical(set_value, "")) {
    expect_gte(as.integer(set_value), 1L)
  }
})

test_that("Poisson covariate models store the log-likelihood of their own betas", {
  # Regression test for the beta1-dropping slice that produced a wrong stored
  # $likelihood and a heap over-read in the C++ kernel (CRAN valgrind, 2.0.3).
  # The Poisson objective must equal the GP objective evaluated at eta = 0.
  length <- 200
  period <- 50
  cov <- cbind(sin = sin(2 * pi / period * (1:length)),
               cos = cos(2 * pi / period * (1:length)))

  set.seed(4711)
  d1 <- suppressWarnings(
    cocoSim(order = 1, type = "Poisson", par = c(0.4, 0.3, -0.2),
            xreg = cov, length = length))
  fit1 <- cocoReg(order = 1, type = "Poisson", data = d1, xreg = cov)
  expect_equal(
    fit1$likelihood,
    -coconots:::likelihoodGP1cov(fit1$par[1], 0, fit1$par[-1], length(d1),
                                 1, d1, cov, fit1$link_function))

  set.seed(4712)
  d2 <- suppressWarnings(
    cocoSim(order = 2, type = "Poisson", par = c(0.2, 0.05, 0.3, 0.3, -0.2),
            xreg = cov, length = length))
  fit2 <- cocoReg(order = 2, type = "Poisson", data = d2, xreg = cov)
  expect_equal(
    fit2$likelihood,
    -coconots:::likelihoodGP2cov(fit2$par[1], fit2$par[2], fit2$par[3], 0,
                                 fit2$par[-(1:3)], length(d2), 1, 2, d2, cov,
                                 fit2$link_function))
})

#-----------------Input validation (counts, short series, xreg rank)--------
test_that("counts above 170 are rejected with an informative error", {
  data <- c(rep(3, 20), 171, rep(4, 20))
  expect_error(cocoReg(type = "Poisson", order = 1, data = data),
               "counts above 170 are not supported .*observed count = 171")
  expect_error(cocoReg(type = "GP", order = 2, data = data, julia = TRUE),
               "counts above 170 are not supported")
  expect_error(cocoSoc(data, print.progress = FALSE), "counts above 170")

  # C++: a factorial above 170 raises the same error instead of returning Inf
  expect_equal(coconots:::fact(170), factorial(170))
  expect_error(coconots:::fact(171), "counts above 170 are not supported")
  expect_error(coconots:::dGP1h_vec(171L, 2L, 1, 0.3, 0.1), "required count = 171")
  expect_error(coconots:::likelihoodGP1(1, 0.3, 0.1, length(data), 1L, as.integer(data)),
               "counts above 170 are not supported")
  expect_error(coconots:::likelihoodGP2(1, 0.2, 0.1, 0.1, 0.1, length(data), 1L, 2L,
                                        as.integer(data), 2L),
               "counts above 170 are not supported")

  set.seed(20261105)
  sim <- cocoSim(order = 1, type = "Poisson", par = c(1, 0.4), length = 100)
  fit <- cocoReg(order = 1, type = "Poisson", data = sim)
  expect_error(cocoScore(fit, max_x = 171), "'max_x' = 171 exceeds the largest supported count 170")
  expect_error(predict(fit, max = 200), "'max' = 200 exceeds the largest supported count 170")
  expect_error(cocoSoc(sim, print.progress = FALSE, max_x_score = 171), "'max_x_score' = 171")
  expect_error(cocoSim(order = 1, type = "Poisson", par = c(1, 0.4), length = 10, init = 171),
               "initial count = 171")
  bad <- fit
  bad$ts[10] <- 200
  expect_error(cocoResid(bad), "observed count = 200")
  expect_error(cocoPit(bad), "observed count = 200")
  expect_error(cocoScore(bad), "observed count = 200")
  expect_error(cocoBoot(bad, rep.Bootstrap = 2), "observed count = 200")
  expect_error(predict(bad), "observed count = 200")

  # the automatic forecast support ends at 170
  fc <- predict(fit)
  expect_lte(max(fc[[1]]$x), 170)
})

test_that("non-finite predictive pmf values raise an error instead of NaN/Inf scores", {
  # GP innovation pmf at x = 170 with eta = 0.95: pow(lambda + 170 eta, 169)
  # overflows in double precision
  coco <- list(ts = c(0, 170, 0, 3), par = c(1, 0.5, 0.95), type = "GP", order = 1,
               cov = NULL, likelihood = -100, link_function = "log")
  class(coco) <- "coco"
  expect_error(cocoScore(coco), "predictive pmf is not finite")
})

test_that("dR2_vec terminates on the finite support of R given (y, z)", {
  # 1 - 1e-17 rounds to 1, which the cumulative mass need not reach exactly;
  # the former loop never terminated then
  for (val in c(1e-14, 1e-17)) {
    p <- coconots:::dR2_vec(5L, 3L, 1, 0.2, 0.1, 0.1, 0.2, val)
    expect_lte(length(p), 5 + 3 + 1)
    expect_equal(sum(p), 1, tolerance = 1e-10)
    expect_true(all(is.finite(p)))
  }
  set.seed(20261106)
  sim <- cocoSim(order = 2, type = "Poisson", par = c(1, 0.2, 0.1, 0.1), length = 100)
  fit <- cocoReg(order = 2, type = "Poisson", data = sim)
  res <- cocoResid(fit, val.num = 1e-17)
  expect_true(all(is.finite(res$fitted)))
  expect_length(res$fitted, 98)
})

test_that("a rank-deficient xreg is rejected naming the aliased columns", {
  n <- 100
  x1 <- sin(2 * pi * (1:n) / 12)
  x2 <- cos(2 * pi * (1:n) / 12)
  set.seed(20261107)
  data <- rpois(n, 2)
  xreg <- cbind(x1 = x1, x2 = x2, x3 = x1 + x2)
  expect_error(cocoReg(type = "Poisson", order = 1, data = data, xreg = xreg),
               "xreg is rank deficient: column\\(s\\) x3 are linear combinations of others")
  expect_error(cocoReg(type = "GP", order = 2, data = data, xreg = cbind(a = x1, b = x1),
                       julia = TRUE),
               "xreg is rank deficient: column\\(s\\) b")
  expect_error(cocoReg(type = "Poisson", order = 1, data = data, xreg = cbind(x1, 0 * x1)),
               "xreg is rank deficient")
  # full rank passes the check
  expect_null(coconots:::coco_check_xreg_rank(cbind(x1, x2)))
})

test_that("too short series raise informative errors", {
  expect_error(cocoReg(type = "Poisson", order = 1, data = c(1, 2)),
               "cocoReg: the time series has 2 observation\\(s\\), but at least 3 are required")
  expect_error(cocoReg(type = "GP", order = 2, data = c(1, 2, 3)),
               "cocoReg: the time series has 3 observation\\(s\\), but at least 4 are required")
  expect_error(cocoReg(type = "GP", order = 2, data = c(1, 2, 3), julia = TRUE),
               "at least 4 are required")
  expect_error(cocoSoc(c(1, 2, 3), print.progress = FALSE),
               "cocoSoc: the time series has 3 observation\\(s\\), but at least 4 are required")

  set.seed(20261108)
  sim <- cocoSim(order = 1, type = "Poisson", par = c(1, 0.4), length = 50)
  fit <- cocoReg(order = 1, type = "Poisson", data = sim)
  expect_error(cocoBoot(fit, numb.lags = 50, rep.Bootstrap = 2),
               "numb.lags = 50 must be smaller than the length of the time series \\(50\\)")
  short <- fit
  short$ts <- short$ts[1]
  expect_error(cocoScore(short), "cocoScore: the time series has 1 observation")
  expect_error(cocoPit(short), "cocoPit: the time series has 1 observation")
  expect_error(cocoResid(short), "cocoResid: the time series has 1 observation")
  expect_error(cocoBoot(short, rep.Bootstrap = 2), "cocoBoot: the time series has 1 observation")
  short$ts <- numeric(0)
  expect_error(predict(short), "predict: the time series has 0 observation")
})
