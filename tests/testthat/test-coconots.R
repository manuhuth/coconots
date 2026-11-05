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

#-----------------PIT time range (t = order + 1, ..., T)-----------------------
# Reference PIT of an order-1 model computed directly from the transition
# densities over t = 2, ..., T (T - 1 time points).
pit_order1_reference <- function(data, lambda, alpha, eta, J) {
  u <- seq(0, 1, length = J + 1)
  rows <- t(vapply(2:length(data), function(t) {
    dvec <- coconots:::dGP1h_vec(data[t], data[t - 1], lambda[t], alpha, eta)
    Px <- sum(dvec[seq_len(data[t] + 1)])
    Px1 <- if (data[t] >= 1) sum(dvec[seq_len(data[t])]) else 0
    pmin(pmax((u - Px1) / (Px - Px1), 0), 1)
  }, numeric(J + 1)))
  diff(colMeans(rows))
}

test_that("cocoPit (R) of order-1 models uses t = 2, ..., T", {
  set.seed(2611)
  data <- cocoSim(order = 1, type = "GP", par = c(1, 0.3, 0.2), length = 80)
  fit <- cocoReg(order = 1, type = "GP", data = data)
  pit <- cocoPit(fit, J = 10)
  expect_equal(pit$`PIT values`,
               pit_order1_reference(data, rep(fit$par[1], length(data)),
                                    fit$par[2], fit$par[3], 10),
               tolerance = 1e-10)
  expect_equal(sum(pit$`PIT values`), 1, tolerance = 1e-12)
  # chi-square statistic on the bin counts of the T - 1 PIT time points
  n_pit <- length(data) - 1
  d <- pit$`PIT values`
  expect_equal(pit$p_value,
               stats::pchisq(sum((d * n_pit - n_pit / 10)^2 / (n_pit / 10)), 9,
                             lower.tail = FALSE))

  # covariate path (per-observation loop)
  n <- 150
  cov <- cbind(sin = sin(2 * pi / 50 * (1:n)), cos = cos(2 * pi / 50 * (1:n)))
  set.seed(2612)
  d1 <- suppressWarnings(
    cocoSim(order = 1, type = "Poisson", par = c(0.4, 0.3, -0.2),
            xreg = cov, length = n))
  fit1 <- cocoReg(order = 1, type = "Poisson", data = d1, xreg = cov)
  lambda1 <- coconots:::applyLinkFunctionVec(fit1$cov %*% fit1$par[-1],
                                             fit1$link_function)
  expect_equal(cocoPit(fit1, J = 10)$`PIT values`,
               pit_order1_reference(d1, lambda1, fit1$par[1], 0, 10),
               tolerance = 1e-10)

  # order 2: the PIT values are the mean over t = 3, ..., T and sum to one
  set.seed(2613)
  data2 <- cocoSim(order = 2, type = "Poisson", par = c(0.5, 0.2, 0.05, 0.3),
                   length = 100)
  fit2 <- cocoReg(order = 2, type = "Poisson", data = data2)
  expect_equal(sum(cocoPit(fit2)$`PIT values`), 1, tolerance = 1e-12)
})

test_that("cocoPit gives the same PIT and p-value with R and Julia", {
  skip_on_cran()
  skip_if_no_julia()
  models <- list(
    list(type = "Poisson", order = 1, par = c(0.8, 0.4)),
    list(type = "GP", order = 1, par = c(0.8, 0.4, 0.2)),
    list(type = "Poisson", order = 2, par = c(0.5, 0.2, 0.05, 0.3)),
    list(type = "GP", order = 2, par = c(0.5, 0.2, 0.05, 0.3, 0.2))
  )
  for (i in seq_along(models)) {
    m <- models[[i]]
    set.seed(3100 + i)
    data <- cocoSim(order = m$order, type = m$type, par = m$par, length = 120)
    fit <- cocoReg(order = m$order, type = m$type, data = data,
                   julia_installed = TRUE)
    pit_r <- cocoPit(fit, J = 10)
    pit_j <- cocoPit(fit, J = 10, julia = TRUE)
    expect_equal(as.numeric(pit_r$`PIT values`), as.numeric(pit_j$`PIT values`),
                 tolerance = 1e-8)
    expect_equal(pit_r$p_value, pit_j$p_value, tolerance = 1e-6)
    expect_equal(pit_r$confidence_bands, pit_j$confidence_bands)
  }
})

test_that("cocoPit confidence bands use conf.alpha, not the model alpha", {
  set.seed(3201)
  data <- cocoSim(order = 1, type = "Poisson", par = c(0.8, 0.4), length = 120)
  fit <- cocoReg(order = 1, type = "Poisson", data = data)
  n_pit <- length(data) - 1
  for (conf in c(0.05, 0.1)) {
    pit <- cocoPit(fit, J = 10, conf.alpha = conf)
    half <- sqrt(stats::qchisq(1 - conf, 9) / 100 * n_pit) / n_pit
    expect_equal(pit$confidence_bands, 0.1 + c(-1, 1) * half)
    expect_equal(pit$alpha, conf)
  }
})

#-----------------Default starting values (no covariates)------------------------
# Default start of cocoReg_base.R at the cocoReg defaults (seasonality c(1, 2),
# replace.start.val 1e-5, iteration.start.val 0.6; GP2 shrinks by 0.99), in the
# R ordering (lambda first) and the Julia ordering (lambda last, order 2 on the
# reparameterized alpha scale).
default_start_reference <- function(type, order, data) {
  rv <- 1e-5
  acf <- forecast::Acf(data, plot = FALSE, lag.max = 3)$acf
  m <- mean(data)
  eta <- 0
  if (type == "GP") {
    eta <- 1 - (m / stats::var(data))^0.5
    if (eta < 0) eta <- rv
    if (eta > 1) eta <- 1 - rv
  }
  eta_par <- if (type == "GP") eta else NULL
  if (order == 1) {
    alpha <- acf[2]
    if (alpha < 0) alpha <- rv
    lambda <- m * (1 - eta) * (1 - alpha)
    if (lambda < 0) lambda <- rv
    return(list(r = c(lambda, alpha, eta_par), julia = c(alpha, eta_par, lambda)))
  }
  idx <- 3:length(data)
  mu111 <- sum((data[idx] - m) * (data[idx - 1] - m) * (data[idx - 2] - m)) /
    length(data)
  a3 <- mu111 * (1 - eta)^4 / (m * (1 + 2 * eta))
  if (a3 < 0) a3 <- rv
  if (a3 > 1) a3 <- 1 - rv
  a2 <- acf[3] - a3
  if (a2 < 0) a2 <- rv
  if (a2 > 1) a2 <- 1 - rv
  a1 <- acf[2] - a3
  if (a1 < 0) a1 <- rv
  shrink <- if (type == "GP") 0.99 else 0.6
  while (a1 + a2 + a3 > 1) {
    a1 <- a1 * shrink
    a2 <- a2 * shrink
    a3 <- a3 * shrink
  }
  while (2 * a1 + a3 > 1) {
    a1 <- a1 * shrink
    a3 <- a3 * shrink
  }
  lambda <- m * (1 - eta) * (1 - a1 - a2 - a3)
  if (lambda < 0) lambda <- rv
  list(r = c(lambda, a1, a2, a3, eta_par),
       julia = c(2 * a1 + a3, a3 / (2 * a1 + a3), a2 / (1 - a1 - a3), eta_par,
                 lambda))
}

start_test_models <- list(
  list(type = "Poisson", order = 1, par = c(0.8, 0.4)),
  list(type = "GP", order = 1, par = c(0.8, 0.4, 0.2)),
  list(type = "Poisson", order = 2, par = c(0.5, 0.2, 0.05, 0.3)),
  list(type = "GP", order = 2, par = c(0.5, 0.2, 0.05, 0.3, 0.2)),
  # strong dependence (alpha1 + alpha2 + alpha3 = 0.9, 2 alpha1 + alpha3 = 0.85)
  list(type = "Poisson", order = 2, par = c(0.2, 0.3, 0.35, 0.25)),
  list(type = "GP", order = 2, par = c(0.2, 0.3, 0.35, 0.25, 0.2))
)

test_that("default_start_reference reproduces the default start of cocoReg (R)", {
  # the same start gives the same deterministic fit
  for (i in seq_along(start_test_models)) {
    m <- start_test_models[[i]]
    set.seed(4200 + i)
    data <- cocoSim(order = m$order, type = m$type, par = m$par, length = 150)
    fit_default <- suppressWarnings(cocoReg(order = m$order, type = m$type,
                                            data = data))
    fit_start <- cocoReg(order = m$order, type = m$type, data = data,
                         start = default_start_reference(m$type, m$order, data)$r)
    expect_equal(fit_start$par, fit_default$par, tolerance = 1e-8)
  }
})

test_that("Coconots.jl computes the same default start as cocoReg (R)", {
  skip_on_cran()
  skip_if_no_julia()
  coconots:::addJuliaFunctions()
  for (i in seq_along(start_test_models)) {
    m <- start_test_models[[i]]
    set.seed(4200 + i)
    data <- cocoSim(order = m$order, type = m$type, par = m$par, length = 150)
    start_julia <- as.numeric(JuliaConnectoR::juliaLet(
      'convert(Vector{Float64}, Coconots.get_starting_values!(type, Int(order), data, nothing, nothing, 0, -Inf))',
      type = m$type, order = m$order, data = data))
    expect_equal(start_julia, default_start_reference(m$type, m$order, data)$julia,
                 tolerance = 1e-10)
  }
})

test_that("PIT chi-square p-value is the upper tail", {
  # a perfectly uniform histogram has statistic 0 and p-value 1
  expect_equal(coconots:::pitChisqPValue(rep(0.1, 10), 99, 10), 1)
  # all mass in one bin: statistic (n - n/J)^2/(n/J) + (J - 1) n/J = n (J - 1)
  d <- c(1, rep(0, 9))
  expect_equal(coconots:::pitChisqPValue(d, 99, 10),
               stats::pchisq(99 * 9, 9, lower.tail = FALSE))
  expect_lt(coconots:::pitChisqPValue(d, 99, 10), 1e-10)
})

test_that("cocoResid uses t = order + 1, ..., T for both orders", {
  set.seed(2711)
  data <- cocoSim(order = 1, type = "GP", par = c(1, 0.3, 0.2), length = 60)
  fit <- cocoReg(order = 1, type = "GP", data = data)
  res <- cocoResid(fit)
  expect_length(res$residuals, length(data) - 1)
  lambda <- fit$par[[1]]; alpha <- fit$par[[2]]; eta <- fit$par[[3]]
  fitted_ref <- alpha * data[1:59] + lambda / (1 - eta)
  expect_equal(unname(res$fitted), fitted_ref, tolerance = 1e-12)
  expect_equal(unname(res$residuals), data[2:60] - fitted_ref, tolerance = 1e-12)

  set.seed(2712)
  data2 <- cocoSim(order = 2, type = "Poisson", par = c(0.5, 0.2, 0.05, 0.3),
                   length = 60)
  fit2 <- cocoReg(order = 2, type = "Poisson", data = data2)
  res2 <- cocoResid(fit2)
  expect_length(res2$residuals, length(data2) - 2)
  expect_equal(unname(res2$residuals), unname(data2[3:60] - res2$fitted))
})

test_that("order-2 conditional pmf and residuals are finite at alphas exactly on the boundary 0", {
  # The Julia backend returns boundary estimates exactly at 0 (e.g. alpha1 = 0
  # for a GP2 fit with covariates); a zero rate beta_i = lambda U alpha_i made
  # the C++ factor beta_i * pow(beta_i, -1) at index 0 NaN (0 * Inf), so
  # cocoResid raised the non-finite pmf error. The limit is the value just
  # inside the boundary.
  for (p in list(c(0, 0.2, 0.3), c(0.2, 0, 0.3), c(0.2, 0.3, 0), c(0.25, 0.1, 0.5))) {
    # last case: 2 alpha1 + alpha3 = 1, i.e. the rate zeta is zero
    for (eta in c(0, 0.15)) {
      for (yz in list(c(0L, 0L), c(3L, 2L), c(1L, 6L))) {
        p_in <- p
        p_in[p_in == 0] <- 1e-13
        if (2 * p[1] + p[3] == 1) p_in[1] <- p[1] - 1e-13
        g0 <- coconots:::dR2_vec(yz[1], yz[2], 1.5, p[1], p[2], p[3], eta, 1e-14)
        g1 <- coconots:::dR2_vec(yz[1], yz[2], 1.5, p_in[1], p_in[2], p_in[3], eta, 1e-14)
        expect_true(all(is.finite(g0)))
        expect_equal(sum(g0), 1, tolerance = 1e-10)
        n <- min(length(g0), length(g1))
        expect_equal(g0[1:n], g1[1:n], tolerance = 1e-8)
      }
    }
  }

  length <- 120
  cov <- cbind(sin(2 * pi / 25 * (1:length)), cos(2 * pi / 25 * (1:length)))
  set.seed(20261107)
  data <- cocoSim(order = 2, type = "GP", par = c(0.15, 0.1, 0.3, 0.1, 0.4, -0.2),
                  xreg = cov, length = length)
  fit <- cocoReg(order = 2, type = "GP", data = data, xreg = cov)
  fit0 <- fit
  fit0$par[1] <- 0
  fit1 <- fit
  fit1$par[1] <- 1e-13
  r0 <- cocoResid(fit0)
  r1 <- cocoResid(fit1)
  expect_true(all(is.finite(r0$residuals)))
  expect_equal(r0$residuals, r1$residuals, tolerance = 1e-8)
  expect_equal(r0$cond.var, r1$cond.var, tolerance = 1e-8)
})

test_that("order-2 likelihood, densities, residuals and PIT are finite at alpha1 = alpha3 = 0", {
  # With alpha1 = alpha3 = 0 the rate lambda U (alpha1 + alpha3) in the
  # bivariate normalizer P(y, z) is zero and the C++ j = 0 summand evaluated
  # (beta1 + beta3) * pow(0, -1) = 0 * Inf = NaN. The limit is the value just
  # inside the boundary. Gradients are compared with one-sided finite
  # differences at the boundary (the AD gradient at alpha = 1e-12 itself loses
  # digits to the cancellation in (beta1 + beta3) / (lambda U (alpha1 + alpha3))).
  fd_grad <- function(f, par, idx, h = 1e-6) {
    # second-order one-sided difference into par[idx] > 0, central otherwise
    sapply(seq_along(par), function(i) {
      e <- replace(numeric(length(par)), i, h)
      if (i %in% idx) {
        (-3 * f(par) + 4 * f(par + e) - f(par + 2 * e)) / (2 * h)
      } else {
        (f(par + e) - f(par - e)) / (2 * h)
      }
    })
  }

  set.seed(20261108)
  data <- cocoSim(order = 2, type = "GP", par = c(1.5, 0.05, 0.3, 0.05, 0.2), length = 150)
  idata <- as.integer(data)
  for (eta in c(0, 0.2)) {
    l0 <- coconots:::likelihoodGP2(1.5, 0, 0.3, 0, eta, length(data), 1L, 2L, idata, 1L)
    l1 <- coconots:::likelihoodGP2(1.5, 1e-12, 0.3, 1e-12, eta, length(data), 1L, 2L, idata, 1L)
    expect_true(is.finite(l0))
    expect_equal(l0, l1, tolerance = 1e-8)
    for (yz in list(c(0L, 0L), c(3L, 2L), c(1L, 6L), c(5L, 5L))) {
      p0 <- coconots:::Pyz(yz[1], yz[2], 1.5, 0, 0.3, 0, eta)
      p1 <- coconots:::Pyz(yz[1], yz[2], 1.5, 1e-12, 0.3, 1e-12, eta)
      expect_true(is.finite(p0) && p0 > 0)
      expect_equal(p0, p1, tolerance = 1e-8)
      d0 <- coconots:::dGP2h_vec(15L, yz[1], yz[2], 1.5, 0, 0.3, 0, eta)
      d1 <- coconots:::dGP2h_vec(15L, yz[1], yz[2], 1.5, 1e-12, 0.3, 1e-12, eta)
      expect_true(all(is.finite(d0)))
      expect_equal(d0, d1, tolerance = 1e-8)
    }
    par0 <- c(1.5, 0, 0.3, 0, eta)
    ad <- coconots:::likelihoodAD(par0, 2L, TRUE, idata, 1L, 2L, NULL, "log", 2L)
    expect_equal(ad$value, l0, tolerance = 1e-10)
    expect_true(all(is.finite(ad$gradient)))
    expect_true(all(is.finite(ad$hessian)))
    f <- function(p) coconots:::likelihoodGP2(p[1], p[2], p[3], p[4], p[5],
                                              length(data), 1L, 2L, idata, 1L)
    idx <- if (eta == 0) c(2, 4, 5) else c(2, 4)
    expect_equal(as.numeric(ad$gradient), fd_grad(f, par0, idx), tolerance = 1e-5)
  }

  # covariates (likelihoodGP2cov, gp2_pyz_cov)
  length <- 150
  cov <- cbind(sin(2 * pi / 25 * (1:length)), cos(2 * pi / 25 * (1:length)))
  set.seed(20261109)
  datac <- cocoSim(order = 2, type = "GP", par = c(0.05, 0.3, 0.05, 0.2, 0.4, -0.2),
                   xreg = cov, length = length)
  idatac <- as.integer(datac)
  for (eta in c(0, 0.2)) {
    lc0 <- coconots:::likelihoodGP2cov(0, 0.3, 0, eta, c(0.4, -0.2), length, 1L, 2L,
                                       idatac, cov, "log", 1L)
    lc1 <- coconots:::likelihoodGP2cov(1e-12, 0.3, 1e-12, eta, c(0.4, -0.2), length, 1L, 2L,
                                       idatac, cov, "log", 1L)
    expect_true(is.finite(lc0))
    expect_equal(lc0, lc1, tolerance = 1e-8)
    parc0 <- c(0, 0.3, 0, eta, 0.4, -0.2)
    adc <- coconots:::likelihoodAD(parc0, 2L, TRUE, idatac, 1L, 2L, cov, "log", 2L)
    expect_equal(adc$value, lc0, tolerance = 1e-10)
    expect_true(all(is.finite(adc$gradient)))
    expect_true(all(is.finite(adc$hessian)))
    fc <- function(p) coconots:::likelihoodGP2cov(p[1], p[2], p[3], p[4], p[5:6], length,
                                                  1L, 2L, idatac, cov, "log", 1L)
    idx <- if (eta == 0) c(1, 3, 4) else c(1, 3)
    expect_equal(as.numeric(adc$gradient), fd_grad(fc, parc0, idx), tolerance = 1e-5)
  }

  # post-processing on fits moved to alpha1 = alpha3 = 0
  fit <- cocoReg(order = 2, type = "GP", data = data)
  fitc <- cocoReg(order = 2, type = "GP", data = datac, xreg = cov)
  for (case in list(list(fit = fit, i = c(2, 4)), list(fit = fitc, i = c(1, 3)))) {
    fit0 <- case$fit
    fit0$par[case$i] <- 0
    fit1 <- case$fit
    fit1$par[case$i] <- 1e-12
    r0 <- cocoResid(fit0)
    r1 <- cocoResid(fit1)
    expect_true(all(is.finite(r0$residuals)))
    expect_equal(r0$residuals, r1$residuals, tolerance = 1e-8)
    expect_equal(r0$cond.var, r1$cond.var, tolerance = 1e-8)
    pit0 <- cocoPit(fit0, J = 10)
    pit1 <- cocoPit(fit1, J = 10)
    expect_true(all(is.finite(pit0$`PIT values`)))
    expect_equal(as.numeric(pit0$`PIT values`), as.numeric(pit1$`PIT values`),
                 tolerance = 1e-8)
    expect_equal(pit0$p_value, pit1$p_value, tolerance = 1e-8)
  }
})

test_that("R and Julia order-2 log-likelihoods agree at alpha1 = alpha3 = 0", {
  skip_on_cran()
  skip_if_no_julia()
  coconots:::addJuliaFunctions()
  set.seed(20261110)
  data <- cocoSim(order = 2, type = "GP", par = c(1.5, 0.05, 0.3, 0.05, 0.2), length = 150)
  length <- 150
  cov <- cbind(sin(2 * pi / 25 * (1:length)), cos(2 * pi / 25 * (1:length)))
  lambdas <- as.numeric(exp(cov %*% c(0.4, -0.2)))
  for (eta in c(0, 0.2)) {
    l_r <- coconots:::likelihoodGP2(1.5, 0, 0.3, 0, eta, length(data), 1L, 2L,
                                    as.integer(data), 1L)
    l_j <- as.numeric(JuliaConnectoR::juliaLet(
      'Coconots.compute_negative_log_likelihood_GP2(1.5, 0.0, 0.3, 0.0, eta, Int.(data))',
      eta = eta, data = as.numeric(data)))
    expect_true(is.finite(l_j))
    expect_equal(l_j, l_r, tolerance = 1e-8)
    lc_r <- coconots:::likelihoodGP2cov(0, 0.3, 0, eta, c(0.4, -0.2), length, 1L, 2L,
                                        as.integer(data), cov, "log", 1L)
    lc_j <- as.numeric(JuliaConnectoR::juliaLet(
      'Coconots.compute_negative_log_likelihood_GP2(lambdas, 0.0, 0.3, 0.0, eta, Int.(data))',
      lambdas = lambdas, eta = eta, data = as.numeric(data)))
    expect_true(is.finite(lc_j))
    expect_equal(lc_j, lc_r, tolerance = 1e-8)
  }
})

test_that("order-1 likelihood, densities, residuals and PIT are finite at alpha = 0 and alpha = 1", {
  # The quasi-binomial thinning term g(r | y) contains
  # alpha * pow(alpha + psi r, r - 1), which is 0 * Inf = NaN at alpha = 0 and
  # r = 0, and (1 - alpha) * pow(1 - alpha + psi (y - r), y - r - 1), which is
  # 0 * Inf = NaN at alpha = 1 (psi = 0) and r = y. The limit is the value just
  # inside the boundary. At alpha = 1, P(x | y) = 0 for x < y, so the
  # likelihood is compared on a nondecreasing series there. Gradients are
  # compared with second-order one-sided finite differences into the
  # parameter space.
  fd_grad <- function(f, par, dir, h = 1e-6) {
    # dir[i] = +1 / -1: one-sided difference upwards / downwards, 0: central
    sapply(seq_along(par), function(i) {
      e <- replace(numeric(length(par)), i, h)
      if (dir[i] != 0) {
        e <- dir[i] * e
        dir[i] * (-3 * f(par) + 4 * f(par + e) - f(par + 2 * e)) / (2 * h)
      } else {
        (f(par + e) - f(par - e)) / (2 * h)
      }
    })
  }

  set.seed(20261111)
  data <- cocoSim(order = 1, type = "GP", par = c(1.5, 0.3, 0.2), length = 150)
  set.seed(20261112)
  data_up <- cumsum(stats::rpois(150, 0.15))
  length <- 150
  cov <- cbind(sin(2 * pi / 25 * (1:length)), cos(2 * pi / 25 * (1:length)))
  for (case in list(list(a = 0, a_in = 1e-12, d = 1, x = data),
                    list(a = 1, a_in = 1 - 1e-12, d = -1, x = data_up))) {
    idata <- as.integer(case$x)
    for (eta in c(0, 0.2)) {
      l0 <- coconots:::likelihoodGP1(1.5, case$a, eta, length, 1L, idata)
      l1 <- coconots:::likelihoodGP1(1.5, case$a_in, eta, length, 1L, idata)
      expect_true(is.finite(l0))
      expect_equal(l0, l1, tolerance = 1e-8)
      for (y in c(0L, 1L, 4L, 9L)) {
        d0 <- coconots:::dGP1h_vec(15L, y, 1.5, case$a, eta)
        d1 <- coconots:::dGP1h_vec(15L, y, 1.5, case$a_in, eta)
        expect_true(all(is.finite(d0)))
        expect_equal(d0, d1, tolerance = 1e-8)
        h0 <- vapply(0:15, function(x) coconots:::dGP1h(x, y, 1.5, case$a, eta), 0)
        expect_equal(h0, d0, tolerance = 1e-12)
      }
      par0 <- c(1.5, case$a, eta)
      ad <- coconots:::likelihoodAD(par0, 1L, TRUE, idata, 1L, 1L, NULL, "log", 2L)
      expect_equal(ad$value, l0, tolerance = 1e-10)
      expect_true(all(is.finite(ad$gradient)))
      expect_true(all(is.finite(ad$hessian)))
      f <- function(p) coconots:::likelihoodGP1(p[1], p[2], p[3], length, 1L, idata)
      dir <- c(0, case$d, if (eta == 0) 1 else 0)
      expect_equal(as.numeric(ad$gradient), fd_grad(f, par0, dir), tolerance = 1e-5)

      # covariates (likelihoodGP1cov)
      lc0 <- coconots:::likelihoodGP1cov(case$a, eta, c(0.4, -0.2), length, 1L, idata,
                                         cov, "log")
      lc1 <- coconots:::likelihoodGP1cov(case$a_in, eta, c(0.4, -0.2), length, 1L, idata,
                                         cov, "log")
      expect_true(is.finite(lc0))
      expect_equal(lc0, lc1, tolerance = 1e-8)
      parc0 <- c(case$a, eta, 0.4, -0.2)
      adc <- coconots:::likelihoodAD(parc0, 1L, TRUE, idata, 1L, 1L, cov, "log", 2L)
      expect_equal(adc$value, lc0, tolerance = 1e-10)
      expect_true(all(is.finite(adc$gradient)))
      expect_true(all(is.finite(adc$hessian)))
      fc <- function(p) coconots:::likelihoodGP1cov(p[1], p[2], p[3:4], length, 1L, idata,
                                                    cov, "log")
      dirc <- c(case$d, if (eta == 0) 1 else 0, 0, 0)
      expect_equal(as.numeric(adc$gradient), fd_grad(fc, parc0, dirc), tolerance = 1e-5)
    }
  }

  # simulation at the boundary: alpha = 0 thins every count to 0, alpha = 1
  # keeps it (X_t = X_{t-1} + innovation)
  set.seed(20261113)
  u <- stats::runif(30)
  inn <- as.integer(stats::rpois(30, 1.5))
  init <- c(5L, integer(29))
  keep <- as.integer(5L + c(0L, cumsum(inn[-1])))
  expect_equal(coconots:::simGP1(1.5, 0, 0.2, 30L, 1L, 1L, init, u, inn)[-1], inn[-1])
  expect_equal(coconots:::simGP1(1.5, 1, 0.2, 30L, 1L, 1L, init, u, inn), keep)
  expect_equal(coconots:::simGP1cov(0, 0.2, c(0.4, -0.2), 30L, 1L, 1L, init, cov[1:29, ],
                                    u, inn[-1], "log")[-1], inn[-1])
  expect_equal(coconots:::simGP1cov(1, 0.2, c(0.4, -0.2), 30L, 1L, 1L, init, cov[1:29, ],
                                    u, inn[-1], "log"), keep)

  # post-processing on fits moved to alpha = 0 and alpha = 1 (for cocoPit at
  # alpha = 1 on a nondecreasing series, where every observed count has
  # positive conditional probability)
  set.seed(20261114)
  datac <- cocoSim(order = 1, type = "GP", par = c(0.3, 0.2, 0.4, -0.2), xreg = cov,
                   length = length)
  fit <- cocoReg(order = 1, type = "GP", data = data)
  fitc <- cocoReg(order = 1, type = "GP", data = datac, xreg = cov)
  for (case in list(list(fit = fit, i = 2), list(fit = fitc, i = 1))) {
    for (b in list(list(a = 0, a_in = 1e-12, ts = NULL),
                   list(a = 1, a_in = 1 - 1e-12, ts = data_up))) {
      fit0 <- case$fit
      fit0$par[case$i] <- b$a
      fit1 <- case$fit
      fit1$par[case$i] <- b$a_in
      r0 <- cocoResid(fit0)
      r1 <- cocoResid(fit1)
      expect_true(all(is.finite(r0$residuals)))
      expect_equal(r0$residuals, r1$residuals, tolerance = 1e-8)
      expect_equal(r0$cond.var, r1$cond.var, tolerance = 1e-8)
      if (!is.null(b$ts)) {
        fit0$ts <- b$ts
        fit1$ts <- b$ts
      }
      pit0 <- cocoPit(fit0, J = 10)
      pit1 <- cocoPit(fit1, J = 10)
      expect_true(all(is.finite(pit0$`PIT values`)))
      expect_equal(as.numeric(pit0$`PIT values`), as.numeric(pit1$`PIT values`),
                   tolerance = 1e-8)
      expect_equal(pit0$p_value, pit1$p_value, tolerance = 1e-8)
    }
  }
})

test_that("R and Julia order-1 log-likelihoods agree at alpha = 0 and alpha = 1", {
  skip_on_cran()
  skip_if_no_julia()
  coconots:::addJuliaFunctions()
  set.seed(20261115)
  data <- cocoSim(order = 1, type = "GP", par = c(1.5, 0.3, 0.2), length = 150)
  set.seed(20261116)
  data_up <- cumsum(stats::rpois(150, 0.15))
  length <- 150
  cov <- cbind(sin(2 * pi / 25 * (1:length)), cos(2 * pi / 25 * (1:length)))
  lambdas <- as.numeric(exp(cov %*% c(0.4, -0.2)))
  for (case in list(list(a = 0, x = data), list(a = 1, x = data_up))) {
    for (eta in c(0, 0.2)) {
      l_r <- coconots:::likelihoodGP1(1.5, case$a, eta, length, 1L, as.integer(case$x))
      l_j <- as.numeric(JuliaConnectoR::juliaLet(
        'Coconots.compute_negative_log_likelihood_GP1(1.5, alpha, eta, Int.(data))',
        alpha = case$a, eta = eta, data = as.numeric(case$x)))
      expect_true(is.finite(l_r))
      expect_true(is.finite(l_j))
      expect_equal(l_j, l_r, tolerance = 1e-8)
      lc_r <- coconots:::likelihoodGP1cov(case$a, eta, c(0.4, -0.2), length, 1L,
                                          as.integer(case$x), cov, "log")
      lc_j <- as.numeric(JuliaConnectoR::juliaLet(
        'Coconots.compute_negative_log_likelihood_GP1(lambdas, alpha, eta, Int.(data))',
        lambdas = lambdas, alpha = case$a, eta = eta, data = as.numeric(case$x)))
      expect_true(is.finite(lc_r))
      expect_true(is.finite(lc_j))
      expect_equal(lc_j, lc_r, tolerance = 1e-8)
    }
  }
})
