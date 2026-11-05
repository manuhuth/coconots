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
