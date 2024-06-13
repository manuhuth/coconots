#-----------------cgeck install julia packages----
test_that("insatll_Julia_packages", {
  skip_on_cran()
  expect_no_error(installJuliaPackages())
})
#-----------------No covariates, constrained--------------------------------
test_that("GP2Works", {
  skip_on_cran()
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
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