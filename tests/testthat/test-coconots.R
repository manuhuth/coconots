#-----------------No covariates, constrained--------------------------------
test_that("GP2Works", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(123499)
  data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit)
  cocoBoot(fit, rep.Bootstrap = 6)
  
})

test_that("Poisson2Works", {
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(12347)
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit)
  cocoBoot(fit, rep.Bootstrap = 6)
})

test_that("GP1Works", {
  length <- 300
  par <- c(0.5, 0.2, 0.2)
  set.seed(12341)
  data.sim <- cocoSim(order = 1, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "GP", data = data)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit)
  cocoBoot(fit, rep.Bootstrap = 6)
})

test_that("Poisson1Works", {
  length <- 300
  par <- c(0.5, 0.2)
  set.seed(12345)
  data.sim <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit)
  cocoBoot(fit, rep.Bootstrap = 6)
})

#-----------------Covariates, constrained--------------------------------
test_that("Poisson1Works_cov", {
 ##Poisson1 model with covariates
 length <- 300
 period <- 50
 sin <- sin(2*pi/period*(1:length))
 cos <- cos(2*pi/period*(1:length))
 cov <- cbind(sin, cos)
 par <- c(0.8, 0.2, -0.2)
 set.seed(1234)
 data.sim <- cocoSim(order = 1, type = "Poisson", par = par,
                     xreg = cov, length = length)
 data <- data.sim$data
 fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov)
 cocoPit(fit)
 cocoResid(fit)
 cocoSummary(fit)
 cocoScore(fit)
 cocoForecast(fit, xcast = cov[1,])
 cocoBoot(fit, rep.Bootstrap = 6)
})

test_that("GP1Works_cov", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.8, 0.2, 0.2, -0.2)
  set.seed(1234567)
  data.sim <- cocoSim(order = 1, type = "GP", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "GP", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit, xcast = cov[1,])
  cocoBoot(fit, rep.Bootstrap = 6)
})

test_that("Poisson2Works_cov", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, -0.2)
  set.seed(12345678)
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit, xcast = cov[1,])
  cocoBoot(fit, rep.Bootstrap = 6)
})

test_that("GP2Works_cov", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, 0.2, -0.2)
  set.seed(123456798)
  data.sim <- cocoSim(order = 2, type = "GP", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit, xcast = cov[1,])
  cocoBoot(fit, rep.Bootstrap = 6)
})

#-----------------No covariates, unconstrained--------------------------------
test_that("GP2Works_constr", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(123499)
  data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data, constrained.optim = FALSE)
  
})

test_that("Poisson2Works_constr", {
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(12347)
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data, constrained.optim = FALSE)
})

test_that("GP1Works_constr", {
  length <- 300
  par <- c(0.5, 0.2, 0.2)
  set.seed(12341)
  data.sim <- cocoSim(order = 1, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "GP", data = data, constrained.optim = FALSE)
})

test_that("Poisson1Works_constr", {
  length <- 300
  par <- c(0.5, 0.2)
  set.seed(12345)
  data.sim <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data, constrained.optim = FALSE)
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
  data.sim <- cocoSim(order = 1, type = "Poisson", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov, constrained.optim = FALSE)
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
  data.sim <- cocoSim(order = 1, type = "GP", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "GP", data = data, xreg = cov, constrained.optim = FALSE)
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
  data.sim <- cocoSim(order = 2, type = "GP", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data, xreg = cov,
                 constrained.optim = FALSE)
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
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data, xreg = cov, constrained.optim = FALSE)
})



#-----------------Simulate initial time series-------------
test_that("initialize_time_series", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(123499)
  data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
  data <- data.sim$data
  data.sim2 <- cocoSim(order = 2, type = "GP", par = par, length = length,
                      init = data)
})
#----------------Custom starting values covariates------------------------------
test_that("custom_start_values_GP2", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(122133499)
  data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data, start = par)
})

test_that("custom_start_values_Poisson2", {
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(13423499)
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data, start = par)
})

test_that("custom_start_values_GP1", {
  length <- 300
  par <- c(0.5,0.2,0.3)
  set.seed(12323499)
  data.sim <- cocoSim(order = 1, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "GP", data = data, start = par)
})

test_that("custom_start_values_Poisson1", {
  length <- 300
  par <- c(0.5,0.2)
  set.seed(12349923)
  data.sim <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data, start = par)
})
#----------------Custom starting values covariates------------------------------
test_that("Poisson1Works_cov", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.8, 0.2, -0.2)
  set.seed(1234)
  data.sim <- cocoSim(order = 1, type = "Poisson", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data, xreg = cov)
})

test_that("GP1Works_cov", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.8, 0.2, 0.2, -0.2)
  set.seed(1234567)
  data.sim <- cocoSim(order = 1, type = "GP", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "GP", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit, xcast = cov[1,])
  cocoBoot(fit, rep.Bootstrap = 6)
})

test_that("Poisson2Works_cov", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, -0.2)
  set.seed(12345678)
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit, xcast = cov[1,])
  cocoBoot(fit, rep.Bootstrap = 6)
})

test_that("GP2Works_cov", {
  ##Poisson1 model with covariates
  length <- 300
  period <- 50
  sin <- sin(2*pi/period*(1:length))
  cos <- cos(2*pi/period*(1:length))
  cov <- cbind(sin, cos)
  par <- c(0.2, 0.05, 0.3, 0.2, 0.2, -0.2)
  set.seed(123456798)
  data.sim <- cocoSim(order = 2, type = "GP", par = par,
                      xreg = cov, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data, xreg = cov)
  cocoPit(fit)
  cocoResid(fit)
  cocoSummary(fit)
  cocoScore(fit)
  cocoForecast(fit, xcast = cov[1,])
  cocoBoot(fit, rep.Bootstrap = 6)
})
#----------------Wrong models no covariates-----------------------------------
#-----------------No covariates, constrained--------------------------------
test_that("wrong_modelsGP2", {
  length <- 300
  par <- c(0.5,0.2,0.05,0.3,0.3)
  set.seed(123499)
  data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data)
  fit <- cocoReg(order = 2, type = "Poisson", data = data)
  fit <- cocoReg(order = 1, type = "GP", data = data)
})

test_that("wrong_modelsPoisson2", {
  length <- 300
  par <- c(0.5, 0.2, 0.05, 0.3)
  set.seed(12347)
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data)
  fit <- cocoReg(order = 1, type = "GP", data = data)
  fit <- cocoReg(order = 2, type = "GP", data = data)
  
})

test_that("wrong_modelsGP1", {
  length <- 300
  par <- c(0.5, 0.2, 0.2)
  set.seed(12341)
  data.sim <- cocoSim(order = 1, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data)
  fit <- cocoReg(order = 2, type = "Poisson", data = data)
  fit <- cocoReg(order = 2, type = "GP", data = data)
})

test_that("wrong_modelsPoisson1", {
  length <- 300
  par <- c(0.5, 0.2)
  set.seed(12345)
  data.sim <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data)
  fit <- cocoReg(order = 1, type = "GP", data = data)
  fit <- cocoReg(order = 2, type = "GP", data = data)

})