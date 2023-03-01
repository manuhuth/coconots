test_that("GP2Works", {
  length <- 1000
  par <- c(0.5,0.2,0.05,0.3,0.3)
  data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data)
})

test_that("Poisson2Works", {
  length <- 1000
  par <- c(0.5, 0.2, 0.05, 0.3)
  data.sim <- cocoSim(order = 2, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "Poisson", data = data)
})

test_that("GP1Works", {
  length <- 1000
  par <- c(0.5, 0.2, 0.2)
  data.sim <- cocoSim(order = 1, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "GP", data = data)
})

test_that("Poisson1Works", {
  length <- 1000
  par <- c(0.5, 0.2)
  data.sim <- cocoSim(order = 1, type = "Poisson", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 1, type = "Poisson", data = data)
})