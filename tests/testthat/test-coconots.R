test_that("GP2Works", {
  length <- 1000
  par <- c(0.5,0.2,0.05,0.3,0.3)
  data.sim <- cocoSim(order = 2, type = "GP", par = par, length = length)
  data <- data.sim$data
  fit <- cocoReg(order = 2, type = "GP", data = data)
})
