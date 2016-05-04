context("Solver")



test_that("test cpp implementation of MCPA, comparison with R code", {

  data("data.for.test", package = "tess3r")

  # compute laplacian
  W <- ComputeHeatKernelWeight(data.for.test$coord, 2.0)
  Lapl <- ComputeGraphLaplacian(W)

  # With R code
  set.seed(0)
  Rres <- SolveTess3Projected(data.for.test$X,
                              data.for.test$K,
                              data.for.test$d,
                              Lapl,
                              lambda = 1.0,
                              max.iteration = 25)

  # cpp code
  set.seed(0)
  # With cpp code
  XBin <- X2XBin(X,data.for.test$d)
  Lapl <- as.matrix(Lapl)
  cppres <- list()
  cppres$G <- matrix(0.0, nrow = (data.for.test$d + 1) * data.for.test$L, ncol = data.for.test$K)
  cppres$Q <- matrix(runif(data.for.test$n * data.for.test$K), data.for.test$n, data.for.test$K)
  cppres$Q <- ProjectQ(cppres$Q)
  ComputeMCPASolution(XBin,
                      data.for.test$K,
                      Lapl,
                      lambdaPrim = 1.0,
                      data.for.test$d + 1,
                      maxIteration = 25, tolerance = 1e-10, Q = cppres$Q, G = cppres$G)

  expect_less_than(ComputeRmseWithBestPermutation(cppres$Q, Rres$Q), 1e-12)
  expect_less_than(ComputeRmseWithBestPermutation(cppres$G, Rres$G), 1e-12)
})
