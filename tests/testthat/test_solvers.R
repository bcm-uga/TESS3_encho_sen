context("Solver")



test_that("test cpp implementation of MCPA, comparison with R code", {

  data("data.for.test", package = "tess3r")

  # compute laplacian
  W <- ComputeHeatKernelWeight(data.for.test$coord, 2.0)
  Lapl <- ComputeGraphLaplacian(W)

  # change type of matrix
  X = matrix(as.integer(data.for.test$X),nrow(data.for.test$X),ncol(data.for.test$X))

  # With R code
  Rres <- SolveTess3Projected(X,
                              data.for.test$K,
                              data.for.test$d,
                              Lapl,
                              lambda = 1.0,
                              max.iteration = 25)

  # With cpp code
  XBin <- ComputeXBin(X,data.for.test$d)
  # cast as double
  XBin <- matrix(as.double(XBin), nrow(XBin), ncol(XBin))
  Lapl <- as.matrix(Lapl)
  cppres <- ComputeMCPASolution(XBin,
                                data.for.test$K,
                                Lapl,
                                lambdaPrim = 1.0,
                                data.for.test$d + 1,
                                maxIteration = 25, tolerance = 1e-10)

  expect_less_than(ComputeRmseWithBestPermutation(cppres$Q, Rres$Q), 1e-2)
  expect_less_than(ComputeRmseWithBestPermutation(cppres$G, Rres$G), 1e-2)
})
