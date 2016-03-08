context("Solver")



test_that("test cpp implementation of MCPA, comparison with R code", {

  data("data.for.test", package = "TESS3enchoSen")

  # compute laplacian
  W <- ComputeHeatKernelWeight(data.for.test$coord, NULL)
  Lapl <- ComputeGraphLaplacian(W)

  # change type of matrix
  X = matrix(as.integer(data.for.test$X),nrow(data.for.test$X),ncol(data.for.test$X))

  # With R code
  Rres <- SolveTess3Projected(X,
                              data.for.test$K,
                              data.for.test$d,
                             Lapl,
                             lambda = 1.0,
                             max.iteration = 20)

  # With cpp code
  XBin <- ComputeXBin(X,data.for.test$d)
  cppres <- ComputeMCPASolution(XBin,
                                data.for.test$K,
                                Lapl,
                                lambda = 1.0,
                                data.for.test$d +1,
                                maxIteration = 20)

  expect_less_than(ComputeRmseWithBestPermutation(cppres$Q, Rres$Q), 1e-2)
  expect_less_than(ComputeRmseWithBestPermutation(cppres$G, Rres$G), 1e-2)
})
