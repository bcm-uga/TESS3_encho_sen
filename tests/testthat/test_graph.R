context("Graph")


test_that("test of graph matrix computation based on variogram", {

  data("data.for.test", package = "tess3r")
  W <- ComputeGraphBasedOnVariogram(coord = data.for.test$coord, X = data.for.test$X, plot = TRUE)
  expect_equal(dim(W), c(data.for.test$n, data.for.test$n))
  expect_equal(isSymmetric(W), TRUE)
  expect_gte(min(W), 0.0)
  expect_lte(max(W), 1.0)
  W <- ComputeGraphBasedOnVariogram(coord = data.for.test$coord, X = data.for.test$X, plot = TRUE, nugget = 1.5e5)
  expect_equal(dim(W), c(data.for.test$n, data.for.test$n))
  expect_equal(isSymmetric(W), TRUE)
  expect_gte(min(W), 0.0)
  expect_lte(max(W), 1.0)
  W <- ComputeGraphBasedOnVariogram(coord = data.for.test$coord, X = data.for.test$X, plot = TRUE, nugget = 1.5e5,
                                    lag = 0.5, lmax = 2)
  expect_equal(dim(W), c(data.for.test$n, data.for.test$n))
  expect_equal(isSymmetric(W), TRUE)
  expect_gte(min(W), 0.0)
  expect_lte(max(W), 1.0)
})


test_that("test of graph matrix computation", {

  data("data.for.test", package = "tess3r")
  W <- ComputeHeatKernelWeight(coord = data.for.test$coord,
                               sigma = 2.0)

  Lapl <- ComputeGraphLaplacian(W)

  coord = data.for.test$coord
})


test_that("test of spectral decompostion of Lapl", {
  skip("debug test")

  data("data.for.test", package = "tess3r")
  W <- ComputeHeatKernelWeight(coord = data.for.test$coord,
                               sigma = 2.0)
  Lapl <- ComputeGraphLaplacian(W)

  res.spectra <- ComputeEigenValuesWithRSpectra(Lapl, nrow(Lapl)-1)
  res.arpack <- ComputeEigenValuesWithIgraphArpack(Lapl, nrow(Lapl)-1)
  res.base <- eigen(Lapl)

  expect_less_than(mean( (sort(head(res.base$values,-1), decreasing = TRUE) - res.spectra$values)^2), 1e-25)
  expect_less_than(mean( (sort(head(res.base$values,-1), decreasing = TRUE) - res.arpack$values)^2), 1e-25)
  expect_less_than(mean( (res.spectra$values - res.arpack$values)^2), 1e-25)

})
