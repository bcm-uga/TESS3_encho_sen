context("Graph")



test_that("test of graph matrix computation", {

  data("data.for.test", package = "TESS3enchoSen")
  W <- ComputeHeatKernelWeight(coord = data.for.test$coord,
                               sigma = 2.0)

  Lapl <- ComputeGraphLaplacian(W)

  coord = data.for.test$coord
})


test_that("test of spectral decompostion of Lapl", {
  skip("debug test")

  data("data.for.test", package = "TESS3enchoSen")
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
