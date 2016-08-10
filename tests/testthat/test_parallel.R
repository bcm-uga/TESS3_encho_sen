context("parallel")

test_that("how fast is the parralel version ", {
  skip("Work only if openMP is installed")

  data("data.at", package = "tess3r")
  set.seed(377645)
  # run tess3
  time.1 <- system.time(tess3.res.1 <- tess3Main(X = data.at$X,
                                   coord = data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20, openMP.core.num = 1))

  time.4 <- system.time(tess3.res.4 <- tess3Main(X = data.at$X,
                                                     coord = data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20, openMP.core.num = 4))

  expect_lt(time.4[3], time.1[3])

  expect_lt(ComputeRmseWithBestPermutation(tess3.res.1$Q, tess3.res.4$Q), 0.02)
  expect_lt(ComputeRmseWithBestPermutation(tess3.res.1$G, tess3.res.4$G), 0.01)

})
