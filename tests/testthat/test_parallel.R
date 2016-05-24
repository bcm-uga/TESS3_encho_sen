context("parallel")

test_that("how fast is the parralel version ", {
  skip_on_cran()
  data("data.at", package = "tess3r")
  set.seed(377645)
  # run tess3
  time.1 <- system.time(tess3.res.1 <- tess3r::tess3(X = data.at$X,
                                   coord = data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20, openMP.core.num = 1))

  time.4 <- system.time(tess3.res.4 <- tess3r::tess3(X = data.at$X,
                                                     coord = data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20, openMP.core.num = 4))

  expect_less_than(time.4[3], time.1[3])

  expect_less_than(ComputeRmseWithBestPermutation(tess3.res.1$Q, tess3.res.4$Q), 0.02)
  expect_less_than(ComputeRmseWithBestPermutation(tess3.res.1$G, tess3.res.4$G), 0.01)

})
