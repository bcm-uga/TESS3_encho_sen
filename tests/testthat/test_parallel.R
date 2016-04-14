context("parallel")

test_that("how fast is the parralel version ", {
  data("data.at", package = "tess3r")

  # run tess3
  system.time(tess3.res.1 <- tess3r::tess3(data.at$X,
                                   data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20, openMP.core.num = 1))

  system.time(tess3.res.8 <- tess3r::tess3(data.at$X,
                                                  data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20, openMP.core.num = 4))

  expect_less_than(ComputeRmseWithBestPermutation(tess3.res.1$Q, tess3.res.8$Q), 0.02)
  expect_less_than(ComputeRmseWithBestPermutation(tess3.res.1$G, tess3.res.8$G), 0.01)

})
