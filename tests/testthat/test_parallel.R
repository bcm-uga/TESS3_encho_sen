context("parallel")

test_that("how fast is the parralel version ", {
  skip_if_not_installed("awsMethods")

  data("data.at", package = "TESS3enchoSen")

  # run tess3
  awsMethods::setCores(1)
  system.time(tess3.res.1 <- TESS3enchoSen::tess3(data.at$X,
                                   data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20))

  awsMethods::setCores(4)
  system.time(tess3.res.8 <- TESS3enchoSen::tess3(data.at$X,
                                                  data.at$coord, K = 3, ploidy = 1, lambda = 1.0, tolerance = 1e-100, max.iteration = 20))

  expect_less_than(ComputeRmseWithBestPermutation(tess3.res.1$Q, tess3.res.8$Q), 0.01)
  expect_less_than(ComputeRmseWithBestPermutation(tess3.res.1$G, tess3.res.8$G), 0.01)

})
