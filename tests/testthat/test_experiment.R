context("experiment")



test_that("tess3 time and it", {
  data("data.for.test", package = "tess3r")
  set.seed(878)
  tim <- system.time(tess3.res <- tess3(X = data.for.test$X,
                                        coord = data.for.test$coord,
                                        K = 6,
                                        ploidy = 1,
                                        lambda = 1.0,
                                        method = "MCPA"))
  expect_equal(tess3.res$it, 34)
  expect_lte(sum(tess3.res$times, na.rm = TRUE), tim[3])


  tim <- system.time(tess3.res <- tess3(X = data.for.test$X,
                                        coord = data.for.test$coord,
                                        K = 2,
                                        ploidy = 1,
                                        lambda = 1.0,
                                        method = "OQA"))
  expect_equal(tess3.res$it, 15)
  expect_lte(sum(tess3.res$times, na.rm = TRUE), tim[3])

})
