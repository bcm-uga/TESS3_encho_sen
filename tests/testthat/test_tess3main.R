context("TESS3 main")



test_that("TESS3 main with method MCPA", {

  data("data.for.test", package = "TESS3enchoSen")
  set.seed(0)
  tess3.res <- TESS3(genotype = data.for.test$X,
                     coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA")

  # consistent error ?
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.2)
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.1)

})

test_that("TESS3 main with method OQA", {

  data("data.for.test", package = "TESS3enchoSen")
  set.seed(0)
  tess3.res <- TESS3(genotype = data.for.test$X,
                     coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "OQA")

  # consistent error ?
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.2)
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.1)
})
