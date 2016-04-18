context("TESS3 main")



test_that("TESS3 main with method MCPA", {

  data("data.for.test", package = "tess3r")
  set.seed(0)
  tess3.res <- tess3(genotype = data.for.test$X,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA")

  # consistent error ?
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.03)
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.06)

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_less_than(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.08)

  # with a Q.init
  set.seed(0)
  K = 3
  Q.init <- matrix(runif(nrow(data.for.test$X) * K), nrow(data.for.test$X), K)
  Q.init <- ProjectQ(Q.init)
  tess3.res <- tess3(genotype = data.for.test$X,
                     geographic.coordinate = data.for.test$coord,
                     K = K,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA",
                     Q.init = Q.init)
  # consistent error ?
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.03)
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.06)
})

test_that("TESS3 main with method OQA", {

  data("data.for.test", package = "tess3r")
  set.seed(0)
  tess3.res <- tess3(genotype = data.for.test$X,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "OQA",
                     tolerance = 0.00001)

  # consistent error ?
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.03)
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.06)

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_less_than(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.08)
})



test_that("TESS3 main check arg", {
  data("data.for.test", package = "tess3r")

  expect_error(tess3.res <- tess3(genotype = data.for.test$X,
                                  geographic.coordinate = data.for.test$coord,
                                  K = 3,
                                  ploidy = 1,
                                  lambda = 1.0,
                                  method = "OA",
                                  tolerance = 0.00001),".*Unknow method name.*")

  W = matrix(2, nrow = 3)
  expect_error(tess3.res <- tess3(genotype = data.for.test$X,
                                  geographic.coordinate = data.for.test$coord,
                                  K = 3,
                                  ploidy = 1,
                                  lambda = 1.0,
                                  W = W,
                                  tolerance = 0.00001),"W must be squared symetric of size.*")
  W = matrix(runif(data.for.test$n ^ 2), data.for.test$n,data.for.test$n)
  expect_error(tess3.res <- tess3(genotype = data.for.test$X,
                                  geographic.coordinate = data.for.test$coord,
                                  K = 3,
                                  ploidy = 1,
                                  lambda = 1.0,
                                  W = W,
                                  tolerance = 0.00001),"W must be squared symetric of size.*")

  expect_error(tess3(genotype = data.for.test$X,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 0,
                     lambda = 1.0,
                     method = "MCPA",
                     tolerance = 0.00001),".*The maximum value of the genotype matrix can not be superior than ploidy \\+ 1.*")

  genotype <- data.for.test$X
  genotype[1,1] <- -9
  expect_error(tess3(genotype = genotype,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA",
                     tolerance = 0.00001),".*Negative values in the genotype matrix are not allowed.*")

  expect_error(tess3(genotype = data.for.test$X,
                     geographic.coordinate = data.for.test$coord[-1,],
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA",
                     tolerance = 0.00001),".*Number of row in the coordinate matrix and the genotype matrix must be the same.*")



})

test_that("TESS3 main with missing value", {
  data("data.for.test", package = "tess3r")

  # mask data
  masked.prop <- 0.1
  masked.X <- data.for.test$X
  masked.X[sample(1:(ncol(masked.X)*nrow(masked.X)), (ncol(masked.X)*nrow(masked.X)) * masked.prop)] <- NA


  # run tess3 with MCPA
  set.seed(0)
  tess3.res <- tess3(genotype = masked.X,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA")

  # consistent error ?
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.05)
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.1)

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_less_than(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.09)

  # run tess3 with OQA
  set.seed(0)
  tess3.res <- tess3(genotype = masked.X,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "OQA")

  # consistent error ?
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.05)
  expect_less_than(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.1)

  # Consistent Fst ?
  expect_less_than(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.09)
})


test_that("TESS3 main Fst", {
  data("data.for.test", package = "tess3r")

  set.seed(0)
  tess3.res <- tess3(genotype = data.for.test$X,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA")

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_less_than(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.08)

  # Debug
  # hist(tess3.res$Fst)
  # hist(tess3.res$p.value)
  # BioCompToolsR::manhattanPlot(tess3.res$p.value)
  # tess3.res$gif


})

