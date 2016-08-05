context("TESS3 main")

test_that("TESS3 copy", {
  set.seed(357467)
  n <- 100
  K <- 3
  ploidy <- 1
  L <- 3001
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)
  data.list$XBin <- matrix(0.0, n, L * (ploidy + 1))
  X2XBin(data.list$X, data.list$ploidy, data.list$XBin)

  expect_error(tess3.res <- tess3Main(X = data.list$X,
                                  XBin = NULL,
                                  coord = data.list$coord,
                                  K,
                                  ploidy,
                                  lambda = 1.0,
                                  W = NULL,
                                  method = "projected.ls",
                                  max.iteration = 200,
                                  tolerance = 1e-5,
                                  openMP.core.num = 4,
                                  Q.init = NULL,
                                  mask = 0.0,
                                  copy = FALSE),
               "To force the function not doing copy of the data, you must set XBin\\.")

  expect_error(tess3.res <- tess3Main(X = NULL,
                                  XBin = data.list$X,
                                  coord = data.list$coord,
                                  K,
                                  ploidy,
                                  lambda = 1.0,
                                  W = NULL,
                                  method = "projected.ls",
                                  max.iteration = 200,
                                  tolerance = 1e-5,
                                  openMP.core.num = 4,
                                  Q.init = NULL,
                                  mask = 0.0,
                                  copy = FALSE),
               "Number of columns of XBin must be a multiple of ploidy \\+ 1")

  tess3.res <- tess3Main(X = NULL,
                     XBin = data.list$XBin,
                     coord = data.list$coord,
                     K,
                     ploidy,
                     lambda = 1.0,
                     W = NULL,
                     method = "projected.ls",
                     max.iteration = 200,
                     tolerance = 1e-5,
                     openMP.core.num = 4,
                     Q.init = NULL,
                     mask = 0.0,
                     copy = FALSE)

})

test_that("TESS3 to dataframe", {
  set.seed(57467)
  n <- 100
  K <- 6
  ploidy <- 4
  L <- 3001
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)
  data.list$XBin <- matrix(0.0, n, L * (ploidy + 1))
  X2XBin(data.list$X, data.list$ploidy, data.list$XBin)

  tess3.res <- tess3Main(X = NULL,
                         XBin = data.list$XBin,
                         coord = data.list$coord,
                         K,
                         ploidy,
                         lambda = 1.0,
                         W = NULL,
                         method = "projected.ls",
                         max.iteration = 200,
                         tolerance = 1e-5,
                         openMP.core.num = 4)

  df <- data.frame(log.pvalue = tess3.res$log.pvalue,
                   fst = tess3.res$Fst,
                   Fscore = tess3.res$Fscore,
                   pvalue = tess3.res$pvalue,
                   index = seq_along(tess3.res$log.pvalue))


  df <- data.frame(tess3.res$Q)
  df <- data.frame(tess3.res$G)

})

test_that("TESS3 main, algo.copy and XBin", {
  set.seed(698)
  n <- 100
  K <- 3
  ploidy <- 2
  L <- 3000
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)
  data.list$XBin <- matrix(0.0, n, L * (ploidy + 1))
  X2XBin(data.list$X, data.list$ploidy, data.list$XBin)
  set.seed(687)
  tess3.res.nocopy <- tess3Main(X = NULL,
                     XBin = data.list$XBin,
                     coord = data.list$coord,
                     K,
                     ploidy,
                     lambda = 1.0,
                     W = NULL,
                     method = "projected.ls",
                     max.iteration = 200,
                     tolerance = 1e-5,
                     openMP.core.num = 4,
                     Q.init = NULL,
                     mask = 0.0,
                     algo.copy = FALSE)

  set.seed(687)
  tess3.res.copy <- tess3Main(X = NULL,
                            XBin = data.list$XBin,
                            coord = data.list$coord,
                            K,
                            ploidy,
                            lambda = 1.0,
                            W = NULL,
                            method = "projected.ls",
                            max.iteration = 200,
                            tolerance = 1e-5,
                            openMP.core.num = 4,
                            Q.init = NULL,
                            mask = 0.0,
                            algo.copy = TRUE)

  set.seed(687)
  tess3.res.copy.X <- tess3Main(X = data.list$X,
                          XBin = NULL,
                          coord = data.list$coord,
                          K,
                          ploidy,
                          lambda = 1.0,
                          W = NULL,
                          method = "projected.ls",
                          max.iteration = 200,
                          tolerance = 1e-5,
                          openMP.core.num = 4,
                          Q.init = NULL,
                          mask = 0.0,
                          algo.copy = TRUE)


  set.seed(687)
  tess3.res.nocopy.X <- tess3Main(X = data.list$X,
                            XBin = NULL,
                            coord = data.list$coord,
                            K,
                            ploidy,
                            lambda = 1.0,
                            W = NULL,
                            method = "projected.ls",
                            max.iteration = 200,
                            tolerance = 1e-5,
                            openMP.core.num = 4,
                            Q.init = NULL,
                            mask = 0.0,
                            algo.copy = FALSE)

  expect_lt(ComputeRmseWithBestPermutation(tess3.res.copy$Q, tess3.res.nocopy$Q), 1e-15)
  expect_lt(ComputeRmseWithBestPermutation(tess3.res.copy$G, tess3.res.nocopy$G), 1e-15)

  expect_lt(ComputeRmseWithBestPermutation(tess3.res.copy.X$G, tess3.res.nocopy$G), 1e-15)
  expect_lt(ComputeRmseWithBestPermutation(tess3.res.copy.X$Q, tess3.res.nocopy$Q), 1e-15)

  expect_lt(ComputeRmseWithBestPermutation(tess3.res.nocopy.X$G, tess3.res.nocopy$G), 1e-200)
  expect_lt(ComputeRmseWithBestPermutation(tess3.res.nocopy.X$Q, tess3.res.nocopy$Q), 1e-200)

})

test_that("TESS3 main with method MCPA", {

  data("data.for.test", package = "tess3r")
  set.seed(0)
  tess3.res <- tess3Main(X = data.for.test$X,
                     coord = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "projected.ls")

  # consistent error ?
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.03)
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.06)

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_lt(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.08)

  # with a Q.init
  set.seed(0)
  K = 3
  Q.init <- matrix(runif(nrow(data.for.test$X) * K), nrow(data.for.test$X), K)
  Q.init <- ProjectQ(Q.init)
  tess3.res <- tess3Main(X = data.for.test$X,
                     coord = data.for.test$coord,
                     K = K,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "projected.ls",
                     Q.init = Q.init)
  # consistent error ?
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.03)
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.06)


  # rmse and cross entropy
  expect_lte(tess3.res$rmse, 0.3712458)
  expect_lte(tess3.res$crossentropy, 0.2111163)
})

test_that("TESS3 main with method OQA", {

  data("data.for.test", package = "tess3r")
  set.seed(0)
  tess3.res <- tess3Main(X = data.for.test$X,
                     coord = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "qp",
                     tolerance = 0.00001)

  # consistent error ?
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.03)
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.06)

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_lt(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.08)

  # rmse and cross entropy
  expect_lte(tess3.res$rmse, 0.3712458)
  expect_lte(tess3.res$crossentropy, 0.2111269)

})



test_that("TESS3 main check arg", {
  data("data.for.test", package = "tess3r")

  expect_error(tess3.res <- tess3Main(X = data.for.test$X,
                                  coord = data.for.test$coord,
                                  K = 3,
                                  ploidy = 1,
                                  lambda = 1.0,
                                  method = "OA",
                                  tolerance = 0.00001),"method must be projected.ls or qp")

  W <- matrix(2, nrow = 3)
  expect_error(tess3.res <- tess3Main(X = data.for.test$X,
                                  coord = data.for.test$coord,
                                  K = 3,
                                  ploidy = 1,
                                  lambda = 1.0,
                                  W = W,
                                  tolerance = 0.00001),"W must be a squared symmetric matrix")
  W = matrix(runif(data.for.test$n ^ 2), data.for.test$n,data.for.test$n)
  expect_error(tess3.res <- tess3Main(X = data.for.test$X,
                                  coord = data.for.test$coord,
                                  K = 3,
                                  ploidy = 1,
                                  lambda = 1.0,
                                  W = W,
                                  tolerance = 0.00001),"W must be a squared symmetric matrix")

  expect_error(tess3Main(X = data.for.test$X,
                     coord = data.for.test$coord,
                     K = 3,
                     ploidy = 0,
                     lambda = 1.0,
                     method = "projected.ls",
                     tolerance = 0.00001),
               "The maximum value of the genotype matrix \\(X\\) cannot be greater than ploidy \\+ 1\\. Missing data must be encoded as NA\\.")

  genotype <- data.for.test$X
  genotype[1,1] <- -9
  expect_error(tess3Main(X = genotype,
                     coord = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "projected.ls",
                     tolerance = 0.00001),
               "Negative values in the genotype matrix \\(X\\) are not allowed\\. Missing data must be encoded as NA\\.")

  expect_error(tess3Main(X = data.for.test$X,
                     coord = data.for.test$coord[-1,],
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "projected.ls",
                     tolerance = 0.00001),"Number of row of coord and X/XBin must be the same")

})

test_that("TESS3 main with missing value", {
  data("data.for.test", package = "tess3r")

  # mask data
  set.seed(0)
  masked.prop <- 0.1
  masked.X <- data.for.test$X
  masked.X[sample(1:(ncol(masked.X)*nrow(masked.X)), (ncol(masked.X)*nrow(masked.X)) * masked.prop)] <- NA


  # run tess3 with MCPA
  tess3.res <- tess3Main(X = masked.X,
                     coord = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "projected.ls")

  # consistent error ?
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.05)
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.1)

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_lt(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.09)

  # rmse and cross entropy
  expect_lte(tess3.res$rmse, 0.3608639)
  expect_lte(tess3.res$crossentropy, 0.2309402)

  # run tess3 with OQA
  set.seed(0)
  tess3.res <- tess3Main(X = masked.X,
                     coord = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "qp")

  # consistent error ?
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$Q, tess3.res$Q), 0.05)
  expect_lt(ComputeRmseWithBestPermutation(data.for.test$G, tess3.res$G), 0.1)

  # Consistent Fst ?
  expect_lt(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.09)

  # rmse and cross entropy
  expect_lte(tess3.res$rmse, 0.3608441)
  expect_lte(tess3.res$crossentropy, 0.2309072)

})


test_that("TESS3 main Fst", {
  data("data.for.test", package = "tess3r")

  set.seed(0)
  tess3.res <- tess3Main(X = data.for.test$X,
                     coord = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "projected.ls")

  # Consistent Fst ?
  Fst <- ComputeFst(data.for.test$Q, data.for.test$G, data.for.test$d + 1)
  expect_lt(sqrt(mean((Fst - tess3.res$Fst) ^ 2)), 0.08)

  # Debug
  # hist(tess3.res$Fst)
  # hist(tess3.res$p.value)
  # BioCompToolsR::manhattanPlot(tess3.res$p.value)
  # tess3.res$gif


})


test_that("TESS3 cross validation", {
  set.seed(354354)
  n <- 100
  K <- 3
  ploidy <- 2
  L <- 1000
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)
  expect_message(tess3.res <- tess3Main(X = data.list$X,
                     coord = data.list$coord,
                     K = K,
                     ploidy = data.list$ploidy,
                     lambda = 1.0,
                     method = "projected.ls",
                     mask = 0.05), "Mask 0.05% of X for cross validation|Missing value detected in genotype")

  expect_lte(ComputeRmseWithBestPermutation(data.list$G, tess3.res$G), 0.0857528)
  expect_lte(ComputeRmseWithBestPermutation(data.list$Q, tess3.res$Q), 0.09725877)

  expect_lte(tess3.res$crossvalid.rmse, 0.3947187)
  expect_lte(tess3.res$rmse, 0.3947187)

  expect_lte(tess3.res$crossentropy, 0.5104649)
  expect_lte(tess3.res$crossvalid.crossentropy, 0.5292129)


  # With already missing values
  set.seed(54354)
  n <- 120
  K <- 4
  ploidy <- 3
  L <- 3000
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)
  data.list$X[sample(seq_along(data.list$X), 0.5 * length(data.list$X))] <- NA
  expect_message(tess3.res <- tess3Main(X = data.list$X,
                                    coord = data.list$coord,
                                    K = K,
                                    ploidy = data.list$ploidy,
                                    lambda = 1.0,
                                    method = "projected.ls",
                                    mask = 0.5), "Mask 0.05% of X for cross validation|Missing value detected in genotype")

  expect_lte(ComputeRmseWithBestPermutation(data.list$G, tess3.res$G), 0.2417362)
  expect_lte(ComputeRmseWithBestPermutation(data.list$Q, tess3.res$Q), 0.2253116)

  expect_lte(tess3.res$crossvalid.rmse, 0.3940670)
  expect_lte(tess3.res$rmse, 0.2698213) # because rmse is computed on naive impuation of X

  expect_lte(tess3.res$crossentropy, 0.8473548)
  expect_lte(tess3.res$crossvalid.crossentropy, 0.9144449)

})
