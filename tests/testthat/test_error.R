context("Error")


test_ComputeAveragedCrossEntropy <- function(P, Q, na.rm = FALSE, rm.logInfandNan = FALSE) {
  if (rm.logInfandNan) {
    Q[Q <= 0] <- 0.000001
  }
  aux <- -P * log(Q)
  return(mean(aux, na.rm = na.rm))
}


test_that("Compute Rmse and AveragedCrossEntropy", {

  Q1 <- matrix(as.integer(0.0), 3, 3)
  Q2 <- matrix(as.integer(1),3,3)
  expect_equal(ComputeRmse(Q1, Q2), sqrt(mean((Q1 - Q2) ^ 2)))
  expect_equal(ComputeRmseWithBestPermutationGreedy(Q1, Q2), ComputeRmse(Q1, Q2))
  expect_equal(ComputeRmseWithBestPermutation(Q1, Q2), ComputeRmse(Q1, Q2))
  expect_equal(ComputeAveragedCrossEntropy(Q1, Q2), test_ComputeAveragedCrossEntropy(Q1, Q2))

  Q1 <- runif(900)
  Q2 <- runif(900)
  expect_equal(ComputeRmse(Q1, Q2), sqrt(mean((Q1 - Q2) ^ 2)))
  expect_equal(ComputeAveragedCrossEntropy(Q1, Q2), test_ComputeAveragedCrossEntropy(Q1, Q2))

  set.seed(102)
  P <- runif(1000,-0.01, 1)
  Q <- runif(1000,-1, 1)
  expect_equal(ComputeAveragedCrossEntropy(P, Q), test_ComputeAveragedCrossEntropy(P, Q, rm.logInfandNan = TRUE))

  # with NA
  P <- runif(1000,-0.01, 1)
  Q <- runif(1000,-1, 1)
  P[sample(1:1000, 30)] <- NA
  Q[sample(1:1000, 30)] <- NA
  expect_equal(ComputeAveragedCrossEntropy(P, Q), test_ComputeAveragedCrossEntropy(P, Q, rm.logInfandNan = TRUE, na.rm = TRUE))
  expect_equal(ComputeRmse(P, Q), sqrt(mean((P - Q) ^ 2, na.rm = TRUE)))

  # try with bad size
  P <- rbinom(100, 1, 0.5)
  Q <- rbinom(1000, 1, 0.5)
  expect_error(ComputeRmse(P, Q), "Different vector size")
  expect_error(ComputeAveragedCrossEntropy(P, Q), "Different vector size")

  # debug
  # Q1 <- matrix(as.numeric(0.0),200,300000)
  # Q2 <- matrix(as.numeric(1.0),200,300000)
  # expect_equal(ComputeRmse(Q1, Q2), 1)
  # expect_equal(ComputeAveragedCrossEntropy(Q1, Q2), 0)
  # rm(Q1)
  # rm(Q2)
})


test_that("Compute Rmse", {

  Q1 <- matrix(as.integer(0.0),3,3)
  Q2 <- matrix(as.integer(1),3,3)
  expect_equal(ComputeRmse(Q1, Q2), sqrt(mean((Q1 - Q2) ^ 2)))
  expect_equal(ComputeRmseWithBestPermutationGreedy(Q1, Q2), ComputeRmse(Q1, Q2))
  expect_equal(ComputeRmseWithBestPermutation(Q1, Q2), ComputeRmse(Q1, Q2))

  Q1 <- matrix(as.numeric(0.0),3,3)
  Q2 <- matrix(as.numeric(1),3,3)
  expect_equal(ComputeRmse(Q1, Q2), sqrt(mean((Q1 - Q2) ^ 2)))
  expect_equal(ComputeRmseWithBestPermutationGreedy(Q1, Q2), ComputeRmse(Q1, Q2))
  expect_equal(ComputeRmseWithBestPermutation(Q1, Q2), ComputeRmse(Q1, Q2))

  Q1 <- matrix(as.numeric(0.0),200,300000)
  Q2 <- matrix(as.numeric(1.0),200,300000)
  expect_equal(ComputeRmse(Q1, Q2), 1)
  rm(Q1)
  rm(Q2)
})

test_that("rmse.tess3", {
  data("data.for.test", package = "tess3r")
  set.seed(878)
  tess3.res <- tess3Main(X = data.for.test$X,
                         coord = data.for.test$coord,
                         K = 6,
                         ploidy = 1,
                         lambda = 1.0,
                         method = "projected.ls")
  expect_lte(rmse.tess3Main(tess3.res, data.for.test$X, 1), 0.37)

  mask <- sample(1:(data.for.test$n * data.for.test$L), data.for.test$n * data.for.test$L * 0.25)
  expect_lte(rmse.tess3Main(tess3.obj = tess3.res, X = data.for.test$X, ploidy = 1, mask = mask),0.368)
})

test_that("Compute spatial reg", {

  data("data.for.test", package = "tess3r")
  set.seed(0)

  # compute input
  W <- ComputeHeatKernelWeight(data.for.test$coord, 2.0)
  Lapl <- as.matrix(ComputeGraphLaplacian(W))
  vpmax <- max(eigen(Lapl)$values)
  lambdaPrim <- 1.0
  X <- data.for.test$X
  XBin <- matrix(0.0, nrow(X), ncol(X) * (data.for.test$d + 1))
  X2XBin(X, data.for.test$d, XBin)
  # cast as double
  XBin <- matrix(as.double(XBin), nrow(XBin), ncol(XBin))

  # with K = 3
  K = 3
  cppres <- list()
  cppres$G <- matrix(0.0, nrow = (data.for.test$d + 1) * data.for.test$L, ncol = K)
  cppres$Q <- matrix(runif(data.for.test$n * K), data.for.test$n, K)
  cppres$Q <- ProjectQ(cppres$Q)
  ComputeSpatialPenalty(cppres$Q, W)

  aux <- capture.output(ComputeMCPASolution(XBin,
                                            K,
                                            Lapl,
                                            lambdaPrim = lambdaPrim,
                                            data.for.test$d + 1,
                                            maxIteration = 200,
                                            tolerance = 1e-10,
                                            Q = cppres$Q,
                                            G = cppres$G,
                                            verbose = FALSE))
  spatial.penalty.k3 <- ComputeSpatialPenalty(cppres$Q, W) / (K * data.for.test$n )

  # with K = 5
  K = 5
  cppres <- list()
  cppres$G <- matrix(0.0, nrow = (data.for.test$d + 1) * data.for.test$L, ncol = K)
  cppres$Q <- matrix(runif(data.for.test$n * K), data.for.test$n, K)
  cppres$Q <- ProjectQ(cppres$Q)
  ComputeSpatialPenalty(cppres$Q, W)

  aux <- capture.output(ComputeMCPASolution(XBin,
                                            K,
                                            Lapl,
                                            lambdaPrim = lambdaPrim,
                                            data.for.test$d + 1,
                                            maxIteration = 200,
                                            tolerance = 1e-10,
                                            Q = cppres$Q,
                                            G = cppres$G,
                                            verbose = FALSE))
  spatial.penalty.k5 <- ComputeSpatialPenalty(cppres$Q, W) / (K * data.for.test$n )

  # must be equal 0.1340827
  expect_lt(abs(spatial.penalty.k5 - spatial.penalty.k3),0.14)

})
