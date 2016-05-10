context("Error")



test_that("rmse.tess3", {
  data("data.for.test", package = "tess3r")
  set.seed(878)
  tess3.res <- tess3(X = data.for.test$X,
                     coord = data.for.test$coord,
                     K = 6,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA")
  expect_lte(rmse.tess3(tess3.res, data.for.test$X, 1), 0.37)

  mask <- sample(1:(data.for.test$n * data.for.test$L), data.for.test$n * data.for.test$L * 0.25)
  expect_lte(rmse.tess3(tess3.obj = tess3.res, X = data.for.test$X, ploidy = 1, mask = mask),0.368)
})

test_that("Compute spatial reg", {

  data("data.for.test", package = "tess3r")
  set.seed(0)

  # compute input
  W <- ComputeHeatKernelWeight(data.for.test$coord, 2.0)
  Lapl <- as.matrix(ComputeGraphLaplacian(W))
  vpmax <- max(eigen(Lapl)$values)
  lambdaPrim <- 1.0
  X = data.for.test$X
  XBin <- X2XBin(X,data.for.test$d)
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
                      maxIteration = 200, tolerance = 1e-10, Q = cppres$Q, G = cppres$G))
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
                      maxIteration = 200, tolerance = 1e-10, Q = cppres$Q, G = cppres$G))
  spatial.penalty.k5 <- ComputeSpatialPenalty(cppres$Q, W) / (K * data.for.test$n )

  # must be equal 0.1340827
  expect_less_than(abs(spatial.penalty.k5 - spatial.penalty.k3),0.14)

})
