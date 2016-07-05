context("genotype data functions")


test_that("test ComputeXBin", {
  set.seed(7456)
  n <- 100
  L <- 2000
  ploidy <- 2
  X <- matrix(as.double(rbinom(n * L, ploidy, 0.5)), n, L)
  XBin <- matrix(0.0, n, L * (ploidy + 1))
  X2XBin(X, ploidy, XBin)
  expect_equal( mean( XBin2X(XBin,ploidy) - X), 0)

  n <- 100
  L <- 2000
  ploidy <- 6
  X <- matrix(as.double(rbinom(n * L, ploidy, 0.5)), n, L)
  XBin <- matrix(0.0, n, L * (ploidy + 1))
  X2XBin(X, ploidy, XBin)
  expect_equal( mean( XBin2X(XBin,ploidy) - X), 0)

  # with missing value
  n <- 5
  L <- 20
  ploidy <- 6
  X <- matrix(as.double(rbinom(n * L, ploidy, 0.5)), n, L)
  mask <- sample(1:length(X), 0.5 * length(X))
  X[mask] <- NA
  XBin <- matrix(0.0, n, L * (ploidy + 1))
  X2XBin(X, ploidy, XBin)
  expect_equal(mean(is.na(XBin)), mean(is.na(X)))
})

test_that("test ComputeFst", {
  D <- 3
  L <- 100
  n <- 30
  K <- 3
  Q <- matrix(runif(n * K), n, K)
  Q <- ProjectQ(Q)
  G <- matrix(runif(L * D * K), L * D, K)
  G <- ProjectG(G, L, D)
  Fst <- ComputeFst(Q, G, D)
  expect_equal(length(Fst[Fst > 1.0]), 0)
  expect_equal(length(Fst[Fst < 0.0]), 0)
  expect_equal(ncol(Fst), 1)
  expect_equal(nrow(Fst), L)
})

test_that("test GtoFreq", {
  set.seed(878)
  n <- 198
  K <- 5
  ploidy <- 5
  L <- 1854
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)
  tess3.res <- tess3Main(X = data.list$X,
                     coord = data.list$coord,
                     K = K,
                     ploidy = ploidy,
                     lambda = 1.0,
                     method = "MCPA")
  Freq <- GtoFreq(tess3.res$G, ploidy)
  expect_equal(dim(Freq),c(L,K))
  expect_gte(min(Freq),0.0)
  expect_lte(max(Freq),1.0)
  expect_lte(abs(mean(data.list$X) / ploidy - mean(tcrossprod(tess3.res$Q, Freq))), 0.0001901845)
})
