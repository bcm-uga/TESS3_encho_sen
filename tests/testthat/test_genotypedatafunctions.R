context("genotype data functions")


test_that("test ComputeXBin", {
  M <- matrix(as.integer(c(1,0,2,0,1,0,1,0,1)),3,3)
  MBin <- ComputeXBin(M, 2)
  expect_equal( mean( ComputeXFromXBin(MBin,2) - M), 0)
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
  expect_equal(length(Fst[Fst>1.0]), 0)
  expect_equal(length(Fst[Fst<0.0]), 0)
  expect_equal(ncol(Fst), 1)
  expect_equal(nrow(Fst), L)
})
