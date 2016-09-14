context("TESS3 on likelihoods genotype")

test_that("TESS3 on likelihoods genotype", {
  set.seed(873)
  n <- 100
  K <- 3
  ploidy <- 3
  L <- 5000
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)

  XProba <- tcrossprod(data.list$Q, data.list$G)

  tess3.res <- tess3Main(X = NULL,
                         XProba = XProba,
                         coord = data.list$coord,
                         K,
                         ploidy,
                         lambda = 0.0)

  expect_lt(ComputeRmseWithBestPermutation(tess3.res$Q, data.list$Q), 2e-4)
  expect_lt(ComputeRmseWithBestPermutation(tess3.res$G, data.list$G), 1e-4)

})
