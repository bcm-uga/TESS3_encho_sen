context("Sampler")

test_that("Sample with ms", {
  set.seed(757575)
  tess3.ms <- "~/BiocompSoftware/msdir/ms"
  n <- 200
  K <- 2
  ploidy <- 1
  data.list <- SampleGenoOFWithMs(n = n,
                                  nsites.neutral = 100000,
                                  nsites.selected = 1000,
                                  crossover.proba = 0.25 * 10 ^ -8,
                                  m.neutral = 0.25 * 10 ^ -6,
                                  m.selected = 0.25 * 10 ^ -7,
                                  mutation.rate.per.site = 0.25 * 10 ^ -8,
                                  N0 = 10 ^ 6,
                                  k = 0.5,
                                  min.maf = 0.05,
                                  plot.debug = FALSE,
                                  tess3.ms = tess3.ms)

  expect_equal(data.list$n,n)
  expect_equal(data.list$K,K)
  expect_equal(data.list$ploidy,ploidy)
  expect_equal(dim(data.list$admixed.genotype),c(n, data.list$L))
  expect_equal(dim(data.list$Q),c(n,K))
  expect_equal(max(data.list$admixed.genotype),ploidy)
  expect_equal(min(data.list$admixed.genotype),0)
  expect_equal(dim(data.list$coord),c(n,2))

  # test Z
  Q <- outer(1:n, 1:K, Vectorize(function(x,y) mean(data.list$Z[x,] == y)))
  expect_lt(ComputeRmseWithBestPermutation(Q, data.list$Q),0.006)

  # with only neutral
  n <- 200
  K <- 2
  ploidy <- 1
  data.list <- SampleGenoOFWithMs(n = n,
                                  nsites.neutral = 100000,
                                  nsites.selected = 0,
                                  crossover.proba = 0.25 * 10 ^ -8,
                                  m.neutral = 0.25 * 10 ^ -6,
                                  m.selected = 0.25 * 10 ^ -7,
                                  mutation.rate.per.site = 0.25 * 10 ^ -8,
                                  N0 = 10 ^ 6,
                                  k = 0.5,
                                  min.maf = 0.05,
                                  plot.debug = TRUE,
                                  tess3.ms = tess3.ms)

  expect_equal(data.list$n,n)
  expect_equal(data.list$K,K)
  expect_equal(data.list$ploidy,ploidy)
  expect_equal(dim(data.list$admixed.genotype),c(n, data.list$L))
  expect_equal(dim(data.list$Q),c(n,K))
  expect_equal(max(data.list$admixed.genotype),ploidy)
  expect_equal(min(data.list$admixed.genotype),0)
  expect_equal(dim(data.list$coord),c(n,2))

  # test Z
  Q <- outer(1:n, 1:K, Vectorize(function(x,y) mean(data.list$Z[x,] == y)))
  expect_lt(ComputeRmseWithBestPermutation(Q, data.list$Q),0.006)

})

test_that("Sample from TESS2.3 generative model", {
  set.seed(1548554)

  ## no trend surface
  K <- 3
  n.by.pop <- c(20,10,30)
  n <- sum(n.by.pop)
  coord <- SampleNormalClusterCoord(n.by.pop, K, 0.5, 0.2)
  # debug
  # plot(coord, col = rep(1:K, times = n.by.pop))
  W <- as.matrix(ComputeHeatKernelWeight(coord, 0.5))
  sigma <- 1.0
  vpmax <-  max(eigen(W)$values)
  rho <- 1 / vpmax / 2.0
  # debug
  # image(W)
  Q <- SampleTESS2.3Q(coord, K, W, sigma, rho, f = function(X) c(1.0,X[1], X[2]), beta = matrix(0.0, 3, K))
  expect_equal(dim(Q),c(n, K))
  expect_equal(apply(Q, 1, sum), rep(1,n))
  expect_gte(min(Q), 0.0)
  expect_lte(max(Q), 1.0)
  # debug
  # plot(Q, coord)

  ## with trend surface
  K <- 2
  n.by.pop <- c(30,30)
  n <- sum(n.by.pop)
  coord <- SampleNormalClusterCoord(n.by.pop, K, 0.5, 0.2)
  coord <- coord[order(coord[,1]),]
  # debug
  # plot(coord, col = rep(1:K, times = n.by.pop))
  W <- as.matrix(ComputeHeatKernelWeight(coord, 0.5))
  sigma <- 0.02
  vpmax <-  max(eigen(W)$values)
  rho <- 1 / vpmax / 2.0
  # debug
  # image(W)
  a <- 1 / (min(coord[,1]) -  max(coord[,1]))
  b <- 1 - a * min(coord[,1])
  # debug :
  # curve( exp( a + b * x ), min(coord[,1]), max(coord[,1]))
  # curve( exp( - a - b * x ), min(coord[,1]), max(coord[,1]))
  Q <- SampleTESS2.3Q(coord, K, W, sigma, rho, f = function(X) c(1.0,X[1], X[2]), beta = matrix(c(a,-b,0.0,-a,b,0.0), 3, K))
  expect_equal(dim(Q),c(n, K))
  expect_equal(apply(Q, 1, sum), rep(1,n))
  expect_gte(min(Q), 0.0)
  expect_lte(max(Q), 1.0)
  # debug
  # plot(Q, coord)
  # plot(Q[,2])

})


test_that("Sample from TESS3 generative model", {
  set.seed(0)
  n <- 100
  K <- 3
  ploidy <- 2
  L <- 1000
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                             Q = SampleUnifQ(n, K),
                                             coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                             ploidy = ploidy)
  expect_equal(data.list$n,n)
  expect_equal(data.list$K,K)
  expect_equal(data.list$ploidy,ploidy)
  expect_equal(dim(data.list$X),c(n, L))
  expect_equal(dim(data.list$Q),c(n,K))
  expect_equal(dim(data.list$G),c(L * (ploidy + 1), K))
  expect_equal(max(data.list$X),ploidy)
  expect_equal(min(data.list$X),0)
  expect_equal(dim(data.list$coord),c(n,2))


  tess3.res <- tess3(genotype = data.list$X,
                     geographic.coordinate = data.list$coord,
                     K = data.list$K,
                     ploidy = data.list$ploidy,
                     lambda = 1.0,
                     method = "MCPA")
  expect_lt(ComputeRmseWithBestPermutation(tess3.res$Q, data.list$Q), 0.082)
  expect_lt(ComputeRmseWithBestPermutation(tess3.res$G, data.list$G), 0.086)
})

test_that("Sample Q", {

  ## SampleUnifQ
  n <- 100
  K <- 3
  Q <- SampleUnifQ(n, K)
  expect_equal(dim(Q),c(n, K))
  expect_equal(apply(Q, 1, sum), rep(1,n))
  expect_gte(min(Q), 0.0)
  expect_lte(max(Q), 1.0)

  ## SampleDistFromCenterQ
  K <- 3
  n.by.pop <- c(20,10,30)
  n <- sum(n.by.pop)
  coord <- SampleNormalClusterCoord(n.by.pop, K, 0.5, 0.2)
  # debug
  # plot(coord, col = rep(1:K, times = n.by.pop))
  Q <- SampleDistFromCenterQ(coord, n.by.pop, K, f = function(D) exp(-D ^ 2 / 2.0))
  expect_equal(dim(Q),c(n, K))
  expect_equal(apply(Q, 1, sum), rep(1,n))
  expect_gte(min(Q), 0.0)
  expect_lte(max(Q), 1.0)
  # debug
  # plot(Q, coord)

  ## SampleDistFromCenterDirichletQ
  K <- 3
  n.by.pop <- c(20,10,30)
  n <- sum(n.by.pop)
  coord <- SampleNormalClusterCoord(n.by.pop, K, 0.5, 0.2)
  # debug
  # plot(coord, col = rep(1:K, times = n.by.pop))
  Q <- SampleDistFromCenterDirichletQ(coord, n.by.pop, K, f = function(D) exp(-D ^ 2 / 2.0))
  expect_equal(dim(Q),c(n, K))
  expect_equal(apply(Q, 1, sum), rep(1,n))
  expect_gte(min(Q), 0.0)
  expect_lte(max(Q), 1.0)
  # debug
  # plot(Q, coord)


  ## SampleFuncQ
  n <- 100
  K <- 2
  coord <- cbind(sort(c(rnorm(n/2, -2, 1), rnorm(n/2, 2, 1))), runif(n))
  Q <- SampleFuncQ(coord)
  expect_equal(dim(Q),c(n, K))
  expect_equal(apply(Q, 1, sum), rep(1,n))
  expect_gte(min(Q), 0.0)
  expect_lte(max(Q), 1.0)
  # debug
  # plot(Q, coord)
  # plot(Q[,1])
})

test_that("Sample G", {

  L <- 1000
  ploidy <- 3
  K <- 3
  G <- SampleUnifDirichletG(L, ploidy, K)
  expect_equal(dim(G),c((ploidy + 1) * L, K))
  G.array <- array(G, dim = c(ploidy + 1, L, K))
  expect_equal(as.vector(apply(G.array, c(2,3), sum)), rep(1, L * K))
  expect_gte(min(G), 0.0)
  expect_lte(max(G), 1.0)


})



test_that("Sample coord", {

  n.by.pop <- 30
  K <- 3
  coord <- SampleNormalClusterCoord(n.by.pop, K, 10, 0.1)
  expect_equal(dim(coord), c(n.by.pop * K, 2))
  # debug
  # plot(coord, col = rep(1:K, times = n.by.pop))

})

