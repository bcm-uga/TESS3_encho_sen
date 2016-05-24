context("Variogram")

test_that("test variogram", {
  data("data.for.test", package = "tess3r")
  em.vario <- CalculateEmpiricalGenSemivariogram(X = data.for.test$X, ploidy = 1, coord = data.for.test$coord, breaks = "FD", na.rm = TRUE)
  # ggplot2::ggplot(em.vario, ggplot2::aes(x = h, y = semi.variance, size = size)) + ggplot2::geom_point()

  em.vario <- CalculateEmpiricalGenSemivariogram(X = data.for.test$X * 2, ploidy = 2, coord = data.for.test$coord, breaks = "FD", na.rm = TRUE)
  # ggplot2::ggplot(em.vario, ggplot2::aes(x = h, y = log(semi.variance), size = size)) + ggplot2::geom_point()

  # calculate range nugget and still
  # em.vario.fit <- FitGeneralSemivariogram(semi.variogram = em.vario, epsilon = 1e-6)
  # expect_lte(abs(em.vario.fit$range - 1.655975), 1e-6 )

  #### test with Q unif
  set.seed(77587)
  set.seed(0)
  n <- 100
  K <- 3
  ploidy <- 2
  L <- 1000
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)

  em.vario <- CalculateEmpiricalGenSemivariogram(X = data.list$X, ploidy = ploidy, coord = data.list$coord, breaks = "FD", na.rm = TRUE)
  # ggplot2::ggplot(em.vario, ggplot2::aes(x = h, y = semi.variance, size = size)) + ggplot2::geom_point()
  # em.vario.fit <- FitGeneralSemivariogram(semi.variogram = em.vario, epsilon = 1e-6)

  #### test with Q as a function
  n <- 100
  K <- 2
  ploidy <- 2
  L <- 1000
  coord <- cbind(sort(c(rnorm(n/2, -2, 1), rnorm(n/2, 2, 1))), runif(n))
  Q <- SampleFuncQ(coord)
  G <- SampleUnifDirichletG(L, ploidy, K)
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = Q,
                                                  coord = coord,
                                                  ploidy = ploidy)
  em.vario <- CalculateEmpiricalGenSemivariogram(X = data.list$X, ploidy = ploidy, coord = data.list$coord, breaks = "FD", na.rm = TRUE)
  # ggplot2::ggplot(em.vario, ggplot2::aes(x = h, y = semi.variance, size = size)) + ggplot2::geom_point()
  # em.vario.fit <- FitGeneralSemivariogram(semi.variogram = em.vario, epsilon = 1e-6)

})


