context("Impute")

test_that("impute", {

  # sample data
  set.seed(2)
  n <- 100
  K <- 3
  ploidy <- 2
  L <- 3001
  data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
                                                  Q = SampleUnifQ(n, K),
                                                  coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
                                                  ploidy = ploidy)

  # mask data
  set.seed(0)
  masked.prop <- 0.1
  masked.X <- data.list$X
  masked.X[sample(1:(ncol(masked.X)*nrow(masked.X)), (ncol(masked.X)*nrow(masked.X)) * masked.prop)] <- NA


  # run tess3
  tess3.obj <- tess3(X = masked.X,
                     coord = data.list$coord,
                     K = 3,
                     ploidy = 2,
                     lambda = 1.0,
                     method = "projected.ls",
                     rep = 2)


  # impute
  expect_error(imputed.X <- ImputeRandom(tess3.obj, masked.X),
               "tess3.res must be a tess3 result of class.*")
  imputed.X.random <- ImputeRandom(Gettess3res(tess3.obj,3), masked.X)
  mean(data.list$X[is.na(masked.X)] != imputed.X.random[is.na(masked.X)]) # pourcentage of error

  imputed.X.round <- ImputeRound(tess3.res = Gettess3res(tess3.obj,3), masked.X)
  mean(data.list$X[is.na(masked.X)] != imputed.X.round[is.na(masked.X)]) # pourcentage of error
})
