context("TESS3 project")



test_that("tess3project constructor", {

  # two rep
  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                            coord = data.for.test$coord,
                            K = 2:3,
                            ploidy = 1,
                            lambda = 1.0,
                            method = "projected.ls",
                            rep = 2,
                            keep = "all")

  expect_equal(class(tess3project.res)[2],"tess3")
  expect_equal(length(tess3project.res),3) # because length for list take last indice...
  expect_equal(length(tess3project.res[[2]]$tess3.run),2)
  summary(tess3project.res)

  # one rep and one K
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "projected.ls",
                                   rep = 1,
                                   keep = "all")
  expect_equal(class(tess3project.res)[2],"tess3Main")
  expect_equal(class(tess3project.res)[3],"tess3")
  summary(tess3project.res)

  # singulare cases
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = NULL,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "projected.ls",
                                   rep = 2,
                                   keep = "best")
  summary(tess3project.res)
  expect_error(tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = NULL,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "projected.ls",
                                   rep = -1,
                                   keep = "best"), "rep must greater than 1")
})


test_that("tess3project plot rmse", {

  # tess3
  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 1:2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "projected.ls",
                                   rep = 3,
                                   keep = "all")
  expect_error(plot(tess3project.res, crossvalid = TRUE),"tess3 was run with mask = 0. Run it with mask > 0.0 to have the cross validation rmse computed")
  plot(tess3project.res)
  plot(tess3project.res, crossentropy = TRUE)

  # tess3Main
  tess3project.res <- tess3(X = data.for.test$X,
                            coord = data.for.test$coord,
                            K = 1,
                            ploidy = 1,
                            lambda = 1.0,
                            method = "projected.ls",
                            rep = 1,
                            keep = "all")
  plot(tess3project.res)
})

test_that("tess3project plot crossvalid.rmse", {

  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 1:2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "projected.ls",
                                   rep = 3,
                                   keep = "all",
                                   mask = 0.2)

  plot(tess3project.res, main = "blabla", type = "l")
  plot(tess3project.res, crossvalid = FALSE, crossentropy = TRUE, main = "blabla", type = "l")
  plot(tess3project.res, crossvalid = TRUE, main = "blabla", type = "l")
  plot(tess3project.res, crossvalid = TRUE, crossentropy = TRUE, main = "blabla", type = "l")
})

test_that("Gettess3res", {

  # K = 2, rep = 1
  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                            coord = data.for.test$coord,
                            K = 2,
                            ploidy = 1,
                            lambda = 1.0,
                            method = "projected.ls",
                            rep = 1,
                            keep = "all")
  expect_false(!is.tess3(tess3project.res))
  expect_null(Gettess3res(tess3project.res, K = 2, rep = 2))
  expect_null(Gettess3res(tess3project.res, K = 1, rep = 1))
  expect_false(!is.tess3Main(Gettess3res(tess3project.res, K = 2, rep = "best")))
  expect_false(!is.tess3Main(Gettess3res(tess3project.res, K = 2, rep = 1)))

  # K = 2:3, rep = 3
  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                            coord = data.for.test$coord,
                            K = 2:3,
                            ploidy = 1,
                            lambda = 1.0,
                            method = "projected.ls",
                            rep = 3,
                            keep = "all")
  expect_false(!is.tess3(tess3project.res))
  tess3Main.res <- Gettess3res(tess3project.res, K = 2, rep = 1)
  expect_false(!is.tess3Main(tess3Main.res))
  expect_equal(ncol(tess3Main.res$Q), 2)
  tess3Main.res <- Gettess3res(tess3project.res, K = 3, rep = "best")
  expect_false(!is.tess3Main(tess3Main.res))
  expect_equal(ncol(tess3Main.res$Q), 3)
  ## null
  tess3Main.res <- Gettess3res(tess3project.res, K = 1, rep = "best")
  expect_null(tess3Main.res)
  tess3Main.res <- Gettess3res(tess3project.res, K = 1, rep = 23)
  expect_null(tess3Main.res)

})

test_that("tess3project bad use", {

  # structure format genotype matrix
  input.file = "http://membres-timc.imag.fr/Olivier.Francois/secondary_contact.str"
  mydata = read.table(input.file)
  mydata[1:10, 1:10]

  expect_error(obj <- tess3(as.matrix(mydata[,-(1:2)]), coord = mydata[,1:2], K = 5, ploidy = 1, openMP.core.num = 4),
               "The maximum value of the genotype matrix \\(X\\) cannot be greater than ploidy \\+ 1. Missing data must be encoded as NA\\.")


  data("data.at", package = "tess3r")
  expect_error(obj <- tess3(as.matrix(mydata[,-(1:2)]), coord = data.at$coord, K = 5, ploidy = 1, openMP.core.num = 4),
               "The maximum value of the genotype matrix \\(X\\) cannot be greater than ploidy \\+ 1. Missing data must be encoded as NA\\.")

  expect_error(obj <- tess3(data.at$X, coord = as.matrix(mydata[,1:2])[1:30,], K = 5, ploidy = 1, openMP.core.num = 4),
               "Number of row of coord and X/XBin must be the same")

})
