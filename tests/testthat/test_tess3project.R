context("TESS3 project")



test_that("tess3project constructor", {

  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                     coord = data.for.test$coord,
                     K = 2:3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA",
                     rep = 2,
                     keep = "all")

  expect_equal(class(tess3project.res)[1],"tess3")
  expect_equal(length(tess3project.res),2)
  expect_equal(length(tess3project.res[[2]]$tess3.run),2)

  summary(tess3project.res)

  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "MCPA",
                                   rep = 2,
                                   keep = "all")
  expect_equal(length(tess3project.res[[1]]$tess3.run),2)


  summary(tess3project.res)

  # singulare cases
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = NULL,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "MCPA",
                                   rep = 2,
                                   keep = "best")
  summary(tess3project.res)
  expect_error(tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = NULL,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "MCPA",
                                   rep = -1,
                                   keep = "best"), "rep must greater than 1")


})


test_that("tess3project plot rmse", {

  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 1:2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "MCPA",
                                   rep = 3,
                                   keep = "all")
  expect_error(plot(tess3project.res, crossvalid = TRUE),"tess3 was run with mask = 0. Run it with mask > 0.0 to have the cross validation rmse computed")
  plot(tess3project.res)
  plot(tess3project.res, crossentropy = TRUE)
})

test_that("tess3project plot crossvalid.rmse", {

  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 1:2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "MCPA",
                                   rep = 3,
                                   keep = "all",
                                   mask = 0.2)

  plot(tess3project.res, main = "blabla", type = "l")
  plot(tess3project.res, crossvalid = FALSE, crossentropy = TRUE, main = "blabla", type = "l")
  plot(tess3project.res, crossvalid = TRUE, main = "blabla", type = "l")
  plot(tess3project.res, crossvalid = TRUE, crossentropy = TRUE, main = "blabla", type = "l")
})

test_that("tess3project is", {
  data("data.for.test", package = "tess3r")
  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 1:2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "MCPA",
                                   rep = 3,
                                   keep = "all")
  expect_false(!is.tess3(tess3project.res))
  expect_false(!is.tess3Main(Gettess3res(tess3project.res, K = 2, rep = 1)))
  expect_false(!is.tess3Main(Gettess3res(tess3project.res, K = 2, rep = "best")))

  tess3project.res <- tess3(X = data.for.test$X,
                                   coord = data.for.test$coord,
                                   K = 1:2,
                                   ploidy = 1,
                                   lambda = 1.0,
                                   method = "MCPA",
                                   rep = 3,
                                   keep = "best")
  expect_false(!is.tess3(tess3project.res))
  expect_false(!is.tess3Main(Gettess3res(tess3project.res, K = 2, rep = 2)))
  expect_false(!is.tess3Main(Gettess3res(tess3project.res, K = 2, rep = "best")))

})
