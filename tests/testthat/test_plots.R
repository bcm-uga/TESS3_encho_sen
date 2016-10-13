context("Plot")

test_that("test of pvalue plot", {
  data("data.for.test", package = "tess3r")
  set.seed(0)
  tess3.res <- tess3Main(X = data.for.test$X,
                         coord = data.for.test$coord,
                         K = 3,
                         ploidy = 1,
                         lambda = 1.0,
                         method = "projected.ls")

  plot(tess3.res$pvalue)
})

test_that("test of Q plots", {

  data("data.at", package = "tess3r")
  set.seed(0)

  # run tess3
  tess3.res = tess3r::tess3(X = data.at$X,
                            coord = data.at$coord, K = 3, ploidy = 1, lambda = 1.0)

  # piechart
  expect_error(plot(tess3.res$Q, data.at$coord, method = "piechart"), "In development")


  # interpolation maps
  coord <- data.at$coord
  resolution <- c(300, 300)
  background <- TRUE
  Q <- tess3.res$Q
  ## window
  window <- ComputeWindow(coord)
  interpolation.model <- FieldsKrigModel()


  plot(tess3.res$Q, data.at$coord)
  plot(tess3.res$Q, data.at$coord, method = "map.all")


})



test_that("test of Q ggplots", {

  data("data.at", package = "tess3r")
  set.seed(0)

  # run tess3
  tess3.res = tess3r::tess3(X = data.at$X,
                            coord = data.at$coord, K = 3, ploidy = 1, lambda = 1.0)

  # interpolation maps
  coord <- data.at$coord
  Q <- tess3.res$Q


  ggtess3Q(Q, coord)

  ggtess3Q(Q, coord, resolution = c(200,200), window = NULL,
           background = TRUE, map.polygon = NULL,
           raster.filename = NULL,
           interpolation.model = FieldsTpsModel(),
           col.palette = CreatePalette())

})
