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
  skip_if_not_installed("mapplots")
  skip_if_not_installed("maps")
  skip_if_not_installed("fields")
  skip_if_not_installed("automap")
  skip_if_not_installed("sp")
  skip_if_not_installed("RColorBrewer")

  data("data.at", package = "tess3r")
  set.seed(0)

  # run tess3
  tess3.res = tess3r::tess3(X = data.at$X,
                            coord = data.at$coord, K = 3, ploidy = 1, lambda = 1.0)

  # piechart
  plot(tess3.res$Q, data.at$coord, plot.type = "piechart", main = "piechart", xlab = "x", ylab = "y", background = FALSE)
  plot(tess3.res$Q, data.at$coord, plot.type = "piechart", window = c(0,20,40,50), main = "piechart", xlab = "x", ylab = "y", col = rainbow(4))
  plot(tess3.res$Q, data.at$coord, plot.type = "piechart", window = c(0,20,40,50), main = "piechart", xlab = "x", ylab = "y", col = rainbow(4), radius = 0.25)

  # interpolation maps
  coord <- data.at$coord
  resolution <- c(300, 300)
  background <- TRUE
  raster.filename <- system.file("extdata/raster","earth.tif",package = "tess3r")
  Q <- tess3.res$Q
  ## window
  window <- ComputeWindow(coord)
  expect_equal(window, c(-8.848306, 38.693732, 37.017768, 63.595338))
  ## grid
  grid <- ComputeGridAndBackground(window, resolution, background, raster.filename)
  image(grid$grid.x, grid$grid.y, grid$background) # must be europe
  ## interpolation.function
  interpolation.function <- idw(1.0)
  list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
  image(grid$grid.x, grid$grid.y,list.grid.z[[2]] * grid$background)
  interpolation.function <- kriging()
  list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
  image(grid$grid.x, grid$grid.y, list.grid.z[[2]] * grid$background)

  plot(tess3.res$Q, data.at$coord, plot.type = "max", main = "max", xlab = "x", ylab = "y",
       resolution = c(300,300),
       window = NULL,
       background = FALSE,
       raster.filename = NULL,
       interpolation.function = idw(2.0), map = FALSE)

  plot(tess3.res$Q, data.at$coord, plot.type = "all", main = "max", xlab = "x", ylab = "y",
       resolution = c(300,300),
       window = NULL,
       background = FALSE,
       raster.filename = NULL,
       interpolation.function = idw(2.0), map = TRUE)

  plot(tess3.res$Q, data.at$coord, plot.type = "max", main = "max", xlab = "x", ylab = "y",
       resolution = c(300,300),
       window = NULL,
       background = TRUE,
       raster.filename = NULL,
       interpolation.function = kriging(), map = TRUE)

  plot(tess3.res$Q, data.at$coord, plot.type = "all", main = "max", xlab = "x", ylab = "y",
       resolution = c(300,300),
       window = NULL,
       background = TRUE,
       raster.filename = NULL,
       interpolation.function = kriging(), map = TRUE)

})

