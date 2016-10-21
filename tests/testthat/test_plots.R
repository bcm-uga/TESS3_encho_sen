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

  ### INTERPOLATION MAPS ###
  coord <- data.at$coord
  resolution <- c(300, 300)
  background <- TRUE
  Q <- tess3.res$Q

  # Testing a .asc rasterfile
  raster.filename <- system.file("extdata/raster","alps.asc",package = "tess3r")
  grid <- ComputeGridAndBackground(window=c(5,15,43,48), resolution, background=T, raster.filename, 10.0)  # long to compute
  image(grid$grid.x, grid$grid.y, grid$background) # must be europe

  raster.filename <- system.file("extdata/raster","earth.tif",package = "tess3r")
  ## window
  window <- ComputeWindow(coord)
  ## grid
  grid <- ComputeGridAndBackground(window, resolution, background, raster.filename)  # long to compute
  image(grid$grid.x, grid$grid.y, grid$background) # must be europe


  ## interpolation.function
  interpolation.function <- idw(1.0)
  list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
  image(grid$grid.x, grid$grid.y,list.grid.z[[2]] * grid$background)
  interpolation.function <- kriging()
  list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
  image(grid$grid.x, grid$grid.y, list.grid.z[[2]] * grid$background)


  plot(tess3.res$Q, data.at$coord, method = "map.max", main = "max", xlab = "x", ylab = "y",
       resolution = c(300,300),
       window = NULL,
       background = TRUE,
       raster.filename = NULL,
       interpolation.function = idw(2.0), map = TRUE)

  # Testing backgound and outline:
  plot(tess3.res$Q, data.at$coord, method = "map.max", main = "max", xlab = "x", ylab = "y",
       resolution = c(50,50), background = FALSE, map = FALSE)
  plot(tess3.res$Q, data.at$coord, method = "map.max", main = "max", xlab = "x", ylab = "y",
       resolution = c(50,50), background = FALSE, map = TRUE)
  plot(tess3.res$Q, data.at$coord, method = "map.max", main = "max", xlab = "x", ylab = "y",
       resolution = c(50,50), background = TRUE, map = FALSE)
  plot(tess3.res$Q, data.at$coord, method = "map.max", main = "max", xlab = "x", ylab = "y",
       resolution = c(50,50), background = TRUE, map = TRUE)


  gr = NULL
  gr = plot(tess3.res$Q, data.at$coord, method = "map.all", main = "all", xlab = "x", ylab = "y",
       resolution = c(100,100),
       window = NULL,
       background = TRUE,
       raster.filename = NULL,
       interpolation.function = idw(2.0), map = TRUE)
  expect_identical(is.null(gr), FALSE)

  # Use the previsouly computed grid + change interpolation
  plot(tess3.res$Q, data.at$coord, method = "map.max", main = "max", xlab = "x", ylab = "y",
       grid = gr,
       resolution = c(100,100), # this has no impact here bc grid is provided
       window = NULL,
       background = TRUE,
       raster.filename = NULL,
       interpolation.function = kriging(), map = TRUE)



  ## Piechart internal: PlotPiechartAncestryCoef ##
  PlotPiechartAncestryCoef(tess3.res$Q, data.at$coord,
                           window=c(0,20,40,50), map=T,
                           col=rainbow(ncol(tess3.res$Q)),
                           add.pie = FALSE, names.pie = NULL)

  PlotPiechartAncestryCoef(tess3.res$Q, data.at$coord,
                           window=c(0,20,40,50), map=F,
                           col=rainbow(ncol(tess3.res$Q)), radius=.05,
                           add.pie = FALSE, names.pie = NULL)

  info.bypop.nocoord = bypop(tess3.res$Q, pop=data.at$countries)
  info.bypop = bypop(tess3.res$Q, data.at$coord, data.at$countries)

  PlotPiechartAncestryCoef(info.bypop$Q, info.bypop$coord,
                           map=T, window=c(0,20,40,50),
                           col=rainbow(ncol(tess3.res$Q)), radius=.01,
                           add.pie = FALSE, names.pie = info.bypop$labels)

  PlotPiechartAncestryCoef(info.bypop$Q, info.bypop$coord,
                           map=T, window=ComputeWindow(coord),
                           col=rainbow(ncol(tess3.res$Q)), radius=.01,
                           add.pie = FALSE, names.pie = info.bypop$labels)

  # pie size proportional to sample size
  PlotPiechartAncestryCoef(info.bypop$Q, info.bypop$coord,
                           map=T, window=c(-10,50,30,70),
                           col=rainbow(ncol(tess3.res$Q)), radius=.02,
                           radius.prop=info.bypop$size,
                           add.pie = FALSE, names.pie = info.bypop$labels,cex=.6)

  # legend bubble depends only on samples inside the window
  PlotPiechartAncestryCoef(info.bypop$Q, info.bypop$coord,
                           map=T, window=c(10,25,45,55),
                           col=rainbow(ncol(tess3.res$Q)), radius=.05,
                           radius.prop=info.bypop$size,
                           add.pie = FALSE, names.pie = info.bypop$labels,cex=.6)

  PlotPiechartAncestryCoef(info.bypop$Q, info.bypop$coord,
                           map=T, window=c(-10,50,30,70),
                           col=rainbow(ncol(tess3.res$Q)), radius=.02,
                           radius.prop=log10(1+info.bypop$size),
                           names.pie = info.bypop$labels,cex=.7)

  ## PIECHARTS external functions ##
  # by individual
  piechartQ(tess3.res$Q, data.at$coord, method = "piechart", main = "piechart", xlab = "x", ylab = "y", map = FALSE)

  plot(tess3.res$Q, data.at$coord, method = "piechart", main = "piechart", xlab = "x", ylab = "y",map = TRUE)

  piechartQ(tess3.res$Q, data.at$coord, method = "piechart", window = c(0,20,40,50), main = "piechart",
           xlab = "x", ylab = "y", col = rainbow(ncol(tess3.res$Q)))

  piechartQ(tess3.res$Q, data.at$coord, method = "piechart", window = c(0,20,40,50),
           main = "piechart", xlab = "x", ylab = "y", col = rainbow(ncol(tess3.res$Q)), radius = .02)
  plot(tess3.res$Q, data.at$coord, method = "piechart", window = c(0,20,40,50),
           main = "piechart", xlab = "x", ylab = "y", col = rainbow(ncol(tess3.res$Q)), radius = .02)
  # by pop
  piechartQ(tess3.res$Q, data.at$coord, method = "piechart.pop",
      pop=data.at$countries,cex=.6)

  plot(tess3.res$Q, data.at$coord, method = "piechart.pop",
       pop=as.factor(data.at$countries), scale=T)

  plot(tess3.res$Q, data.at$coord, method = "piechart.pop",
       pop=as.factor(data.at$countries), scale=T, legend=F)

  plot(tess3.res$Q, data.at$coord, method = "piechart.pop",
       pop=data.at$countries,
       radius = 0.08, scale=T,
       window = c(0,20,40,50))

  plot(tess3.res$Q, data.at$coord, method = "piechart.pop",
       pop=data.at$countries,
       radius = .08, scale=T,
       window = c(0,20,40,55), label.distx=.5, leg.bubble=NULL)

  ### Interpolation maps + pie charts on top ###
  # pie charts by individual
  gr = plot(tess3.res$Q, data.at$coord, method = "map.max",
       resolution=c(100,100),
       interpol=kriging(10))
  plot(tess3.res$Q, data.at$coord, method = "piechart",add.pie=T)

  # by groups
  poprnd = sample(1:5,replace=T,nrow(tess3.res$Q))
  poprnd[1] = 6 # force a group to have only one element (this should not result in an error)
  plot(tess3.res$Q, data.at$coord, method = "map.max",grid=gr)
  plot(tess3.res$Q, data.at$coord, method = "piechart.pop",
       pop=poprnd, add=T)

  ## map with legend + graphics.reset + piechart
  plot(tess3.res$Q, data.at$coord, method = "map.max",grid=gr, legend=T, graphics.reset=TRUE)
  plot(tess3.res$Q, data.at$coord, method = "piechart",add.pie=F)

  plot(tess3.res$Q, data.at$coord, method = "map.max",grid=gr, legend=T, graphics.reset=FALSE)
  plot(tess3.res$Q, data.at$coord, method = "piechart",add.pie=T)
  dev.off()

  plot(tess3.res$Q, data.at$coord, method = "map.all",grid=gr, legend=T,graphics.reset=FALSE)
  plot(tess3.res$Q, data.at$coord, method = "piechart",add.pie=T)
  dev.off()

  plot(tess3.res$Q, data.at$coord, method = "map.all",grid=gr, legend=T,graphics.reset=TRUE)
  plot(tess3.res$Q, data.at$coord, method = "piechart",add.pie=F)

  # Generqates a warning
  plot(tess3.res$Q, data.at$coord, method = "map.all",grid=gr, legend=T,graphics.reset=TRUE,legend.width=.5,legend.lab="test")


  plot(tess3.res$Q, data.at$coord, method = "map.all",grid=gr,
       legend=T,graphics.reset=TRUE,legend.width=.5,leg.extra.args=list(legend.lab="test"))

})


