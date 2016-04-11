context("Plot")

test_that("test of pvalue plot", {
  data("data.for.test", package = "TESS3enchoSen")
  set.seed(0)
  tess3.res <- tess3(genotype = data.for.test$X,
                     geographic.coordinate = data.for.test$coord,
                     K = 3,
                     ploidy = 1,
                     lambda = 1.0,
                     method = "MCPA")

  p <- plot(tess3.res$pvalue)


})

test_that("test of Q plots ", {
  skip_if_not_installed("mapplots")
  skip_if_not_installed("maps")
  skip_if_not_installed("raster")
  skip_if_not_installed("gstat")
  skip_if_not_installed("automap")
  skip_if_not_installed("sp")
  skip_if_not_installed("grid")
  skip_if_not_installed("ggplot2")

  data("data.at", package = "TESS3enchoSen")
  set.seed(0)

  # run tess3
  tess3.res = TESS3enchoSen::tess3(data.at$X,
                                           data.at$coord, K = 4, ploidy = 1, lambda = 1.0)

  # piechart
  plot(tess3.res$Q, data.at$coord, plot.type = "piechart")
  plot(tess3.res$Q, data.at$coord, plot.type = "piechart", window = c(0,20,40,50))

  # florastyle
  plot(tess3.res$Q, data.at$coord, plot.type = "florastyle", resolution = c(300,300), window = NULL, background = TRUE, raster.filename = NULL, interpolation.function = NULL)
  plot(tess3.res$Q, data.at$coord, plot.type = "florastyle", resolution = c(100,100), window = c(0,30,40,60), background = FALSE, raster.filename = NULL, interpolation.function = NULL)

  # ggplot
  plot(tess3.res$Q, data.at$coord, plot.type = "ggplot", resolution = c(100,100), window = NULL, background = TRUE, raster.filename = NULL, interpolation.function = kriging())
  plot(tess3.res$Q, data.at$coord, plot.type = "ggplot", resolution = c(200,200), window = c(-15,50,30,70), background = TRUE, raster.filename = NULL, interpolation.function = idw(3))
  plot(Q = tess3.res$Q, coord = data.at$coord, plot.type = "ggplot", resolution = c(200,200), window = c(-15,50,30,70), background = FALSE, raster.filename = NULL, interpolation.function = idw(3))

})

test_that("test of plot with raster pkg", {

  skip("debug test")

  # retrieve package file
  raster.filename <- system.file("extdata/raster","earth.tif",package = "TESS3enchoSen")

  # data
  K = 4
  at.geno = "~/PatatorHomeDir/Data/At/little_sample/Athaliana.geno"
  at.coord = "~/PatatorHomeDir/Data/At/little_sample/Athaliana.coord"
  data.at = list()
  data.at$X = LEA::read.geno(at.geno)
  data.at$coord = as.matrix(read.table(at.coord))

  # run tess3
  tess3enchosen.obj = TESS3enchoSen::tess3(data.at$X,
                                           data.at$coord, K = K, ploidy = 1, lambda = 1.0)

  # interpolate with idw
  ancestry.raster <- RasterAncestryCoef(ancestry.coef = tess3enchosen.obj$Q[,2],
                                        coord = data.at$coord,
                                        interpolation.function = idw(),
                                        resolution = list(ncol = 400,nrow = 400))

  # plot with ggplot and rasterVis
  rasterVis::gplot(ancestry.raster) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient(low = 'grey', high = 'red', na.value = "white") +
    ggplot2::coord_equal() +
    ggplot2::geom_point(data = as.data.frame(data.at$coord),ggplot2::aes(x = V1, y = V2))

  # interpolate with kriging
  ancestry.raster <- RasterAncestryCoef(tess3enchosen.obj$Q[, 1], data.at$coord, interpolation.function = kriging(),
                                        resolution = list(ncol = 100, nrow = 100))

  # plot with raster pkg
  plot(ancestry.raster)

  # plot all ancestral pop rgb
  ancestry.raster.rgb <- RasterAncestryCoefRgb(Q = tess3enchosen.obj$Q,
                                               coord = data.at$coord,
                                               interpolation.function = idw(3),
                                              resolution = list(ncol = 100, nrow = 100))

  raster::plotRGB(ancestry.raster.rgb)

  # plot with ggplot
  ancestry.raster.stack <- RasterAncestryCoefStack(Q = tess3enchosen.obj$Q,
                                               coord = data.at$coord,
                                               interpolation.function = idw(3),
                                               resolution = list(ncol = 100, nrow = 100))

})

