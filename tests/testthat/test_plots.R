context("Plot")



test_that("test of PlotAncestryCoefWithIdw", {

  skip_on_cran()
  # retrieve package file
  raster.filename <- system.file("extdata/raster","earth.tif",package = "TESS3enchoSen")

  # data
  K = 3
  at.geno = "~/PatatorHomeDir/Data/At/little_sample/Athaliana.geno"
  at.coord = "~/PatatorHomeDir/Data/At/little_sample/Athaliana.coord"
  data.at = list()
  data.at$X = LEA::read.geno(at.geno)
  data.at$coord = read.table(at.coord)

  # run tess3
  tess3enchosen.obj = TESS3enchoSen::TESS3(data.at$X,
                                           data.at$coord, K = K, ploidy = 1, lambda = 1.0)

  # interpolate with idw
  ancestry.raster <-RasterAncestryCoef(ancestry.coef = tess3enchosen.obj$Q[,2],
                                      coord = data.at$coord,
                                      interpolation.function = kriging(),
                                      resolution = list(ncol = 400,nrow=400))

    rasterVis::gplot(ancestry.raster) +
      ggplot2::geom_raster(ggplot2::aes(fill = value)) +
      ggplot2::scale_fill_gradient(low = 'grey', high = 'red', na.value = "white") +
      ggplot2::coord_equal() +
      ggplot2::geom_point(data = as.data.frame(data.at$coord),ggplot2::aes(x=V1,y=V2))

  # interpolate with kriging
  ancestry.raster <- RasterAncestryCoef(tess3enchosen.obj$Q[,1],data.at$coord, interpolation.function = kriging(),
                   resolution = list(ncol = 100,nrow=100))


  # plot all ancestral pop
  ancestry.raster.rgb <- RasterAncestryCoefRGB(Q = tess3enchosen.obj$Q,
                                               coord = data.at$coord,
                                               interpolation.function = idw(3),
                                              resolution = list(ncol = 100,nrow=100))
  breaks = seq(0,255,30)
  ancestry.raster.rgb.quant <- raster::calc( raster::cut(ancestry.raster.rgb,breaks = breaks), function(x) breaks[x])
  raster::plotRGB(ancestry.raster.rgb)


})
