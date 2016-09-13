context("TESS3 on true data")


test_that("TESS3 on arabidopsis thaliana", {
  if (Sys.info()["nodename"] != "timc-bcm-15.imag.fr") {
    skip("test only on timc-bcm-15.imag.fr")
  }
  skip("Too long")
  data.pipe <- pipe('ssh cayek@patator.imag.fr "cat /home/cayek/Data/At/At.RData"')
  load(data.pipe)
  At$XBin <- matrix(0.0, nrow(At$X), ncol(At$X) * (2 + 1))
  X2XBin(At$X, 2, At$XBin)

  tess3.res <- tess3Main(X = NULL,
                     XProba = At$XBin,
                     coord = At$coord,
                     K = 3,
                     ploidy = 2,
                     lambda = 1.0,
                     W = NULL,
                     method = "projected.ls",
                     max.iteration = 200,
                     tolerance = 1e-5,
                     openMP.core.num = 4,
                     Q.init = NULL,
                     mask = 0.0)
  gc()
  tess3project.res <- tess3(X = NULL,
                                   XProba = At$XBin,
                                   coord = At$coord,
                                   K = 3,
                                   ploidy = 2,
                                   lambda = 1.0,
                                   rep = 1,
                                   W = NULL,
                                   method = "projected.ls",
                                   max.iteration = 200,
                                   tolerance = 1e-5,
                                   openMP.core.num = 4,
                                   Q.init = NULL,
                                   mask = 0.0,
                                   keep = "all")
  gc()

  palette.step = 9
  col.palette = list(
    c(RColorBrewer::brewer.pal(palette.step,"Reds"))[4:9],
    c(RColorBrewer::brewer.pal(palette.step,"Greens"))[4:9],
    c(RColorBrewer::brewer.pal(palette.step,"Blues"))[4:9]
  )
  plot(tess3.res$Q, data.at$coord, plot.type = "max", main = "max", xlab = "x", ylab = "y",
       resolution = c(300,300),
       window = c(-15,45,30,65),
       background = TRUE,
       raster.filename = NULL,
       interpolation.function = kriging(), map = TRUE,
       col.palette = col.palette)

})
