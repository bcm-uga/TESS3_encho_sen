###################################################
#######################Methods#####################
###################################################

#' Plot ancestry coeficient on a map.
#'
#' @param Q
#' @param coord
#' @param plot.type
#'
#' @return
#' @export
#'
#' @examples
plot.tess3Q <- function(Q, coord, plot.type = "piechart", resolution = c(300,300), window = NULL, background = TRUE, raster.filename = NULL, interpolation.function = idw(), col = NULL, col.palette = NULL, map = TRUE, palette.step = 9, ...) {

  # parameters

  ## col
  if (is.null(col)) {
    col <- topo.colors(ncol(Q), alpha = 0.6)
  }

  ## col.palette
  if (is.null(col.palette)) {
    TestRequiredPkg("RColorBrewer")
    col.palette = list(
      c(RColorBrewer::brewer.pal(palette.step,"Reds")),
      c(RColorBrewer::brewer.pal(palette.step,"Greens")),
      c(RColorBrewer::brewer.pal(palette.step,"Blues")),
      c(RColorBrewer::brewer.pal(palette.step,"YlOrBr")),
      c(RColorBrewer::brewer.pal(palette.step,"RdPu")),
      c(RColorBrewer::brewer.pal(palette.step,"Greys")),
      c(RColorBrewer::brewer.pal(palette.step,"Purples")),
      c(RColorBrewer::brewer.pal(palette.step,"Oranges"))
    )
  }
  if (length(col.palette) < ncol(Q) & (plot.type == "max" | plot.type == "all")) {
    stop("col.palette must of lenght ncol(Q)")
  }


  if (length(col) < ncol(Q)) {
    stop("col must be of length ncol(Q)")
  }

  ## compute of window if window is NULL
  if (is.null(window)) {
    window <- ComputeWindow(coord)
  }

  ##  raster filename
  if (is.null(raster.filename)) {
    raster.filename <- system.file("extdata/raster","earth.tif",package = "tess3r")
  }

  if (plot.type == "piechart") {
    PlotPiechartAncestryCoef(Q, coord, window, background, col, ...)
  } else if (plot.type == "max") {
    message("Compute grid and background.")
    grid <- ComputeGridAndBackground(window, resolution, background, raster.filename)
    message("Interpolate value on grid.")
    list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
    message("Plot.")
    PlotInterpotationMax(coord, list.grid.z, grid$grid.x, grid$grid.y, grid$background, col.palette, map, ...)
  } else if (plot.type == "all") {
    message("Compute grid and background.")
    grid <- ComputeGridAndBackground(window, resolution, background, raster.filename)
    message("Interpolate value on grid.")
    list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
    message("Plot.")
    PlotInterpotationAll(coord, list.grid.z, grid$grid.x, grid$grid.y, grid$background, col.palette, map, ...)
  }
}
