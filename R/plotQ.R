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
plot.tess3Q <- function(Q, coord, plot.type = "piechart", resolution = c(300,300), window = NULL, background = TRUE, raster.filename = NULL, interpolation.function = idw(), col = NULL, col.palette = NULL, map = TRUE, ...) {

  # parameters

  ## col
  if (is.null(col)) {
    col <- topo.colors(ncol(Q), alpha = 0.6)
  }

  ## col.palette
  if (is.null(col.palette)) {
    TestRequiredPkg("RColorBrewer")
    col.palette = list(
      c("gray95",RColorBrewer::brewer.pal(9,"Reds")),
      c("gray95",RColorBrewer::brewer.pal(9,"Greens")),
      c("gray95",RColorBrewer::brewer.pal(9,"Blues")),
      c("gray95",RColorBrewer::brewer.pal(9,"YlOrBr")),
      c("gray95",RColorBrewer::brewer.pal(9,"RdPu")),
      c("gray95",RColorBrewer::brewer.pal(9,"Greys")),
      c("gray95",RColorBrewer::brewer.pal(9,"Purples")),
      c("gray95",RColorBrewer::brewer.pal(9,"Oranges"))
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
