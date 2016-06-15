###################################################
#######################Methods#####################
###################################################
#' This function displays a Structure-like barplot representation of the
#' ancestry coefficient matrix. It includes a sort-by-Q option.
#' @title Barplot representation of a Q-matrix
#' @author Olivier François
#'
#' @param Q an object of class \code{tess3Q} containing a matrix of ancestry coefficients computed from \code{tess3}.
#' @param sort.by.Q a Boolean value indicating whether the individuals should be sorted by their ancestry level or not.
#' @param col.palette is a list of color palettes. If \code{NULL}, a default list with 8 color palettes is used.
#' @param palette.length an integer value for the number of colors in each element of the palette list.
#' @return None
#' @examples
#' data(data.at)
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5, ploidy = 1, openMP.core.num = 4)
#' Qmatrix <- obj$Q
#' barplot(Qmatrix, border = NA, space = 0, xlab = "individuals", ylab = "Ancestry proportions", main = "Ancestry matrix")
#' @export
barplot.tess3Q = function(Q, sort.by.Q = TRUE, col.palette = NULL, palette.length = 9,...){
  if (class(Q) != "tess3Q") stop("Object Q not of class tess3Q.")
  ## color palette
  if (is.null(col.palette)) {
    if (ncol(Q) > 8)  stop("The default color palette contains 8 colors, and expects less than 9 clusters.")
    TestRequiredPkg("RColorBrewer")
    col.palette = list(
      c(RColorBrewer::brewer.pal(palette.length,"Reds")),
      c(RColorBrewer::brewer.pal(palette.length,"Greens")),
      c(RColorBrewer::brewer.pal(palette.length,"Blues")),
      c(RColorBrewer::brewer.pal(palette.length,"YlOrBr")),
      c(RColorBrewer::brewer.pal(palette.length,"RdPu")),
      c(RColorBrewer::brewer.pal(palette.length,"Greys")),
      c(RColorBrewer::brewer.pal(palette.length,"Purples")),
      c(RColorBrewer::brewer.pal(palette.length,"Oranges"))
    )
  }
  ## colors
  mycol = sapply(col.palette, FUN = function(x) x[palette.length/2])

  if (sort.by.Q){
    gr = apply(Q, MARGIN = 1, which.max)
    gm = max(gr)
    gr.o = order(sapply(1:gm, FUN = function(g) mean(Q[,g])))
    gr = sapply(gr, FUN = function(i) gr.o[i])
    or = order(gr)
    Qm = t(Q[or,])
    class(Qm) = "matrix"
    barplot(Qm, col = mycol, ...)
    }
  else {
    Qm = t(Q)
    class(Qm) = "matrix"
    barplot(Qm, col = mycol, ...)
  }
}


#' Plot individual ancestry coefficients on a map.
#'
#' @param Q an object of class \code{tess3Q} containing a matrix of ancestry coefficients computed from \code{tess3}.
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
    stop("col.palette must of length ncol(Q)")
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
