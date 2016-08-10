###################################################
#######################Methods#####################
###################################################
#' Displays a barplot representation of the
#' ancestry coefficient matrix. Includes a sort-by-Q option.
#' @title Barplot representation of a Q-matrix
#' @author Kevin Caye, Olivier François
#' @param Q an object of class \code{tess3Q} (Q matrix) containing a matrix of ancestry coefficients computed from \code{tess3} or converted from other program formats.
#' @param sort.by.Q a Boolean value indicating whether individuals should be sorted by their ancestry level or not.
#' @param col.palette is a list of color palettes. If \code{NULL}, a default list with 8 color palettes is used.
#' @param palette.length an integer value for the number of colors in each element of the palette list.
#' @return Generates a graphical output.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{as.qmatrix}} \code{\link{CreatePalette}}
#' @examples
#' data(data.at)
#' obj <- tess3Main(data.at$X, coord = data.at$coord, K = 5, ploidy = 1, openMP.core.num = 4)
#' Q.matrix <- obj$Q
#' barplot(Q.matrix, border = NA, space = 0, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
#' axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .6)
#' @export
barplot.tess3Q = function(Q, sort.by.Q = TRUE, col.palette = NULL, palette.length = 9, lab = FALSE, ...){
  if (class(Q)[1] != "tess3Q") {warning("Object Q not of class tess3Q.")}
  ## defines a default color palette (8 colors)
  if (is.null(col.palette)) {
    cat("Use CreatePalette() to define color palettes.\n")
    if (ncol(Q) > 8)
      stop("The default color palette expects less than 9 clusters.")
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
  if (palette.length > 3){
  colpal = sapply(col.palette, FUN = function(x) x[palette.length/2])}
  else {
  colpal = sapply(col.palette, FUN = function(x) x[palette.length])
  }

  if (sort.by.Q){
    gr = apply(Q, MARGIN = 1, which.max)
    gm = max(gr)
    gr.o = order(sapply(1:gm, FUN = function(g) mean(Q[,g])))
    gr = sapply(gr, FUN = function(i) gr.o[i])
    or = order(gr)
    Qm = t(Q[or,])
    class(Qm) = "matrix"
    barplot(Qm, col = colpal,...)
    return(list(order = or))
    }
  else {
    Qm = t(Q)
    class(Qm) = "matrix"
    barplot(Qm, col = colpal, ...)
  }
}


#' @title Displays geographic maps for ancestry coefficients
#' @author Kevin Caye, Flora Jay, Olivier François
#'
#' @param Q an object of class \code{tess3Q} containing a matrix of ancestry coefficients.
#' @param
#' @param col.palette a list of color palettes. If \code{NULL}, a default list with 8 color palettes is used.
#' @param palette.length an integer value for the number of colors in each element of the palette list.
#' @details details on interpolation methods
#' \itemize{
#' \item \code{kriging()} \code{\link[fields]{Krig}}
#' \item
#' }
#' @return None
#' @examples
#' data(data.at)
#' obj <- tess3Main(data.at$X, coord = data.at$coord, K = 5, ploidy = 1, openMP.core.num = 4)
#' Qmatrix <- obj$Q
#' plot(obj$Q, data.at$coord, method = "map.max", interpol = kriging(10), cex = .4, xlab = "Longitude", ylab= "Latitude", main = "Ancestry coefficients")
#' @export
plot.tess3Q <- function(Q, coord, method = "map.max", resolution = c(300,300), window = NULL,
                        background = TRUE, raster.filename = NULL, interpolation.function = kriging(10),
                        col = NULL, col.palette = NULL, map = TRUE, palette.length = 9, ...) {

  # parameters

  ## col
  if (is.null(col)) {
    col <- topo.colors(ncol(Q), alpha = 0.6)
  }

  ## col.palette
  if (is.null(col.palette)) {
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
  if (length(col.palette) < ncol(Q) & (method == "map.max" | method == "map.all")) {
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

  if (method == "piechart") {
    PlotPiechartAncestryCoef(Q, coord, window, background, col, ...)
  } else if (method == "map.max") {
    message("Computing grid and background...")
    grid <- ComputeGridAndBackground(window, resolution, background, raster.filename)
    message("Interpolating maximum ancestry coefficients on a grid...")
    list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
    message("Plotting the results...")
    PlotInterpotationMax(coord, list.grid.z, grid$grid.x, grid$grid.y, grid$background, col.palette, map, ...)
  } else if (method == "map.all") {
    message("Computing grid and background...")
    grid <- ComputeGridAndBackground(window, resolution, background, raster.filename)
    message("Interpolating all ancestry coefficients on a grid...")
    list.grid.z <- interpolation.function(Q, coord, grid$grid.x, grid$grid.y)
    message("Plotting the results...")
    PlotInterpotationAll(coord, list.grid.z, grid$grid.x, grid$grid.y, grid$background, col.palette, map, ...)
  }
}
