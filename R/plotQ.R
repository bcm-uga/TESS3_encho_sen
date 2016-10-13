###################################################
#######################helpers#####################
###################################################

CheckPlotParam <- function(x, coord,
                           method,
                           resolution, window,
                           background, map.polygon,
                           raster.filename,
                           interpolation.model,
                           col.palette) {
  ## Params checks
  if (length(col.palette) < ncol(x) & (method == "map.max" | method == "map.all")) {
    stop("col.palette must of length ncol(Q)")
  }

}

ComputeInterpolStack <- function(Q, coord,
                                 interpolation.model, resolution, window,
                                 background, map.polygon, raster.filename) {
  ## Compute of window if window is NULL
  if (is.null(window)) {
    window <- ComputeWindow(coord)
  }

  ## Make grid
  raster.grid <- MakeRasterGrid(window, resolution)

  ## Map polygon
  if (background && is.null(map.polygon) && is.null(raster.filename)) {
    TestRequiredPkg("rworldmap")
    map.polygon <- rworldmap::getMap(resolution = "low")
  }

  ## Interpolate
  interpol.stack <- InterpolRaster(coord, Q, raster.grid, interpolation.model)


  ## Crop
  if (background) {
    interpol.stack <- Crop(interpol.stack, map.polygon, raster.filename)
  }
  return(interpol.stack)
}

RasterStackToMatrix <- function(interpol.stack) {
  res <- list()
  for (i in seq_along(names(interpol.stack))) {
    res[[i]] <- raster::as.matrix(interpol.stack[[i]])
    res[[i]] <- t(apply(t(res[[i]]), 1, rev))# transpose and rev col to plot ...
  }
  return(res)
}

###################################################
#######################Utils#######################
###################################################

#' This function creates a list of color palettes for the plot and barplot
#' functions
#' @title Create a list of palettes
#' @author Kevin Caye, Olivier François
#' @param color.vector a vector of R colors.
#' @param palette.length an integer number of colors in each palette.
#' @return An object of class \code{list} containing a list of color palettes.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{barplot.tess3Q}}
#' @examples
#' library(tess3r)
#'
#' ## an A. thaliana example
#' data(data.at)
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5, ploidy = 1,
#'              openMP.core.num = 4)
#' Qmatrix <- qmatrix(obj,K=5)
#' my.colors <- c("tomato", "yellow", "blue", "wheat","olivedrab")
#' my.palette <- CreatePalette(my.colors, 9)
#' plot(obj$Q, data.at$coord, method = "mapping.max", col.palette = my.palette,
#'      cex = .4, xlab = "Longitude", ylab= "Latitude",
#'      main = "Ancestry coefficients")
#' @export
CreatePalette <- function(color.vector = c("tomato", "chartreuse", "gold", "blue", "violet", "wheat","olivedrab"), palette.length = 9){
  ll = NULL
  for (i in 1:length(color.vector)) {
    ll[[i]] = colorRampPalette(c("grey96", color.vector[i]))(palette.length)
  }
  ll
}

###################################################
#######################Methods#####################
###################################################
#' Displays a barplot representation of the
#' ancestry coefficient matrix. Includes a sort-by-Q option.
#' @title Barplot representation of a Q-matrix
#' @author Kevin Caye, Olivier François
#'
#' @param height an object of class \code{tess3Q} (Q matrix) containing a matrix of ancestry coefficients computed from \code{tess3} or converted from other program formats.
#' @param sort.by.Q a Boolean value indicating whether individuals should be sorted by their ancestry level or not.
#' @param col.palette is a list of color palettes. If \code{NULL}, a default list with 8 color palettes is used.
#' @param lab a list of individual labels.
#' @param ... other parameters of the function \code{\link{barplot.default}}.
#' @param palette.length an integer value for the number of colors in each element of the palette list.
#'
#' @return Generates a graphical output.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{as.qmatrix}} \code{\link{CreatePalette}}
#' @examples
#' library(tess3r)
#'
#' # Retrieve a dataset
#' data(data.at)
#'
#' # Run of TESS3
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5,
#'                  ploidy = 1, method = "projected.ls", openMP.core.num = 4)
#'
#' # Get the ancestry matrix
#' Q.matrix <- qmatrix(obj, K = 5)
#'
#' # Plot the barplot
#' barplot(Q.matrix, border = NA, space = 0, xlab = "Individuals",
#'         ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
#' axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .4)
#' @export
barplot.tess3Q = function(height, sort.by.Q = TRUE, col.palette = NULL, palette.length = 9, lab = FALSE, ...){
  # Because prototype of this S3 method must be the same that barplot.default
  Q = height
  rm(height)

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
  if (palette.length > 3) {
  colpal = sapply(col.palette, FUN = function(x) x[palette.length/2])}
  else {
  colpal = sapply(col.palette, FUN = function(x) x[palette.length])
  }

  if (sort.by.Q) {
    gr = apply(Q, MARGIN = 1, which.max)
    gm = max(gr)
    gr.o = order(sapply(1:gm, FUN = function(g) mean(Q[,g])))
    gr = sapply(gr, FUN = function(i) gr.o[i])
    or = order(gr)
    Qm = t(Q[or,])
    class(Qm) = "matrix"
    graphics::barplot(Qm, col = colpal,...)
    return(list(order = or))
    }
  else {
    Qm = t(Q)
    class(Qm) = "matrix"
    graphics::barplot(Qm, col = colpal, ...)
  }
}


#' Displays geographic maps for ancestry coefficients.
#' @title Displays geographic maps for ancestry coefficients
#' @author Kevin Caye, Flora Jay, Olivier François
#'
#' @param coord The numeric matrix of size \eqn{n \times 2} where \eqn{n} is the
#' number of individuals.
#' @param resolution An integer vector of the resolution of the grid used to
#' computed the interpolating surface
#' @param window The window size, such that \code{window = c(xmin, xmax, ymin, ymax)}
#' contains the window mina nd max coordinates.
#' @param background If TRUE the raster file is used as a stencil to render only
#' raster pixel on earth.
#' @param map.polygon The \code{sp::SpatialPolygonsDataFrame} object used to crop the interpolating surfaces.
#' If \code{NULL}, the function \code{\link[rworldmap]{getMap}} is used.
#' @param raster.filename The raster file name used to compute the background stencil.
#' This is an alternative method to crop the interpolating surfaces. The default method uses \code{map.polygon}.
#' @param col.palette A list of color palettes. You can use the function \code{\link{CreatePalette}}.
#' @param interpolation.model The interpolation model used to compute compute
#' the interpolating surface. You can use functions \code{\link{FieldsTpsModel}},
#' \code{\link{FieldsKrigModel}}.
#' @param x An object of class \code{tess3Q}. The ancestry coefficient matrix.
#' @param method \code{"map.max"} or \code{"map.all"}. If \code{"map.all"}, a
#' interpolating surface of the ancestry coefficient is plotted for each ancestral
#' population. If \code{"map.max"} only the maximum union of the interpolating surfaces is plotted.
#' @param ... \code{\link{plot.default}} other parameters.
#'
#' @return None
#' @examples
#' library(tess3r)
#'
#' # Retrieve a dataset
#' data(data.at)
#'
#' # Run of TESS3
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5,
#'                  ploidy = 1, method = "projected.ls", openMP.core.num = 4)
#'
#' # Get the ancestry matrix
#' Q.matrix <- qmatrix(obj, K = 5)
#'
#' # Plot the spatial interpolation of the ancestry matrix
#' plot(Q.matrix, data.at$coord, method = "map.max",
#'      resolution = c(400,400),
#'      interpolation.model = FieldsKrigModel(10), cex = .4,
#'      xlab = "Longitude", ylab= "Latitude", main = "Ancestry coefficients")
#' @export
plot.tess3Q <- function(x, coord,
                        method = "map.max",
                        resolution = c(300,300), window = NULL,
                        background = TRUE, map.polygon = NULL,
                        raster.filename = NULL,
                        interpolation.model = FieldsKrigModel(10),
                        col.palette = CreatePalette(), ...) {

  CheckPlotParam(x, coord,
                 method,
                 resolution, window,
                 background, map.polygon,
                 raster.filename,
                 interpolation.model,
                 col.palette)


  if (method == "piechart") {
    stop("In development.")
    # PlotPiechartAncestryCoef(x, coord, window, background, col, ...)
  } else {
    interpol.stack <- ComputeInterpolStack(x, coord,
                                           interpolation.model, resolution, window,
                                           background, map.polygon, raster.filename)
    list.grid.z <- RasterStackToMatrix(interpol.stack)
    if (method == "map.max") {

      PlotInterpotationMax(coord, list.grid.z,
                           grid.x = raster::xFromCol(interpol.stack),
                           grid.y = rev(raster::yFromRow(interpol.stack)),
                           col.palette, map = TRUE, ...)
    } else if (method == "map.all") {

      PlotInterpotationAll(coord, list.grid.z,
                           grid.x = raster::xFromCol(interpol.stack),
                           grid.y = rev(raster::yFromRow(interpol.stack)),
                           col.palette, map = TRUE, ...)
    }
  }
}


#' Displays geographic maps for ancestry coefficients.
#' @title Displays geographic maps for ancestry coefficients
#' @author Kevin Caye, Flora Jay, Olivier François
#'
#' @param Q An object of class \code{tess3Q}. The ancestry coefficient matrix.
#' @param coord The numeric matrix of size \eqn{n \times 2} where \eqn{n} is the
#' number of individuals.
#' @param resolution An integer vector of the resolution of the grid used to
#' computed the interpolating surface
#' @param window The window size, such that \code{window = c(xmin, xmax, ymin, ymax)}
#' contains the window mina nd max coordinates.
#' @param background If TRUE the raster file is used as a stencil to render only
#' raster pixel on earth.
#' @param map.polygon The \code{sp::SpatialPolygonsDataFrame} object used to crop the interpolating surfaces.
#' If \code{NULL}, the function \code{\link[rworldmap]{getMap}} is used.
#' @param raster.filename The raster file name used to compute the background stencil.
#' This is an alternative method to crop the interpolating surfaces. The default method uses \code{map.polygon}.
#' @param col.palette A list of color palettes. You can use the function \code{\link{CreatePalette}}.
#' @param interpolation.model The interpolation model used to compute compute
#' the interpolating surface. You can use functions \code{\link{FieldsTpsModel}},
#' \code{\link{FieldsKrigModel}}.
#'
#' @return None
#' @export
#'
#' @examples
#' library(tess3r)
#'
#' # Retrieve a dataset
#' data(data.at)
#'
#' # Run of TESS3
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5,
#'                  ploidy = 1, method = "projected.ls", openMP.core.num = 4)
#'
#' # Get the ancestry matrix
#' Q.matrix <- qmatrix(obj, K = 5)
#'
#' # Plot the spatial interpolation of the ancestry matrix
#' ggtess3Q(Q.matrix, data.at$coord)
ggtess3Q <- function(Q, coord, resolution = c(300,300), window = NULL,
                     background = TRUE, map.polygon = NULL,
                     raster.filename = NULL,
                     interpolation.model = FieldsKrigModel(10),
                     col.palette = CreatePalette()) {


  CheckPlotParam(Q, coord,
                 method = "map.max",
                 resolution, window,
                 background, map.polygon,
                 raster.filename,
                 interpolation.model,
                 col.palette)

  interpol.stack <- ComputeInterpolStack(Q, coord,
                                         interpolation.model, resolution, window,
                                         background, map.polygon, raster.filename)


  toplot <- data.frame(raster::rasterToPoints(interpol.stack))
  ## compute breaks
  col.breaks <- apply(toplot[1:ncol(Q) + 2], 2,
                      function(c) seq(min(c),
                                      max(c),
                                      length.out = length(col.palette[[1]]) + 1))
  ## compute color for each tile
  color <- function(coef, col.palette, col.breaks) {
    max.i <- which.max(coef)
    c <- max(which(col.breaks[,max.i] - as.numeric(coef[max.i]) >= 0)[1] - 1,1)
    return(col.palette[[max.i]][c])
    # return(c)
  }
  toplot$color <- apply(toplot[1:ncol(Q) + 2], 1,
                        function(r) color(r, col.palette, col.breaks))

  mappl <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = toplot, ggplot2::aes(x = x, y = y, fill = color)) +
    ggplot2::scale_fill_identity()

  return(mappl)
}


