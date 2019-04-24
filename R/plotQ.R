###################################################
#######################helpers#####################
###################################################

CheckPlotParam <- function(x, coord,
                           method,
                           resolution, window,
                           background, map.polygon,
                           raster.filename,
                           interpolation.model,
                           col.palette,
                           add.map) {
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
#' @return an object of class \code{list} containing a list of color palettes.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{barplot.tess3Q}}
#' @examples
#' library(tess3r)
#'
#' ## Load A. thaliana example
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
#' This function displays a barplot representation of the
#' ancestry coefficient matrix. It includes a sort-by-Q option.
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
#' @return A permutation of individual labels used in the sort.by.Q option (order). Displays the
#' Q matrix.
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
#' # Display a barplot for the Q matrix
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
    return(list(order = 1:nrow(Q)))
      }
}


#' This function displays interpolated values of ancestry coefficients on geographic maps.
#' @title Display geographic maps of ancestry coefficients
#' @author Kevin Caye, Flora Jay, Olivier François
#'
#' @param coord a numeric matrix of dimension \eqn{n} times 2 where \eqn{n} is the
#' number of individuals. The matrix must contain (Longitude, Latitude) coordinates for all individuals.
#' @param resolution an integer vector \code{resolution = c(rx,ry)} for the resolution of the grid used to
#' computed the interpolating surface. \code{rx} and \code{ry} are resolution numbers for the x-axis and y-axis respectively.
#' @param window an integer vector for the plotting window, such that \code{window = c(xmin, xmax, ymin, ymax)}
#' contains the window's min and max coordinates.
#' @param background if \code{TRUE} a raster file is used as a stencil to render only raster pixel on earth.
#' @param map.polygon an object of class \code{sp::SpatialPolygonsDataFrame} used to crop the interpolating surface.
#' If \code{NULL}, the function \code{\link[rworldmap]{getMap}} is used.
#' @param raster.filename a raster file name used to compute the background stencil.
#' This is an alternative method to crop the interpolating surfaces. The default method uses \code{map.polygon}.
#' @param col.palette a list of color palettes. Color palettes can be defined by using the function \code{\link{CreatePalette}}.
#' @param interpolation.model an interpolation model used to compute the interpolating surface. Interpolation models can use the
#' functions \code{\link{FieldsTpsModel}} or \code{\link{FieldsKrigModel}}.
#' @param x an object of class \code{tess3Q} containing an ancestry coefficient matrix computed from \code{\link{tess3}}.
#' @param method a character string \code{"map.max"} or \code{"map.all"}. If \code{"map.all"}, interpolating surfaces are displayed
#' for the ancestry coefficients of all populations. If \code{"map.max"} the union of interpolating surfaces is displayed.
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
                        col.palette = CreatePalette(),
                        add.map = TRUE,
                        ...) {

  CheckPlotParam(x, coord,
                 method,
                 resolution, window,
                 background, map.polygon,
                 raster.filename,
                 interpolation.model,
                 col.palette,
                 add.map)


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
                           col.palette, map = add.map,...)
    } else if (method == "map.all") {

      PlotInterpotationAll(coord, list.grid.z,
                           grid.x = raster::xFromCol(interpol.stack),
                           grid.y = rev(raster::yFromRow(interpol.stack)),
                           col.palette, map = add.map, ...)
    }
  }
}


#' This function displays geographic maps of ancestry coefficients using the ggplot syntax.
#' @title Display geographic maps of ancestry coefficients using the ggplot grammar
#' @author Kevin Caye, Flora Jay, Olivier François
#'
#' @param Q an object of class \code{tess3Q} corresponding to an ancestry coefficient matrix obtained \code{tess3}.
#' @param coord a numeric matrix of dimension \eqn{n} times 2 where \eqn{n} is the
#' number of individuals. The matrix must contain (Longitude, Latitude) coordinates for all individuals.
#' @param resolution an integer vector \code{resolution = c(rx,ry)} for the resolution of the grid used to
#' computed the interpolating surface. \code{rx} and \code{ry} are resolution numbers for the x-axis and y-axis respectively.
#' @param window an integer vector for the plotting window, such that \code{window = c(xmin, xmax, ymin, ymax)}
#' contains the window's min and max coordinates.
#' @param background if \code{TRUE} a raster file is used as a stencil to render only raster pixel on earth.
#' @param map.polygon an object of class \code{sp::SpatialPolygonsDataFrame} used to crop the interpolating surface.
#' If \code{NULL}, the function \code{\link[rworldmap]{getMap}} is used.
#' @param raster.filename a raster file name used to compute the background stencil.
#' This is an alternative method to crop interpolating surfaces. The default method uses \code{map.polygon}.
#' @param col.palette a list of color palettes. Color palettes can be defined by using the function \code{\link{CreatePalette}}.
#' @param interpolation.model an interpolation model used to compute the interpolating surface. Interpolation models can use the
#' functions \code{\link{FieldsTpsModel}} or \code{\link{FieldsKrigModel}}.
#' @return None
#' @export
#'
#' @examples
#' library(tess3r)
#'
#' # Load Arabidopsis data
#' data(data.at)
#'
#' # Run of TESS3
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5,
#'                  ploidy = 1, method = "projected.ls", openMP.core.num = 4)
#'
#' # Get the ancestry matrix
#' Q.matrix <- qmatrix(obj, K = 5)
#'
#' # Show a spatial interpolation of the ancestry matrix
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


