###################################################
################# Methods #########################
###################################################

###################################################
################# Barplot #########################
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
#' @param palette.length An integer value in [3,9] for the number of colors in
#'   each element of the palette list. Not used if a palette is provided with
#'   \code{col.palette}.
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

  if (!inherits(Q,"tess3Q")) {warning("Object Q not of class tess3Q.")}
  ## defines a default color palette ( colors)
  if (is.null(col.palette)) {
    message("You can use CreatePalette() to define new color palettes.\n")
    if (ncol(Q) > 8)  stop("The default color palette expects less than 8 clusters.")
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
  if (length(col.palette) < ncol(Q))  stop("col.palette must of length at least ncol(Q)")
  # Pick one color for each cluster (in the gradient range)
  # Works if col.palette is a list of vectors (pick one color for each first ncol(Q) vectors of the list)
  # Works if col.palette is a vector (pick the first ncol(Q) elements)
  colpal <- unlist(lapply(col.palette[1:ncol(Q)], FUN=function(gradient) gradient[ ceiling(length(gradient)*2/3) ] ))


  #if (palette.length > 3){
  #colpal = sapply(col.palette, FUN = function(x) x[palette.length/2])}
  #else {
  #colpal = sapply(col.palette, FUN = function(x) x[palette.length])
  #}

  if (sort.by.Q){
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

###################################################
############ Interpolation map ####################
###################################################

#' Displays geographic maps for ancestry coefficients (interpolation).
#' @title Displays geographic maps of ancestry coefficients.
#' @description Ancestry coefficients are interpolated on a map and color
#'   gradients indicate their value
#' @author Kevin Caye, Flora Jay, Olivier François
#'
#' @param x an object of class \code{tess3Q}. The ancestry coefficient matrix.
#' @param col.palette a list of color palettes. If \code{NULL}, a default list
#'   with 8 color palettes is used.
#' @param coord The numeric matrix of size \eqn{n \times 2} where \eqn{n} is the
#'   number of individuals.
#' @param method \code{"map.max"} or \code{"map.all"}. If \code{"map.all"}, a
#'   interpolating surface of the ancestry coefficient is plotted for each
#'   ancestral population. If \code{"map.max"} only the maximum union of the
#'   interpolating surfaces is plotted.
#' @param resolution An integer vector of the resolution of the grid used to
#'   computed the interpolating surface.
#' @param window Vector defining boundaries of the map eg
#'   c(long_min,long_max,lat_min,lat_max)
#' @param background if TRUE compute a background stencil.
#' @param raster.filename Name of a raster file.
#' @param interpolation.function Function used to interpolate ancestry
#'   coefficients on map. Eg. \code{kriging(5), idw(1.0) }.
#' @param map If \code{TRUE} plot the continental outlines.
#' @param grid Use a previously saved grid. This saves computational time if you
#'   plot several times similar maps.
#' @param palette.length An integer value in [3,9] for the number of colors in
#'   each element of the palette list. Not used if a palette is provided with
#'   \code{col.palette}.
#' @param ... Extra arguments given to \code{PlotInterpotationMax}, \code{PlotInterpotationAll}, \code{\link[graphics]{image}} and \code{\link[graphics]{points}}.

#'
#' @details Details on interpolation methods: \code{\link{kriging}},
#'   \code{\link{idw}}, \code{\link[fields]{Krig}}, \code{\link[gstat]{krige}}
#' @seealso \code{\link[tess3r]{piechart.tess3Q}} to display piecharts on the
#'   map.
#'
#' @return invisible(grid) The newly built grid if it was not already provided as an input.
#'   This grid can be used in a later call of \code{plot.tess3Q} unless you want
#'   to change the resolution or the input raster file. Note that the grid is invisible unless stored in an object.
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
#' savedgrid <- plot(Q.matrix, data.at$coord, method = "map.max",
#'      resolution = c(400,400),
#'      interpol = kriging(10), cex = .4,
#'      xlab = "Longitude", ylab= "Latitude", main = "Ancestry coefficients")
#'
#'  # Unless you want to change the resolution, the input raster file or the argument \code{background},
#'  # you can re-use the saved grid for all your next plots, eg:
#'  plot(Q.matrix, data.at$coord, method = "map.max",
#'      grid = savedgrid,
#'      interpol = idw(), cex = .4, col.main="red",
#'      xlab = "Longitude", ylab= "Latitude", main = "My new title")
#'
#'  # Display legend
#'  plot(Q.matrix, data.at$coord, method = "map.max",grid=savedgrid, legend=T)
#'
#'  plot(tess3.res$Q, data.at$coord, method = "map.max",grid=savedgrid, legend=T,
#'        horizontal=T,legend.lab="ancestry coef",legend.cex=.75)
#'
#'  # Display piecharts on top of the map, use  graphics.reset=FALSE
#'  # See ?piechart for more examples
#'  plot(Q.matrix, data.at$coord, method = "map.max", grid=savedgrid, legend=T, graphics.reset=FALSE)
#'  plot(Q.matrix, data.at$coord, method = "piechart", add.pie=T, radius=.006)
#' @export
map.tess3Q <- function(x, coord, method = "map.max", resolution = c(300,300), window = NULL,
                        background = TRUE, raster.filename = NULL, interpolation.function = kriging(10),
                        col.palette = NULL, map = TRUE, palette.length = 9, legend=FALSE, grid=NULL, ...) {


  if (is.null(coord)) stop("Argument coord mandatory for map")
  if (!inherits(x,"tess3Q")) {warning("Object x not of class tess3Q.")}

  ## col.palette
  if (is.null(col.palette)) {
    message("You can use CreatePalette() to define new color palettes.\n")
    if (ncol(x) > 8)  stop("The default color palette expects less than 9 clusters.")
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
  if (length(col.palette) < ncol(x) & (method == "map.max" | method == "map.all")) {
    stop("col.palette must of length ncol(Q)")
  }


  ## compute of window if window is NULL
  if (is.null(window)) {
    window <- ComputeWindow(coord)
  }

  ##  raster filename
  if (is.null(raster.filename)) {
    raster.filename <- system.file("extdata/raster","earth.tif",package = "tess3r")
  }

  if (!method %in% c("map.max","map.all"))  stop("method=",method, "is not a valid option")

  if (is.null(grid)) {
    message("Computing grid and background...")
    grid <- ComputeGridAndBackground(window, resolution, background, raster.filename)
    is.grid.new <- TRUE
  } else {
    is.grid.new <- FALSE
    message("Grid was given as input - Resolution, rasterfile  and background arguments are NOT used.")
  }

  if (method == "map.max") {
    message("Interpolating maximum ancestry coefficients on a grid...")
    list.grid.z <- interpolation.function(x, coord, grid$grid.x, grid$grid.y)
    message("Plotting the results...")
    PlotInterpotationMax(coord, list.grid.z, grid$grid.x, grid$grid.y, grid$background, col.palette, map, legend = legend, ...)
  } else if (method == "map.all") {
    message("Interpolating all ancestry coefficients on a grid...")
    list.grid.z <- interpolation.function(x, coord, grid$grid.x, grid$grid.y)
    message("Plotting the results...")
    PlotInterpotationAll(coord, list.grid.z, grid$grid.x, grid$grid.y, grid$background, col.palette, map, legend = legend, ...)
  }

    if (is.grid.new) invisible(grid)
}




###################################################
################# Piechart map ###################
###################################################

#'Displays piecharts of ancestry coefficients on geographic maps.
#'@title Displays piecharts of ancestry coefficients on geographic maps.
#'@description Piecharts are displayed at each individual location or at the
#'  average location of all samples belonging to a same group. Piechart color
#'  slices represent ancestry coefficients for all clusters.
#'@author Kevin Caye, Flora Jay, Olivier François
#'
#'@param x an object of class \code{tess3Q}. The ancestry coefficient matrix.
#'@param coord The numeric matrix of size \eqn{n \times 2} where \eqn{n} is the
#'  number of individuals.
#'@param method \code{"piechart"} or \code{"piechart.pop"}. If
#'  \code{"piechart"}, one piechart is plotted for each individual at its
#'  sampled location. If \code{"piechart.pop"} one piechart is plotted for each
#'  population at the location averaged over all individuals belonging to this
#'  population. Individual population labels are required, see \code{pop}.
#'@param window Vector defining boundaries of the map eg \code{c(long_min,
#'  long_max,l at_min, lat_max)}
#'@param col.palette The list of color palettes used for interpolation maps. If
#'  \code{NULL}, a default list with 8 color palettes is used.
#'@param col A vector of colors for piecharts of length ncol(x). By default
#'  colors are chosen to match col.palette.
#'@param map If \code{TRUE} plot the continental outlines.
#'@param pop Required only when using \code{method="piechart.pop"}. Population
#'  labels (for each individual). Vector of \code{n} elements of class
#'  \code{numeric}, \code{character} or \code{factor}
#'@param add.pie If \code{TRUE} pie charts are displayed on top of current plot.
#'  If \code{FALSE} creating a new plot.
#'@param names.pie Labels for pie charts (one per pie). If \code{NULL} and
#'  \code{method="piechart.pop"}, labels are set up automatically from
#'  \code{pop}. If \code{""}: no labels. Else vector of npop elements
#'  (\code{method="piechart.pop"}) or nindiv elements
#'  (\code{method="piechart"}).
#'@param radius A single value defining the radius = percentage of the total
#'  window occupied by the largest pie
#'@param scale If TRUE pie areas will be scaled to the sampling size of each
#'  group (defined by \code{pop}). Only relevant when
#'  \code{method="piechart.pop"}.
#'@param palette.length An integer value in [3,9] for the number of colors in
#'  each element of the palette list. Not used if a palette is provided with
#'  \code{col.palette}.
#'@param ... Additional parameters will be passed to
#'  \code{\link[graphics]{plot}} and \code{\link[mapplots]{add.pie}}.
#'@return None
#'@details \code{piechart(q, coord, method="piechart", ...)} # equivalent to:
#'  \code{plot(q, coord, method="piechart", ...)}
#'@seealso  \code{\link{map.tess3Q}} to plot piecharts on top of interpolation
#'  maps \code{\link[mapplots]{legend.bubble}}
#' @examples
#' library(tess3r)
#' data(data.at)
#' obj <- tess3(data.at$X, coord = data.at$coord,
#'                  K = 5, ploidy = 1, openMP.core.num = 4)
#' Q.matrix <- qmatrix(obj, K = 5)
#' piechart(Q.matrix, data.at$coord, method="piechart")
#'
#' # Display piecharts for groups of individuals (here group=country)
#' # And change the default location of the pie labels with label.distx/y
#' piechart(Q.matrix,data.at$coord, method="piechart.pop",
#'    pop = data.at$countries, scale=F,
#'    window=c(2,20,45,55),
#'    radius = .02, cex =.6,
#'    label.distx=1.3,label.disty=-.4)
#'
#'  # Scale the pie to the group sampling size
#'  piechart(Q.matrix,data.at$coord, method="piechart.pop",
#'    col=topo.colors(ncol(Q.matrix), alpha = 0.6),
#'    pop = substr(data.at$countries,1,2),
#'    window=c(2,20,45,55),
#'    scale=T, radius=.06, cex = 1.2)
#'  text(3,53,"Sampling size",cex=.7)
#'
#'  # Change legend using leg.bubble.par a list of arguments passed to legend.bubble{mapplots}
#'  piechart(Q.matrix,data.at$coord, method="piechart.pop",
#'    pop = substr(data.at$countries,1,2),
#'    window=c(2,20,45,55),
#'    scale=T, radius=.06, cex = 1.2,
#'    leg.bubble.par=list(x="bottomright",bty="n",lwd=4))
#'
#'
#'  # Display piecharts on top of a interpolated ancestry coefficient map
#'  plot(Q.matrix, data.at$coord, method = "map.max",
#'      resolution = c(300,300),
#'      interpol = kriging(), cex = .4,
#'      xlab = "Longitude", ylab= "Latitude", main = "Ancestry coefficients")
#'  piechart(Q.matrix, data.at$coord, method = "piechart.pop",
#'      pop = substr(data.at$countries,1,2),
#'      scale=F, add.pie = T,radius=.015, cex = 1.2)
#'
#'
#'@export
piechart.tess3Q <- function(x, coord, method = "piechart", window = NULL,
                        col.palette = NULL, map = TRUE, palette.length = 9, grid=NULL,
                        add.pie=FALSE, pop=NULL, names.pie=NULL, radius=0.01, scale= FALSE , ...) {

  if (is.null(coord)) stop("Argument coord mandatory for piechart")
  if (!inherits(x,"tess3Q")) {warning("Object x not of class tess3Q.")}
  # Colors
  if (is.null(col.palette)) {
    message("Use CreatePalette() to define new color palettes.\n")
    if (ncol(x) > 8)  stop("The default color palette expects less than 9 clusters.")
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
  if (length(col.palette) < ncol(x))  stop("col.palette must have a length larger than the number of clusters")
  # Pick one color for each cluster (in the gradient range)
  # Works if col.palette is a list of vectors (pick one color for each first ncol(Q) vectors of the list)
  # Works if col.palette is a vector (pick the first ncol(Q) elements)
  col <- unlist(lapply(col.palette[1:ncol(Q)], FUN=function(gradient) gradient[ ceiling(length(gradient)*2/3) ] ))


  ## compute of window if window is NULL
  if (is.null(window)) {
    window <- ComputeWindow(coord)
  }

  if (method == "piechart") {
    PlotPiechartAncestryCoef(x, coord, window, col=col, map=map, add.pie=add.pie, names.pie=names.pie, radius=radius,...)

  } else if (method == "piechart.pop") {
    if (is.null(pop)) stop("pop argument required when using method=\"piechart.pop\" ")
    message("Grouping results by population")
    info.bypop = bypop(x, coord, pop)
    if (is.null(names.pie)) names.pie = info.bypop$labels
    if (scale) { # pie area proportional to group size
      PlotPiechartAncestryCoef(info.bypop$Q, info.bypop$coord, window, col=col, map=map,
                               add.pie=add.pie, names.pie=names.pie, radius=radius, radius.prop = info.bypop$size, ...)
    } else {
      PlotPiechartAncestryCoef(info.bypop$Q, info.bypop$coord, window, col=col, map=map,
                               add.pie=add.pie, names.pie=names.pie, radius=radius, ...)
    }
  }
}



#' @title Plotting ancestry coefficients
#' @description Wrapper for all the different ways of displaying ancestry
#'   coefficients (barplot, piecharts, gradients).
#' @author Kevin Caye, Flora Jay, Olivier François
#'
#' @param Q An object of class \code{tess3Q}. The ancestry coefficient matrix.
#' @param method \itemize{ \item{ \code{"barplot"}. Barplot of ancestry
#'   coefficients. See \code{\link[tess3r]{barplot}}} \item{ \code{"map.max"} or
#'   \code{"map.all"}. Gradients of ancestry coefficients are displayed on a
#'   map. See  \code{\link[tess3r]{map}}} \item{ \code{"piechart"} or
#'   \code{"piechart.pop"}. Piecharts of ancestry coefficients for each
#'   individual or each population are displayed at sampled locations on a map.
#'   See  \code{\link[tess3r]{piechart}} } }
#' @param ... Additionnal parameters. Check corresponding functions
#'   \code{\link[tess3r]{barplot}}, \link[tess3r]{map} and
#'   \code{\link[tess3r]{piechart}} about parameter requirements.
#' @export
plot.tess3Q <- function(Q, coord=NULL, method="map.max", ...){
  if (method %in% c("map.all","map.max")) return(map.tess3Q(Q, coord = coord, method = method,...))
  if (method %in% c("piechart","piechart.pop")) return(piechart.tess3Q(Q, coord = coord, method = method,...))
  if (method %in% c("barplot")) return(barplot.tess3Q(Q, ...))

}

