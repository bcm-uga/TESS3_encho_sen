###################################################
#######################helpers#####################
###################################################

#' Compute the graph window.
#'
#' @param coord Coordinate matrix.
ComputeWindow <- function(coord) {
  sd.x <- sd(coord[,1])
  sd.y <- sd(coord[,2])
  window <- c(min(coord[,1]) - 0.4 * sd.x, max(coord[,1]) + 0.4 * sd.x, min(coord[,2]) - 0.4 * sd.y,
              max(coord[,2]) + 0.4 * sd.y)
  return(window)
}


SanitizeRaster <- function(raster,threshold=0.0) {
  raster[ raster < threshold] <- NA
  raster[ raster >= threshold] <- 1  # line order matters
  return(raster)
}

SanitizeRasterCalc <- function(raster,threshold=0.0) {
  fun <- function(x) {
    x[x < threshold] <- NA
    x[x >= threshold] <- 1  # line order matters
    return(x)
  }
  raster <- raster::calc(raster, fun) #~4 times slower when not using calc
  return(raster)
}


###################################################
######### pie chart internal function #############
###################################################

#' Plot pie chart of ancestry coefficient.
#' @title Plot pie chart of ancestry coefficient.
#' @description  This function is used by \code{\link{piechartQ}}. It can be used
#'   directly but it is advisable to use \code{\link{piechartQ}} instead.
#' @param Q Ancestry coefficient matrix.
#' @param coord Coordinate matrix.
#' @param window Window for the plot corresponding to
#'   c(min_xlim,max_xlim,min_ylim,max_ylim)
#' @param map If \code{TRUE} plot the continental outlines.
#' @param col Vector of ncol(Q) colors
#' @param radius A single value defining the radius = percentage of the total window occupied by the largest pie
#' @param radius.prop . Vector of nrow(Q) integers giving the number of sampled
#'   individuals for each group. Each pie area will be scaled to this number.
#' @param add.pie If TRUE pie chart are added to current plot, if FALSE a new
#'   plot is created.
#' @param label.distx For pie label location: distance to pie center on the
#'   xaxis.
#' @param label.disty For pie label location: distance to the pie top edge on
#'   the yaxis.
#' @param leg.bubble.args A list of arguments to pass to
#'   \code{\link[mapplots]{legend.bubble}}. Argument names must
#'   match. If \code{leg.bubble.args = NULL} no legend is shown.
#' @param legend Boolean. If \code{legend=FALSE} no legend is shown.
#' @param ... Additional parameters will be passed to \code{plot} and
#'   \code{\link[mapplots]{add.pie}}.
#'
PlotPiechartAncestryCoef <- function(Q, coord, window, col, map=T,
                                     radius = 0.01,
                                     radius.prop = NULL,
                                     add.pie = FALSE, names.pie = NULL,
                                     label.distx = 0, label.disty = (window[4] - window[3])*0.02 ,
                                     leg.bubble.args = list(x="topleft",y=NULL) , legend=T,
                                     ...) {

  if (is.null(radius)) {
    radius <- sqrt((window[2] - window[1]) ^ 2 + (window[4] - window[3]) ^ 2) * 0.01
  } else {
    radius <- sqrt((window[2] - window[1]) ^ 2 + (window[4] - window[3]) ^ 2) * radius
  }

  TestRequiredPkg("mapplots")
  if (!add.pie) {
    #plot(coord[(coord[,1] >= window[1]) & (coord[,1] <= window[2]) & (coord[,2] >= window[3]) & (coord[,2] <= window[4]),(window[2] - window[1])],type = "n", ...)
    plot(window[1:2],window[3:4], type = "n", ...)

    if (map) {
      TestRequiredPkg("maps")
      require("maps")
      message("This function required to attach maps namespace.")
      maps::map(add = TRUE, col = "grey90", fill = TRUE)
    }
  }

  isInWindow <- unlist(lapply(1:nrow(coord),
                              function(i) (coord[i,1] >= window[1] && coord[i,1] <= window[2] && coord[i,2] >= window[3] && coord[i,2] <= window[4]) ))

  if (is.null(radius.prop)) {
    scaling <- rep(1, nrow(Q))
  } else {
    scaling <- sqrt(radius.prop)/sqrt(max(radius.prop[isInWindow],na.rm=T))
  }
  radius <- radius * scaling


  if (is.null(names.pie) || names.pie=="") names.pie=rep("",nrow(Q))
  if (length(names.pie) != nrow(Q) ) stop("Argument names.pie should be NULL or \"\" or have the same length as rows in Q  (number of individuals or number of populations)")

    for (i in which(isInWindow) ) {
      mapplots::add.pie(z = Q[i,],
                     x = coord[i,1],
                     y = coord[i,2],
                     radius = radius[i],
                     col = col, labels=NA, ...)

      text(coord[i,1]+label.distx,coord[i,2]+radius[i]+label.disty,labels=names.pie[i],...) #+radius[i]
    }

    if (!is.null(radius.prop) && !is.null(leg.bubble.args) && legend) {
      leg.bubble.args$z = max(radius.prop[isInWindow])
      leg.bubble.args$maxradius = max(radius[isInWindow])
      do.call( mapplots::legend.bubble, leg.bubble.args)
    }

}


###################################################
########### internal interpolated maps ############
###################################################

#' Plot map with the maximum values of the interpolation surfaces.
#' @description  This function is used by \code{\link{plot.tess3Q}}. It can be
#'   used directly but it is advisable to use \code{\link{plot.tess3Q}} instead.
#' @param coord Coordinate matrix.
#' @param list.grid.z List of interpolation surface matrices.
#' @param grid.x x-coordinates of grid cell centers.
#' @param grid.y y-coordinates of grid cell centers.
#' @param background Boolean marix.
#' @param col.palette List of color palette.
#' @param map If true map function of maps package is call to plot polygon from
#'   map database.
#' @param legend Boolean. Plot the color gradient legend?
#' @param horizontal Boolean. Horizontal position for legend key?
#' @param graphics.reset If \code{FALSE} the plotting parameters will not be
#'   reset and one can add more information onto the image plot. Use false if
#'   you want to display piecharts on top of the map.
#' @param layout.nkeys Integer. Number of color keys by column (resp. row) in
#'   legend when \code{horizontal = FALSE} (resp. \code{TRUE})
#' @param legend.width Width of each column (rep row) for the vertical (resp
#'   horizontal) legend (default is 1). The total size allocated to the whole
#'   plot (map+legend) is set to 10.
#' @param leg.extra.args List of extra arguments given to \code{\link{image.plot.legend}}
#' @param ... Extra arguments given to \code{\link[graphics]{image}} and
#'   \code{\link[graphics]{points}}.
#'
#'  @examples
#'  ## DO NOT RUN
#'  PlotInterpotationMax(coord, list.grid.z, grid.x, grid.y, background = T, col.palette, map = T,
#'    legend = T, horizontal = FALSE, graphics.reset = TRUE, layout.nkeys = 2, legend.width = 1.2,
#'    leg.extra.args = list(legend.lab="Ancestry Coef.",legend.cex=1.5))
#'
PlotInterpotationMax <- function(coord, list.grid.z, grid.x, grid.y, background, col.palette, map,
                                 legend, horizontal = FALSE, graphics.reset = TRUE,
                                 layout.nkeys = 3, legend.width = 1, leg.extra.args = list(), ...) {

  # rmk : bag data structure for list.grid.z ...

  # which is max
  for (i in 1:length(grid.x)) {
    for (j in 1:length(grid.y)) {
      aux <- sapply(1:length(list.grid.z), function(k) list.grid.z[[k]][i, j])
      k.max <- which.max(aux)
      for (k in 1:length(list.grid.z)) {
        if (k != k.max) {
          list.grid.z[[k]][i,j] <- NA
        }
      }
    }
  }

  if (legend) {
    layout.width <-10
    old.par <- par(no.readonly = TRUE)
    K <- length(list.grid.z)
    layout.nrow <- layout.nkeys
    layout.ncol <- ceiling(K/layout.nrow) + 1
    if (!horizontal) {
      layout(matrix(c(rep(1,layout.nrow),2:((layout.ncol-1)*layout.nrow+1)),ncol=layout.ncol),
             width=c(layout.width - legend.width*(layout.ncol-1),rep(legend.width,layout.ncol-1)))
    } else {
      layout(t(matrix(c(rep(1,layout.nrow),2:((layout.ncol-1)*layout.nrow+1)),ncol=layout.ncol)),
             height=c(layout.width - legend.width*(layout.ncol-1),rep(legend.width,layout.ncol-1)))
    }
  }
  for (k in 1:length(list.grid.z)) {
    plotting <- try(image(grid.x, grid.y,
          list.grid.z[[k]] * background,
          col = col.palette[[k]],
          add = (k > 1),
          ...)
    )
    if (class(plotting) == "try-error") {
      par(old.par)
      stop("Error: \n", plotting, "\nSetting graphics parameters back to previous.
           If error persists try to enlarge the plot area, or reset graphic device by calling dev.off() until all X devices are shut down
           [warning: this will delete all your current plots].")
    }
  }
  map.par <- par(no.readonly = TRUE)
  points(coord, pch = 19, ...)
  if (map) {
    TestRequiredPkg("maps")
    require("maps")
    message("This function required to attach maps namespace.")
    maps::map(add = TRUE, interior = FALSE)
  }

  # Potting gradient legend keys
  if (legend) {
    for (k in 1:length(list.grid.z)) {
      plotting <- try( do.call(image.plot.legend, c(list(x = grid.x, y = grid.y, z = list.grid.z[[k]] * background,
                                                       col = col.palette[[k]], add = F,
                                                       horizontal = horizontal, main.par = old.par),
                                                    leg.extra.args))
      )
      if (class(plotting) == "try-error") {
        par(old.par)
        stop("Error, cannot display legend: \n", plotting, "\nSetting graphics parameters back to previous.
           If error persists try to enlarge the plot area, or reset graphic device by calling dev.off() until all X devices are shut down
           [warning: this will delete all your current plots].")
      }
    }
    if (!graphics.reset) {
      par(mfg=c(1,1))
      par(plt=map.par$plt)
      par(usr=map.par$usr)
      set.par <- par(no.readonly = TRUE)
    } else {
      par(old.par)
    }
  }
}


#' Plot all the interpolation surfaces on several graphs.
#' @description  This function is used by \code{\link{plot.tess3Q}}. It can be
#'   used directly but it is advisable to use \code{\link{plot.tess3Q}} instead.
#' @param coord Coordinate matrix.
#' @param list.grid.z List of interpolation surface matrices.
#' @param grid.x x-coordinates of grid cell centers.
#' @param grid.y y-coordinates of grid cell centers.
#' @param background Boolean marix.
#' @param col.palette List of color palette.
#' @param map If true map function of maps package is call to plot polygon from
#'   map database.
#' @param legend Boolean. Plot the color gradient legend?
#' @param horizontal Boolean. Horizontal position for legend key?
#' @param graphics.reset If \code{FALSE} the plotting parameters will not be
#'   reset and one can add more information onto the image plot. Use FALSE if
#'   you want to display piecharts on top of the map.
#' @param legend.width Width in characters of the legend strip.
#' @param leg.extra.args List of extra arguments given to \code{\link[fields]{image.plot}}
#' @param ... Extra arguments given to \code{\link[graphics]{image}},
#'   \code{\link[fields]{image.plot}} and \code{\link[graphics]{points}}.
#'
#'  @examples
#'  ## DO NOT RUN
#'  PlotInterpotationAll(coord, list.grid.z, grid.x, grid.y, background = T, col.palette, map = T,
#'    legend = T, horizontal = FALSE, graphics.reset = TRUE, legend.width = 1,
#'    leg.extra.args = list(legend.lab="Ancestry Coef.",legend.cex=1.5))

#'
PlotInterpotationAll <- function(coord, list.grid.z, grid.x, grid.y, background, col.palette, map, legend,
                                 horizontal = FALSE, graphics.reset = TRUE, legend.width = 1, leg.extra.args = list(), ...) {

  if (legend) {
    TestRequiredPkg("fields")
    require("fields")
  }
  old.par <- par(no.readonly = TRUE)
  for (k in 1:length(list.grid.z)) {
      if (legend) {
         do.call( fields::image.plot, c( list( x = grid.x, y = grid.y, z = list.grid.z[[k]] * background,
                                               col = col.palette[[k]], legend.width = legend.width, horizontal = horizontal),
                                         leg.extra.args, list(...)) )
      } else {
        image(grid.x, grid.y,
          list.grid.z[[k]] * background,
          col = col.palette[[k]],
          ...)
      }
    points(coord, pch = 19, ...)
    if (map) {
      TestRequiredPkg("maps")
      requireNamespace("maps")
      message("This function required to attach maps namespace.")
      maps::map(add = TRUE, interior = FALSE)
    }
  }
  if (graphics.reset) par(old.par)

}


#' Compute a grid for interpolation
#' @description For internal usage, this function is called by \code{\link{plot.tess3Q}}.
#' @param resolution An integer vector of the resolution of the grid used to
#'   computed the interpolating surface.
#' @param window Vector defining boundaries of the map eg
#'   c(long_min,long_max,lat_min,lat_max)
#' @param background if TRUE compute a background stencil.
#' @param raster.filename Name of a raster file.
#' @param threshold A \code{numeric}. If \code{background=TRUE}, interpolation will be performed only at cells having values above \code{threshold}
#'
#' @return a grid of the area on which ancestry coefficients will be interpolated
#'
#' @export
#'
ComputeGridAndBackground <- function(window, resolution, background, raster.filename,threshold=0.0) {
  TestRequiredPkg("raster")
  if (background) {
    imported.raster <- raster::raster(raster.filename)
    imported.raster <- raster::resample(imported.raster,raster::raster(raster::extent(window), ncol = resolution[1], nrow = resolution[2]))
    imported.raster <- SanitizeRasterCalc(imported.raster,threshold) #FJ: moving to third step because matrix is smaller after resampling

  } else {
    imported.raster <- raster::raster(raster::extent(window), ncol = resolution[1], nrow = resolution[2], vals = 1)
  }
  res <- list(grid.x = raster::xFromCol(imported.raster),
              grid.y = rev(raster::yFromRow(imported.raster)),
              background = raster::as.matrix(imported.raster))
  res$background <- t(apply(t(res$background), 1, rev)) # transpose and rev col
  return(res)
}

#' Return an interpolation function. Return a wrapper function of the function idw of
#' the package gstat. See \code{\link[gstat]{idw}}
#'
#'
#' @export
#'
#' @param idp numeric; specify the inverse distance weighting power.
idw <- function(idp=1.0){
  TestRequiredPkg("gstat")
  TestRequiredPkg("sp")

  return(function(z, coord, grid.x, grid.y) {

    # grid
    grd <- expand.grid(X = grid.x, Y = grid.y)
    sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
    sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object

    res <- list()
    for (k in 1:ncol(z)) {
      # data to interpolate
      dat <- data.frame(X = coord[,1], Y = coord[,2], Z = z[,k])
      # Force data frame object into a SpatialPointsDataFrame object
      sp::coordinates(dat) <- c("X","Y")
      res[[k]] <- t(apply(matrix(gstat::idw(Z ~ 1, dat, newdata = grd, idp = idp)$var1.pred, length(grid.x), length(grid.y)), 1, rev))
    }
    return(res)
  })
}

#' Return an interpolation function. Return a wrapper function of the function autoKrige of
#' the package automap. See \code{\link[automap]{autoKrige}}
#'
#'
#' @param formula Formula that defines the dependent variable as a linear model of independent variables.
#'
#' @export
universalkriging <- function(formula = as.formula(Z ~ X + Y)){
  TestRequiredPkg("automap")
  TestRequiredPkg("sp")

  return(function(z, coord, grid.x, grid.y) {
    # grid
    grd <- expand.grid(X = grid.x, Y = grid.y)
    sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
    sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object

    res <- list()
    for (k in 1:ncol(z)) {
      # data to interpolate
      dat <- data.frame(X = coord[,1], Y = coord[,2], Z = z[,k])
      # Force data frame object into a SpatialPointsDataFrame object
      sp::coordinates(dat) <- c("X","Y")
      res[[k]] <- t(apply(matrix(automap::autoKrige(formula, dat, grd)$krige_output$var1.pred, length(grid.x), length(grid.y)), 1, rev))
    }
    return(res)
  })
}

#' Return an interpolation function. Return a wrapper function of the function Krig of
#' the package fields. See \code{\link[fields]{Krig}}.
#'
#' @param theta Numeric.
#'
#' @export
kriging <- function(theta = 10){
  TestRequiredPkg("fields")

  return(function(z, coord, grid.x, grid.y) {
    # grid
    grd <- fields::make.surface.grid(list(grid.x, grid.y))

    res <- list()
    for (k in 1:ncol(z)) {
      clust <- fields::Krig(coord, z[,k], theta = theta)
      look <- predict(clust, grd) # evaluate on a grid of points
      out <- fields::as.surface(grd, look)
      res[[k]] = out[[8]]
    }
    return(res)
  })
}
