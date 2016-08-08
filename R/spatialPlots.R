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


###################################################
##################pie chart #######################
###################################################

#' Plot pie chart of ancestry coefficient.
#'
#' @param Q Ancestry coefficient matrix.
#' @param coord Coordinate matrix.
#' @param window Window of the graph.
#' @param background if TRUE compute a background stencil.
#' @param col colors
#' @param radius radius of pies.
#' @param ... TODOC
#'
PlotPiechartAncestryCoef <- function(Q, coord, window, background, col, radius = sqrt((window[2] - window[1]) ^ 2 + (window[4] - window[3]) ^ 2) * 0.01, ...) {
  TestRequiredPkg("mapplots")

  plot(coord[(coord[,1] >= window[1]) & (coord[,1] <= window[2]) & (coord[,2] >= window[3]) & (coord[,2] <= window[4]),],
       type = "n", ...)

  if (background) {
    TestRequiredPkg("maps")
    require("maps")
    message("This function required to attach maps namespace.")
    maps::map(add = TRUE, col = "grey90", fill = TRUE)
  }

  for (i in 1:nrow(Q)) {
   mapplots::add.pie(z = Q[i,],
                     x = coord[i,1],
                     y = coord[i,2],
                     radius = radius,
                     col = col,
                     labels = "")
  }
}


###################################################
##################interpolated maps ###############
###################################################

#' Plot map with the maximum values of the interpolation surfaces.
#'
#' @param coord Coordinate matrix.
#' @param list.grid.z List of interpolation surface matrices.
#' @param grid.x TODOC
#' @param grid.y TODOC
#' @param background Boolean marix.
#' @param col.palette List of color palette.
#' @param map If true map function of maps package is call to plot polygon from
#' map database.
#' @param ... TODOC
PlotInterpotationMax <- function(coord, list.grid.z, grid.x, grid.y, background, col.palette, map,...) {

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

  for (k in 1:length(list.grid.z)) {
    image(grid.x, grid.y,
          list.grid.z[[k]] * background,
          col = col.palette[[k]],
          add = (k > 1),
          ...)
  }
  points(coord, pch = 19, ...)
  if (map) {
    TestRequiredPkg("maps")
    require("maps")
    message("This function required to attach maps namespace.")
    maps::map(add = TRUE, interior = FALSE)
  }
}

#' Plot all the interpolation surfaces on several graphs.
#'
#' @param coord Coordinate matrix.
#' @param list.grid.z List of interpolation surface matrices.
#' @param grid.x TODOC
#' @param grid.y TODOC
#' @param background Boolean marix.
#' @param col.palette List of color palette.
#' @param map If true map function of maps package is call to plot polygon from
#' map database.
#' @param ... TODOC
#'
PlotInterpotationAll <- function(coord, list.grid.z, grid.x, grid.y, background, col.palette, map,...) {
  for (k in 1:length(list.grid.z)) {
    image(grid.x, grid.y,
          list.grid.z[[k]] * background,
          col = col.palette[[k]],
          ...)
    points(coord, pch = 19, ...)
    if (map) {
      TestRequiredPkg("maps")
      require("maps")
      message("This function required to attach maps namespace.")
      maps::map(add = TRUE, interior = FALSE)
    }
  }
}

ComputeGridAndBackground <- function(window, resolution, background, raster.filename) {
  TestRequiredPkg("raster")
  if (background) {
    imported.raster <- raster::raster(raster.filename)
    imported.raster <- raster::resample(imported.raster,raster::raster(raster::extent(window), ncol = resolution[1], nrow = resolution[2]))
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
#' the package gstat.
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
#' the package automap.
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
#' the package fields
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
