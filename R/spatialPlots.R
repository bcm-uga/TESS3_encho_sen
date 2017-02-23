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


SanitizeRaster <- function(raster) {
  raster[raster >= 0.0] <- 1
  raster[raster < 0.0] <- NA
  return(raster)
}


MakeRasterGrid <- function(window, resolution) {
  raster::raster(raster::extent(window), ncol = resolution[1], nrow = resolution[2], vals = 1)
}


InterpolRaster <- function(coord, Q, raster.grid, interpolation.model) {
  interpol.stack <- raster::stack()
  for (j in seq_along(Q[1,])) {
    model <- interpolation.model(coord, Q[,j])
    interpol.stack <- raster::stack(interpol.stack, raster::interpolate(raster.grid, model))
  }
  return(interpol.stack)
}


Crop <- function(interpol.stack, map.polygons, raster.filename) {
  if (!is.null(map.polygons)) {
    return(raster::mask(interpol.stack, map.polygons))
  } else if (!is.null(raster.filename)) {
    imported.raster <- raster::raster(raster.filename)
    imported.raster <- raster::resample(imported.raster, interpol.stack)
    imported.raster <- SanitizeRaster(imported.raster)
    return(interpol.stack * imported.raster)
  } else {
    return(interpol.stack)
  }
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
#' @param col.palette List of color palette.
#' @param map If true map function of maps package is call to plot polygon from
#' map database.
#' @param ... TODOC
PlotInterpotationMax <- function(coord, list.grid.z, grid.x, grid.y, col.palette, map,...) {

  ## UGLY function !!!!

  # rmk : bag data structure for list.grid.z ...

  # which is max
  for (i in 1:length(grid.x)) {
    for (j in 1:length(grid.y)) {
      aux <- sapply(1:length(list.grid.z), function(k) list.grid.z[[k]][i, j])
      k.max <- which.max(aux)
      if (length(k.max) != 0) {
        for (k in 1:length(list.grid.z)) {
          if (k != k.max) {
            list.grid.z[[k]][i,j] <- NA
          }
        }
      }
    }
  }

  for (k in 1:length(list.grid.z)) {
    image(grid.x, grid.y,
          list.grid.z[[k]],
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
#' @param col.palette List of color palette.
#' @param map If true map function of maps package is call to plot polygon from
#' map database.
#' @param ... TODOC
#'
PlotInterpotationAll <- function(coord, list.grid.z, grid.x, grid.y, col.palette, map,...) {


  for (k in 1:length(list.grid.z)) {
    image(grid.x, grid.y,
          list.grid.z[[k]],
          col = col.palette[[k]],
          ...)
    points(coord, pch = 19, ...)
    if (map) {
      TestRequiredPkg("maps")
      requireNamespace("maps")
      message("This function required to attach maps namespace.")
      maps::map(add = TRUE, interior = FALSE)
    }
  }
}



###################################################
##################model functions##################
###################################################


#' Return an interpolation function. Return a wrapper function of the function Tps of
#' the package fields. See \code{\link[fields]{Tps}}.
#'
#' @param ... Parameters of \code{\link[fields]{Tps}}.
#'
#' @export
FieldsTpsModel <- function(...){
  TestRequiredPkg("fields")

  return(function(coord, z) {
    return(fields::Tps(coord, z, ...))
  })
}

#' Return an interpolation function. Return a wrapper function of the function Krig of
#' the package fields. See \code{\link[fields]{Krig}}.
#'
#' @param theta Numeric.
#' @param ... Parameters of \code{\link[fields]{Krig}}.
#'
#' @export
FieldsKrigModel <- function(theta = 10, ...){
  TestRequiredPkg("fields")

  return(function(coord, z) {
    return(fields::Krig(coord, z, theta = theta, ...))
  })
}
