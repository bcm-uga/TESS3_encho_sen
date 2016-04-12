###################################################
#######################helpers#####################
###################################################

#' Title
#'
#' @param coord
#'
#' @return
#'
#' @examples
ComputeWindow <- function(coord) {
  sd.x <- sd(coord[,1])
  sd.y <- sd(coord[,2])
  window <- c(min(coord[,1]) - 0.05 * sd.x, max(coord[,1]) + 0.05 * sd.x, min(coord[,2]) - 0.05 * sd.y, max(coord[,2]) + 0.05 * sd.y)
  return(window)
}


###################################################
##################pie chart #######################
###################################################

#' Title
#'
#' @param Q
#' @param coord
#' @param window
#'
#' @return
#' @export
#'
#' @examples
PlotPiechartAncestryCoef <- function(Q, coord, window, background, col, ...) {
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
                     radius = sqrt((window[2] - window[1]) ^ 2 + (window[4] - window[3]) ^ 2) * 0.01,
                     col = col,
                     labels = "")
  }
}


###################################################
##################interpolated maps ###############
###################################################

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
  points(coord, pch = 19)
  if (map) {
    TestRequiredPkg("maps")
    require("maps")
    message("This function required to attach maps namespace.")
    maps::map(add = TRUE, interior = FALSE)
  }
}

PlotInterpotationAll <- function(coord, list.grid.z, grid.x, grid.y, background, col.palette, map,...) {
  for (k in 1:length(list.grid.z)) {
    image(grid.x, grid.y,
          list.grid.z[[k]] * background,
          col = col.palette[[k]],
          ...)
    points(coord, pch = 19)
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

#' TODO
#'
#'
#' @export
#' @param idp
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

#' TODO
#'
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

#' TODO
#'
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

