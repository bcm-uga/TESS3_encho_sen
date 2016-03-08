#' TODO
#'
#'
TestRequiredPkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(pkg,"needed for this function to work. Please install it."),
         call. = FALSE)
  }
}

#' TODO
#'
#'
idw <- function(idp=1.0){
  TestRequiredPkg("gstat")

  return(function(dat,grd) {
    return(gstat::idw(Z~1,dat,newdata=grd,idp=idp))
  })

}

#' TODO
#'
#'
kriging <- function(formula = as.formula(Z ~ X + Y)){
  TestRequiredPkg("automap")

  return(function(dat,grd) {
    kr = automap::autoKrige(formula, dat, grd)
    return(kr$krige_output)
  })

}

#' TODO
#'
#'
ComputeBackgroundandGrid <- function(dat, raster.filename, resolution) {
  imported.raster <- raster::raster(raster.filename)

  # crop to dat
  imported.raster <- raster::resample(imported.raster,raster::raster(raster::extent(dat), ncol = resolution$ncol, nrow = resolution$nrow))
  imported.raster <- raster::calc(imported.raster, function(x) {ifelse(x>0,1,0)})
  # plot(imported.raster)

  # interpolate
  # Create an empty grid where n is the total number of cells
  grd              <- as.data.frame(raster::rasterToPoints(imported.raster)[,1:2])
  names(grd)       <- c("X", "Y")
  sp::coordinates(grd) <- c("X", "Y")
  sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
  sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object

  return(list(grd = grd, background = imported.raster))

}


#' TODO
#'
#'
RasterAncestryCoef <- function(ancestry.coef, coord, raster.filename = system.file("extdata/raster","earth.tif",package = "TESS3enchoSen"), interpolation.function = idw(), resolution = list(ncol = 300, nrow = 300)) {

  # test installed pkg
  TestRequiredPkg("sp")
  TestRequiredPkg("raster")
  TestRequiredPkg("gstat")
  # TestRequiredPkg("rasterVis")
  # TestRequiredPkg("ggplot2")

  # data to interpolate
  dat <- data.frame(X = coord[,1], Y = coord[,2], Z = ancestry.coef)
  # Force data frame object into a SpatialPointsDataFrame object
  sp::coordinates(dat) <- c("X","Y")

  grid <- ComputeBackgroundandGrid(dat, raster.filename, resolution)

  # Interpolate the surface using a power value of 2 (idp=2.0)
  dat.interpolated <- interpolation.function(dat,grid$grd)

  # mask
  ancestry.raster <- grid$background * raster::raster(dat.interpolated)
  # plot
#   ancestry.plot <- rasterVis::gplot(ancestry.raster) +
#     ggplot2::geom_raster(ggplot2::aes(fill = value)) +
#     ggplot2::scale_fill_gradient(low = 'green', high = 'red', na.value = "white") +
#     ggplot2::coord_equal() +
#     ggplot2::geom_point(data = as.data.frame(dat),ggplot2::aes(x=X,y=Y))

  return(ancestry.raster)
}


#' TODO
#'
#'
RasterAncestryCoefRGB <- function(Q, coord, raster.filename = system.file("extdata/raster","earth.tif",package = "TESS3enchoSen"), interpolation.function = idw(), resolution = list(ncol = 300, nrow = 300)) {

  # compute grid and background
  dat <- data.frame(X = coord[,1], Y = coord[,2], Z = Q[,1])
  # Force data frame object into a SpatialPointsDataFrame object
  sp::coordinates(dat) <- c("X","Y")
  grid <- ComputeBackgroundandGrid(dat, raster.filename, resolution)

  # compute raster for each ancestral population
  ancestry.raster.list = NULL
  for (i in 1:ncol(Q)){
    dat$Z = Q[,i]
    # Interpolate the surface using a power value of 2 (idp=2.0)
    dat.interpolated <- interpolation.function(dat,grid$grd)

    # plot
    ancestry.raster.list <- c(ancestry.raster.list, grid$background * raster::raster(dat.interpolated))

  }

  # compute color
  cols = col2rgb(rainbow(ncol(Q)))

  # overlay raster
  GetCoefFunction <- function(cols) {
    return(function(s) {
      idmax <- which.max(s)
      if(length(idmax) == 0) {
        return(c(NA,NA,NA))
      } else {
        return(c(cols[1,idmax],cols[2,idmax],cols[3,idmax]) * s[idmax])
      }
    })
  }

#   GetCoefFunction <- function(cols) {
#     return(function(s) {
#
#         return(cols %*% s)
#
#     })
#   }

  ancestry.raster.rgb <- raster::overlay(raster::stack(ancestry.raster.list), fun=GetCoefFunction(cols))



#   raster::plotRGB(ancestry.raster.rgb)
#   sp::plot(dat, add=TRUE, pch=16, cex=0.5)

  return(ancestry.raster.rgb)

}


