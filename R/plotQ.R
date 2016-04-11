###################################################
##################pie chart functions##############
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
PlotPiechartAncestryCoef <- function(Q, coord, window) {
  TestRequiredPkg("maps")
  TestRequiredPkg("mapplots")

  if (is.null(window)) {
    sd.x <- sd(coord[,1])
    sd.y <- sd(coord[,1])
    window <- c(min(coord[,1]) - 0.05 * sd.x, max(coord[,1]) + 0.05 * sd.x, min(coord[,2]) - 0.05 * sd.y, max(coord[,2]) + 0.05 * sd.y)
  }

  require("maps")
  message("This function required to load maps package.")

  maps::map(database = "world",
               ylim = window[3:4],
               xlim = window[1:2],
               col = "grey80",
               fill = TRUE)

  col <- topo.colors(ncol(Q), alpha = 0.6)
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
##################with raster functions############
###################################################


#' TODO
#'
#'
#' @export
#' @param idp
idw <- function(idp=1.0){
  TestRequiredPkg("gstat")

  return(function(dat,grd) {
    return(gstat::idw(Z~1,dat,newdata=grd,idp=idp))
  })

}

#' TODO
#'
#'
#' @export
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
ComputeBackgroundandGrid <- function(dat, raster.filename, resolution, window = NULL) {
  imported.raster <- raster::raster(raster.filename)

  if (is.null(window)) {
    # window from coord
    window.extent <- raster::extent(dat)
  } else {
    window.extent <- raster::extent(window)
  }

  # crop to dat
  imported.raster <- raster::resample(imported.raster,raster::raster(window.extent, ncol = resolution$ncol, nrow = resolution$nrow))
  # imported.raster <- raster::calc(imported.raster, function(x) {ifelse(x>0,1,0)})
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


RasterAncestryCoefStack <- function(Q, coord, raster.filename = NULL, interpolation.function = idw(), resolution = c(300, 300), background = TRUE, window = NULL) {

  if (is.null(raster.filename)) {
    raster.filename <- system.file("extdata/raster","earth.tif",package = "TESS3enchoSen")
  }

  # compute grid and background
  dat <- data.frame(X = coord[,1], Y = coord[,2], Z = Q[,1])
  # Force data frame object into a SpatialPointsDataFrame object
  sp::coordinates(dat) <- c("X","Y")
  resolution = list(ncol = resolution[1], nrow = resolution[2])
  grid <- ComputeBackgroundandGrid(dat, raster.filename, resolution, window = window)

  # compute raster for each ancestral population
  ancestry.raster.list = NULL
  for (i in 1:ncol(Q)) {
    dat$Z = Q[,i]
    # Interpolate the surface using a power value of 2 (idp=2.0)
    dat.interpolated <- interpolation.function(dat, grid$grd)

    if (background) {
      ancestry.raster.list <- c(ancestry.raster.list, grid$background * raster::raster(dat.interpolated))
    } else {
      ancestry.raster.list <- c(ancestry.raster.list, raster::raster(dat.interpolated))
    }

  }
  return(raster::stack(ancestry.raster.list))
}


#' TODO
#'
#'
RasterAncestryCoefRgb <- function(Q, coord, raster.filename = system.file("extdata/raster","earth.tif",package = "TESS3enchoSen"), interpolation.function = idw(), resolution = list(ncol = 300, nrow = 300)) {

  raster.stack <- RasterAncestryCoefStack(Q, coord, raster.filename, interpolation.function, resolution)

  # compute color
  #palette <- c("red4","green4","blue4","orange4")
  cols = col2rgb(rainbow(ncol(Q)))
  #cols = col2rgb(palette[1:ncol(Q)])
  # overlay raster
  GetCoefFunction <- function(cols) {
    return(function(s) {
      idmax <- which.max(s)
      if (length(idmax) == 0) {
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

  ancestry.raster.rgb <- raster::overlay(raster.stack, fun = GetCoefFunction(cols))



  #   raster::plotRGB(ancestry.raster.rgb)
  #   sp::plot(dat, add=TRUE, pch=16, cex=0.5)

  return(ancestry.raster.rgb)

}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  TestRequiredPkg("grid")

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

PlotRasterStack <- function(raster, coord) {

  TestRequiredPkg("ggplot2")
  TestRequiredPkg("raster")

  K = dim(raster)[3]
  df <- raster::as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  # rename column (col 1 and 2 are x y)
  names(df)[-c(1,2)] <- sapply(1:K, function(k) paste0("layer.",k))
  toplot <- list()
  for (k in 1:K) {
    toplot[[k]] <-  ggplot2::ggplot(df) +
      ggplot2::geom_raster(ggplot2::aes_string(x = "x", y = "y", fill = paste0("layer.",k))) +
      ggplot2::scale_fill_gradient(low = "grey", high = "red4", name = "ancestry coef") +
      ggplot2::ggtitle(paste0("Ancestral Population ",k)) + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_point(data = as.data.frame(coord), ggplot2::aes(x = V1, y = V2))
  }
  return(multiplot(plotlist = toplot, cols = ceiling(sqrt(K))))
}

###################################################
##################flora's functions################
###################################################

lColorGradients = list(
  c("gray95",RColorBrewer::brewer.pal(9,"Reds")),
  c("gray95",RColorBrewer::brewer.pal(9,"Greens")),
  c("gray95",RColorBrewer::brewer.pal(9,"Blues")),
  c("gray95",RColorBrewer::brewer.pal(9,"YlOrBr")),
  c("gray95",RColorBrewer::brewer.pal(9,"RdPu")),
  c("gray95",RColorBrewer::brewer.pal(9,"Greys"))
)

#' TODO
#'
PlotAncestryCoef <- function(Q,
                             coord,
                             resolution = c(300,300),
                             window = NULL,
                             background = TRUE,
                             raster.filename = NULL) {
  # read raster file
  if (is.null(raster.filename)) {
    raster.filename <- system.file("extdata/raster","earth.tif",package = "TESS3enchoSen")
  }
  grid <- CreateGridAndBackgroundFromRasterFile(raster.filename, coord, resolution, window)

  K = ncol(Q)

  if(!background) {
    grid$background[] <- 1
  }

  maps(matrix = Q,
       coord = coord,
       grid=grid$fields.grid,
       constraints=grid$background,
       method="max",
       main=paste("ancestry coefficient with K =",K))


}


#' TODO
#'
CreateGridAndBackgroundFromRasterFile <- function(file, coord, resolution, window)
{

  raster.filename <- system.file("extdata/raster","earth.tif",package = "TESS3enchoSen")
  imported.raster <- raster::raster(raster.filename)

  # resize raster
  if (is.null(window)) {
    # window from coord
    names(coord) <- c("x","y")
    window.extent <- raster::extent(coord)
  } else {
    window.extent <- raster::extent(window)
  }
  imported.raster <- raster::resample(imported.raster,raster::raster(window.extent,
                                                                     ncol = resolution[1],
                                                                     nrow = resolution[2]))
  # plot(imported.raster)

  # background matrix
  background <- raster::as.matrix(imported.raster)
  background[is.na(background)] <- 0
  background <- t(apply(t(background), 1, rev)) # transpose and rev col
  # image(background)


  fields.grid <- fields::make.surface.grid(list( raster::xFromCol(imported.raster),rev(raster::yFromRow(imported.raster))))

  return(list(background = background, fields.grid = fields.grid))
}

#' TODO
#'
maps <- function(matrix, coord, grid, constraints=NULL, method="treshold", colorGradientsList = lColorGradients, onemap = T, onepage = T, ...)
{
  if ( (method != "treshold") & (method != "max")) {stop(paste("Unknown method",method))}
  if (class(constraints)!= "NULL") {
    if ( nrow(grid) != nrow(constraints)*ncol(constraints) ) {
      stop(paste("Argument grid assumes", nrow(grid), "pixels, but argument constaints assumes", nrow(constraints)*ncol(constraints),"pixels"))
    }
  }

  if (onemap & method=="max") {
    mapsMethodMax(matrix=matrix,coord=coord,grid=grid,constraints=constraints,colorGradientsList=colorGradientsList,...)

  } else {
    K=ncol(matrix)
    if (length(colorGradientsList)<K)
    {
      stop(paste(K,"clusters detected but only",length(colorGradientsList),"color gradient(s) defined.",
                 "You should complete colorGradientsList to have as many gradients as clusters."))
    }

    if (!onemap & onepage) {KK = as.integer(K/2); par(mfrow = c(2,KK+1)) }

    for (k in 1:K)
    {
      clust=NULL
      clust= fields::Krig(coord, matrix[,k], theta = 10)
      look<- predict(clust,grid) # evaluate on a grid of points
      out<- fields::as.surface( grid, look)

      if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }

      ncolors=length(colorGradientsList[[k]])
      if (onemap)
      {
        out[[8]][ out[[8]] < .5 ] = NA
        graphics::image(out,add=(k>1),col=colorGradientsList[[k]][(ncolors-4):ncolors],breaks=c(seq(.5,.9,.1),+200),...)
      } else {
        graphics::image(out,col=colorGradientsList[[k]][(ncolors-9):ncolors],breaks=c(-200,.1,seq(.2,.9,.1),+200),...)
        graphics::points(coord,pch=19)
      }
    }

    if (onemap) { graphics::points(coord,pch=19) }
  }

}

mapsMethodMax <- function(matrix,coord,grid,constraints,colorGradientsList,...)

{
  K=ncol(matrix)
  if (length(colorGradientsList)<K)
  {
    stop(paste(K,"clusters detected but only",length(colorGradientsList),"color gradient(s) defined.",
               "You should complete colorGradientsList to have as many gradients as clusters."))
  }

  listOutClusters=NULL
  matrixOfVectors =NULL
  for (k in 1:K)
  {
    clust=NULL
    clust= fields::Krig(coord, matrix[,k], theta = 10)
    look<- predict(clust,grid) # evaluate on a grid of points
    out<- fields::as.surface( grid, look)
    listOutClusters[[k]] = out[[8]]
    matrixOfVectors = cbind(matrixOfVectors,c(out[[8]]))
  }
  long = out[[1]]
  lat = out[[2]]

  whichmax = matrix(apply(matrixOfVectors ,MARGIN=1,FUN=which.max),nrow=length(long))

  for (k in 1:K)
  {
    ncolors=length(colorGradientsList[[k]])
    if (class(constraints)!= "NULL") { listOutClusters[[k]][ !constraints ] = NA }
    listOutClusters[[k]][ whichmax != k ] = NA
    image(long,lat,listOutClusters[[k]],add=(k>1),col=colorGradientsList[[k]][(ncolors-9):ncolors],breaks=c(-200,.1,seq(.2,.9,.1),+200),...)
  }
  points(coord,pch=19)
}


###################################################
#######################Methods#####################
###################################################

#' Plot ancestry coeficient on a map.
#'
#' @param Q
#' @param coord
#' @param plot.type
#'
#' @return
#' @export
#'
#' @examples
plot.tess3Q <- function(Q, coord, plot.type = "florastyle", resolution = c(300,300), window = NULL, background = TRUE, raster.filename = NULL, interpolation.function = idw()) {

  if (plot.type == "piechart") {
    PlotPiechartAncestryCoef(Q, coord, window)
    if (is.null(window))
      PlotPiechartAncestryCoef(Q, coord, window)
  } else if (plot.type == "florastyle") {
    PlotAncestryCoef(Q = Q, coord = coord, resolution, window, background, raster.filename)
  } else if (plot.type == "ggplot") {
    raster.stack <- RasterAncestryCoefStack(Q = Q,
                            coord = coord,
                            raster.filename = raster.filename,
                            interpolation.function = interpolation.function,
                            resolution = resolution,
                            background = background,
                            window = window)
    PlotRasterStack(raster = raster.stack, coord = coord)
  }
}
