#################################
##### Useful extra functions ####
#################################
#' Group by population
#' @title Group by population
#' @description Computes ancestry coefficients and sampled coordinates averaged
#'   over specified populations.
#' @author Flora Jay, Kevin Caye, Olivier François
#' @param Q An object of class \code{tess3Q} or \code{matrix}. The ancestry
#'   coefficient matrix.
#' @param coord The numeric matrix of size \eqn{n \times 2} where n is the
#'   number of individuals. If \code{coord=NULL}, only the ancestry coefficients
#'   are averaged.
#' @param pop Population labels. Vector of \code{n} elements of class
#'   \code{numeric}, \code{character} or \code{factor}
#'
#' @return An object of class \code{tess3bypop}, which is a list of:
#'
#'\describe{ \item{Q}{\code{tess3Q} object. Average admixture coefficients for
#'each population sample of dimensions \code{npopulations x K}}
#'
#'\item{coord}{\code{matrix} Average location for each population sample of
#'dimensions \code{npopulations x 2}}
#'
#'\item{labels}{Vector of unique population names} }
#'
#' @examples
#' data(data.for.test)
#' pop = rep(c("popA","popB","popC"),each=50)
#' # or
#' pop = rep(1:3,each=50)
#' # or (only if individuals are ordered by population)
#' pop = unlist(lapply(1:data.for.test$K,FUN= function(i) rep(paste("pop",i,sep=""),data.for.test$n.by.pop[i])))
#'
#' info.bypop = bypop(data.for.test$Q, data.for.test$coord, pop)
#' barplot(info.bypop$Q)
#'
#' @export
bypop = function(Q,coord=NULL,pop) {

  ulabels=unique(pop)
  nlabels=length(ulabels)
  qpop = matrix(NA, ncol = ncol(Q), nrow = nlabels)
  rownames(qpop) = ulabels
  if (is.null(coord)) {
    coord.pop = NULL
  } else {
    coord.pop = matrix(NA, ncol = 2, nrow = nlabels)
    rownames(coord.pop) = ulabels
  }
  pop.size = NULL
  for (i in 1:nlabels){
    qpop[i,] = apply(Q[pop == ulabels[i],,drop=FALSE], 2, mean)
    if (!is.null(coord)) {
      coord.pop[i,] = apply(coord[pop == ulabels[i],,drop=FALSE], 2, mean)
    }
    pop.size = c(pop.size, sum(pop == ulabels[i]))
  }
  class(qpop) = "tess3Q"
#  names(pop.size) = ulabels
  res = list(Q=qpop, coord=coord.pop, labels=ulabels, size=pop.size)
  class(res)="tess3bypop"
  return(res)
}



#' Internal function for plotting gradient legend keys
#' @description Slightly adapted from \link[fields]{image.plot}. Plots only legend key, no image.
#' To make it compatible with multiple key plotting using layout
#'
#' @details See \code{\link[fields]{image.plot}} for further description
#'
image.plot.legend = function (..., add = FALSE, breaks = NULL, nlevel = 64, col = NULL,
                              horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2,
                              legend.mar = NULL, legend.lab = NULL,
                              legend.line = 2, lab.breaks = NULL,
                              axis.args = NULL, legend.args = NULL, legend.cex = 1, midpoint = FALSE,
                              border = NA, lwd = 1, verbose = FALSE, main.par=NULL)
{
  old.par <- par(no.readonly = TRUE)
  if (is.null(col)) {
    col <- tim.colors(nlevel)
  }
  else {
    nlevel <- length(col)
  }
  info <- fields::imagePlotInfo(..., breaks = breaks, nlevel = nlevel)
  breaks <- info$breaks
  if (verbose) {
    print(info)
  }
  #breaks[1]=min(0,info$zlim[1])
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  if (verbose) {
    print(breaks)
    print(midpoints)
    print(ix)
    print(iy)
    print(iz)
    print(col)
  }
  if (is.null(legend.mar)) {
    if (is.null(main.par)) {
        legend.mar <- old.par$mar
      } else {
        legend.mar <- main.par$mar
      }
    if (horizontal) {
      legend.mar[3] <-  0
      legend.mar[1] <- 5.1
    } else  {
      legend.mar[2] <-  0
      legend.mar[4] <- 4.1
    }
  }
  par(mar=legend.mar)
  if (!horizontal) {
    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
          ylab = "", col = col, breaks = breaks, add = add)
  }
  else {
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
          ylab = "", col = col, breaks = breaks, add = add)
  }
  if (!is.null(lab.breaks)) {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
                   axis.args)
  }
  do.call("axis", axis.args)
  box()
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal,
                                                         1, 4), line = legend.line, cex = legend.cex)
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }

  #par(old.par)
}

