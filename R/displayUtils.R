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
#' pop = unlist(lapply(1:data.for.test$K,
#'                    FUN= function(i) rep(paste("pop",i,sep=""),
#'                    data.for.test$n.by.pop[i])))
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
#' @description Slightly adapted from \link[fields]{image.plot}. Plots only
#'   legend key, no image. To make it compatible with multiple key plotting
#'   using layout. Extra parameters for tuning legend configuration (eg spacing
#'   between keys)
#' @param legend.mar A vector of size 4 used for the number of lines of extra
#'   margins around color keys (see ?par mar for details)
#' @details See \code{\link[fields]{image.plot}} for description of other
#'   parameters.
#'
img.plot.legend = function (..., add = FALSE, breaks = NULL, nlevel = 64, col = NULL,
                              horizontal = FALSE, legend.shrink = NA, legend.width = NA,
                              legend.mar = rep(0,4),
                              legend.lab = NULL, legend.line = 2, lab.breaks = NULL,
                              axis.args = NULL, legend.args = NULL, legend.cex = 1, midpoint = FALSE,
                              border = NA, lwd = 1, verbose = FALSE)
{
  #  old.par <- par(no.readonly = TRUE)
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

  if (!is.null(legend.mar))  par(mar = legend.mar)
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
  graphics::box()
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal,
                                                         1, 4), line = legend.line, cex = legend.cex)
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }

}


#' Setup correct layout for multiple key legend
#' @description Internal function called by \link{PlotInterpotationMax}
#' @param K Number of clusters (ie number of color keys)
#' @param legend.ncol Number of columns (resp. rows)  in legend when
#'   \code{horizontal = FALSE} (resp. \code{TRUE}).
#' @param legend.space vector of size two  (x-s, y-axis) for the space size
#'   (number of lines) between color keys.
#' @param legend.width Width (number of lines) of each column (resp. row) for
#'   the vertical (resp. horizontal) legend.
#' @param horizontal Boolean. Horizontal position for legend key?
#' @return the number of figures in the layout (invisibly)
#' @seealso \code{\link[graphics]{layout}}
#' @export
#' @examples
#' require(graphics)
#' require(tess3r)
#' im <- img.plot.legend.setup(K=5, legend.ncol=2,
#'            legend.space=c(3,4), legend.width=3, horizontal=F)
#' # Layout on which the map and keys will be displayed
#' layout.show(im$nn)
#' # Reset graphic device for future plots
#' dev.off()
#'
img.plot.legend.setup = function(K, legend.ncol, legend.space, legend.width, horizontal) {
  op <- par()
  column.nkey <- ceiling(K/legend.ncol) #umber of color keys per layout column
  layout.ncol <- legend.ncol*2+2 # adding 1 column for main plot and space columns
  layout.nrow <- column.nkey*2 + 1 # counting empty columns used for spacing + right margin

  if (horizontal) {
    nlines <- op$cra[1]/op$cin[1]
    key.size <- (nlines - op$mar[4] - op$mar[2] - (column.nkey-1)*legend.space[2])/column.nkey
  } else {
    nlines <- op$cra[2]/op$cin[2]
    key.size <- (nlines - op$mar[3] - op$mar[1] - (column.nkey-1)*legend.space[1])/column.nkey
  }
  if (key.size <0) stop("Error: Figure too large, decrease legend.space or increase legend.ncol")
  #
  #     if (!horizontal) {
  #       layout(matrix(c(rep(1,layout.nrow),2:((layout.ncol-1)*layout.nrow+1)),ncol=layout.ncol),
  #              width=c(layout.width - legend.width*(layout.ncol-1),rep(legend.width,layout.ncol-1)))
  #     } else {
  #       layout(t(matrix(c(rep(1,layout.nrow),2:((layout.ncol-1)*layout.nrow+1)),ncol=layout.ncol)),
  #              height=c(layout.width - legend.width*(layout.ncol-1),rep(legend.width,layout.ncol-1)))
  #     }

  refvec = c(1:layout.nrow,rep(layout.nrow+1,layout.nrow))+1
  layout.indexes = c(sapply(1:legend.ncol, FUN=function(col) refvec+(layout.nrow+1)*(col-1)))
  layout.indexes = c(layout.indexes, rep(max(layout.indexes)+1,layout.nrow))

  if (!horizontal) {
    nn <- graphics::layout(matrix(c(rep(1,layout.nrow),layout.indexes),ncol=layout.ncol),
                 width = c(nlines - legend.width*(layout.ncol-1),
                           rep(c(legend.width,legend.space[1]),legend.ncol), op$mar[4]),
                 height = c(op$mar[3],rep(c(legend.space[2],key.size),column.nkey)[-1],op$mar[1])
    )
  } else {
    nn <- graphics::layout(t(matrix(c(rep(1,layout.nrow),layout.indexes),ncol=layout.ncol)),
                 height = c(nlines - legend.width*(layout.ncol-1),
                            rep(c(legend.width,legend.space[2]),legend.ncol), op$mar[1]),
                 width = c(op$mar[2],rep(c(legend.space[1],key.size),column.nkey)[-1],op$mar[4])
    )
  }

  res <- list(nn = nn, column.nkey = column.nkey)
  return(res)
}
