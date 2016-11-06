#' Get tess3 run result
#'
#' Get the \code{\link{tess3Main}} object returned by the TESS3 algorithm.
#'
#' @param tess3 A tess3 object.
#' @param K Number of ancestral populations.
#' @param rep Repetition index.
#'
#' @return An object of class \code{\link{tess3Main}}.
#'
#'
#' @seealso \code{\link{tess3}}
#'
#' @export
Gettess3res <- function(tess3, K, rep = "best") {
  if (is.tess3Main(tess3)) {
    if (rep != 1 & rep != "best") {
      # stop("Tess3 algorithm was run only one time.")
      return(NULL)
    }
    if (K != ncol(tess3$Q)) {
      # stop("This value of K is not available.")
      return(NULL)
    }
    return(tess3)
  }
  if (!is.tess3(tess3)) {
    stop("tess3 must of class tess3")
  }
  if (rep == "best") {
    best.rep <- min(which.min(tess3[[K]]$rmse)[1],length(tess3[[K]]$tess3.run))
  } else {
    best.rep <- min(as.numeric(rep),length(tess3[[K]]$tess3.run))
  }
  return(tess3[[K]]$tess3.run[[best.rep]])
}

#' Get Q matrix estimate
#'
#' Get the \code{Q} ancestry coefficient matrix computed by the tess3 function.
#'
#' @param tess3 a tess3 object.
#' @param K number of ancestral populations.
#' @param rep an integer indicating which run to display (default = lowest error run).
#'
#' @seealso \code{\link{tess3}} \code{\link{barplot.tess3Q}} \code{\link{plot.tess3Q}} \code{\link{CreatePalette}}
#'
#' @return a \code{Q} matrix of ancestry coefficients from a specified run.
#'
#' @examples
#' library(tess3r)
#'
#' # load Arabidopsis data
#' data(data.at)
#'
#' # Running tess3
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
qmatrix <- function(tess3, K, rep = "best") {
  return(Gettess3res(tess3, K, rep)$Q)
}


#' Get p-values from a genome scan for selection
#'
#' Get the \code{pvalue} vector returned by the tess3 algorithm.
#'
#' @param tess3 an object of class \code{tess3}.
#' @param K a number of ancestral populations.
#' @param rep an integer correponding to a run number (default is 'best' run).
#'
#' @return Corrected pvalues for all loci to be used in a scan for adaptive alleles.
#'
#' @seealso \code{\link{tess3}}
#'
#' @export
pvalue <- function(tess3, K, rep = "best") {
  return(Gettess3res(tess3, K, rep)$pvalue)
}
