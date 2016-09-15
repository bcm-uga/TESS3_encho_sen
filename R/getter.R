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

#' Get Q matrix estimate.
#'
#' Get the \code{Q} ancestry coefficient matrix returned by the TESS3 algorithm.
#'
#' @param tess3 A tess3 object.
#' @param K Number of ancestral populations.
#' @param rep Repetition index.
#'
#' @seealso \code{\link{tess3}}
#'
#' @return The \code{Q} ancestry coefficient matrix.
#'
#' @export
qmatrix <- function(tess3, K, rep = "best") {
  return(Gettess3res(tess3, K, rep)$Q)
}


#' Get pvalue vector.
#'
#' Get the \code{pvalue} vector returned by the TESS3 algorithm.
#'
#' @param tess3 A tess3 object.
#' @param K Number of ancestral populations.
#' @param rep Repetition index.
#'
#' @return The \code{pvalue} ancestry coefficient matrix.
#'
#' @seealso \code{\link{tess3}}
#'
#' @export
pvalue <- function(tess3, K, rep = "best") {
  return(Gettess3res(tess3, K, rep)$pvalue)
}
