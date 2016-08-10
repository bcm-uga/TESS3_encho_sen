#' Get tess3 run result
#'
#' @param tess3 A tess3 object.
#' @param K Number of ancestral populations.
#' @param rep Repetition index.
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
#' @param tess3 A tess3 object.
#' @param K Number of ancestral populations.
#' @param rep Repetition index.
#'
#' @export
qmatrix <- function(tess3, K, rep = "best") {
  return(Gettess3res(tess3, K, rep)$Q)
}


#' Get Q matrix estimate.
#'
#' @param tess3 A tess3 object.
#' @param K Number of ancestral populations.
#' @param rep Repetition index.
#'
#' @export
pvalue <- function(tess3, K, rep = "best") {
  return(Gettess3res(tess3, K, rep)$pvalue)
}
