#' Title
#'
#' @param tess3.res
#' @param masked.X
#'
#' @return
#' @export
#'
#' @examples
ImputeRandom <- function(tess3.res, masked.X) {
  if (!is.tess3Main(tess3.res)) {
    stop("tess3.res must be a tess3 result of class tess3Main. You can use function tess3Main or Gettess3res on a tess3 object.")
  }
  ploidy <- max(masked.X, na.rm = TRUE)
  Prob <- tess3.res$Q %*% t(tess3.res$G)
  Prob <- array(Prob, c(nrow(Prob), ploidy + 1, ncol(Prob) / (ploidy + 1)))

  for (i in 1:nrow(masked.X)) {
    for (j in 1:ncol(masked.X)) {
      if (is.na(masked.X[i, j])) {
        masked.X[i, j] <- sample(0:ploidy, 1, replace = FALSE, prob = Prob[i,,j])
      }
    }
  }
  return(masked.X)
}

ImputeRound <- function(tess3.res, masked.X) {
  if (!is.tess3Main(tess3.res)) {
    stop("tess3.res must be a tess3 result of class tess3Main. You can use function tess3Main or Gettess3res on a tess3 object.")
  }
  ploidy <- max(masked.X, na.rm = TRUE)
  P <- tess3.res$Q %*% t(GtoFreq(tess3.res$G, ploidy))
  masked.id <- which(is.na(masked.X))
  masked.X[masked.id] <- round(P[masked.id])
  return(masked.X)
}
