#' Calculate empirical semi-variance.
#'
#' @param Dz Variable distance matrix.
#' @param Dx Spatial distance matrix.
#' @param breaks Same parameter that in hist R base function.
#' @param na.rm A logical indicating whether missing values should be removed.
#'
#' @return Semi-variance.
#' @export
CalculateEmpiricalSemivariogram <- function(Dz, Dx, breaks = "FD", na.rm = TRUE) {

  message("Computing semi variance")
  breaks <- graphics::hist(Dx, breaks = breaks, plot = FALSE)$breaks
  cuts <- base::cut(Dx, breaks = length(breaks))
  levs <- base::levels(cuts)
  vario.res <- 1:length(levs)
  vario.size <- 1:length(levs)
  h <- seq(0,max(Dx) - min(Dx), length.out = length(levs))
  for (k in seq_along(levs)) {
    vario.size[k] <- length(which(!is.na(Dz[cuts == levs[k]])))
    if (vario.size[k] > 0) {
      vario.res[k] <- mean(Dz[cuts == levs[k]], na.rm = na.rm)
    } else {
      vario.res[k] <- NA
    }

  }
  return(data.frame(h = h, interval = levs, size = vario.size, semi.variance = vario.res))
}



#' Calculate empirical semi variance from genotype matrix.
#'
#' @param X Genotype matrix.
#' @param ploidy The number of chromosome.
#' @param coord Coordinate matrix.
#' @param breaks Same parameter that in hist R base function.
#' @param na.rm A logical indicating whether missing values should be removed.
#'
#' @return Semi-variance
#' @export
#'
#' @examples
#' library(tess3r)
#'
#' data("data.for.test", package = "tess3r")
#' em.vario <- CalculateEmpiricalGenSemivariogram(X = data.for.test$X,
#'                                 ploidy = 1,
#'                                 coord = data.for.test$coord,
#'                                 breaks = "FD",
#'                                 na.rm = TRUE)
#' ggplot2::ggplot(em.vario, ggplot2::aes(x = h, y = semi.variance, size = size)) +
#'                ggplot2::geom_point()
CalculateEmpiricalGenSemivariogram <- function(X, ploidy, coord, breaks = "FD", na.rm = TRUE) {
  # ensure type of X
  X <- matrix(as.double(X), nrow(X), ncol(X))
  XBin <- matrix(as.double(X), nrow(X), ncol(X) * (ploidy + 1))
  CheckXCoord(X, coord)
  X2XBin(X, ploidy, XBin)
  rm(X)
  message("Computing distance matrices")
  dx <- stats::dist(XBin, method = "manhattan") / ncol(XBin)
  dgeo <- stats::dist(coord)

  return(CalculateEmpiricalSemivariogram(dx, dgeo, breaks = breaks, na.rm = na.rm))
}
