#' Title
#'
#' @param Dz
#' @param Dx
#' @param breaks
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
CalculateEmpiricalSemivariogram <- function(Dz, Dx, breaks = "FD", na.rm = TRUE) {

  message("Computing semi variance")
  breaks <- hist(Dx, breaks = breaks, plot = FALSE)$breaks
  cuts <- cut(Dx, breaks = length(breaks))
  levs <- levels(cuts)
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



#' Title
#'
#' @param X
#' @param ploidy
#' @param coord
#' @param breaks
#' @param na.rm
#'
#' @return
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
  CheckXCoord(X, ploidy, coord)
  X2XBin(X, ploidy, XBin)
  rm(X)
  message("Computing distance matrices")
  dx <- dist(XBin, method = "manhattan") / ncol(XBin)
  dgeo <- dist(coord)

  return(CalculateEmpiricalSemivariogram(dx, dgeo, breaks = breaks, na.rm = na.rm))
}


#' todo
#'
#' @param semi.variogram
#' @param epsilon
#'
#' @return
#'
#' @examples
FitGeneralSemivariogram <- function(semi.variogram, epsilon = 1e-6) {
  message("DO NOT work")
  res <- list()
  res$nugget <- semi.variogram$semi.variance[1]
  # find k when the semi variance reach asymptote
  # aux.mean <- sapply(seq_along(semi.variogram$semi.variance), function(k) mean(semi.variogram$semi.variance[-(1:k)]))
  aux.var <- sapply(seq_along(semi.variogram$semi.variance), function(k) var(semi.variogram$semi.variance[-(1:k)]))
  # aux.mean.delta <- sapply(seq_along(head(aux.mean,-3)), function(k) abs(aux.mean[k] - aux.mean[k + 1]))
  aux.var.delta <- sapply(seq_along(head(aux.var,-3)), function(k) abs(aux.var[k] - aux.var[k + 1]))
  k <- min(which(aux.var.delta < epsilon))
  res$still <- mean(semi.variogram$semi.variance[-(1:(k-1))])
  res$range <- semi.variogram$h[k]
 return(res)
}
