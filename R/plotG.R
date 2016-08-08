#' Plot the manhattan plot of -log(p.value).
#'
#' @param x tess3pvalue object.
#' @param ... TODOC
#'
#' @return TODOC
#' @export
plot.tess3pvalue <- function(x, ...) {
  loci.index <- 1:length(x)
  plot(loci.index, -log(x), ...)
}
