#' Manhattan plot for tess3 significance values.
#'
#' @param x tess3pvalue object.
#' @param ... TODOC
#'
#' @return TODOC
#' @export
plot.tess3pvalue <- function(x, ...) {
  locus.index <- 1:length(x)
  plot(locus.index, -log10(x),...)
}
