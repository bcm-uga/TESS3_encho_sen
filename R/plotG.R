#' Manhattan plot for tess3 significance values.
#'
#' @param x tess3pvalue object.
#' @param ... \code{\link{plot.default}} other parameters.
#'
#' @return none
#' @seealso See \code{\link{tess3Main}} more explainations and examples.
#' @export
plot.tess3pvalue <- function(x, ...) {
  locus.index <- 1:length(x)
  plot(locus.index, -log10(x),...)
}
