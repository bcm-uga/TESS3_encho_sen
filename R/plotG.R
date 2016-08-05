#' Plot the manhattan plot of -log(p.value).
#'
#' @param pvalue tess3pvalue object.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.tess3pvalue <- function(pvalue, ...) {
  loci.index <- 1:length(pvalue)
  plot(loci.index, -log(pvalue), ...)
}

