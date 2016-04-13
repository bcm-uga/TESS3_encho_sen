#' Manhattan plot
#'
#' @param pvalue
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
