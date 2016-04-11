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
  TestRequiredPkg("ggplot2")
  return(ggplot2::ggplot(data.frame( pvalue = as.vector(pvalue), index = 1:length(pvalue))) + ggplot2::geom_point(ggplot2::aes(x = index, y = -log(pvalue))))
}
