###################################################
#' This function tests Q matrix objects, and converts \code{matrix} objects into valid Q matrices
#' @title Converts Q matrix
#' @author Kevin Caye, Olivier François
#' @param Q an object of class \code{matrix} containing a matrix of ancestry coefficients.
#' @return An object of class \code{tess3Q}.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{barplot.tess3Q}}
#' @examples
#' ## an example with 3 individuals and 2 clusters
#' Qmatrix <- matrix(c(0.4,0.6,0.3,0.7, 0.2, 0.8), byrow = T, nrow = 3)
#' Qmatrix <- as.qmatrix(Qmatrix)
#' barplot(Qmatrix, space = 0, xlab = "individuals", ylab = "Ancestry proportions", main = "Ancestry matrix")
#' @export
as.qmatrix <- function(Q){
  if (class(Q) != "matrix") stop("Input matrix is not an ancestry matrix.")
  if (min(Q) < 0) stop("Q contains negative elements.")
  sumofq <- apply(Q, MARGIN = 1, sum)
  if ( sum(sumofq) != nrow(Q)) stop("Input matrix is not an ancestry matrix.")
  class(Q) = "tess3Q"
  return(Q)
}

#' This function creates a list of color palettes for the plot and barplot functions
#' @title Create a list of palettes
#' @author Kevin Caye, Olivier François
#' @param color.vector a vector of R colors.
#' @param palette.length an integer number of colors in each palette.
#' @return An object of class \code{list} containing a list of color palettes.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{barplot.tess3Q}}
#' @examples
#' ## an A. thaliana example
#' data(data.at)
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5, ploidy = 1, openMP.core.num = 4)
#' Qmatrix <- obj$Q
#' my.colors <- c("tomato", "yellow", "blue", "wheat","olivedrab")
#' my.palette <- CreatePalette(my.colors, 9)
#' plot(obj$Q, data.at$coord, method = "mapping.max", col.palette = my.palette, interpol = kriging(10), cex = .4, xlab = "Longitude", ylab= "Latitude", main = "Ancestry coefficients")
#' @export
CreatePalette <- function(color.vector = c("tomato", "chartreuse", "gold", "blue", "violet", "wheat","olivedrab"), palette.length = 9){
  ll = NULL
   for (i in 1:length(color.vector)){
     ll[[i]] = colorRampPalette(c("grey96", color.vector[i]))(palette.length)
   }
 ll
}

#' Convert Fst into t-score and compute p value
#'
#'
#'
ComputeTscoreAndPvalue <- function(Fst, K, n) {
  res <- list()
  res$Fscore = Fst / (1 - Fst) * (n - K) / (K - 1)
  # Pvalue, we assume F.score ~ gif * F(df1 = K -1, df2 = n - K)
  res$gif = median(res$Fscore, na.rm = TRUE) / qf(0.5, df1 = K - 1, df2 = n - K)
  res$pvalue <- pf(res$Fscore / res$gif, df1 = K - 1, df2 = n - K, lower.tail = FALSE)
  return(res)
}

#' Convert Fst into chi2 and compute p value
#'
#'
#'
ComputeChi2AndPvalue <- function(Fst, K, n) {
  res <- list()
  # Convert Fst into chi 2
  res$chi2 = Fst * (n - K)/(1 - Fst)
  # compute the gif
  res$gif = median(res$chi2) / qchisq(1 / 2, df = K - 1)
  # compute adjusted p-values from the combined z-scores
  res$pvalue = as.numeric(pchisq(res$chi2 / res$gif, df = K - 1, lower.tail = FALSE))
  return(res)
}
