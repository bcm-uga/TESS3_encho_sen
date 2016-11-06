#' Compute allele frequency matrix from genotype frequency matrix.
#'
#' @param G genotype frequency matrix.
#' @param ploidy the number of chromosomes.
#'
#' @return a matrix of allele frequencies for each locus.
#' @export
GtoFreq <- function(G, ploidy) {
  G.array <- array(G, c(ploidy + 1, nrow(G) / (ploidy + 1), ncol(G)))
  Freq <- apply(G.array, c(2,3), function(g) sum(0:ploidy * g / ploidy))
  return(Freq)
}
