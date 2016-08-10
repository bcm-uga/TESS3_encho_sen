#' Compute allele frquencies matrix from the genotype frequencie matrix.
#'
#' @param G Genotype frequencie matrix.
#' @param ploidy The number of chromosome.
#'
#' @return TODOC
#' @export
GtoFreq <- function(G, ploidy) {
  G.array <- array(G, c(ploidy + 1, nrow(G) / (ploidy + 1), ncol(G)))
  Freq <- apply(G.array, c(2,3), function(g) sum(0:ploidy * g / ploidy))
  return(Freq)
}
