#' Compute the probability matrix
#'
#'
#'
#' @export
ComputeProba <- function(genotype, ploidy) {
  genotype[is.na(genotype)] <- as.integer(-1)
  # To be sure we have an integer genotype matrix
  genotype <- matrix(as.integer(genotype), nrow(genotype), ncol(genotype))
  proba <- ComputeXBin(genotype, ploidy)
  proba[proba < 0] <- NA
  return(proba)
}
