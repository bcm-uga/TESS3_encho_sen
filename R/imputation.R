#' Impute missing values
#'
#' \code{ImputeRandom} imputes missing values with the Q ancestry coeficient matrix and G ancestral
#' frequencies estimates. Missing values are sampled with the estimated genotype
#' frequencies (P = Q G^T).
#'
#' @param tess3.res tess3Main object with Q and G estimates.
#' @param masked.X Genotype matrix with the missing values (NA values).
#'
#' @return the imputed genotype matrix.
#' @export
#'
#' @examples
#' library(tess3r)
#'
#  #' sample data
#' n <- 100
#' K <- 3
#' ploidy <- 2
#' L <- 3001
#' data.list <- SampleGenoFromGenerativeModelTESS3(G = SampleUnifDirichletG(L, ploidy, K),
#'                                                 Q = SampleUnifQ(n, K),
#'                                                 coord = SampleNormalClusterCoord(n.by.pop = n, K = 1),
#'                                                 ploidy = ploidy)

#' # mask data
#' set.seed(0)
#' masked.prop <- 0.1
#' masked.X <- data.list$X
#' masked.X[sample(1:(ncol(masked.X)*nrow(masked.X)), (ncol(masked.X)*nrow(masked.X)) * masked.prop)] <- NA
#'
#'
#' # run tess3
#' tess3.obj <- tess3(X = masked.X,
#'                    coord = data.list$coord,
#'                    K = 3,
#'                    ploidy = 2,
#'                    lambda = 1.0,
#'                    method = "projected.ls",
#'                    rep = 2)
#'
#' imputed.X.random <- ImputeRandom(Gettess3res(tess3.obj,3), masked.X)
#' mean(data.list$X[is.na(masked.X)] != imputed.X.random[is.na(masked.X)]) # pourcentage of error
#'
#' imputed.X.round <- ImputeRound(tess3.res = Gettess3res(tess3.obj,3), masked.X)
#' mean(data.list$X[is.na(masked.X)] != imputed.X.round[is.na(masked.X)]) # pourcentage of error
#'
#'
ImputeRandom <- function(tess3.res, masked.X) {
  if (!is.tess3Main(tess3.res)) {
    stop("tess3.res must be a tess3 result of class tess3Main. You can use function tess3Main or Gettess3res on a tess3 object.")
  }
  ploidy <- max(masked.X, na.rm = TRUE)
  Prob <- tess3.res$Q %*% t(tess3.res$G)
  Prob <- array(Prob, c(nrow(Prob), ploidy + 1, ncol(Prob) / (ploidy + 1)))

  for (i in 1:nrow(masked.X)) {
    for (j in 1:ncol(masked.X)) {
      if (is.na(masked.X[i, j])) {
        masked.X[i, j] <- sample(0:ploidy, 1, replace = FALSE, prob = Prob[i,,j])
      }
    }
  }
  return(masked.X)
}

#' \code{ImputeRound} imputes missing values with the Q ancestry coeficient matrix and G ancestral
#' frequencies estimates. Missing values are computed as a round of estimated genotype
#' frequencies (P = Q G^T).
#'
#'
#' @export
#' @rdname ImputeRandom
ImputeRound <- function(tess3.res, masked.X) {
  if (!is.tess3Main(tess3.res)) {
    stop("tess3.res must be a tess3 result of class tess3Main. You can use function tess3Main or Gettess3res on a tess3 object.")
  }
  ploidy <- max(masked.X, na.rm = TRUE)
  P <- tess3.res$Q %*% t(GtoFreq(tess3.res$G, ploidy))
  masked.id <- which(is.na(masked.X))
  masked.X[masked.id] <- round(P[masked.id])
  return(masked.X)
}
