#' Compute the value spatial penalty used in TESS3 objectif function.
#'
#' @param Q
#' @param W
#'
#' @return
#' @export
#'
#' @examples
ComputeSpatialPenalty <- function(Q, W) {
  Lapl <- as.matrix(ComputeGraphLaplacian(W))
  aux <- crossprod(Lapl, Q)
  aux <- crossprod(aux, Q)
  return(sum(diag(aux)))
}


#' Compute the root mean square error between matrix by permuting
#' matrix column such that we have the best rmse.
#'
#'
#'
#' @param Q1
#' @param Q2
#'
#' @export
ComputeRmseWithBestPermutation <- function(Q1, Q2) {
  TestRequiredPkg("permute")

  aux = ComputeRmse(Q1,Q2)

  K = dim(Q1)[2]
  #Because of R !!!!
  if ( K == 2 ) {
    perms = matrix(c(2,1),nrow = 1,ncol = 2)
  } else {
    perms = permute::allPerms(K)
  }
  for (i in 1:(dim(perms)[1])) {

    aux1 = ComputeRmse(Q1,Q2[,perms[i,]])
    if (aux1 < aux) {
      aux = aux1
    }

  }
  return(aux)
}

#' Compute the root mean square error between matrix by permuting
#' matrix column such that we have the best rmse. Warning : Greedy algorithm !
#'
#'
#'
#' @export
ComputeRmseWithBestPermutationGreedy <- function(Q1, Q2) {

  # find perm
  perm = 1:ncol(Q1)

  taken = c()

  for (i in 1:ncol(Q1)) {
    min = .Machine$double.xmax
    if (i > 1) {
      taken = c(taken,perm[i - 1])
    }
    for (j in 1:ncol(Q1)) {
      aux = ComputeRmse(Q1[,i],Q2[,j])
      if (aux < min && !(j %in% taken)) {
        perm[i] = j
        min = aux
      }
    }
  }

  if (sum(duplicated(perm)) != 0) {
    stop("this is not a permutation")
  }
  return(ComputeRmse(Q1,Q2[,perm]))
}




