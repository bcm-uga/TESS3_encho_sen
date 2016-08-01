ComputeGraphLaplacian <- function(W) {
  D <- diag(apply(W,1,sum))
  return(D - W)
}

ComputeEigenValuesWithIgraphArpack <- function(Lapl, k) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph needed for ComputeEigenValuesWithIgraph function to work. Please install it.",
         call. = FALSE)
  }

  f2 <- function(x, extra=NULL) { cat("."); as.vector(Lapl %*% x) }
  baev <- igraph::arpack(f2, sym=TRUE, options=list(n=nrow(Lapl), nev = k, ncv = min(nrow(Lapl), max(2*k+1, 20)),
                                            which="LM", maxiter=200))

  return(baev)
}

ComputeEigenValuesWithRSpectra <- function(Lapl, k) {
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("RSpectra needed for ComputeEigenValuesWithRSpectra function to work. Please install it.",
         call. = FALSE)
  }

  res <- RSpectra::eigs_sym(Lapl, k, which = "LM")

  return(res)
}

ComputeMeanDist <- function(coord) {
  W <- matrix(0,nrow(coord),nrow(coord))
  for (i in 1:nrow(coord)) {
    for (j in 1:nrow(coord)){
      W[i,j] <- sqrt(sum((coord[i,]-coord[j,])^2))
    }
  }
  return(mean(W))
}
