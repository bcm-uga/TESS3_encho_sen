ComputeHeatKernelWeight <- function(coord,sigma) {
  W <- matrix(0,nrow(coord),nrow(coord))
  for(i in 1:nrow(coord)) {
    for(j in 1:nrow(coord)){
      W[i,j] <- sqrt(sum((coord[i,]-coord[j,])^2))
    }
  }
  if(is.null(sigma)) {
    sigma <- 0.05*mean(W)
    cat("sigma = ",sigma)
  }
  W <- apply(W, c(1,2), function(d){exp(-d^2/sigma^2)})
  return(W)
}

ComputeGraphLaplacian <- function(W) {
  D <- diag(apply(W,1,sum))
  return(D - W)
}
