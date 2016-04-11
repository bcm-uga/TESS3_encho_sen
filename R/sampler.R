#' sample locus position such lp ~ unif(0,1)
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleUnifLocusPosition <- function() {
  function(n, n.by.pop, L, d, K) {
    lp = runif(L)
    return(sort(lp))
  }
}
#' sample Q such as alpha ~ unif
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleUnifQ <- function() {
  function(coord, n, n.by.pop, L, d, K) {
    alpha = matrix(runif(n*K,min = 0, max = 1),n,K)
    Q = t(apply(alpha, 1, function(r){gtools::rdirichlet(1,r)}))
    return(Q)
  }
}

#' sample coord such as a mixture of K cluster distributed with gaussian law
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleNormalClusterCoord <- function() {
  function(n, n.by.pop, L, d, K) {
    # sample pop center
    mu = MASS::mvrnorm(K, c(0,0), diag(c(1,1)))

    # sample coord
    if(length(n.by.pop) == 1) {
      n.by.pop = rep(n.by.pop,K)
    }
    assertthat::assert_that(length(n.by.pop) == K)
    coord=c()
    for(k in 1:K) {
      coord = rbind(coord, MASS::mvrnorm(n.by.pop[k], mu[k,], diag(c(0.2,0.2)) ))
    }
    # debug
    # plot(coord, col = rep(2:(K+1),times = n.by.pop))
    return(coord)
  }
}

#' sample G such as uniform diriclet
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleDirichletG <- function() {
  function(n, n.by.pop, L, d, K) {
    # sample G
    G = c()
    for( l in 1:L) {
      G = rbind(G,t(gtools::rdirichlet(K,rep(1.0/(d+1),(d+1)))))
    }
    return(G)
  }
}


#' exp(-||X-Y||^2 / 2 * sigma^2)
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
CovFunctionSquaredExp <- function(sigma) {
  function(X,Y){
    return(exp( - sum((X-Y)^2) / (2 * sigma) ))
  }
}

#' sample Q such Q = logit(a)
#' a ~ N(0,W) W_ij = cov.function(X_i,X_j)
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleGaussianProcessQ <- function(cov.function) {
  function(coord, n, n.by.pop, L, d, K) {
    # compute cov
    W = matrix(0,nrow(coord),nrow(coord))
    for(i in 1:nrow(coord)) {
      for(j in 1:nrow(coord)){
        W[i,j] = cov.function(coord[i,], coord[j,])
      }
    }

    alpha = t(MASS::mvrnorm(n = K, rep(0,n),W))
    alpha = boot::inv.logit(alpha)
    Q = t(apply(alpha, 1, function(r){gtools::rdirichlet(1,r)}))
    return(Q)
  }
}

#' gaussian function
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleDistDFromCenterQ <- function(sigma) {
  function(coord, n, n.by.pop, L, d, K) {

    # compute center of each cluster
    mu <- t(sapply(split(1:n, ceiling(seq_along(1:n)/n.by.pop)), function(l) {apply(coord[l,],2,mean)}))

    # sample Q
    dist.from.center.aux <- function(c) {
      apply(mu, 1, function(r){return(sqrt(sum((r-c)^2)))})
    }
    dist.from.center = t(apply(coord, 1, dist.from.center.aux))
    dist.from.center = exp(-dist.from.center/sigma)
    Q = dist.from.center
    for(i in 1:nrow(Q)){
      s = sum(Q[i,])
      if(s != 0.0) {
        Q[i,] = Q[i,] / s
      }
    }
    return(Q)
  }
}

#' sample pop center and compute dist.from.center as dirichlet parameter
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleDistFromCenterDirichletQ <- function() {
  function(coord, n, n.by.pop, L, d, K) {

    # compute center of each cluster
    mu <- t(sapply(split(1:n, ceiling(seq_along(1:n)/n.by.pop)), function(l) {apply(coord[l,],2,mean)}))

    # sample Q
    dist.from.center.aux <- function(c) {
      apply(mu, 1, function(r){return(sqrt(sum((r-c)^2)))})
    }
    dist.from.center = t(apply(coord, 1, dist.from.center.aux))
    dist.from.center = exp(-dist.from.center)
    #   W = matrix(0,nrow(coord),nrow(coord))
    #   for(i in 1:nrow(coord)) {
    #     for(j in 1:nrow(coord)){
    #       W[i,j] = sqrt(sum((coord[i,]-coord[j,])^2))
    #     }
    #   }
    #   theta = mean(W)
    #   W = apply(W, c(1,2), function(d){exp(-d/theta)})
    #   D = diag(apply(W,1,sum))
    #   Lap = D - W
    #   alpha = apply(dist.from.center, 2, function(c){MASS::mvrnorm(1,c,Lap)})
    #   alpha = apply(alpha,1:2, function(e){ e-min(alpha)})
    Q = t(apply(dist.from.center, 1, function(r){gtools::rdirichlet(1,r)}))
    # debug
    #   library(tess3r)
    #   asciiFile=system.file("extdata/","lowResEurope.asc",package = "tess3r")
    #   lat.pix=seq(from=-3,by=0.01,length=600)
    #   long.pix=seq(from=-3,by=0.01,length=600)
    #   grid=fields::make.surface.grid( list( long.pix,lat.pix))
    #   constraints=matrix(TRUE, 600,600)
    #   maps(matrix = Q,
    #        coord = coord,
    #        grid=grid,constraints=constraints,method="max",main="Ancestry Coefficient with K = 3")

    return(Q)
  }
}



#' Sample X such that P(X_i_dl + j) = Sum(Q_ik G_kj).
#' Coord ...
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleGenoFromGenerativeModel <- function(n.by.pop, L, d, K,
                                          G.sampler = SampleDirichletG(),
                                          Q.sampler = SampleUnifQ(),
                                          coord.sampler = SampleNormalClusterCoord(),
                                          locus.position.sampler = SampleUnifLocusPosition()) {
  n = K * n.by.pop

  # sample locus position
  locus.position = locus.position.sampler(n, n.by.pop, L, d, K)

  # sample coord
  coord = coord.sampler(n, n.by.pop, L, d, K)

  # sample Q
  Q = Q.sampler(coord, n, n.by.pop, L, d, K)

  # sample G
  G = G.sampler(n, n.by.pop, L, d, K)

  # P
  P = tcrossprod(Q,G)

  # X
  X = matrix(0,nrow(Q),L)
  allele = 0:(d)
  for (i in 1:nrow(Q)) {
    for( j in 1:L) {
      X[i,j] = allele %*% rmultinom(1,1, P[i,1 + ((j-1)*(d+1)):((j)*(d+1)-1)])
    }
  }

  return( list(X = X,G = G,Q = Q, coord = coord,n = nrow(X),L = L,n.by.pop = n.by.pop,K = K, d = d, locus.position = locus.position))
}
