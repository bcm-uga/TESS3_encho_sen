#' Main function.
#'
#' @param K
#' @param ploidy
#' @param lambda
#' @param W
#' @param method
#' @param max.iteration
#' @param tolerance
#' @param X
#' @param openMP.core.num
#' @param Q.init
#' @param coord
#' @param mask
#' @param XBin
#' @param algo.copy if TRUE data will be copy to speed the algorithm
#'
#' @return
#' @export
#'
#' @examples
tess3Main <- function(X,
                      XBin = NULL,
                      coord,
                      K,
                      ploidy,
                      lambda = 1.0,
                      W = NULL,
                      method = "MCPA",
                      max.iteration = 200,
                      tolerance = 1e-5,
                      openMP.core.num = 1,
                      Q.init = NULL,
                      mask = 0.0,
                      copy = TRUE,
                      algo.copy = TRUE,
                      verbose = FALSE)
{
  # mem <- c()
  # mem <- c(mem, pryr::mem_used())
  res = list()

  ################################################
  # Init openMP
  InitOpenMP(openMP.core.num)

  ################################################

  # mem <- c(mem, pryr::mem_used())

  ################################################
  # copy X
  if (copy & !is.null(X)) {
    X <- matrix(as.double(X), nrow(X), ncol(X))
    CheckX(X, ploidy)
  } else if (!copy & is.null(XBin)) {
    stop("To force the function not doing copy of the data, you must set XBin.")
  }
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # check type of input

  ## geographic.coordinate
  if (is.null(coord) && is.null(W)) {
    stop("If no graph-weight matrix is specified, geographic coordinates must be provided as coord parameter")
  }

   ## K
  if (!is.numeric(K)) {
    stop("K must be numeric")
  }
  ## ploidy
  if (!is.numeric(ploidy)) {
    stop("ploidy must be numeric")
  }
  ## Q.init
  if (!is.null(Q.init)) {
    if (!is.matrix(Q.init) ) {
      stop("Q.init must be a matrix")
    }
  }
  ## method
  if (!is.character(method)) {
    stop("method must be the name of the method to use")
  }
  ## max.iteration
  if (!is.numeric(max.iteration)) {
    stop("max.iteration must be numeric")
  }
  ## tolerance
  if (!is.numeric(tolerance)) {
    stop("tolerance must be numeric")
  }

  if (mask > 1.0 | mask < 0.0) {
    stop("mask is the proportion of the genotype masked for cross validation and must between 0 and 1")
  }

  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # Check consistence of input

  # Compute W
  if (is.null(W)) {
    W <- ComputeHeatKernelWeight(coord, ComputeMeanDist(coord) * 0.05)
  }


  # Q.init
  if (!is.null(Q.init)) {
    if (nrow(Q.init) != nrow(X) | ncol(Q.init) != K) {
      stop("Q.init must be of size n * K")
    }
  }

  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # Compute parameters
  ## Laplacian
  CheckW(W)
  Lapl <- ComputeGraphLaplacian(W)

  ## Compute number of loci and indiv
  if (!is.null(X)) {
    res$L <- ncol(X)
  } else if (!is.null(XBin)) {
    if (ncol(XBin) %% (ploidy + 1) != 0) {
      stop("Number of columns of XBin must be a multiple of ploidy + 1.")
    }
    res$L <- ncol(XBin) %/% (ploidy + 1)
  } else {
    stop("Both X and XBin can not be null")
  }

  res$n <- nrow(coord)
  res$ploidy <- ploidy
  res$K <- K

  if (is.null(XBin)) {
    XBin <- matrix(0.0, res$n, res$L * (res$ploidy + 1))
    X2XBin(X, ploidy, XBin)
    rm(X)
  }
  CheckXBinWCoord(XBin, ploidy, W, coord)
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # mask if asked
  if (mask != 0.0) {
    message("Mask ", mask, "% of X for cross validation")
    missing.index.X <- sample(1:(length(XBin)), length(XBin) * mask)
    masked.X.value <- XBin[missing.index.X]
    XBin[missing.index.X] <- NA
  }
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # check if there is missing data and compute the binary representation
  if (any(is.na(XBin))) {
    message("Missing value detected in genotype")
    # Naive imputation by mean
    geno.freq <- apply(XBin, 2, function(x) mean(x,na.rm = TRUE))
    geno.freq <- matrix(geno.freq,1)[rep(1, res$n),]
    XBin[is.na(XBin)] <- geno.freq[is.na(XBin)]
    rm(geno.freq)
  }
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # compute Q and G matrix
  if (method == "OQA") {
    res <- c(res, SolveTess3QP(XBin,
                               K,
                               ploidy,
                               Lapl,
                               lambda,
                               max.iteration = max.iteration,
                               tolerance = tolerance,
                               verbose = verbose))
  } else if (method == "MCPA") {
    Lapl <- as.matrix(Lapl)
    # Q and G
    res$G <- matrix(0.0, nrow = (res$ploidy + 1) * res$L, ncol = res$K)
    if (!is.null(Q.init)) {
      res$Q <- Q.init
    } else {
      res$Q <- matrix(runif(res$n * K), res$n, res$K)
      res$Q <- ProjectQ(res$Q)
    }

    # mem <- c(mem,pryr::mem_used())
    if (!algo.copy) {
      ComputeMCPASolutionNoCopyX(X = XBin,
                                 K = K,
                                 Lapl = Lapl,
                                 lambdaPrim = lambda,
                                 D = ploidy + 1,
                                 maxIteration = max.iteration,
                                 tolerance = tolerance,
                                 Q = res$Q,
                                 G = res$G,
                                 verbose = verbose)
    } else {
      ComputeMCPASolution(X = XBin,
                          K = K,
                          Lapl = Lapl,
                          lambdaPrim = lambda,
                          D = ploidy + 1,
                          maxIteration = max.iteration,
                          tolerance = tolerance,
                          Q = res$Q,
                          G = res$G,
                          verbose = verbose)
    }

  } else {
    stop("Unknow method name")
  }
  class(res$Q) <- c("tess3Q", class(res$Q))
  class(res$G) <- c("tess3G", class(res$G))

  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # Selection scan
  # Compute stat for test
  if (K > 1) {
    res$Fst <- ComputeFst(res$Q, res$G, ploidy + 1)

    # To avoid numerical issues
    res$Fst[res$Fst < 0.0] = 0.0

    # Convert Fst into t score
    res <- c(res, ComputeTscoreAndPvalue(res$Fst, K, res$n))
    class(res$pvalue) <- c("tess3pvalue", class(res$pvalue))
  }

  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # Compute rmse
  QtG <- tcrossprod(res$Q, res$G)
  res$rmse <- ComputeRmse(XBin, QtG)
  res$crossentropy <- ploidy * ComputeAveragedCrossEntropy(XBin, QtG) #because this function compute mean also by allele
  if (mask > 0.0) {
    res$crossvalid.rmse <- ComputeRmse(masked.X.value, QtG[missing.index.X])
    res$crossvalid.crossentropy <- ploidy * ComputeAveragedCrossEntropy(masked.X.value, QtG[missing.index.X])
  }
  ################################################
  class(res) <- c("tess3Main", class(res))

  # mem <- c(mem,pryr::mem_used())

  return(res)
}


#' Summary of tess3 object.
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.tess3 <- function(object, ...) {
  cat(paste("=== Object of class tess3 ===\n"))
  cat(paste("Number of individuals n:", object$n,"\n"))
  cat(paste("Number of loci L:", object$L,"\n"))
  cat(paste("Ploidy:", object$ploidy,"\n"))
  cat(paste("Number of ancestral populations K:", object$K,"\n"))
  cat(paste("RMSE(genotype, Q * G^T):", object$rmse,"\n"))
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.tess3Main <- function(x) {
  inherits(x, "tess3Main")
}

#' Title
#'
#' @param tess3.obj
#' @param genotype
#' @param ploidy
#'
#' @return
#' @export
#'
#' @examples
rmse.tess3Main <- function(tess3.obj, X, ploidy, mask = NULL) {
  if (!is.tess3Main(tess3.obj)) {
    stop("tess3.obj must of class tess3Main")
  }
  CheckX(X, ploidy)
  if (!is.null(mask)) {
    X[-mask] <- NA
  }
  XBin <- matrix(0.0, nrow(X), ncol(X) * (ploidy + 1))
  X2XBin(X, ploidy, XBin)
  return(ComputeRmse(XBin, tcrossprod(tess3.obj$Q, tess3.obj$G)))
}

#' tess3r : estimation of spatial population structure
#'
#' This R package implements the TESS3 method and tools useful to plot program outputs.
#'
#' @docType package
#'
#' @name tess3r
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' @useDynLib tess3r
NULL
