#' Main function.
#'
#' @param genotype
#' @param geographic.coordinate
#' @param K
#' @param ploidy
#' @param lambda
#' @param W
#' @param method
#' @param max.iteration
#' @param tolerance
#'
#' @return
#' @export
#'
#' @examples
tess3 <- function(genotype,
                  geographic.coordinate,
                  K,
                  ploidy,
                  lambda,
                  W = NULL,
                  method = "MCPA",
                  max.iteration = 200,
                  tolerance = 1e-5,
                  openMP.core.num = 1,
                  Q.init = NULL)
{
  res = list()

  ################################################
  # Init openMP
  InitOpenMP(openMP.core.num)

  ################################################

  ################################################
  # check type of the input
  ## genotype
  if (!is.matrix(genotype)) {
    stop("genotype must be a matrix")
  }
  ## geographic.coordinate
  if (is.null(geographic.coordinate) && is.null(W)) {
    stop("If no W graph weight matrix is specified, geographic.coordinate is mandatory")
  }
  if (!is.matrix(geographic.coordinate)) {
    stop("geographic.coordinate must be a matrix")
  }

   ## K
  if (!is.numeric(K)) {
    stop("K must be numeric")
  }
  ## ploidy
  if (!is.numeric(ploidy)) {
    stop("ploidy must be numeric")
  }
  ## W
  if (!is.null(W)) {
    if (!is.matrix(W) && !(attr(class(W),"pacakge") == "Matrix")) {
      stop("W must be a matrix. Use R base matrix type or package Matrix")
    }
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

  ################################################


  ################################################
  # Check consistence of input

  # To be sure we have an double coordinate matrix
  geographic.coordinate <- matrix(as.double(geographic.coordinate),
                                  nrow(geographic.coordinate),
                                  ncol(geographic.coordinate))

  genotype <- matrix(as.integer(genotype),
                     nrow(genotype),
                     ncol(genotype))

  ## check dim of geno, coord,  are consistent
  if (nrow(genotype) != nrow(geographic.coordinate)) {
    stop("Number of row in the coordinate matrix and the genotype matrix must be the same")
  }
  ## check if rang of genotype matrix
  if (max(genotype, na.rm = TRUE) > ploidy) {
    stop("The maximum value of the genotype matrix can not be superior than ploidy + 1")
  }
  if (min(genotype, na.rm = TRUE) < 0) {
    stop("Negative values in the genotype matrix are not allowed")
  }

  # Compute W
  if (is.null(W)) {
    W <- ComputeHeatKernelWeight(geographic.coordinate, ComputeMeanDist(geographic.coordinate) * 0.05)
  }
  if (is.matrix(W)) {
    # convert into sparse matrix
    W <- Matrix::Matrix(W, sparse = TRUE)
  }

  # check W
  if (nrow(W) != ncol(W) || nrow(W) != nrow(genotype) || !Matrix::isSymmetric(W)) {
    stop("W must be squared symetric of size nrow(genotype) x nrow(genotype) ")
  }

  # Q.init
  if (!is.null(Q.init)) {
    if (nrow(Q.init) != nrow(genotype) | ncol(Q.init) != K) {
      stop("Q.init must be of size nbIndiv * K")
    }
  }

  ################################################

  ################################################
  # Compute parameters
  ## Laplacian
  Lapl <- ComputeGraphLaplacian(W)

  ## Compute number of loci and indiv
  res$L <- ncol(genotype)
  res$n <- nrow(genotype)
  res$ploidy <- ploidy
  res$K <- K

  ################################################
  # check if there is missing data and compute the binary representation
  missing <- FALSE
  if (any(is.na(genotype))) {
    message("Missing value detected in genotype")
    missing <- TRUE
    genotype[is.na(genotype)] <- as.integer(-1)
  }
  genotype <- ComputeXBin(genotype, ploidy)
  # To be sure we have an double genotype matrix
  genotype <- matrix(as.double(genotype), nrow(genotype), ncol(genotype))
  if (missing) {
    # Naive imputation by mean
    genotype[genotype < 0] <- NA
    geno.freq <- apply(genotype, 2, function(x) mean(x,na.rm = TRUE))
    geno.freq <- matrix(geno.freq,1)[rep(1, res$n),]
    genotype[is.na(genotype)] <- geno.freq[is.na(genotype)]
  }
  ################################################

  ################################################
  # compute Q and G matrix
  if (method == "OQA") {
    res <- c(res, SolveTess3QP(genotype,
                        K,
                        ploidy,
                        Lapl,
                        lambda,
                        max.iteration = max.iteration,
                        tolerance = tolerance))
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

    ComputeMCPASolution(X = genotype,
                        K = K,
                        Lapl = Lapl,
                        lambdaPrim = lambda,
                        D = ploidy + 1,
                        maxIteration = max.iteration,
                        tolerance = tolerance,
                        Q = res$Q,
                        G = res$G )
  } else {
    stop("Unknow method name")
  }
  class(res$Q) <- "tess3Q"
  class(res$G) <- "tess3G"

  ################################################

  ################################################
  # Selection scan
  # Compute stat for test
  if (K > 1) {
    res$Fst <- ComputeFst(res$Q, res$G, ploidy + 1)

    # To avoid numerical issues
    res$Fst[res$Fst < 0.0] = 0.0

    # Convert Fst into t score
    res <- c(res, ComputeTscoreAndPvalue(res$Fst, K, res$n))
    class(res$pvalue) <- "tess3pvalue"
  }

  ################################################

  ################################################
  # Compute rmse
  res$rmse <- ComputeRmse(genotype, res$Q %*% t(res$G))
  ################################################
  class(res) <- "tess3"
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
is.tess3 <- function(x) {
  inherits(x, "tess3")
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
rmse.tess3 <- function(tess3.obj, genotype, ploidy) {
  if (!is.tess3(tess3.obj)) {
    stop("tess3.obj must of class tess3")
  }
  if (typeof(genotype) != "integer") {
    genotype <- matrix(as.integer(genotype), nrow(genotype), ncol(genotype))
  }
  genotype <- ComputeXBin(genotype, ploidy)
  return(ComputeRmse(genotype, tcrossprod(tess3.obj$Q, tess3.obj$G)))
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
