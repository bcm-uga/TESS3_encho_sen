#' Estimates spatial population structure
#'
#' \code{tess3Main} estimates spatial population structure using a graph based non
#' negative matrix factorization. After estimating the population structure is used to
#' compute a Fst statistic for each locus. See references for more details.
#'
#' @param K An integer which corresponds to
#' the number of ancestral populations.
#' @param ploidy An integer which corresponds to
#' the number of set of chromosomes.
#' @param lambda A numeric which corresponds to the
#' spatial regularization parameter.
#' @param W A numeric matrix which corresponds to the graph weiht matrix.
#' If NULL, W it is computed as
#' \code{W[i,j] = exp( - (coord[i] - coord[j])^2 / sigma^2)}.
#' Where \code{coord[i]} is
#' the geographic coordinate for the individual i and
#' \code{sigma} equals 5 percent of the average geographic distance between individuals.
#' @param method \code{"projected.ls"} or \code{"qp"}. If \code{"projected.ls"},
#' an aleternated projected least squares algorithm is used. If \code{"qp"},
#' an alternated quadratic programing algorithm is used. See references for more
#' details
#' @param max.iteration The max number of
#' iteration of the optimization algorithm.
#' @param tolerance A numeric which corresponds to the tolerance of the
#' stopping criteria of the optimization algorithm.
#' @param X A numeric matrix which corresponds to the genotype matrix. This matrix
#' must be of size \eqn{n \times L} where \eqn{n} is the number of individual and
#' \eqn{L} is the number of loci. Values of this matrix are integer corresponding
#' to the number of variant alleles observed at a locus. If \code{NULL}, \code{XProba}
#' is used.
#' @param openMP.core.num If openMP is available on your computer, it is the
#' number of core used by the algorithm.
#' @param Q.init A numeric matrix which corresponds to the initial value of
#' \code{Q} for the algorithm.
#' @param coord The numeric matrix of size \eqn{n \times 2} where \eqn{n} is the
#' number of individuals.
#' @param mask If not \code{NULL} this the proportion of the data matrix which
#' is masked to compute the cross validation criteria.
#' @param algo.copy if TRUE data will be copy 1 time to speed the algorithm.
#' @param copy if TRUE data will be copy 1 time.
#' @param verbose If \code{TRUE} more information are printed.
#' @param XProba A numeric matrix which correspond to the probability for the genotype.
#' This matrix must be of size \eqn{n \times (ploidy + 1)L} where
#' \eqn{n} is the number of individual, \eqn{L} is the number of loci. Values of
#' this matrix are numeric between 0 and 1 corresponding
#' to the genome probability. It is the matrix used in graph based non negative
#' factorization matrix. If \code{NULL}, it is computed from the genotype matrix \code{X}.
#' See reference for more details.
#'
#' @return An object of class tess3Main which is a list with components:
#' \describe{
#'    \item{L}{The number of loci.}
#'    \item{n}{The number of individuals.}
#'    \item{ploidy}{The number of set of chromosomes.}
#'    \item{K}{The number of ancestral population.}
#'    \item{G}{The ancestral genotype frequency matrix.}
#'    \item{Q}{The ancestry coefficient matrix.}
#'    \item{Fst}{The Fst statistic computed for each locus.}
#'    \item{Fscore}{The Fscore computed from Fst.}
#'    \item{pvalue}{The pvalues computed from Fscore.}
#'    \item{log.pvalue}{The \eqn{log(pvalue)}.}
#'    \item{rmse}{The root square mean error between \code{XProba} and
#'    \code{tcrossprod(Q, G)}.}
#'    \item{crossentropy}{The cross entropy error between \code{XProba} and
#'    \code{tcrossprod(Q, G)}.}
#'   \item{crossvalid.rmse}{If masked not \code{NULL}.
#'   The root square mean error between \code{XProba[masked]} and
#'   \code{tcrossprod(Q, G)[masked]}.}
#'    \item{crossvalid.crossentropy}{If masked not \code{NULL}.
#'    The cross entropy error between \code{XProba[masked]} and
#'    \code{tcrossprod(Q, G)[masked]}.}
#' }
#'
#'
#' @export
#' @examples
#' library(tess3r)
#'
#' # Arabidopsis thaliana data set
#' data(data.at)
#' genotype <- data.at$X
#' coordinates <- data.at$coord
#'
#' # Run of tess3 main algorithm
#' tess3.obj <- tess3Main(X = genotype,
#'                       coord = coordinates,
#'                       K = 3,
#'                       method = "projected.ls",
#'                       ploidy = 1)
#'
#' # Run of tess3 main algorithm with cross validation errors computation.
#' tess3.obj <- tess3Main(X = genotype,
#'                       coord = coordinates,
#'                       K = 3,
#'                       method = "projected.ls",
#'                       ploidy = 1,
#'                       mask = 0.05)
#'
#'
#' @references \url{http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12471/full}
#' @seealso \code{\link{tess3}}
tess3Main <- function(X,
                      XProba = NULL,
                      coord,
                      K,
                      ploidy,
                      lambda = 1.0,
                      W = NULL,
                      method = "projected.ls",
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
  } else if (!copy & is.null(XProba)) {
    stop("To force the function not doing copy of the data, you must set XProba")
  }
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # check type of input

  ## geographic.coordinate
  if (is.null(coord) && is.null(W)) {
    stop("If no graph matrix is specified, geographic coordinates must be provided as coord parameter")
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
    stop("method must be a character string for the method to use")
  }
  if (method != "projected.ls" & method != "qp") {
    stop("method must be projected.ls or qp")
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

  # check coord
  CheckCoord(coord)

  # Compute W
  if (is.null(W)) {
    W <- ComputeHeatKernelWeight(coord, ComputeMeanDist(coord) * 0.05)
  }


  # Q.init
  if (!is.null(Q.init)) {
    if (nrow(Q.init) != nrow(X) | ncol(Q.init) != K) {
      stop("Q.init must be of dimensions (number of individuals) x K.")
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
    res$n <- nrow(X)
  } else if (!is.null(XProba)) {
    if (ncol(XProba) %% (ploidy + 1) != 0) {
      stop("Number of columns of XProba must be a multiple of ploidy + 1.")
    }
    res$L <- ncol(XProba) %/% (ploidy + 1)
    res$n <- nrow(coord)
  } else {
    stop("X or XProba must be non-null")
  }

  res$ploidy <- ploidy
  res$K <- K

  if (is.null(XProba)) {
    XProba <- matrix(0.0, res$n, res$L * (res$ploidy + 1))
    X2XBin(X, ploidy, XProba)
    rm(X)
  }
  CheckXBin(XProba, ploidy)
  CheckCoord(coord)
  CheckXBinCoord(XProba, coord)
  CheckW(W)
  CheckWCoord(W, coord)
  CheckXBinW(XProba, W)
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # mask if asked
  if (mask != 0.0) {
    message("Mask ", mask, "% of genotypes for cross validation")
    missing.index.X <- sample(1:(length(XProba)), length(XProba) * mask)
    masked.X.value <- XProba[missing.index.X]
    XProba[missing.index.X] <- NA
  }
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # check if there is missing data and compute the binary representation
  if (any(is.na(XProba))) {
    message("Missing value detected in genotype")
    # Naive imputation by mean
    geno.freq <- apply(XProba, 2, function(x) mean(x,na.rm = TRUE))
    geno.freq <- matrix(geno.freq,1)[rep(1, res$n),]
    XProba[is.na(XProba)] <- geno.freq[is.na(XProba)]
    rm(geno.freq)
  }
  ################################################

  # mem <- c(mem,pryr::mem_used())

  ################################################
  # compute Q and G matrix
  if (method == "qp") {
    res <- c(res, SolveTess3QP(XProba,
                               K,
                               ploidy,
                               Lapl,
                               lambda,
                               max.iteration = max.iteration,
                               tolerance = tolerance,
                               verbose = verbose))
  } else if (method == "projected.ls") {
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
      ComputeMCPASolutionNoCopyX(X = XProba,
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
      ComputeMCPASolution(X = XProba,
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
  res$rmse <- ComputeRmse(XProba, QtG)
  res$crossentropy <- ploidy * ComputeAveragedCrossEntropy(XProba, QtG) #because this function compute mean also by allele
  if (mask > 0.0) {
    res$crossvalid.rmse <- ComputeRmse(masked.X.value, QtG[missing.index.X])
    res$crossvalid.crossentropy <- ploidy * ComputeAveragedCrossEntropy(masked.X.value, QtG[missing.index.X])
  }
  ################################################
  class(res) <- c(class(res), "tess3Main")

  # mem <- c(mem,pryr::mem_used())

  return(res)
}


#' Summary of tess3 object.
#'
#' @param object TODOC
#' @param ... TODOC
#'
#' @return TODOCC
#' @export
summary.tess3Main <- function(object, ...) {
  cat(paste("=== Object of class tess3Main ===\n"))
  cat(paste("Number of individuals (n):", object$n,"\n"))
  cat(paste("Number of loci (L):", object$L,"\n"))
  cat(paste("Ploidy:", object$ploidy,"\n"))
  cat(paste("Number of ancestral populations (K):", object$K,"\n"))
  cat(paste("Residual error:", object$rmse,"\n"))
}


#' Title
#'
#' @param x TODOC
  #'
#' @return TODOC
#' @export
is.tess3Main <- function(x) {
  inherits(x, "tess3Main")
}

#' Title
#'
#' @param tess3.obj TODOCC
#' @param X TODOC
#' @param ploidy TODOC
#' @param mask TODOC
#'
#' @export
rmse.tess3Main <- function(tess3.obj, X, ploidy, mask = NULL) {
  if (!is.tess3Main(tess3.obj)) {
    stop("tess3.obj must of class tess3Main.")
  }
  CheckX(X, ploidy)
  if (!is.null(mask)) {
    X[-mask] <- NA
  }
  XBin <- matrix(0.0, nrow(X), ncol(X) * (ploidy + 1))
  X2XBin(X, ploidy, XBin)
  return(ComputeRmse(XBin, tcrossprod(tess3.obj$Q, tess3.obj$G)))
}

#' Title
#'
#' @param x tess3Main object.
#' @param ... TODOC
#'
#' @export
plot.tess3Main <- function(x, ...) {
  message("Nothing to plot")
}

#' tess3r : estimation of spatial population structure
#'
#' This R package implements the TESS3 method and tools useful to plot program outputs.
#'
#' @docType package
#'
#' @name tess3r
#' @importFrom Rcpp evalCpp
#' @importFrom graphics barplot hist image par plot points segments
#' @importFrom grDevices colorRampPalette rainbow topo.colors
#' @importFrom stats as.formula dist median pchisq pf predict qchisq qf rnorm runif sd
#' @importFrom utils capture.output
#' @import RcppEigen
#' @useDynLib tess3r
NULL

#' TODOC
#'
#' @name data.at
#' @docType data
#' @keywords data
NULL

#' TODOC
#'
#' @name data.for.test
#' @docType data
#' @keywords data
NULL

#' TODOC
#'
#' @name durand09
#' @docType data
#' @keywords data
NULL
