#' Estimate ancestry coefficients and run genome scans for selection
#'
#' \code{tess3Main} estimates spatial population structure using a graph based non
#' negative matrix factorization. After estimating the population structure is used to
#' compute a Fst statistic for each locus. See references for more details.
#'
#' @param K an integer corresponding to
#' the number of ancestral populations.
#' @param ploidy an integer which corresponds to the ploidy of the number of copy of chromosomes.
#' @param lambda a nonnegative numeric which corresponds to the
#' spatial regularization parameter.
#' @param W a numeric matrix which corresponds to the graph weight matrix.
#' If \code{NULL}, W is computed as
#' \code{W[i,j] = exp(-(coord[i] - coord[j])^2 / sigma^2)},
#' where \code{coord[i]} is the set of geographic coordinates for individual i and
#' \code{sigma} equals 5 percent of the average geographic distance between individuals.
#' @param method \code{"projected.ls"} or \code{"qp"}. If \code{"projected.ls"},
#' an alternating projected least squares algorithm is used. If \code{"qp"},
#' an alternating quadratic programing algorithm is used. See references for details.
#' @param max.iteration the maximum number of
#' iterations in the optimization algorithm.
#' @param tolerance a numeric value which corresponds to the tolerance paramter in the
#' stopping criterion of the optimization algorithm.
#' @param X a numeric matrix which corresponds to the genotype matrix. This matrix
#' must be of size \eqn{n \times L} where \eqn{n} is the number of individuals and
#' \eqn{L} is the number of loci. Values of this matrix are integers corresponding
#' to the number of variant alleles observed at a locus. If \code{NULL}, \code{XProba}
#' is used.
#' @param openMP.core.num number of core used by the algorithm. It requires that openMP is
#' installed.
#' @param Q.init a numeric matrix which corresponds to the initial value of
#' \code{Q} for the algorithm.
#' @param coord a numeric matrix of size \eqn{n \times 2} where \eqn{n} is the
#' number of individuals. It contains the geographic coordinates (Longitude, Latitude) of
#' all sampled individuals.
#' @param mask if not \code{NULL}, this numeric value is the proportion of genotypic matrix entries
#' which are masked when computing the cross validation criterion.
#' @param algo.copy if TRUE, data will be copied in order to speed the algorithm.
#' @param copy if TRUE data will be copied once.
#' @param verbose If \code{TRUE} run information is printed.
#' @param XProba a numeric matrix which corresponds to genotype likelihoods (probabilities).
#' This matrix must be of size \eqn{n \times (ploidy + 1)L} where
#' \eqn{n} is the number of individuals and \eqn{L} is the number of loci. Entries of
#' this matrix are numeric values between 0 and 1 corresponding
#' to genotype probability. If \code{NULL}, this matrix is computed from the genotype matrix \code{X}.
#' See reference for details.
#'
#' @return An object of class tess3Main which is a list with the following attributes:
#' \describe{
#'    \item{L}{the number of loci}
#'    \item{n}{the number of individuals}
#'    \item{ploidy}{the number of copies of chromosomes.}
#'    \item{K}{the number of ancestral populations.}
#'    \item{G}{the ancestral genotype frequency matrix.}
#'    \item{Q}{the ancestry coefficient matrix.}
#'    \item{Fst}{Fst statistic computed at each locus.}
#'    \item{Fscore}{Fscores computed from the Fst statistics.}
#'    \item{pvalue}{pvalues computed from the Fscores.}
#'    \item{log.pvalue}{The \eqn{log(pvalue)}.}
#'    \item{rmse}{root square mean error between \code{XProba} and
#'    \code{tcrossprod(Q, G)}.}
#'    \item{crossentropy}{cross-entropy error between \code{XProba} and
#'    \code{tcrossprod(Q, G)}.}
#'   \item{crossvalid.rmse}{if masked is not \code{NULL}, root square mean error
#'   between \code{XProba[masked]} and \code{tcrossprod(Q, G)[masked]}.}
#'    \item{crossvalid.crossentropy}{if masked not \code{NULL},
#'    the cross-entropy error between \code{XProba[masked]} and
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
#' @references \url{https://hal.archives-ouvertes.fr/hal-01222555/} \url{http://biorxiv.org/content/early/2016/10/12/080291}
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
