#' Estimate ancestry coefficients and run genome scans for selection
#'
#' \code{tess3} is the main function of the \code{tess3r} package. It runs
#' a graph-based nonnegative matrix factorization algorithm that includes
#' geographic data in the estimation of spatial population structure.
#' The function requires individual genotypes, geographic coordinates, and
#' it can be run for multiple values of the number of ancestral populations.
#' In addition, the function uses the estimates of ancestry coefficients to
#' compute an Fst statistic for each locus, and to return test significance
#' values for a null hypothesis of selective neutrality. See the references
#' for more details.
#'
#' @param K an integer or a vector of integers corresponding to
#' the number of ancestral populations.
#' @param ploidy an integer corresponding to ploidy of the studied
#' organism. Haploids have ploidy = 1, diploids have ploidy = 2, etc.
#' @param lambda a numeric value for the spatial regularization parameter.
#' The default value lambda = 1 attributes equal weights to the loss function
#' and to the penalty function.
#' @param W a matrix which corresponds to the graph weightings.
#' If NULL, W is computed as
#' \code{W[i,j] = exp(-(coord[i] - coord[j])^2 / sigma^2)},
#' where \code{coord[i]} represents
#' the geographic coordinates for individual i, and where
#' \code{sigma} is equal to 5 percent of the average geographic distance between individuals.
#' @param method a character string \code{"projected.ls"} or \code{"qp"}. If \code{"projected.ls"},
#' an alternating projected least squares algorithm is used. If \code{"qp"},
#' an alternating quadratic programing algorithm is used. See references for more
#' details
#' @param max.iteration the maximum number of iterations of the optimization algorithm.
#' @param tolerance a numeric value corresponding to the
#' stopping criteria of the optimization algorithm.
#' @param X a matrix of individual genotypes. This matrix must
#' have \eqn{n} rows and  \eqn{L} columns where \eqn{n} is the number of individuals and
#' \eqn{L} is the number of loci. The entries of this matrix are integers between 0 and
#' ploidy, that correspond to the number of derived/reference alleles observed at each locus.
#' If \code{NULL}, the matrix of genotype likelihood, \code{XProba}, is used.
#' @param openMP.core.num integer representing the number of cores used by the optimization
#' algorithm. It requires that the openMP library is installed in your OS (default for macOS is no).
#' @param Q.init a matrix for initial values of ancestry coefficients for the algorithm. The
#'  default value is a random matrix.
#' @param coord a matrix of size \eqn{n \times 2} where \eqn{n} is the
#' number of individuals encoding the longitude and latitude of each individual (numeric system).
#' @param mask If not \code{NULL}, this numeric value is the proportion of masked data
#' when computing the cross-validation criterion.
#' @param algo.copy boolean. If TRUE data is copied to speed up the algorithm.
#' @param verbose If \code{TRUE} more information is printed.
#' @param XProba A matrix which contains individual genotype likelihoods (probabilities) for each
#' locus. This matrix must contain \eqn{n} rows and \eqn{(ploidy + 1)L} columns where
#' \eqn{n} is the number of individuals, and \eqn{L} is the number of loci. The entries of
#' this matrix are numeric values ranging between 0 and 1, and corresponding to genotype
#' probabilities for each locus. If \code{NULL}, deterministic values are computed from the
#' genotypic matrix \code{X}. See the references for more details.
#' @param rep integer. The number of time the algorithm will be repeated for each value of
#' \code{K}.
#' @param keep If \code{"best"}, only the result with the lowest \code{rmse} score will be kept
#' for each value of \code{K}. If \code{"all"}, all results will be kept
#' and returned for each value of \code{K}. The second option uses more space in memory.
#'
#' @return An object of class tess3 which corresponds to a list of length \code{length(K)}.
#' Each element of this list has the following attributes
#' \describe{
#'    \item{K}{the number of ancestral populations}
#'    \item{tess3.run}{if \code{keep = "best"}, the \code{\link{tess3Main}} result
#'    with the lowest value of the \code{rmse} (loss) function. If \code{keep = "all"},
#'    a list of \code{\link{tess3Main}} results for each repetition}
#'    \item{rmse}{root mean squared error between the genotypic matrix \code{XProba} and the
#'    fitted matrix for each program repetition}
#'    \item{crossentropy}{cross-entropy between the genotypic matrix \code{XProba} and the
#'    fitted matrix for each program repetition}
#'    \item{crossvalid.rmse}{root square mean error between the masked values of genotypic matrix
#'     \code{XProba[masked]} and their fitted values for each repetition. If mask is FALSE, then \code{NULL}.}
#'    \item{crossvalid.crossentropy}{cross-entropy between the masked values of genotypic matrix
#'     \code{XProba[masked]} and their fitted values for each repetition. If mask is FALSE, then \code{NULL}.}
#' }
#' Methods available for this class:
#' \itemize{
#'   \item \code{\link{plot.tess3}}
#'   \item \code{\link{summary.tess3}}
#'   \item \code{\link{is.tess3}}
#'   \item \code{\link{Gettess3res}}
#'   \item \code{\link{qmatrix}}
#'   \item \code{\link{pvalue}}
#' }
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
#' # Running the tess3 function
#' tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:4,
#'                    method = "projected.ls",
#'                    ploidy = 1)
#'
#' # Plot error
#' plot(tess3.obj, pch = 19, col = "blue",
#'      xlab = "Number of ancestral populations",
#'      ylab = "Cross-validation score")
#'
#' # Retrieve the Q-matrix for K = 3 clusters
#' q.matrix <- qmatrix(tess3.obj, K = 3)
#'
#' ## STRUCTURE-like barplot for the Q-matrix
#' barplot(q.matrix, border = NA, space = 0,
#'        xlab = "Individuals", ylab = "Ancestry proportions",
#'        main = "Ancestry matrix") -> bp
#' axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)
#'
#' ## Spatial interpolation of ancestry coefficient
#' my.colors <- c("tomato", "orange", "lightblue")
#' my.palette <- CreatePalette(my.colors, 9)
#' plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),
#'      main = "Ancestry coefficients",
#'      xlab = "Longitude", ylab = "Latitude",
#'      resolution = c(500,500), cex = .4,
#'      col.palette = my.palette)
#'
#' ## Genome scan p-values for K = 3
#' p.values <- pvalue(tess3.obj, K = 3)
#' hist(p.values, col = "lightblue")
#'
#' ## Manhatan plot
#' plot(p.values, main = "Manhattan plot",
#'     xlab = "Locus id",
#'     ylab = "-log10(P-values)",
#'     cex = .3, col = "grey")
#'
#' @references
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12471/full}
#' Caye, Kevin et al. (2016) Fast Inference of Individual Admixture Coefficients Using Geographic Data. bioRxiv
#' doi:10.1101/080291. \url{http://biorxiv.org/content/early/2016/10/12/080291}
#'
#' @seealso \code{\link{tess3Main}}, \code{\link{plot.tess3Q}},
#' \code{\link{barplot.tess3Q}}
tess3 <- function(X,
                  XProba = NULL,
                  coord,
                  K,
                  ploidy,
                  lambda = 1.0,
                  rep = 1,
                  W = NULL,
                  method = "projected.ls",
                  max.iteration = 200,
                  tolerance = 1e-5,
                  openMP.core.num = 1,
                  Q.init = NULL,
                  mask = 0.0,
                  algo.copy = TRUE,
                  keep = "best",
                  verbose = FALSE)
{
  # test param :
  if (rep < 1) {
    stop("rep must greater than 1")
  }

  copy = TRUE
  if (copy & !is.null(X)) {
    X <- as.matrix(X) # to handle type conversion
    X <- matrix(as.double(X), nrow(X), ncol(X)) # to ensure we have a matrix of double
    CheckX(X, ploidy)
  } else if (!copy & is.null(XProba)) {
    stop("To force the function not doing copy of the data, you must set XProba")
  }

  # Compute XBin
  if (!is.null(X)) {
    XProba <- matrix(0.0, nrow(X), ncol(X) * (ploidy + 1))
    X2XBin(X, ploidy, XProba)
    rm(X)
  }


  # if user want only 1 run of tess3 we return a list of result
  if (length(K) == 1 & rep == 1) {
    res <- tess3Main(X = NULL,
                     XProba = XProba,
                     coord = coord,
                     K = K,
                     ploidy = ploidy,
                     lambda = lambda,
                     W = W,
                     method = method,
                     max.iteration = max.iteration,
                     tolerance = tolerance,
                     openMP.core.num = openMP.core.num,
                     Q.init = Q.init,
                     mask = mask,
                     copy = TRUE,
                     algo.copy = algo.copy,
                     verbose = verbose)
    class(res) <- c(class(res), "tess3")
    return(res)
  }

  res <- list()
  for (i in seq_along(K)) {
    rmse.max = Inf
    rmse = 1:rep
    crossentropy = 1:rep
    crossvalid.rmse = 1:rep
    crossvalid.crossentropy = 1:rep
    tess3.run <- list()
    for (r in 1:rep) {
      tess3.aux <- tess3Main(X = NULL,
                             XProba = XProba,
                             coord = coord,
                             K = K[i],
                             ploidy = ploidy,
                             lambda = lambda,
                             W = W,
                             method = method,
                             max.iteration = max.iteration,
                             tolerance = tolerance,
                             openMP.core.num = openMP.core.num,
                             Q.init = Q.init,
                             mask = mask,
                             copy = copy,
                             algo.copy = algo.copy,
                             verbose = verbose)
      rmse[r] <- tess3.aux$rmse
      crossentropy[r] <- tess3.aux$crossentropy
      crossvalid.rmse[r] <- ifelse(!is.null(tess3.aux$crossvalid.rmse),
                                   tess3.aux$crossvalid.rmse, -1)
      crossvalid.crossentropy[r] <-
        ifelse(!is.null(tess3.aux$crossvalid.crossentropy),
               tess3.aux$crossvalid.crossentropy, -1)
      if (keep == "best") {
        if (rmse[r] < rmse.max) {
          tess3.run[[1]] <- tess3.aux
          rmse.max <- rmse[r]
        }
      } else {
        tess3.run[[r]] <- tess3.aux
      }
    }
    res[[K[i]]] <- list(K = K[i], tess3.run = tess3.run,
                        rmse = rmse,
                        crossentropy = crossentropy,
                        crossvalid.rmse = crossvalid.rmse,
                        crossvalid.crossentropy = crossvalid.crossentropy)
  }
  class(res) <- c(class(res), "tess3")
  return(res)
}


#' Title
#'
#' @param object tess3 object.
#' @param ... TODOC
#'
#' @export
#'
summary.tess3 <- function(object, ...) {
  cat(paste("=== Object of class tess3 ===\n"))
  if (length(object) > 0) {
    cat(paste("Number of individuals n:", object[[1]]$tess3.run[[1]]$n,"\n"))
    cat(paste("Number of loci L:", object[[1]]$tess3.run[[1]]$L,"\n"))
    cat(paste("Ploidy:", object[[1]]$tess3.run[[1]]$ploidy,"\n"))
    K = paste0(object[[1]]$K)
    for (i in seq_along(object)[-1]) {
      K <- paste0(K, ", ", object[[i]]$K)
    }
    cat(paste("Number of ancestral populations K:", K,"\n"))
  }
}

#' Plot cross-validation errors for all values of number of ancestral populations
#'
#' @param x a tess3 object.
#' @param crossvalid if TRUE, errors are evaluated on masked data.
#' @param crossentropy If TRUE, the cross-entropy error is used. If FALSE,
#' the root mean square error is used.
#' @param ... other graphic parameters.
#'
#' @export
#'
plot.tess3 <- function(x, crossvalid = FALSE, crossentropy = FALSE, ...) {
  if (length(x) > 0) {
    if (crossvalid) {
      # test if cross valid rmse is not null
      if (x[[1]]$crossvalid.rmse[1] == -1)
        stop("tess3 was run with mask = 0. Run it with mask > 0.0 to have the cross validation rmse computed")
    }
    med <- seq_along(x)
    min <- seq_along(x)
    max <- seq_along(x)
    K <- seq_along(x)
    for (i in seq_along(x)) {
      K[i] <- x[[i]]$K
      if (!crossentropy) {
        if (!crossvalid) {
          med[i] <- median(x[[i]]$rmse)
          min[i] <- min(x[[i]]$rmse)
          max[i] <- max(x[[i]]$rmse)
        } else {
          med[i] <- median(x[[i]]$crossvalid.rmse)
          min[i] <- min(x[[i]]$crossvalid.rmse)
          max[i] <- max(x[[i]]$crossvalid.rmse)
        }
      } else {
        if (!crossvalid) {
          med[i] <- median(x[[i]]$crossentropy)
          min[i] <- min(x[[i]]$crossentropy)
          max[i] <- max(x[[i]]$crossentropy)
        } else {
          med[i] <- median(x[[i]]$crossvalid.crossentropy)
          min[i] <- min(x[[i]]$crossvalid.crossentropy)
          max[i] <- max(x[[i]]$crossvalid.crossentropy)
        }
      }
    }

    plot(K, med, ...)
    epsilon = 0.02
    segments(K, min , K, max)
    segments(K - epsilon, min , K + epsilon, min)
    segments(K - epsilon, max , K + epsilon, max)
  }
}

#' Test if x is a tess3project object
#'
#' @param x an object.
#'
#' @return TRUE if x is a tess3 object.
#' @export
is.tess3 <- function(x) {
  inherits(x, "tess3")
}
