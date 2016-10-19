#' Estimates spatial population structure
#'
#' \code{tess3} estimates spatial population structure using a graph based non
#' negative matrix factorization for several value of ancestry population number.
#' After estimating the population structure is used to
#' compute a Fst statistic for each locus. See references for more details.
#'
#' @param K An integer vector which corresponds to
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
#' @param verbose If \code{TRUE} more information are printed.
#' @param XProba A numeric matrix which correspond to the probability for the genotype.
#' This matrix must be of size \eqn{n \times (ploidy + 1)L} where
#' \eqn{n} is the number of individual, \eqn{L} is the number of loci. Values of
#' this matrix are numeric between 0 and 1 corresponding
#' to the genome probability. It is the matrix used in graph based non negative
#' factorization matrix. If \code{NULL}, it is computed from the genotype matrix \code{X}.
#' See reference for more details.
#' @param rep The number of time the algorithm will be repeated for each value of
#' \code{K}.
#' @param keep If \code{"best"}, for each value of \code{K} only result with best \code{rmse} will be keep. If \code{"all"}, for each value of \code{K} all results will be keep
#' and returned. This second option take more room in memory.
#'
#' @return An object of class tess3 which is a list of size \code{length(K)}.
#' Each element of this list is a list with:
#' \describe{
#'    \item{K}{The number of ancestral population.}
#'    \item{tess3.run}{If \code{keep = "best"}, the \code{\link{tess3Main}} result
#'    with the best \code{rmse}. If \code{keep = "all"}, a list of \code{\link{tess3Main}} result
#'    for each repetition.}
#'    \item{rmse}{The list of the root square mean error between \code{XProba} and
#'                \code{tcrossprod(Q, G)} for each repetition.}
#'    \item{crossentropy}{The list of the cross entropy error between \code{XProba} and
#'                        \code{tcrossprod(Q, G)} for each repetition.}
#'    \item{crossvalid.rmse}{If masked not \code{NULL}.
#'       The list of the root square mean error between \code{XProba[masked]} and
#'        \code{tcrossprod(Q, G)[masked]} for each repetition.}
#'    \item{crossvalid.crossentropy}{If masked not \code{NULL}. The list of
#         the cross entropy error between \code{XProba[masked]} and
#'        \code{tcrossprod(Q, G)[masked]} for each repetition.}
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
#' # Run of tess3 algorithm
#' tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:4,
#'                    method = "projected.ls",
#'                    ploidy = 1)
#'
#' # Plot error
#' plot(tess3.obj, pch = 19, col = "blue",
#'      xlab = "Number of ancestral populations",
#'      ylab = "Cross-validation score")
#'
#' # Retrieve tess3 Q matrix for K = 3 clusters
#' q.matrix <- qmatrix(tess3.obj, K = 3)
#' ## STRUCTURE-like barplot for the Q-matrix
#' barplot(q.matrix, border = NA, space = 0,
#'        xlab = "Individuals", ylab = "Ancestry proportions",
#'        main = "Ancestry matrix") -> bp
#' axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)
#' ## Spatial interpolation of ancestry coefficient
#' my.colors <- c("tomato", "orange", "lightblue")
#' my.palette <- CreatePalette(my.colors, 9)
#' plot(q.matrix, coordinates, method = "map.max", interpol = kriging(10),
#'      main = "Ancestry coefficients",
#'      xlab = "Longitude", ylab = "Latitude",
#'      resolution = c(500,500), cex = .4,
#'      col.palette = my.palette)
#'
#' # Retrieve tess3 results for K = 3
#' p.values <- pvalue(tess3.obj, K = 3)
#' hist(p.values, col = "lightblue")
#' ## Manhatan plot
#' plot(p.values, main = "Manhattan plot",
#'     xlab = "Locus id",
#'     ylab = "-log10(P-values)",
#'     cex = .3, col = "grey")
#'
#' @references \url{http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12471/full}
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
    X <- matrix(as.double(X), nrow(X), ncol(X))
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

#' Plot RMSE(X, Q * t(G)) for all K number of ancestral population with error bars.
#'
#' @param x A tess3 object.
#' @param crossvalid If TRUE error is computed on masked data.
#' @param crossentropy If TRUE the cross entropy error metric is used. If FALSE
#' the RMSE is used.
#' @param ... TODOC
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
#' @param x An object.
#'
#' @return TRUE if x is an tess3 object.
#' @export
is.tess3 <- function(x) {
  inherits(x, "tess3")
}
