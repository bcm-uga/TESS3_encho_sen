#' Run tess3 for multiple value of K and with repetition.
#'
#' @param genotype
#' @param geographic.coordinate
#' @param K
#' @param ploidy
#' @param lambda
#' @param rep
#' @param W
#' @param method
#' @param max.iteration
#' @param tolerance
#' @param keep
#'
#' @return
#' @export
#'
#' @examples
tess3project <- function(genotype,
                         geographic.coordinate,
                         K,
                         ploidy,
                         lambda,
                         rep = 1,
                         W = NULL,
                         method = "MCPA",
                         max.iteration = 50,
                         tolerance = 1e-5,
                         openMP.core.num = 4,
                         keep = "all")
{
  # test param :
  if (rep < 1) {
    stop("rep must greater than 1")
  }

  res <- list()
  for (i in seq_along(K)) {
    rmse.max = Inf
    rmse = 1:rep
    tess3.run <- list()
    for (r in 1:rep) {
      tess3.aux <- tess3(genotype = genotype,
                         geographic.coordinate = geographic.coordinate,
                         K = K[i],
                         ploidy = ploidy,
                         lambda = lambda,
                         W = W,
                         method = method,
                         max.iteration = max.iteration,
                         tolerance = tolerance,
                         openMP.core.num = openMP.core.num)
      rmse[r] <- tess3.aux$rmse
      if (keep == "best") {
        if (rmse[r] < rmse.max) {
          tess3.run[[1]] <- tess3.aux
          rmse.max <- rmse[r]
        }
      } else {
        tess3.run[[r]] <- tess3.aux
      }

    }
    res[[i]] <- list(K = K[i], tess3.run = tess3.run, rmse = rmse)
  }
  class(res) <- "tess3project"
  return(res)
}


#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.tess3project <- function(object, ...) {
  cat(paste("=== Object of class tess3project ===\n"))
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
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.tess3project <- function(object, ...) {
  if (length(object) > 0) {
    rmse.med <- seq_along(object)
    rmse.min <- seq_along(object)
    rmse.max <- seq_along(object)
    K <- seq_along(object)
    for (i in seq_along(object)) {
      K[i] <- object[[i]]$K
      rmse.med[i] <- median(object[[i]]$rmse)
      rmse.min[i] <- min(object[[i]]$rmse)
      rmse.max[i] <- max(object[[i]]$rmse)
    }

    plot(K, rmse.med, ...)
    epsilon = 0.02
    segments(K, rmse.min , K, rmse.max)
    segments(K - epsilon, rmse.min , K + epsilon, rmse.min)
    segments(K - epsilon, rmse.max , K + epsilon, rmse.max)
  }

}

#' Test if x is a tess3project object
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.tess3project <- function(x) {
  inherits(x, "tess3project")
}

#' Get tess3 run result
#'
#' @param x
#' @param K
#' @param rep
#'
#' @return
#' @export
#'
#' @examples
Gettess3res <- function(tess3project, K, rep = "best") {
  if (!is.tess3project(tess3project)) {
    stop("tess3project must of class tess3project")
  }
  if (rep == "best") {
    best.rep <- min(which.min(tess3project[[K]]$rmse)[1],length(tess3project[[K]]$tess3.run))
  } else {
    best.rep <- min(as.numeric(rep),length(tess3project[[K]]$tess3.run))
  }
  return(tess3project[[K]]$tess3.run[[best.rep]])
}

