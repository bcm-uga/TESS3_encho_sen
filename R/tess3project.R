#' Run tess3 for multiple value of K and with repetition.
#'
#' @param ploidy TODOC
#' @param lambda TODOC
#' @param rep TODOC
#' @param W TODOC
#' @param method TODOC
#' @param max.iteration TODOC
#' @param tolerance TODOC
#' @param keep TODOC
#' @param X TODOC
#' @param coord TODOC
#' @param openMP.core.num TODOC
#' @param mask TODOC
#' @param XBin TODOC
#' @param K TODOC
#' @param copy TODOC
#' @param algo.copy If TRUE data will be copy to speed the algorithm.
#' @param verbose TODOC
#' @param Q.init TODOC
#'
#' @return An object of class tess3.
#' @export
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
