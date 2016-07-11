#' Run tess3 for multiple value of K and with repetition.
#'
#' @param ploidy
#' @param lambda
#' @param rep
#' @param W
#' @param method
#' @param max.iteration
#' @param tolerance
#' @param keep
#' @param X
#' @param coord
#' @param openMP.core.num
#' @param mask
#' @param XBin
#' @param K
#' @param Q.init
#' @param no.copy
#'
#' @return
#' @export
#'
#' @examples
tess3 <- function(X,
                  XBin = NULL,
                  coord,
                  K,
                  ploidy,
                  lambda = 1.0,
                  rep = 1,
                  W = NULL,
                  method = "MCPA",
                  max.iteration = 200,
                  tolerance = 1e-5,
                  openMP.core.num = 1,
                  Q.init = NULL,
                  mask = 0.0,
                  copy = TRUE,
                  algo.copy = TRUE,
                  keep = "best",
                  verbose = FALSE)
{
  # test param :
  if (rep < 1) {
    stop("rep must greater than 1")
  }

  if (!copy & is.null(XBin)) {
    stop("To force the function not doing copy of the data, you must set XBin.")
  }
  # Compute XBin
  if (!is.null(X)) {
    XBin <- matrix(0.0, nrow(X), ncol(X) * (ploidy + 1))
    X2XBin(X, ploidy, XBin)
    rm(X)
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
                         XBin = XBin,
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
      crossvalid.rmse[r] <- ifelse(!is.null(tess3.aux$crossvalid.rmse),tess3.aux$crossvalid.rmse, -1)
      crossvalid.crossentropy[r] <- ifelse(!is.null(tess3.aux$crossvalid.crossentropy),tess3.aux$crossvalid.crossentropy, -1)
      if (keep == "best") {
        if (rmse[r] < rmse.max) {
          tess3.run[[1]] <- tess3.aux
          rmse.max <- rmse[r]
        }
      } else {
        tess3.run[[r]] <- tess3.aux
      }
    }
    res[[i]] <- list(K = K[i], tess3.run = tess3.run, rmse = rmse, crossentropy = crossentropy, crossvalid.rmse = crossvalid.rmse, crossvalid.crossentropy = crossvalid.crossentropy)
  }
  class(res) <- c("tess3",class(res))
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
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.tess3 <- function(object, crossvalid = FALSE, crossentropy = FALSE, ...) {
  if (length(object) > 0) {
    if (crossvalid) {
      # test if cross valid rmse is not null
      if (object[[1]]$crossvalid.rmse[1] == -1)
        stop("tess3 was run with mask = 0. Run it with mask > 0.0 to have the cross validation rmse computed")
    }
    med <- seq_along(object)
    min <- seq_along(object)
    max <- seq_along(object)
    K <- seq_along(object)
    for (i in seq_along(object)) {
      K[i] <- object[[i]]$K
      if (!crossentropy) {
        if (!crossvalid) {
          med[i] <- median(object[[i]]$rmse)
          min[i] <- min(object[[i]]$rmse)
          max[i] <- max(object[[i]]$rmse)
        } else {
          med[i] <- median(object[[i]]$crossvalid.rmse)
          min[i] <- min(object[[i]]$crossvalid.rmse)
          max[i] <- max(object[[i]]$crossvalid.rmse)
        }
      } else {
        if (!crossvalid) {
          med[i] <- median(object[[i]]$crossentropy)
          min[i] <- min(object[[i]]$crossentropy)
          max[i] <- max(object[[i]]$crossentropy)
        } else {
          med[i] <- median(object[[i]]$crossvalid.crossentropy)
          min[i] <- min(object[[i]]$crossvalid.crossentropy)
          max[i] <- max(object[[i]]$crossvalid.crossentropy)
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
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.tess3 <- function(x) {
  inherits(x, "tess3")
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

### cette fonction a un bug

Gettess3res <- function(tess3, K, rep = "best") {
  if (!is.tess3(tess3)) {
    stop("tess3 must of class tess3")
  }
  if (rep == "best") {
    best.rep <- min(which.min(tess3[[K]]$rmse)[1],length(tess3[[K]]$tess3.run))
  } else {
    best.rep <- min(as.numeric(rep),length(tess3[[K]]$tess3.run))
  }
  return(tess3[[K]]$tess3.run[[best.rep]])
}

