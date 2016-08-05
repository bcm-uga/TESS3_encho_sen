# checks input files, coord, X, W
#


#' Title
#'
#' @param X
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckX <- function(X, ploidy) {
  #cat("Checking genotype input file.\n")
  if (!is.matrix(X)) {
    stop("X must be a matrix.")
  }
  if (!is.double(X)) {
    stop("Elements of X must be of type double or integer.")
  }
  ## check if rang of genotype matrix
  if (max(X, na.rm = TRUE) > ploidy) {
    stop("The maximum value of the genotype matrix (X) cannot be greater than ploidy + 1. Missing data must be encoded as NA.")
  }
  if (min(X, na.rm = TRUE) < 0) {
    stop("Negative values in the genotype matrix (X) are not allowed. Missing data must be encoded as NA.")
  }
}

#' Title
#'
#' @param XBin
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckXBin <- function(XBin, ploidy) {
  if (!is.matrix(XBin)) {
    stop("XBin must be a matrix")
  }
  if (!is.double(XBin)) {
    stop("Elements of XBin must be of type double")
  }
  if (ncol(XBin) %% (ploidy + 1) != 0) {
    stop("The number of columns in Xbin must be a multiple of (ploidy + 1)")
  }
  ## check if rang of genotype matrix
  if (max(XBin, na.rm = TRUE) != 1) {
    stop("The maximum value of the XBin matrix must be 1")
  }
  if (min(XBin, na.rm = TRUE) < 0) {
    stop("Negative values in the XBin matrix are not allowed")
  }
}


#' Title
#'
#' @param W
#'
#' @return
#'
#' @examples
CheckW <- function(W) {
  if (!is.matrix(W)) {
    stop("W must be a matrix")
  }
  if (!is.double(W)) {
    stop("Elements of W must be of type double")
  }
  if (nrow(W) != ncol(W) | !Matrix::isSymmetric(W)) {
    stop("W must be a squared symmetric matrix")
  }
}

#' Title
#'
#' @param Coord
#'
#' @return
#'
#' @examples
CheckCoord <- function(Coord) {
  if (!is.matrix(Coord)) {
    stop("Coord must be a matrix")
  }
  if (!is.double(Coord)) {
    stop("Elements of Coord must be of type double")
  }
}


#' Title
#'
#' @param X
#' @param W
#'
#' @return
#'
#' @examples
CheckXW <- function(X, W) {
  if (nrow(W) != nrow(X)) {
    stop("W must be of size nrow(X) times nrow(X)")
  }
}

#' Title
#'
#' @param XBin
#' @param W
#'
#' @return
#'
#' @examples
CheckXBinW <- function(XBin, W) {
  if (nrow(W) != nrow(XBin)) {
    stop("W must be of size nrow(X) times nrow(X)")
  }
}

#' Title
#'
#' @param X
#' @param coord
#'
#' @return
#'
#' @examples
CheckXCoord <- function(X, coord) {
  ## check dim of geno, coord,  are consistent
  if (nrow(X) != nrow(coord)) {
    stop("Number of row of coord and X must be the same")
  }
}

#' Title
#'
#' @param XBin
#' @param coord
#'
#' @return
#'
#' @examples
CheckXBinCoord <- function(XBin, coord) {
  ## check dim of geno, coord,  are consistent
  if (nrow(XBin) != nrow(coord)) {
    stop("Number of row of coord and X/XBin must be the same")
  }
}

CheckWCoord <- function(W, coord) {
  if ((nrow(W) != nrow(coord)) & (ncol(W) != nrow(coord))) {
    stop("Number of row of coord and W must be the same")
  }
}

