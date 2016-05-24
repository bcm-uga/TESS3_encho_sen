#' Title
#'
#' @param X
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckX <- function(X, ploidy) {
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  if (!is.double(X)) {
    stop("Element of X must of type double or integer")
  }
  ## check if rang of genotype matrix
  if (max(X, na.rm = TRUE) > ploidy) {
    stop("The maximum value of the X matrix can not be superior than ploidy + 1")
  }
  if (min(X, na.rm = TRUE) < 0) {
    stop("Negative values in the X matrix are not allowed")
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
    stop("Element of XBin must of type double")
  }
  if (ncol(XBin) %% (ploidy + 1) != 0) {
    stop("The number of column must be a multiple of (ploidy + 1)")
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
    stop("Element of W must of type double")
  }
  if (nrow(W) != ncol(W) | !Matrix::isSymmetric(W)) {
    stop("W must be squared symetric matrix")
  }
}

#' Title
#'
#' @param W
#'
#' @return
#'
#' @examples
CheckCoord <- function(Coord) {
  if (!is.matrix(Coord)) {
    stop("Coord must be a matrix")
  }
  if (!is.double(Coord)) {
    stop("Element of Coord must of type double")
  }
}


#' Title
#'
#' @param X
#' @param W
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckXW <- function(X, ploidy, W) {
  CheckW(W)
  CheckX(X, ploidy)
  if (nrow(W) != nrow(X)) {
    stop("W must be of size nrow(X) x nrow(X)")
  }
}

#' Title
#'
#' @param X
#' @param W
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckXBinW <- function(XBin, ploidy, W) {
  CheckW(W)
  CheckXBin(XBin, ploidy)
  if (nrow(W) != nrow(XBin)) {
    stop("W must be of size nrow(X) x nrow(X)")
  }
}

#' Title
#'
#' @param X
#' @param coord
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckXCoord <- function(X, ploidy, coord) {
  CheckCoord(coord)
  CheckX(X, ploidy)
  ## check dim of geno, coord,  are consistent
  if (nrow(X) != nrow(coord)) {
    stop("Number of row of coord and X must be the same")
  }
}

#' Title
#'
#' @param X
#' @param W
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckXWCoord <- function(X, ploidy, W, coord) {
  CheckXW(X, ploidy, W)
  CheckCoord(coord)

  ## check dim of geno, coord,  are consistent
  if (nrow(X) != nrow(coord)) {
    stop("Number of row of coord and X must be the same")
  }
}


#' Title
#'
#' @param X
#' @param W
#' @param ploidy
#'
#' @return
#'
#' @examples
CheckXBinWCoord <- function(XBin, ploidy, W, coord) {
  CheckXBinW(XBin, ploidy, W)
  CheckCoord(coord)

  ## check dim of geno, coord,  are consistent
  if (nrow(XBin) != nrow(coord)) {
    stop("Number of row of coord and X must be the same")
  }
}
