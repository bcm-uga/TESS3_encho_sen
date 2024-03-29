% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tess2tess3}
\alias{tess2tess3}
\title{Import input files from the STRUCTURE and TESS formats}
\usage{
tess2tess3(dataframe = NULL, TESS = TRUE, diploid = TRUE, FORMAT = 1,
  extra.row = 0, extra.column = 0)
}
\arguments{
\item{dataframe}{a data frame read from a STRUCTURE or a TESS input file. Missing data must be encoded by "-9" or by any negative value.}

\item{TESS}{a boolean value set to \code{TRUE} if the TESS format is used, \code{FALSE} if the STRUCTURE format is used. If \code{TRUE}, the
geographic coordinates (Longitude, Latitude) must be binded left to the matrix of genetic markers.}

\item{diploid}{a boolean value set to \code{TRUE} for diploids and \code{FALSE} for haploids.}

\item{FORMAT}{an integer value equal to 1 for markers encoded using one row of data for each individual, and
2 for markers encoded using two rows of data for each individual.}

\item{extra.row}{an integer value indicating the number of extra rows in the header of the input file (marker ids).}

\item{extra.column}{an integer value indicating the number of extra columns in the input file. Extra columns can include individual ids, pop ids,
phenotypes, and they come before the geographic coordinates in TESS input files. Geographic coordinates must be considered as extra columns if the
flag \code{TESS} is set to \code{TESS = FALSE}.}
}
\value{
An object of class \code{list} containing a genotype matrix (X) and individual geographic coordinates (coord).

X a numeric matrix of genotypes with values 0,1,2 or NA.

coord a numeric matrix of geographic coordinates.
}
\description{
This function converts data imported from the STRUCTURE format or from the TESS 2.3 format to the tess3 matrix format.
}
\examples{
library(tess3r)
data(durand09)
d09tess3 <- tess2tess3(durand09, FORMAT = 2, extra.column = 1)
obj <- tess3(X = d09tess3$X, coord = d09tess3$coord,
             K = 1:3, ploidy = 2, openMP.core.num = 4)
Qmatrix <- qmatrix(obj, K = 3)
barplot(Qmatrix, sort.by.Q = FALSE, border = NA,
        space = 0, xlab = "Individuals", ylab = "Ancestry coefficients") -> bp
axis(1, at = 1:nrow(Qmatrix), labels = bp$order, las = 3, cex.axis = .2)
}
\seealso{
\code{\link{tess3}}
}
\author{
Kevin Caye, Flora Jay, Olivier François
}
