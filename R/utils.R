
#' This function converts data sets imported from the STRUCTURE and the TESS 2.3 input file formats to the tess3 formats.
#' @title Import input files from STRUCTURE and TESS 2.3
#' @author Kevin Caye, Flora Jay, Olivier François
#' @param dataframe a data frame read from a STRUCTURE or a TESS file.
#' @param TESS a boolean value set to \code{TRUE} if the TESS 2.3 format is used. If \code{TRUE}, the
#' geographic coordinates (Longitude, Latitude) must be binded left to the genetic marker matrix.
#' @param diploid a boolean value set to \code{TRUE} for diploids and \code{FALSE} for haploids.
#' @param FORMAT an integer value equal to 1 for markers encoded using one row of data for each individual, and
#' \code{FORMAT = 2} for markers encoded using two rows of data for each individual.
#' @param extra.row an integer value indicating the number of extra rows in the header of the input file (marker ids).
#' @param extra.column an integer value indicating the number of extra columns in the input file. Extra columns come before the geographic coordinates
#' in TESS 2.3 input files, and include individual ids, pop ids, phenotypes, etc. Geographic coordinates must be considered as extra columns if the flag \code{TESS}
#' is set to \code{TESS = FALSE}. Missing data are encoded as "-9" or any negative values
#' @return An object of class \code{list} containing a genotype matrix and individual geographic coordinates.
#' @return X a matrix of genotypes with values 0,1,2 or NA.
#' @return coord a matrix of geographic coordinates.
#' @seealso \code{\link{tess3}}
#' @examples
#' library(tess3r)
#' data(durand09)
#' d09tess3 <- tess2tess3(durand09, FORMAT = 2, extra.column = 1)
#' obj <- tess3(X = d09tess3$X, coord = d09tess3$coord,
#'              K = 1:3, ploidy = 2, openMP.core.num = 4)
#' Qmatrix <- Gettess3res(obj, K = 3)$Q
#' barplot(Qmatrix, sort.by.Q = FALSE, border = NA,
#'         space = 0, xlab = "Individuals", ylab = "Ancestry coefficients")
#' @export
tess2tess3 <- function(dataframe = NULL, TESS = TRUE, diploid = TRUE, FORMAT = 1, extra.row = 0, extra.column = 0){

  if (!diploid & FORMAT == 2) stop("FORMAT = 2 is for diploids only.")

  if (is.null(dataframe)) stop("dataframe cannot be NULL.")

  if (!is.data.frame(dataframe)) stop("dataframe must be a data.frame object.")

  dat = dataframe

  if (TESS == FALSE){
    if (extra.row > 0) dat = dat[-(1:extra.row),]
    if (extra.column > 0) dat = dat[,-(1:extra.column)]
    n = dim(dat)[1]
    L = dim(dat)[2]
    if (FORMAT == 1 & diploid == FALSE) {n.ind = n; n.loc = L}
    if (FORMAT == 1 & diploid == TRUE) {n.ind = n; n.loc = L/2}
    if (FORMAT == 2 & diploid == TRUE) {n.ind = n/2; n.loc = L}
    cat("Input file in the STRUCTURE format. The genotypic matrix has", n.ind, "individuals and", n.loc,
        "markers.","\n")
    cat("The number of extra rows is", extra.row,"and the number of extra columns is",extra.column,".\n")
  }


  if (TESS == TRUE){
    if (extra.row > 0) dat = dat[-(1:extra.row),]
    if (extra.column > 0) dat = dat[,-(1:extra.column)]

    n = dim(dat)[1]
    L = dim(dat)[2]

    if (FORMAT == 1 & diploid == FALSE) {
      n.ind = n; n.loc = L-2;
      coord = dat[,1:2]
      dat = dat[,-(1:2)]
    }
    if (FORMAT == 1 & diploid == TRUE) {
      n.ind = n; n.loc = L/2 - 1;
      coord = dat[,1:2]
      dat = dat[,-(1:2)]
    }
    if (FORMAT == 2 & diploid == TRUE) {
      n.ind = n/2; n.loc = L - 2;
      coord = dat[seq(1,n,by = 2), 1:2]
      dat = dat[,-(1:2)]
    }
    cat("Input file in the TESS format. The genotypic matrix has", n.ind, "individuals and", n.loc,
        "markers.","\n")
    cat("The number of extra rows is", extra.row,
        "and the number of extra columns is",extra.column,".\n")
  }

  dat = as.matrix(dat)

  unique.dat = unique(as.numeric(dat))
  missing.dat = unique.dat[unique.dat < 0]

  if (length(missing.dat) == 0)  cat("The input file contains no missing genotypes.","\n")
  if (length(missing.dat) == 1)  cat("Missing alleles are encoded as",missing.dat,".\n")
  if (length(missing.dat) > 1) stop("Multiple values for missing data.","\n")



  # Convert allelic data into absence/presence data at each locus
  # Results are stored in the "dat.binary" object

  L = dim(dat)[2]

  if (FORMAT == 1 & diploid == FALSE) {

    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]== i )
      LL = dim(dat.binary)[2]
      ind = which(allele < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}

  if (FORMAT == 1 & diploid == TRUE) {
    dat.2 = matrix(NA, ncol = L/2, nrow = 2*n)
    for (ii in 1:n){
      dat.2[2*ii-1,] = dat[ii, seq(1,L,by = 2)]
      dat.2[2*ii,] = dat[ii,seq(2,L,by = 2)]
    }
    L = dim(dat.2)[2]

    dat.binary = NULL

    for (j in 1:L){
      allele = sort(unique(dat.2[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat.2[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(allele < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}



  if (FORMAT == 2 & diploid == TRUE) {
    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(allele < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}

  # Compute a genotype count for each allele (0,1,2 or 9 for a missing value)
  # The results are stored in 'genotype'

  n = dim(dat.binary)[1]

  if (diploid == TRUE){
    n = n/2
    genotype = matrix(NA,nrow=n,ncol=dim(dat.binary)[2])
    for(i in 1:n){
      genotype[i,]= dat.binary[2*i-1,]+dat.binary[2*i,]
      genotype[i, (genotype[i,] < 0)] = NA
    }}

  if (FORMAT == 1 & diploid == FALSE){
    genotype = dat.binary
    for(i in 1:n){
      genotype[i, (genotype[i,] < 0)] = NA
    }}

  return(list(X = as.matrix(genotype), coord = as.matrix(coord)))
}


###################################################
#' This function tests Q matrix objects, and converts \code{matrix} objects into valid Q matrices
#' @title Converts Q matrix
#' @author Kevin Caye, Olivier François
#' @param Q an object of class \code{matrix} containing a matrix of ancestry coefficients.
#' @return An object of class \code{tess3Q}.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{barplot.tess3Q}}
#' @examples
#' library(tess3r)
#' ## an example with 3 individuals and 2 clusters
#' Qmatrix <- matrix(c(0.4,0.6,0.3,0.7, 0.2, 0.8), byrow = TRUE, nrow = 3)
#' Qmatrix <- as.qmatrix(Qmatrix)
#' barplot(Qmatrix, space = 0, xlab = "individuals",
#'         ylab = "Ancestry proportions", main = "Ancestry matrix")
#' @export
as.qmatrix <- function(Q){
  if (class(Q) != "matrix") stop("Input matrix is not an ancestry matrix.")
  if (min(Q) < 0) stop("Q contains negative elements.")
  sumofq <- apply(Q, MARGIN = 1, sum)
  if ( sum(sumofq) != nrow(Q)) stop("Input matrix is not an ancestry matrix.")
  class(Q) = "tess3Q"
  return(Q)
}

#' This function creates a list of color palettes for the plot and barplot
#' functions
#' @title Create a list of palettes
#' @author Kevin Caye, Olivier François
#' @param color.vector a vector of R colors.
#' @param palette.length an integer number of colors in each palette.
#' @return An object of class \code{list} containing a list of color palettes.
#' @seealso \code{\link{plot.tess3Q}} \code{\link{barplot.tess3Q}}
#' @examples
#' ## an A. thaliana example
#' data(data.at)
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5, ploidy = 1,
#'              openMP.core.num = 4)
#' Qmatrix <- obj$Q
#' my.colors <- c("tomato", "yellow", "blue", "wheat","olivedrab")
#' my.palette <- CreatePalette(my.colors, 9)
#' plot(obj$Q, data.at$coord, method = "mapping.max", col.palette = my.palette,
#'      interpol = kriging(10), cex = .4, xlab = "Longitude", ylab= "Latitude",
#'      main = "Ancestry coefficients")
#' @export
CreatePalette <- function(color.vector = c("tomato", "chartreuse", "gold", "blue", "violet", "wheat","olivedrab"), palette.length = 9){
  ll = NULL
   for (i in 1:length(color.vector)){
     ll[[i]] = colorRampPalette(c("grey96", color.vector[i]))(palette.length)
   }
 ll
}

#' Convert Fst into t-score and compute p value
#'
#'
#'
#' @param Fst Matrix of Fst.
#' @param K Number of ancestral populations.
#' @param n Number of individuals.
ComputeTscoreAndPvalue <- function(Fst, K, n) {
  res <- list()
  res$Fscore = Fst / (1 - Fst) * (n - K) / (K - 1)
  # Pvalue, we assume F.score ~ gif * F(df1 = K -1, df2 = n - K)
  res$gif = median(res$Fscore, na.rm = TRUE) / qf(0.5, df1 = K - 1, df2 = n - K)
  res$pvalue <- pf(res$Fscore / res$gif, df1 = K - 1, df2 = n - K,
                   lower.tail = FALSE)
  res$log.pvalue <- pf(res$Fscore / res$gif, df1 = K - 1, df2 = n - K,
                   lower.tail = FALSE, log = TRUE)
  return(res)
}

#' Convert Fst into chi2 and compute p value
#'
#'
#'
#' @param Fst Matrix of Fst.
#' @param K Number of ancestral populations.
#' @param n Number of individuals.
ComputeChi2AndPvalue <- function(Fst, K, n) {
  res <- list()
  # Convert Fst into chi 2
  res$chi2 = Fst * (n - K)/(1 - Fst)
  # compute the gif
  res$gif = median(res$chi2) / qchisq(1 / 2, df = K - 1)
  # compute adjusted p-values from the combined z-scores
  res$pvalue = as.numeric(pchisq(res$chi2 / res$gif, df = K - 1, lower.tail = FALSE))
  return(res)
}
