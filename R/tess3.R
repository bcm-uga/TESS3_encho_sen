#' Estimates individual ancestry coefficients, ancestral allele frequencies and an ancestral allele frequency differentiation statistic.
#'
#' \code{\link{TESS3}} estimates admixture coefficients using a graph based Non-Negative
#' Matrix Factorization algorithms, and provide STRUCTURE-like outputs. \code{\link{TESS3}} also computes
#' an ancestral allele frequency differentiation statistic
#'
#' @param input.file A character string containing a the path to the input genotype file,
#' a genotypic matrix in the \code{\link[LEA]{geno}} format.
#' @param input.coord A character string containing a the path to the input coordinate file,
#' a coordinate matrix in the \code{\link{coord}} format.
#' @param K An integer vector corresponding to the number of ancestral populations for
#' which the TESS3 algorithm estimates have to be calculated.
#' @param project A character string among "continue", "new", and "force". If "continue",
#' the results are stored in the current project. If "new", the current
#' project is removed and a new one is created to store the result. If
#' "force", the results are stored in the current project even if the input
#' file has been modified since the creation of the project.
#' @param repetitions An integer corresponding with the number of repetitions for each value of \code{K}.
#' @param alpha A numeric value corresponding to the TESS3 regularization parameter.
#' The results depend on the value of this parameter.
#' @param tolerance A numeric value for the tolerance error.
#' @param entropy A boolean value. If true, the cross-entropy criterion is calculated.
#' @param percentage A numeric value between 0 and 1 containing the percentage of
#' masked genotypes when computing the cross-entropy
#' criterion. This option applies only if \code{entropy == TRUE}.
#' @param I The number of SNPs to initialize the algorithm. It starts the algorithm
#' with a run of snmf using a subset of nb.SNPs random SNPs. If this option
#' is set with nb.SNPs, the number of randomly chosen SNPs is the minimum
#' between 10000 and 10 \% of all SNPs. This option can considerably speeds
#' up snmf estimation for very large data sets.
#' @param iterations An integer for the maximum number of iterations in algorithm.
#' @param ploidy 1 if haploid, 2 if diploid, n if n-ploid.
#' @param seed A seed to initialize the random number generator. By default, the seed is randomly chosen.
#' @param CPU A number of CPUs to run the parallel version of the algorithm. By default, the number of CPUs is 1.
#' @param Q.input.file A character string containing a path to an initialization file for Q, the individual admixture coefficient matrix.
#'
#' @return \code{TESS3} returns an object of class \code{tess3Project}.
#' The following methods can be applied to the object of class {tess3Project}:
#' \item{plot}{
#'   Plot the minimal cross-entropy in function of K.
#' }
#' \item{show}{
#'   Display information about the analyses.
#' }
#' \item{summary}{
#'   Summarize the analyses.
#' }
#' \item{Q}{
#'   Return the admixture coefficient matrix for the chosen run with K
#'   ancestral populations.
#' }
#' \item{G}{
#'   Return the ancestral allele frequency matrix for the chosen run with K
#'   ancestral populations.
#' }
#' \item{FST}{
#'   Return ancestral allele frequency differentiation statistic matrix for the chosen run with K
#'   ancestral populations.
#' }
#' \item{cross.entropy}{
#'   Return the cross-entropy criterion for the chosen runs with K
#'   ancestral populations.
#' }
#' \item{load.snmfProject(file.tess3Project)}{
#'   Load the file containing an tess3Project objet and return the tess3Project
#'   object.
#' }
#' \item{remove.snmfProject(file.tess3Project)}{
#'   Erase a \code{tess3Project} object. Caution: All the files associated with
#'   the object will be removed.
#' }
#'
#' @examples
#' ### Example of analyses using snmf ###
#' # dataset simulated from the plant species Arabidopsis thaliana
#' # It contains 26943 SNPs for 170 individuals.
#' athaliana.genofile <- system.file("extdata/Athaliana","Athaliana.geno",package = "tess3r")
#' athaliana.coord <- system.file("extdata/Athaliana","Athaliana.coord",package = "tess3r")
#'
#' #################
#' # runs of TESS3 #
#' #################
#'
#' # main options, K: (the number of ancestral populations),
#' #        entropy: calculate the cross-entropy criterion,
#'
#' # Runs with K between 1 and 5 with cross-entropy and 2 repetitions.
#' project <- TESS3( athaliana.genofile, athaliana.coord, K=1:5, entropy = TRUE, repetitions = 2,
#'                   project = "new")
#'
#' # plot cross-entropy criterion of all runs of the project
#' plot(project, lwd = 5, col = "red", pch=1)
#'
#' # get the cross-entropy of each run for K = 4
#' ce = cross.entropy(project, K = 4)
#'
#' # select the run with the lowest cross-entropy
#' best = which.min(ce)
#'
#' # plot the best run for K = 3 (ancestry coefficients).
#' barplot(t(Q(project, K = 3, run = best)), col = c(2:4) )
#'
#'
#' @aliases Q G FST cross.entropy load.snmfProject remove.snmfProject
#'
#' @export
TESS3 <- function(genotype,
                  coordinate,
                  K,
                  ploidy,
                  lambda,
                  method = "MCPA",
                  max.iteration = 20)
{
  # TODO : check args

  # compute laplacian
  W <- ComputeHeatKernelWeight(coordinate, NULL)
  Lapl <- ComputeGraphLaplacian(W)

  # compute Q and G matrix
  if (method == "MCPA") {
    res <- SolveTess3Projected(genotype,
                               K,
                               ploidy,
                               Lapl,
                               lambda,
                               max.iteration = max.iteration)
  } else if (method == "OQA") {
    res <- SolveTess3QP(genotype,
                        K,
                        ploidy,
                        Lapl,
                        lambda,
                        max.iteration = max.iteration)
  }
  return(res)

}

#' tess3r : TESS3 R Package
#'
#' This R package implements the TESS3 method and tools useful to plot program outputs.
#'
#' @docType package
#'
#' @name tess3r
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' @useDynLib TESS3enchoSen
NULL
