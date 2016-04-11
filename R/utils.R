#' Convert Fst into t-score and compute p value
#'
#'
#'
ComputeTscoreAndPvalue <- function(Fst, K, n) {
  res <- list()
  res$Fscore = Fst / (1 - Fst) * (n - K) / (K - 1)
  # Pvalue, we assume F.score ~ gif * F(df1 = K -1, df2 = n - K)
  res$gif = median(res$Fscore, na.rm = TRUE) / qf(0.5, df1 = K - 1, df2 = n - K)
  res$pvalue <- pf(res$Fscore / res$gif, df1 = K - 1, df2 = n - K, lower.tail = FALSE)
  return(res)
}

#' Convert Fst into chi2 and compute p value
#'
#'
#'
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
