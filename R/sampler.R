###################################################
#######################Q sampler###################
###################################################

#' sample Q = f(coord)
#'
#' @param n
#' @param K
#'
#' @return
#' @export
#'
#' @examples
SampleFuncQ <- function(coord, f = function(X) c(1 / (1 + exp(-0.5 * X[1])), 1 - 1 / (1 + exp(-0.5 * X[1])))) {
  Q = t(apply(coord, 1, f))
  class(Q) <- "tess3Q"
  return(Q)
}

#' sample Q such as alpha ~ unif(0,1) and Q ~ Dirichlet(alpha)
#'
#' @param n
#' @param K
#'
#' @return
#' @export
#'
#' @examples
SampleUnifQ <- function(n, K) {
  TestRequiredPkg("gtools")
  alpha = matrix(runif(n*K,min = 0, max = 1),n,K)
  Q = t(apply(alpha, 1, function(r){gtools::rdirichlet(1,r)}))
  class(Q) <- "tess3Q"
  return(Q)
}

#' For population the center mu_k is computed. Then D_i_k = ||coord_i - mu_k|| is computed. finally Q = f(D) and Q is projected to statisfy constraints.
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleDistFromCenterQ <- function(coord, n.by.pop, K, f = function(D) exp(-D / 0.2)) {

  if (length(n.by.pop) == 1) {
    n.by.pop = rep(n.by.pop,K)
  }
  if (length(n.by.pop) != K) {
    stop("n.by.pop must be a vector of size K or 1 if each pop have the same effective")
  }
  n = nrow(coord)

  # compute center of each cluster
  mu <- t(sapply(split(1:n, rep(1:K, times = n.by.pop)), function(l) {apply(coord[l,],2,mean)}))

  # sample Q
  dist.from.center.aux <- function(c) {
    apply(mu, 1, function(r){return(sqrt(sum((r - c) ^ 2)))})
  }
  dist.from.center = t(apply(coord, 1, dist.from.center.aux))
  Q = f(dist.from.center)
  Q = ProjectQ(Q)
  class(Q) <- "tess3Q"
  return(Q)
}


#'  For population the center mu_k is computed. Then D_i_k = ||coord_i - mu_k|| is computed. finally alpha = f(D) and Q ~ Dirichlet(alpha).
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleDistFromCenterDirichletQ <- function(coord, n.by.pop, K, f = function(D) exp(-D / 0.2)) {


  if (length(n.by.pop) == 1) {
    n.by.pop = rep(n.by.pop,K)
  }
  if (length(n.by.pop) != K) {
    stop("n.by.pop must be a vector of size K or 1 if each pop have the same effective")
  }
  n = nrow(coord)

  # compute center of each cluster
  mu <- t(sapply(split(1:n, rep(1:K, times = n.by.pop)), function(l) {apply(coord[l,],2,mean)}))

  # sample Q
  dist.from.center.aux <- function(c) {
    apply(mu, 1, function(r){return(sqrt(sum((r - c) ^ 2)))})
  }
  dist.from.center = t(apply(coord, 1, dist.from.center.aux))
  alpha = f(dist.from.center)
  Q = t(apply(alpha, 1, function(r){gtools::rdirichlet(1,r)}))
  class(Q) <- "tess3Q"
  return(Q)
}

#' sample Q as describe in "Spatial Inference of Admixture Proportions and Segondary Contact Zones" E. Durand et al.
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleTESS2.3Q <- function(coord, K, W, sigma, rho, f = function(X) c(1.0,X[1], X[2]), beta = matrix(0.0, 3, K)) {
  TestRequiredPkg("gtools")
  TestRequiredPkg("MASS")
  vpmax <-  max(eigen(W)$values)
  if (rho < 0.0 | rho > 1.0 / vpmax) {
    stop("rho must be in (0, 1 / vpmax)")
  }
  if (sigma < 0.0) {
    stop("sigma must be positive")
  }
  n <- nrow(coord)
  Y <- t(MASS::mvrnorm(K, mu = rep(0.0, n), Sigma = sigma ^ 2 * solve(diag(1, nrow = n, ncol = n) - rho * W)))
  log.alpha <- t(apply(coord, 1, function(r) f(r) %*% beta)) + Y
  alpha <- exp(log.alpha)
  # plot(alpha[,1]) # debug
  Q = t(apply(alpha, 1, function(r){gtools::rdirichlet(1,r)}))
  # plot(Q[,1]) # debug
  # plot(Q[,2]) # debug
  class(Q) <- "tess3Q"
  return(Q)
}

###################################################
#######################G sampler###################
###################################################


#' sample G such as G_dl+._k ~ Dirichlet(1/(ploidy + 1))
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @param L
#'
#' @param ploidy
#' @param K
#'
#' @export
SampleUnifDirichletG <- function(L, ploidy, K) {
    # sample G
    G = c()
    for (l in 1:L) {
      G = rbind(G,t(gtools::rdirichlet(K, rep(1.0/(ploidy + 1),(ploidy + 1)))))
    }
    class(G) <- "tess3G"
    return(G)
}



###################################################
###################coord sampler###################
###################################################

#' sample coord such as a mixture of K cluster distributed with gaussian law
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleNormalClusterCoord <- function(n.by.pop, K, sigma1 = 1.0, sigma2 = 0.2 ) {

  TestRequiredPkg("MASS")
  # sample pop center
  mu = matrix(MASS::mvrnorm(K, c(0,0), sigma1 * diag(c(1,1))),K,2)

  # sample coord
  if (length(n.by.pop) == 1) {
    n.by.pop = rep(n.by.pop,K)
  }
  if (length(n.by.pop) != K) {
    stop("n.by.pop must be a vector of size K or 1 if each pop have the same effective")
  }
  coord = c()
  for (k in 1:K) {
    coord = rbind(coord, MASS::mvrnorm(n.by.pop[k], mu[k,], sigma2 * diag(c(1.0,1.0)) ))
  }
  # debug
  # plot(coord, col = rep(2:(K+1),times = n.by.pop))
  return(coord)
}


###################################################
###################TESS3 sampler###################
###################################################


#' Sample X such that P(X_i_dl + j) = Sum(Q_i_k G_k_dl + j).
#'
#'
#' TODO
#'
#' @return TODO
#'
#' @examples
#' TODO
#
#' @export
SampleGenoFromGenerativeModelTESS3 <- function(Q, G, coord, ploidy) {
  res <- list()

  res$n <- nrow(Q)
  res$K <- ncol(Q)
  res$ploidy <- ploidy
  res$L <- nrow(G) / (ploidy + 1)
  res$Q <- Q
  res$G <- G
  res$coord <- coord
  res$X <- matrix(0, res$n ,res$L)

  P = tcrossprod(res$Q, res$G)

  allele = 0:(res$ploidy)
  for (i in 1:res$n) {
    for (j in 1:res$L) {
      res$X[i,j] = allele %*% rmultinom(1, 1, P[i,1 + ((j - 1)*(ploidy + 1)):((j)*(ploidy + 1) - 1)])
    }
  }

  return(res)
}


###################################################
#####################ms sampler####################
###################################################

run.ms <- function(ms.file, nsam, nreps, theta, rho, nsites, M) {
  res <- list()
  tmp.file <- paste0(tempfile(),".geno")
  ms.command <- paste(ms.file, nsam, nreps,
                      "-t", format(theta, scientific = FALSE),
                      "-r", format(rho, scientific = FALSE), format(nsites, scientific = FALSE),
                      "-I 2", nsam / 2 , nsam / 2, M, " >",tmp.file, sep = " ")
  message(paste0("ms command : ", ms.command))
  system(ms.command)

  # read locus position
  ms.res <- readLines(tmp.file)
  res$locus.pos <- scan(text = strsplit(ms.res[6], split = ":")[[1]][2], dec = ".")
  res$L <- length(res$locus.pos)

  # write only geno
  writeLines(ms.res[-(1:6)],tmp.file)
  rm(ms.res)

  # read geno with LEA
  res$X <- t(LEA::read.geno(tmp.file))
  return(res)
}


#' Title
#'
#' @param k
#' @param min.maf the locus with a maf less than this parameter are removed
#' @param plot.debug if TRUE plot at different stage of the simulation
#' @param n number of indivudual to sample
#' @param nsites.neutral number of site between which recombination occur for neutral loci
#' @param nsites.selected number of site between which recombination occur for selected loci
#' @param m.neutral migration rate for neutral loci
#' @param m.selected migration rate for selected loci
#' @param crossover.proba corss-over probability between adjacent site per generation
#' @param mutation.rate.per.site mutation rate per site
#' @param N0 population size
#'
#' @return
#' @export
#'
#' @examples
SampleGenoOFWithMs <- function(n, nsites.neutral, nsites.selected, crossover.proba, m.neutral, m.selected, mutation.rate.per.site, N0 = 10 ^ 6, k = 0.5, min.maf = 0.05, plot.debug = FALSE) {

  #######################
  #########Init##########
  #######################
  res <- list()
  res$n <- n
  res$ploidy <- 1
  res$N0 <- N0
  res$nsites.neutral <- nsites.neutral
  res$nsites.selected <- nsites.selected
  res$crossover.proba <- crossover.proba
  res$m.neutral <- m.neutral
  res$m.selected <- m.selected
  res$mutation.rate.per.site <- mutation.rate.per.site
  res$min.maf <- min.maf
  res$k <- k

  maf <- function(x) {
    n = length(x)
    min(sum(x) / n,1 - sum(x) / n)
  }


  TestRequiredPkg("LEA")
  TestRequiredPkg("foreach")
  require("foreach")
  message("This function required to attach maps namespace.")

  if (plot.debug) {

    # define function used for plot debug
    fst <- function(project,run = 1, K, ploidy = 2){
      library(LEA)
      ll = dim(G(project, K = K, run = run))[1]
      if (ploidy == 2) {freq = G(project, K = K, run = run)[seq(2,ll,by = 3),]/2 + G(project, K = K, run = run)[seq(3,ll,by = 3),] }
      else {freq = G(project, K = K, run = run)[seq(2,ll,by = 2),]}
      q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
      H.s = apply(freq*(1 - freq), MARGIN = 1, FUN = function(x) sum(q*x) )
      P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x) )
      return(1 - H.s/P.t/(1 - P.t))
    }

  }

  # find ms
  if (is.null(.GlobalEnv$tess3.ms) || file.access(.GlobalEnv$tess3.ms, mode = 1)) {
    stop("You need to define the variable tess3.ms as the path to ms executable")
  }

  #######################
  ###simulate neutral####
  #######################
  message("# Simulate neutral locus")
  ## ms parameter rho = 4 * N0 * r
  r <- res$crossover.proba * (res$nsites.neutral - 1)
  rho <- 4 * N0 * r
  # ms parameter theta = 4 * N0 * mu
  mu <- res$mutation.rate.per.site * res$nsites.neutral
  theta <- 4 * N0 * mu
  # ms parameter nsites
  nsites <- res$nsites.neutral
  # ms parameter nsam = 2 * n
  nsam <- 2 * res$n
  # ms parameter nreps = 1
  nreps <- 1
  # ms parameter M = 4 * N0 * m
  M <- 4 * N0 * res$m.neutral
  ms.res <- run.ms(ms.file = .GlobalEnv$tess3.ms,
                   nsam = nsam,
                   nreps = nreps,
                   theta = theta,
                   rho = rho,
                   M = M,
                   nsites = nsites)
  neutral.L <- ms.res$L
  res$neutral.locus.pos <- ms.res$locus.pos
  neutral.X <- ms.res$X
  rm(ms.res)

  # filter minor allele frequencies
  spectrum <- apply(neutral.X, MARGIN = 2, FUN = maf)
  if (plot.debug) {
    par(mfrow = c(2, 1))
    plot(table(spectrum), main = paste("Folded spectrum of neutral genotype before cutting maf less than",min.maf))
  }
  neutral.X <- neutral.X[, spectrum > min.maf]
  neutral.L <- ncol(neutral.X)
  if (plot.debug) {
    spectrum <- apply(neutral.X, MARGIN = 2, FUN = maf)
    plot(table(spectrum), main = paste("Folded spectrum of neutral genotype after cutting maf less than",min.maf))
    par(mfrow = c(1, 1))
  }

  if (plot.debug) {
    message("## debug plots on neutral data")
    K = 2
    tmp.file <- paste0(tempfile(),".geno")
    LEA::write.geno(neutral.X, tmp.file)
    capture.output(obj <- LEA::snmf(tmp.file, K = K, entropy = FALSE, ploidy = 1, project = "new", alpha = 100), file = "/dev/null")
    q.K = LEA::Q(obj, K = K, run = 1)
    barplot(t(q.K), col = rainbow(2), main = "snmf Q computed from neutral dataset")

    fst.values = fst(obj, K = K, ploidy = 1)

    n = dim(q.K)[1]
    fst.values[fst.values < 0] = 0.000001
    z.scores = sqrt(fst.values * (n - K) / (1 - fst.values))

    lambda = median(z.scores ^ 2) / qchisq(1 / 2, df = K - 1)
    message("lambda = ",lambda)

    adj.p.values = pchisq(z.scores ^ 2 / lambda, df = K - 1, lower.tail = FALSE)

    hist(adj.p.values, col = "green", main = "Fst computed from neutral dataset")
    par(mfrow = c(2, 1))
    plot(-log10(adj.p.values), main = "Manhattan plot of p.value computed from neutral dataset", xlab = "Locus", cex = .5, pch = 19, col = "green4")
    plot(-log10(1 - fst.values), main = "Manhattan plot of 1-fst computed from neutral dataset", xlab = "Locus", cex = .5, pch = 19, col = "orange")
    par(mfrow = c(1, 1))
  }

  #######################
  ###simulate selected###
  #######################
  message("# Simulate selected locus")
  ## ms parameter rho = 4 * N0 * r
  r <- res$crossover.proba * (res$nsites.selected - 1)
  rho <- 4 * N0 * r
  # ms parameter theta = 4 * N0 * mu
  mu <- res$mutation.rate.per.site * res$nsites.selected
  theta <- 4 * N0 * mu
  # ms parameter nsites
  nsites <- res$nsites.selected
  # ms parameter nsam = 2 * n
  nsam <- 2 * res$n
  # ms parameter nreps = 1
  nreps <- 1
  # ms parameter M = 4 * N0 * m
  M <- 4 * N0 * res$m.selected
  ms.res <- run.ms(ms.file = .GlobalEnv$tess3.ms,
                   nsam = nsam,
                   nreps = nreps,
                   theta = theta,
                   rho = rho,
                   M = M,
                   nsites = nsites)
  res$selected.locus.pos <- ms.res$locus.pos
  selected.X <- ms.res$X
  rm(ms.res)

  # filter minor allele frequencies
  spectrum = apply(selected.X, MARGIN = 2, FUN = maf)
  if (plot.debug) {
    par(mfrow = c(2, 1))
    plot(table(spectrum), main = paste("Folded spectrum of selected genotype before cutting maf less than",min.maf))
  }
  selected.X <- selected.X[, spectrum > min.maf]
  selected.L <- ncol(selected.X)
  if (plot.debug) {
    spectrum = apply(selected.X, MARGIN = 2, FUN = maf)
    plot(table(spectrum), main = paste("Folded spectrum of selected genotype after cutting maf less than",min.maf))
    par(mfrow = c(1, 1))
  }

  if (plot.debug) {
    message("## debug plots on selected data")
    K = 2
    tmp.file <- paste0(tempfile(),".geno")
    LEA::write.geno(selected.X, tmp.file)
    capture.output(obj <- LEA::snmf(tmp.file, K = K, entropy = T, ploidy = 1, project = "new", alpha = 100), file = "/dev/null")
    q.K = LEA::Q(obj, K = K, run = 1)
    barplot(t(q.K), col = rainbow(2), main = "snmf Q computed from selected dataset")

    fst.values = fst(obj, K = K, ploidy = 1)


    n = dim(q.K)[1]
    fst.values[fst.values < 0] = 0.000001
    z.scores = sqrt(fst.values * (n - K) / (1 - fst.values))

    lambda = median(z.scores ^ 2)/qchisq(1 / 2, df = K - 1)
    message("lambda = ",lambda)

    adj.p.values = pchisq(z.scores ^ 2 / lambda, df = K - 1, lower.tail = FALSE)

    hist(adj.p.values, col = "green", main = "Fst computed from selected dataset")
    par(mfrow = c(2, 1))
    plot(-log10(adj.p.values), main = "Manhattan plot of p.value computed from selected dataset", xlab = "Locus", cex = .5, pch = 19, col = "green4")
    plot(-log10(1 - fst.values), main = "Manhattan plot of 1-fst computed from selected dataset", xlab = "Locus", cex = .5, pch = 19, col = "orange")
    par(mfrow = c(1, 1))
  }

  #######################
  ###selected + neutral##
  #######################
  X <- cbind(neutral.X, selected.X)

  if (plot.debug) {
    message("## debug plots on neutral + selected data")
    K = 2
    tmp.file <- paste0(tempfile(),".geno")
    LEA::write.geno(X, tmp.file)
    capture.output(obj <- LEA::snmf(tmp.file, K = K, entropy = T, ploidy = 1, project = "new", alpha = 100), file = "/dev/null")
    q.K = LEA::Q(obj, K = K, run = 1)
    barplot(t(q.K), col = rainbow(2), main = "snmf Q computed from neutral + selected dataset")

    fst.values = fst(obj, K = K, ploidy = 1)


    n = dim(q.K)[1]
    fst.values[fst.values < 0] = 0.000001
    z.scores = sqrt(fst.values * (n - K) / (1 - fst.values))

    lambda = median(z.scores ^ 2)/qchisq(1 / 2, df = K - 1)
    message("lambda = ",lambda)

    adj.p.values = pchisq(z.scores ^ 2 / lambda, df = K - 1, lower.tail = FALSE)

    hist(adj.p.values, col = "green", main = "Fst computed from neutral + selected dataset")
    par(mfrow = c(2, 1))
    plot(-log10(adj.p.values), main = "Manhattan plot of p.value computed from neutral + selected dataset", xlab = "Locus", cex = .5, pch = 19, col = "green4")
    plot(-log10(1 - fst.values), main = "Manhattan plot of 1-fst computed from neutral + selected dataset", xlab = "Locus", cex = .5, pch = 19, col = "orange")
    par(mfrow = c(1, 1))
  }

  # keep not admixed genotype
  res$selected.locus.index <- (neutral.L + 1):(neutral.L + selected.L)
  res$neutral.locus.index <- 1:neutral.L
  res$L <- neutral.L + selected.L
  res$not.admixed.genotype <- array(NA, dim = c(res$n, res$L, 2))
  res$not.admixed.genotype[,,1] <- X[1:res$n,]
  res$not.admixed.genotype[,,2] <- X[(res$n + 1):(2 * res$n),]
  res$K <- 2
  res$G <- apply(res$not.admixed.genotype, c(2,3), mean)
  # compute fst before convolution
  auxFst <- function(f) {
    sigma2s <- sum(0.5 * f * (1 - f))
    sigma2T <- sum(0.5 * f) * (1 - sum(0.5 * f) )
    return( 1 - sigma2s / sigma2T)
  }
  res$Fst <- apply(res$G, 1, auxFst)

  #######################
  #######admixture#######
  #######################
  message("# Admixe genotype")

  res$coord <- cbind(sort(c(rnorm(res$n / 2, -2, 1), rnorm(res$n / 2, 2, 1))), runif(res$n))
  sigmoid <- function(x){ 1/(1 + exp(-x)) }
  res$Q <- SampleFuncQ(res$coord, f = function(X) c(sigmoid(k * X[1]), 1 - sigmoid(k * X[1])))

  if (plot.debug) {
    barplot(t(res$Q), col = rainbow(2), main = "admixture Q")
  }

  # compute latent factor matrix
  # res$Z <- outer(1:res$n, 1:res$L, Vectorize(function(x,y) sample(c(1,2),1,prob = res$Q[x,]) )) # Too slow !
  res$Z <- ComputeZHelper(res$Q, res$n, res$L)

  # compute admixed matrix
  # res$admixed.genotype <- outer(1:res$n, 1:res$L, Vectorize(function(i,j) res$not.admixed.genotype[i,j,res$Z[i,j]])) # too slow
  res$admixed.genotype <- ComputeAdmixtedGeno(res$not.admixed.genotype,res$Z, res$n, res$L)

  if (plot.debug) {
    message("## debug plots on admixed data")
    K = 2
    tmp.file <- paste0(tempfile(),".geno")
    LEA::write.geno(res$admixed.genotype, tmp.file)
    capture.output(obj <- LEA::snmf(tmp.file, K = K, entropy = T, ploidy = 1, project = "new", alpha = 100), file = "/dev/null")
    q.K = LEA::Q(obj, K = K, run = 1)
    barplot(t(q.K), col = rainbow(2), main = "snmf Q computed from admixed dataset")

    fst.values = fst(obj, K = K, ploidy = 1)


    n = dim(q.K)[1]
    fst.values[fst.values < 0] = 0.000001
    z.scores = sqrt(fst.values * (n - K) / (1 - fst.values))

    lambda = median(z.scores ^ 2)/qchisq(1 / 2, df = K - 1)
    message("lambda = ",lambda)

    adj.p.values = pchisq(z.scores ^ 2 / lambda, df = K - 1, lower.tail = FALSE)

    hist(adj.p.values, col = "green", main = "Fst computed from admixed dataset")
    par(mfrow = c(2, 1))
    plot(-log10(adj.p.values), main = "Manhattan plot of p.value computed from admixed dataset", xlab = "Locus", cex = .5, pch = 19, col = "green4")
    points(neutral.L, max(-log10(adj.p.values)) , type = "h")
    plot(-log10(1 - fst.values), main = "Manhattan plot of 1-fst computed from admixed dataset", xlab = "Locus", cex = .5, pch = 19, col = "orange")
    par(mfrow = c(1, 1))


    message("## debug plots : comparaison of fst")
    par(mfrow = c(2, 1))
    plot(-log10(adj.p.values), main = "Manhattan plot of p.value computed after admixture", xlab = "Locus", cex = .5, pch = 19, col = "green4")
    points(neutral.L, max(-log10(adj.p.values)) , type = "h")

    fst.values <- res$Fst
    fst.values[fst.values < 0] = 0.000001
    z.scores = sqrt(fst.values * (n - K) / (1 - fst.values))
    lambda = median(z.scores ^ 2)/qchisq(1 / 2, df = K - 1)
    adj.p.values = pchisq(z.scores ^ 2 / lambda, df = K - 1, lower.tail = FALSE)
    plot(-log10(adj.p.values), main = "Manhattan plot of p.value computed before admixture", xlab = "Locus", cex = .5, pch = 19, col = "green4")
    points(neutral.L, max(-log10(adj.p.values)) , type = "h")
    par(mfrow = c(1, 1))
  }
  return(res)

}

q <- function() {
  res <- list()
  res$C <- matrix(0,8,2)
  x <- foreach(i = 1:8, .combine='rbind') %:%
    foreach(j = 1:2, .combine='c') %dopar% {
      l <- runif(1, i, 100)
      print(res$C[i,j])
      i + j + l

    }
}
