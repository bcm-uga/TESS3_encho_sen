SolveLeastSquare <- function(A,b) {
  return(solve( crossprod(A,A) ,crossprod(A,b)))
}

ProjectQ <- function(Q) {
  Q <- apply(Q, 1:2, function(e){max(e,0)})
  for(i in 1:nrow(Q)){
    s <- sum(Q[i,])
    if(s != 0.0) {
      Q[i,] <- Q[i,] / s
    }
  }
  return(Q)
}

ProjectG <- function(G,L,D) {
  G <- apply(G, 1:2, function(e){max(e,0)})
  for(k in 1:ncol(G)) {
    for(l in 1:L) {
      s <- sum(G[(1+D*(l-1)):(D*l) ,k])
      if(s != 0.0) {
        G[(1+D*(l-1)):(D*l) ,k] <- G[(1+D*(l-1)):(D*l) ,k] / s
      }
    }
  }
  return(G)
}


SolveTess3Projected <- function(X, K, d, Lapl, lambda, max.iteration) {

  if(is.null(lambda)) {
    Lapl <- diag(1,nrow(X),nrow(X))
    lambda <- 0.0
  }

  # init
  n <- nrow(X)
  L <- ncol(X)
  X <- ComputeXBin(X,d)
  D <- d+1
  G <- matrix(0, nrow = D * L, ncol = K)
  Q <- matrix(runif(n*K),n,K)
  Q <- ProjectQ(Q)
  normilized.residual.error <- rep(0.0,max.iteration)
  proj.time <- rep(0.0,max.iteration)
  X.norm <- norm(X,"F")

  # compute Lapl diag
  ei <- eigen(Lapl,symmetric = TRUE)
  R <- t(ei$vectors)
  # rotate X
  RX <- R %*% X

  # lamba : reg parameter
  #lambda <- lambda * (D * L * n) / (K * sum(diag(Lapl)))
  #lambda <- lambda * (D * L * n) / (K * n)
  lambda <- lambda * (D * L * n) / (K * n* max(ei$values))
  cat("lambda = ", lambda,"\n")

  # constant
  Ik <- diag(1,K,K)

  # algo
  for(it in 1:max.iteration) {
    ptm <- proc.time()
    # update G
    G <- t(SolveLeastSquare(Q,X))
    G <- ProjectG(G,L,D)

    # update Q
    RQ <- R %*% Q
    for (i in 1:n) {
      RQ[i,] <- t( solve( crossprod(G,G) + lambda * ei$values[i] * Ik ,crossprod(G,t(RX[i, , drop = FALSE]))) )
    }
    Q <- t(R) %*% RQ
    Q <- ProjectQ(Q)

    # compute residual error
    normilized.residual.error[it] <- norm(X - tcrossprod(Q,G),type = "F") / X.norm
    t <- proc.time() - ptm
    proj.time[it] <- t[1] + t[2] + t[4] + t[5]

    cat("iteration : ",it,"& error : ",normilized.residual.error[it],"in time = ",proj.time[it], "\n")

  }

  cat("mean time", mean(proj.time),"\n")

  return(list(Q = Q, G = G, normilized.residual.error = normilized.residual.error, time = mean(proj.time)))

}



SolveTess3QP <- function(X, K, d, Lapl, lambda, max.iteration, tolerance) {

  if(is.null(lambda)) {
    Lapl <- diag(1,nrow(X),nrow(X))
    lambda <- 0.0
  }

  # init
  n <- nrow(X)
  # X <- ComputeXBin(X,d) # done before this function
  D <- d+1
  L <- ncol(X) / D
  G <- matrix(0, nrow = D * L, ncol = K)
  Q <- matrix(runif(n*K),n,K)
  Q <- ProjectQ(Q)
  normilized.residual.error <- rep(0.0,max.iteration)
  QP.time <- rep(0.0,max.iteration)
  X.norm <- norm(X,"F")

  ei <- eigen(Lapl,symmetric = TRUE)
  # lamba : reg parameter
  #lambda = lambda * (D * L * n) / (K * sum(diag(Lapl)))
  # lambda = lambda * (D * L * n) / (K * n)
  lambda <- lambda * (D * L * n) / (K * n* max(ei$values))
  # cat("lambda = ", lambda,"\n") # for debug

  # constant
  Ik <- diag(1,K,K)

  # G QP constant
  QP.G.aux <- list()
  QP.G.aux$A <- kronecker(matrix(1,1,D),diag(1,K,K))
  QP.G.aux$A <- rbind(QP.G.aux$A,diag(1,K*D,K*D))
  QP.G.aux$A <- t(QP.G.aux$A) # same convention that in quadprog solver
  QP.G.aux$b_0 <- matrix(c(rep(1,K),rep(0,K * D)),K*D + K,1)
  QP.G.aux$Id_D <- diag(1,D)
  QP.G.aux$Id_K <- diag(1,K)

  # Q QP constant
  QP.Q.aux <- list()
  QP.Q.aux$A <- kronecker(diag(1,n,n),matrix(1,1,K))
  QP.Q.aux$A <- rbind(QP.Q.aux$A,diag(1,K*n,K*n))
  QP.Q.aux$A <- t(QP.Q.aux$A) # same convention that in quadprog solver
  QP.Q.aux$b_0 <- matrix(c(rep(1,n),rep(0,K * n)),K*n + n,1)
  QP.Q.aux$Id_n <- diag(1,n)
  QP.Q.aux$L_Idk <- lambda * kronecker(Lapl,diag(1,K,K))

  # algo
  it <- 1
  converg = FALSE
  err = -10
  errAux = 0.0
  while (!converg && it <= max.iteration) {
    ptm <- proc.time()
    # update G with QP
    for(l in 1:L) {
      dvec <- crossprod(Q,X[,(1+D*(l-1)):(D*l) ])
      dim(dvec) <- c(K*D,1) # Vec operation
      Dmat <- kronecker(QP.G.aux$Id_D,crossprod(Q,Q))
      aux <- quadprog::solve.QP(Dmat, dvec, QP.G.aux$A, QP.G.aux$b_0, meq=K)
      G[(1+D*(l-1)):(D*l),] <- t(matrix(aux$solution,K,D))
    }

    # update Q with QP
    dvec <- t(X %*% G)
    dim(dvec) <- c(K*n,1) # Vec operation
    Dmat <- kronecker(QP.Q.aux$Id_n,crossprod(G,G)) + QP.Q.aux$L_Idk
    # try to compute chol and ind matrix
    # CDmat <- chol(Dmat)
    # CDmat.inv <- solve(CDmat)
    aux <- quadprog::solve.QP(Dmat, dvec, QP.Q.aux$A, QP.Q.aux$b_0, meq = n, factorized = FALSE)
    Q <- t(matrix(aux$solution,K,n))

    t <- proc.time() - ptm
    QP.time[it] <- t[1] + t[2] + t[4] + t[5]

    # filter NaN (0/0)
    if (length(which(is.na(Q))) > 0) {
      warning("NaN detected in Q, try again !", immediate. = TRUE)
      return(list(Q = NULL, G = NULL, normilized.residual.error = NULL, err = TRUE))
    } else {

      # compute residual error
      normilized.residual.error[it] <- norm(X - tcrossprod(Q,G),type = "F") / X.norm
      errAux <- normilized.residual.error[it]

      # cat("iteration : ",it,"& error : ",normilized.residual.error[it],"in time = ",QP.time[it], "\n") # For debug
      cat("---iteration: ",it,"\n")
      # Test the convergence
      converg = (abs(errAux - err) < tolerance)
      err = errAux
      it <- it + 1
    }
  }

  # cat("mean time", mean(QP.time),"\n") # for debug

  return(list(Q = Q, G = G, normilized.residual.error = normilized.residual.error, err = FALSE, time = mean(QP.time)))

}

