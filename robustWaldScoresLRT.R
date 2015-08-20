robustWaldScoresLRT <- function(y=y, x=x, beta0=beta0, betaA=betaA, sigma=sigma, 
                                Amat=Amat, bvec=bvec, meq=meq) { 
  
  
  if (qr(Amat)$rank < nrow(Amat)) {
    stop("Constraints matrix must have full row-rank.")
  }
  if (is.null(Amat)) {
    stop("No constraints matrix has been specified.")
  }
  if (meq == nrow(Amat)) {
    stop("Test not applicable with equality restrictions only.")
  }
  if (is.null(bvec)) {
    bvec <- rep(0, nrow(Amat))
  }
  if (!is.vector(bvec)) {
    stop("bvec must be a vector.")
  }
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
   
  #tukey bisquare tuning constant
  c = 4.685061
  
  #Calculate M, Q, V 
  res0 <- y - x%*%beta0
  res2 <- y - x%*%betaA
  
  rstar0 <- res0/sigma                                                         
  rstar2 <- res2/sigma
  
  psi0   <- tukeyChi(rstar0, cc=c, deriv=1)  
  psi1   <- tukeyChi(rstar2, cc=c, deriv=1) 
  psideriv0 <- tukeyChi(rstar0, cc=c, deriv=2) 
  psideriv1 <- tukeyChi(rstar2, cc=c, deriv=2) 
  
  #compute M 
  weightsM <- psideriv1 / sigma             
  WM <- weightsM %*% rep(1, p)
  xwM <- x * WM
  M <- t(x) %*% xwM / n
  
  #compute Q 
  weightsQ <- psi1^2
  WQ <- weightsQ %*% rep(1, p)
  xwQ <- x * WQ
  Q <- t(x) %*% xwQ / n
  
  
  idx1 <- which(colSums(abs(Amat)) > 0L)
  idx0 <- which(colSums(abs(Amat)) == 0L)
  
  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)

# V[(p0+1):pa,(p0+1):pa]
  
  V22 <- V[idx1,idx1]
  
  #inverse information matrix
  #V.inv <- solve(V)
  #submatrix of V
  V22.inv <- solve(V22)
  
  #Schur complement?
  M221 <- M[idx1,idx1] - M[idx1,idx0] %*% solve(M[idx0,idx0,drop=FALSE], 
                                                M[idx0,idx1,drop=FALSE])

#  M221 <- M[(p0+1):p,(p0+1):p] - M[(p0+1):p,1:p0] %*% solve(M[1:p0,1:p0,drop=FALSE], M[1:p0,(p0+1):p,drop=FALSE])
 
  weightsZ <- psi0
  Z <- t(x) %*% weightsZ / n  
  
  
  #Wald-type test statistic, Silvapulle (1996, eq. 2.6)
  #FIXME Amat%*%betaA, e.g., c(0,-1,0,0) %*% c(b1,b2,b3,b4) produces a negative value for b2.
  Tn <- sqrt(n)*(Amat%*%betaA)
  Dn <- Tn
  Dmat <- solve(Amat[,idx1,drop=FALSE]%*%V22%*%t(Amat[,idx1,drop=FALSE]))
  dvec <- t(Dn)%*%Dmat
  out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq) 
  b <- out$solution
  TS_W  <- (t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b)) 
  #unconstrained
  #n * betaA[idx1] %*% solve(V22, betaA[idx1]) 

  
  #Score-type test (Silvapulle, 1996, eq. 2.6)
  An <- sqrt(n) * Amat[,idx1,drop=FALSE]%*%(solve(M221) %*% Z[idx1])
  Dn <- An
  Dmat <- solve(Amat[,idx1,drop=FALSE]%*%V22%*%t(Amat[,idx1,drop=FALSE]))
  dvec <- t(Dn)%*%Dmat
  out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq) 
  b <- out$solution
  TS_S <- (t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b))

  #unconstrained
  #result.C <- M221 %*% V22 %*% t(M221)
  #result.R <- n * t(Z[idx1]) %*% solve(result.C, Z[idx1])

  return(list(TS_W = TS_W,          
              #TS_W2 = TS_W2,
              #result_W = result_W,
              TS_S = TS_S,
              #TS_MH = TS_MH,
              #TS_S2 = TS_S2,
              V22 = Dmat
  )
  )
  
}

