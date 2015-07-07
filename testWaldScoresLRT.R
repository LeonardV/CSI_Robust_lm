#Parts of this R code are taken from the Robustbase package.
testWaldScoresLRT <- function(y, x, beta0, beta2, sigma, Amat, meq) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  idx0 <- which(colSums(Amat) == 0)
  idx1 <- which(colSums(Amat) == 1)
  
  #tukey bisquare tuning constant
  c2 = 4.685061
  
  #Calculate M, Q, V 
  res0 <- y - x%*%beta0
  res2 <- y - x%*%beta2
    
  rstar0 <- res0/sigma                                                         
  rstar2 <- res2/sigma
  
  psi0.c2   <- tukeyChi(rstar0, cc=c2, deriv=1)  
  psi1.c2   <- tukeyChi(rstar2, cc=c2, deriv=1) 
  psideriv0.c2 <- tukeyChi(rstar0, cc=c2, deriv=2) 
  psideriv1.c2 <- tukeyChi(rstar2, cc=c2, deriv=2) 
  
  #compute M 
  weightsM <- psideriv1.c2 / sigma             
  WM <- weightsM %*% rep(1, p)
  xwM <- x * WM
  M <- t(x) %*% xwM / n
  
  #compute Q 
  weightsQ <- psi1.c2^2
  WQ <- weightsQ %*% rep(1, p)
  xwQ <- x * WQ
  Q <- t(x) %*% xwQ / n
  
  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)
  V22 <- Amat %*% V %*% t(Amat)
  #inverse information matrix
  #V.inv <- solve(V)
  #submatrix of V
  V22.inv <- solve(V22)
#  V22.inv <- Amat %*% V.inv %*% t(Amat)
  
  M221 <- M[idx1,idx1] - M[idx1,idx0] %*% solve(M[idx0,idx0,drop=FALSE], 
                                                M[idx0,idx1,drop=FALSE])

#  p0 <- length(idx0)
#  M[(p0+1):p,(p0+1):p] - M[(p0+1):p,1:p0] %*% solve(M[1:p0,1:p0,drop=FALSE], 
#                                                    M[1:p0,(p0+1):p,drop=FALSE])

  weightsZ <- psi0.c2 
  Z <- t(x) %*% weightsZ / n  

  #Wald-type test statistic, Silvapulle (1996, eq. 2.6)
  Tn <- sqrt(n)*beta2[idx1]
  Dn <- Tn
  Dmat <- V22.inv
  dvec <- t(Dn)%*%V22.inv
  #constrained optimization
  out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat%*%t(Amat), bvec=bvec, meq=meq) 
  b <- out$solution
  TS_W2 <- (t(Dn)%*%Dmat%*%Dn)   
  TS_W  <- (t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b)) 
  result_W <- n * beta2[idx1] %*% solve(V22, beta2[idx1])
  
  #Score-type test (Silvapulle, 1996, eq. 2.6)
  An <- sqrt(n) * solve(M221) %*% Z[idx1]     
  Dn <- An
  dvec <- t(Dn)%*%V22.inv
  out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat%*%t(Amat), bvec=bvec, meq=meq)   
  b <- out$solution
  TS_S2 <- (t(Dn)%*%Dmat%*%Dn)
  TS_S <- (t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b))

  #Score-type test, extention Markautou and Hettmansperger (1990) (see Silvapulle, 1996)
  #result.u <- M221 %*% V[idx1,idx1] %*% t(M221)
  result.u <- Q[idx1,idx1]-M[idx1,idx0]%*% solve(M[idx0,idx0,drop=F], Q[idx0,idx1,drop=F]) - 
              Q[idx1,idx0]%*% solve(M[idx0,idx0,drop=F], M[idx0,idx1,drop=F])+
              M[idx1,idx0]%*% solve(M[idx0,idx0,drop=F], Q[idx0,idx0,drop=F])%*%
               solve(M[idx0,idx0,drop=F], M[idx0,idx1,drop=F])
  
  U <- sqrt(n) * solve(M221) %*% Z[idx1]     
  A.tilde <- solve(M221) %*% result.u %*% t(solve(M221))
  Dmat <- A.tilde
  dvec <- t(U)%*%A.tilde
  out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat%*%t(Amat), bvec=bvec, meq=meq)   
  b <- out$solution
  TS_MH <- (t(U)%*%Dmat%*%U) - (t(U-b)%*%Dmat%*%(U-b))
  
  return(list(TS_W = TS_W,          
              TS_W2 = TS_W2,
              result_W = result_W,
              TS_S = TS_S,
              TS_MH = TS_MH,
              TS_S2 = TS_S2,
              V22 = V22
             )
        )
  
}

