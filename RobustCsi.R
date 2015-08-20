robustCsi <- function(X, y, Amat, bvec, meq) { 
  
  pvalues_chi <- pvalues_f <- statistics <- rep(as.numeric(NA), 6)
  convergence <- vector("numeric", 3)
  names(pvalues_chi) <- names(pvalues_f) <- names(statistics) <- c("W1", "W2", "S", "Fm", "MH", "Fbar")
  df_error <- nrow(X) - ncol(X)
  df1a <- 0:(nrow(Amat) - meq) 
  df1b <- meq:nrow(Amat)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  #fit linear model baed on OLS estimates
  fit.lm <- lm(y~X[,-1])
  
  #Fbar test-statistic
  fit.clm <- csi.lm(fit.lm, Amat=Amat, bvec=bvec, meq=meq, overall=FALSE,
                    pvalue=TRUE, mix.weights="mvtnorm") 
  
  pvalues_f[6] <- pvalues_chi[6] <- fit.clm$p.value[1]
  statistics[6] <- fit.clm$T.obs[1]
  
  #MM-estimator
  rfit0 <- crlm(x=X, y=y, Amat=Amat, meq=nrow(Amat), bvec=bvec, method="MM", 
                maxit=5000, acc = 1e-9) 
  rfit1 <- crlm(x=X, y=y, Amat=Amat, meq=meq, bvec=bvec, method="MM", 
                maxit=5000, acc = 1e-9)
  rfit2 <- crlm(x=X, y=y, method="MM", maxit=5000, acc = 1e-9)
  
  #check for convergence
  convergence[1] <- rfit0$converged
  convergence[2] <- rfit1$converged 
  convergence[3] <- rfit2$converged
  
  #regression coefficients
  beta0 <- as.vector(round(rfit0$coefficients, 4))
  beta1 <- as.vector(round(rfit1$coefficients, 4))
  beta2 <- as.vector(round(rfit2$coefficients, 4))
  
  #scale parameter
  sigma <- rfit2$s
  #standard deviation
  tau.hat <- summary(rfit2)$stddev  
  
  #R-square and R-square adjusted (see lmrob.R)
  #    fit.rob <- lmrob(y~X[,-1])
  #    summary(fit.rob)
  
  #compute std.error
  #    se <- sqrt(diag(vcov.MM(rfit.h1)))
  
  
  ## robust Wald (W1, W2) and score test-statistics ##
  TS <- robustWaldScoresLRT(y=y, x=X, beta0, beta2, sigma, Amat, meq=meq)
  statistics[c(2,3,5)] <- unlist(TS)[c("TS_W", "TS_S", "TS_MH")]
  
  cov <- TS$V22
  
  wt <- function(cov, meq) {
    if (meq == 0L) {
      wt.bar <- ic.infer:::ic.weights(cov)
    }
    else if (meq > 0) {
      wt.bar <- ic.infer:::ic.weights(solve(solve(cov)[-(1:meq), -(1:meq)]))
    }
    wt.bar
  }
  
  wt_bar <- wt(cov, meq=meq)
  #null distribution mixture of chi^2 distributions
  pvalues_chi[c(2,3,5)] <- sapply(statistics[c(2,3,5)], function(x) 1-pchibarA(x, df1a=df1a, wt=rev(wt_bar))) 
  #null distribution mixture of F distributions
  pvalues_f[c(2,3,5)] <- sapply(statistics[c(2,3,5)], function(x) 1-pfbarA(x, df1a=df1a, df2a=df_error, wt=rev(wt_bar)))  
  
  
  ## robust Fm test statistic ##
  #tuning constant bisquare weight
  cc <- 4.685061
  #compute residuals under h0 and h1
  resid0 <- y - X%*%beta0
  resid1 <- y - X%*%beta1
  
  #residuals / scale
  rstar0 <- c(resid0/sigma)                                               
  rstar1 <- c(resid1/sigma)
  
  L0 <- sum(tukeyChi(rstar0, cc, deriv = 0))
  L1 <- sum(tukeyChi(rstar1, cc, deriv = 0))
  
  #first derivative psi function
  psi.prime.h1 <- tukeyChi(rstar1, cc, deriv=1)                     
  #second derivative psi function
  psi.prime2.h1 <- tukeyChi(rstar1, cc, deriv=2)                    
  
  #asymptotic covariance matrix standardizing constant
  l.h1 <- ( 0.5 * (1/(n - p)) * sum(psi.prime.h1^2) ) / ( (1/n)*sum(psi.prime2.h1) )
  statistics[4] <- 1/l.h1 * (L0 - L1) 
  
  wt_bar <- wt(cov=Amat%*%(solve(t(X)%*%X))%*%t(Amat), meq=meq)
  pvalues_chi[4] <- 1-pchibarA(x=statistics[4], df1a=df1a, wt=rev(wt_bar))
  pvalues_f[4] <- 1-pfbarA(x=statistics[4], df1a=df1a, df2a=df_error, wt=rev(wt_bar))  
  
  ## robust Wald statistic, Silvapulle (1992) ##
  niX <- X
  statistics[1] <- ( (t(beta2-beta0)%*%(t(niX)%*%niX)%*%(beta2-beta0)) - 
                       (t(beta2-beta1)%*%(t(niX)%*%niX)%*%(beta2-beta1)) ) / tau.hat^2 
  
  pvalues_chi[1] <- 1-pchibarA(x=statistics[1], df1a=df1a, wt=rev(wt_bar))
  pvalues_f[1]   <- 1-pfbarA(x=statistics[1], df1a=df1a, df2a=df_error, wt=rev(wt_bar))
  
  OUT <- data.frame(statistics=statistics, pvalues_chi=pvalues_chi, pvalues_f=pvalues_f)
  
  OUT
}
