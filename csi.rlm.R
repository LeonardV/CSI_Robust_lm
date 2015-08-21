csi.rlm <- function(object=NULL, Amat=NULL, bvec=NULL, meq=0) { 
  
  #  if(class(object)[1] != "rlm") {
  #    stop("It only works for rlm()")
  #  }
  
  if (is.null(Amat)) {
    stop("no constraints matrix has been specified.")
  }
  if (meq == nrow(Amat)) {
    stop("test not applicable with equality restrictions only.")
  }
  if (is.null(bvec)) {
    bvec <- rep(0, nrow(Amat))
  }
  if (!is.vector(bvec)) {
    stop("bvec must be a vector.")
  }
  
  
  pvalues_chi <- pvalues_f <- statistics <- rep(as.numeric(NA), 4)
  convergence <- vector("numeric", 3)
  names(pvalues_chi) <- names(pvalues_f) <- names(statistics) <- c("W1", "S", "Fm", "W2")
  
  mfit <- object$model
  Y <- model.response(mfit)
  X <- model.matrix(object)[,,drop = FALSE]
  
  n = nrow(X)
  b <- coef(object)
  p <- length(b)
  
  
  #Fbar test-statistic
  #  fit.clm <- csi.lm(fit.lm, Amat=Amat, bvec=bvec, meq=meq, overall=FALSE,
  #                    pvalue=TRUE, mix.weights="mvtnorm") 
  
  #  pvalues_f[6] <- pvalues_chi[6] <- fit.clm$p.value[1]
  #  statistics[6] <- fit.clm$T.obs[1]
  
  weights <- object$call$weights
  if(!length(weights)) weights <- rep(1, nrow(X))
  w <- object$call$w
  if(!length(w)) w <- rep(1, nrow(X))
  init <- object$call$init
  if(is.null(init)) init <- "ls"
  psi <- MASS:::psi.bisquare
  scale.est <- object$call$scale.est
  k2 <- object$call$k2
  if(is.null(k2)) k2 <- 1.345
  method <- object$call$method
  wt.method <- object$call$wt.method
  if(is.null(wt.method)) wt.method <- "inv.var"
  maxit <- object$call$maxit
  if(is.null(maxit)) maxit <- 5000
  acc <- object$call$acc
  if(is.null(acc)) acc <- 1e-09
  test.vec <- object$call$test.vec
  if(is.null(test.vec)) test.vec <- "resid"
  lqs.control <- object$call$lqs.control
  
  #MM-estimator
  
  rfit0 <- crlm(x=X, y=Y, Amat=Amat, meq=nrow(Amat), bvec=bvec, weights=weights, w=w, 
                init=init, psi=psi, scale.est=scale.est, k2=k2, method=method, wt.method=wt.method,
                maxit=maxit, acc=acc, test.vec=test.vec, lqs.control=lqs.control) 
    
  rfit1 <- crlm(x=X, y=Y, Amat=Amat, meq=meq, bvec=bvec, weights=weights, w=w, 
                init=init, psi=psi, scale.est=scale.est, k2=k2, method=method, wt.method=wt.method,
                maxit=maxit, acc=acc, test.vec=test.vec, lqs.control=lqs.control)
  
  rfit2 <- crlm(x=X, y=Y, weights=weights, w=w, init=init, psi=psi, scale.est=scale.est, k2=k2, 
                method=method, wt.method=wt.method, maxit=maxit, acc=acc, test.vec=test.vec, 
                lqs.control=lqs.control)
  
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
  TS <- robustWaldScoresLRT(y=Y, x=X, beta0=beta0, betaA=beta2, sigma=sigma, 
                            Amat=Amat, bvec=bvec, meq=meq)
  statistics[c(1,2)] <- unlist(TS)[c("TS_W", "TS_S")]
  
  cov <- as.matrix(TS$V22)
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
  df1 <- 0:(nrow(Amat) - meq)
  df.error <- nrow(X)-ncol(X)
  #df2 <- meq:nrow(Amat)
  
  pvalues_chi[c(1,2)] <- sapply(statistics[c(1,2)], function(x) 1-pchibar(x, df1=df1, wt=rev(wt_bar))) 
  #null distribution mixture of F distributions
  pvalues_f[c(1,2)] <- sapply(statistics[c(1,2)], function(x) 1-pfbar(x, df1=df1, df2=df.error, wt=rev(wt_bar)))  
  
  ## robust Fm test statistic ##
  #tuning constant bisquare weight
  cc <- 4.685061
  #compute residuals under h0 and h1
  resid0 <- Y - X%*%beta0
  resid1 <- Y - X%*%beta1
  
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
  statistics[3] <- 1/l.h1 * (L0 - L1) 
  wt_bar <- wt(cov=Amat%*%(solve(t(X)%*%X))%*%t(Amat), meq=meq)
  pvalues_chi[3] <- 1-pchibar(x=statistics[3], df1=df1, wt=rev(wt_bar))
  pvalues_f[3] <- 1-pfbar(x=statistics[3], df1=df1, df2=df.error, wt=rev(wt_bar))  
  
  ## robust Wald statistic, Silvapulle (1992) ##
  niX <- X
  statistics[4] <- ( (t(beta2-beta0)%*%(t(niX)%*%niX)%*%(beta2-beta0)) - 
                       (t(beta2-beta1)%*%(t(niX)%*%niX)%*%(beta2-beta1)) ) / tau.hat^2 
  pvalues_chi[4] <- 1-pchibar(x=statistics[4], df1=df1, wt=rev(wt_bar))
  pvalues_f[4]   <- 1-pfbar(x=statistics[4], df1=df1, df2=df.error, wt=rev(wt_bar))
  
  OUT <- data.frame(statistics=statistics, pvalues_chi=pvalues_chi, pvalues_f=pvalues_f)
  
  OUT
}
