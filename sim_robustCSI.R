sim_robustCsi <- function(N=c(15,25,50), d=0, eta.vec=12, cont=0.10, cont.x=2, 
                          Amat, bvec, meq, R=9, ncpus=1, seed=123, verbose=TRUE) {

  #checks
  if(length(N) > 1) { 
    stopifnot(length(d) == 1 && length(eta.vec) == 1) 
    S <- length(N)
  } 
  else if(length(d) > 1) { 
    stopifnot(length(N) == 1 && length(eta.vec) == 1) 
    S <- length(d)
  }
  else if(length(eta.vec) > 1) { 
    stopifnot(length(d) == 1 && length(N) == 1) 
    S <- length(eta.vec)
  }
   
  result <- list()
  
  for(j in 1:S) {
    if(length(N) > 1) { n <- N[j] } else { n <- N }
    if(length(d) > 1) { betas <- c(1,1,d[j],d[j],d[j]) } else { betas <- c(1,1,d,d,d)}
    if(length(eta.vec) > 1) { eta <- eta.vec[j] } else { eta <- eta.vec }
        
    statistics <- pvalues_chi <- pvalues_f <- matrix(NA, R, 6)
      colnames(statistics) <- colnames(pvalues_chi) <- colnames(pvalues_f) <- c("W1", "W2", "S", "Fm", "MH", "Fbar")
    
    for(i in 1:R) {
      set.seed(seed+i)
      if(verbose) { cat("iteration:", j, "...", i, "\n") }
      X <- cbind(1, mvtnorm:::rmvnorm(n, mean=rep(0,4), sigma=diag(4)))
      y <- cbind(c(betas%*%t(X))) + rnorm(n)   
      
      idx <- sample(1:nrow(y), n*cont, replace=FALSE)                             
      X[idx,cont.x] <- rnorm(length(idx), 5, 0.1)  
      y[idx] <- rnorm(length(idx), eta, 0.1)  
      
      out <- robustCsi(X, y, Amat, bvec, meq)
      statistics[i,1:6]  <- out$statistics
      pvalues_chi[i,1:6] <- out$pvalues_chi
      pvalues_f[i,1:6]   <- out$pvalues_f
    }
    
    result[[j]] <- list(statistics=statistics, pvalues_chi=pvalues_chi, pvalues_f=pvalues_f)
        
  }
  
  result  
}


