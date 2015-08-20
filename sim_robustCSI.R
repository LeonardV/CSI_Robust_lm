sim_robustCsi <- function(constraints=NULL, 
                          N=c(15,25,50), d=0, eta.vec=12, cont=0.10, cont.x=2, 
                          R=999, ncpus=1, seed=123, verbose=TRUE) {
  
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
    
    statistics <- pvalues_chi <- pvalues_f <- matrix(NA, R, 4)
    colnames(statistics) <- colnames(pvalues_chi) <- colnames(pvalues_f) <- c("W1", "S", "W2", "Fm")
    
    for(i in 1:R) {
      set.seed(seed+i)
      if(verbose) { cat("iteration:", j, "...", i, "\n") }
      X <- cbind(1, mvtnorm:::rmvnorm(n, mean=rep(0,4), sigma=diag(4)))
      y <- cbind(c(betas%*%t(X))) + rnorm(n)   
      
      idx <- sample(1:nrow(y), n*cont, replace=FALSE)                             
      X[idx,cont.x] <- rnorm(length(idx), 5, 0.1)  
      y[idx] <- rnorm(length(idx), eta, 0.1)  

      Data <- data.frame(y=y, X=X[,-1])
        colnames(Data) <- c("y","x1","x2","x3","x4")
            
      object <- MASS:::rlm(y~x1+x2+x3+x4, data=Data, method="MM", maxit=5000)
      
      # build a bare-bones parameter table for this object
      lavpartable <- lav_partable(object, est = TRUE, label = TRUE)
      
      # parse the constraints
      CON <- lavaan:::lav_constraints_parse(constraints = constraints,
                                            partable = lavpartable,
                                            debug = FALSE)
      
      Amat <- res_constraints_con_amat(object, constraints = constraints)
      meq <- nrow(res_constraints_ceq_amat(object, constraints = constraints))
      bvec <- res_constraints_rhs_bvec(object, constraints=constraints)
      
      
      out <- restriktor(object, constraints=constraints, control=list(maxit=5000))
      statistics[i,1:4]  <- out$statistics
      pvalues_chi[i,1:4] <- out$pvalues_chi
      pvalues_f[i,1:4]   <- out$pvalues_f
    }
    
    result[[j]] <- list(statistics=statistics, pvalues_chi=pvalues_chi, pvalues_f=pvalues_f)
    
  }
  
  result  
}

