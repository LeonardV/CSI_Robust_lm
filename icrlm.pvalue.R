#compute p value icrlm.test
icrlm.pvalue <- function(object, ...) {

  if (!"icrlm.test" %in% class(object)) {
    stop("object must be of class icrlm.test.")
  }
  
  type <- object$type
  Ts.org <- object$Ts
  Ts.name <- names(Ts.org)
  Amat <- object$Amat
  meq <- object$meq
  meq.alt <- object$meq.alt
  mm <- model.matrix(object$object.rlm)[,,drop=FALSE]
  df.residual <- nrow(mm)-ncol(mm) 
  cov <- object$cov
  
  if (type != "C") {
    if ((qr(Amat)$rank < nrow(Amat))) {
      stop(paste("constraint matrix must have full row-rank"))
    }
      wt.bar <- 
        if (Ts.name %in% c("Wald2", "Score")) {
          wt(Amat%*%cov%*%t(Amat), meq=meq)
        } 
        else if (Ts.name %in% c("Wald1", "Fm")) {
          wt(Amat%*%solve(cov)%*%t(Amat), meq=meq) 
        }
  }
    
  if (type %in% c("A","B")) {
    if(type == "A") {
      # compute df
      df.bar <- 0:(nrow(Amat) - meq)
      # p value based on F-distribution
      pf.value <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual, 
                            wt = rev(wt.bar))
      # p value based on chi^2-distribution
      pchi.value <- 1 - pchibar(Ts.org, df1 = df.bar, wt = rev(wt.bar))
    }
    else if (type == "B") {
        # compute df
        df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
        # p value based on F-distribution
        pf.value <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual, 
                              wt = wt.bar)
        # p value based on chi^2-distribution
        pchi.value <- 1 - pchibar(Ts.org, df1 = df.bar, wt = wt.bar)
    }
  }
  else if (type == "C") {
    # t-distribution
    pt.value <- 1 - pt(Ts.org, df.residual)
  }

  OUT <- list(Ts.org = Ts.org, df.residual = df.residual)
  
  if (type %in% c("A","B")) {
    OUT$df.bar <- df.bar 
    OUT$wt.bar <- wt.bar
    OUT$pf.value <- pf.value 
    OUT$pchi.value <- pchi.value
  } else if (type == "C") {
    OUT$pt.value <- pt.value
  }
              
  return(OUT)
}
