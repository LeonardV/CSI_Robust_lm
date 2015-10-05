#compute p value icrlm.test
icrlm.pvalue <- function(object, constraints = constraints, 
                         mix.weights = c("mvtnorm", "none"),
                         R = 9999, parallel = c("no", "multicore", "snow"), 
                         ncpus = 1L, cl = NULL, seed = NULL, verbose = TRUE) {
  
  if (!"icrlm.test" %in% class(object)) {
    stop("object must be of class icrlm.test.")
  }
  
  if (missing(mix.weights)) 
    mix.weights <- getOption("mix.weights", "mvtnorm")
  mix.weights <- match.arg(mix.weights)
  
  type <- object$type
  Ts.org <- object$Ts
  Ts.name <- names(Ts.org)
  Amat <- object$Amat
  meq <- object$meq
  meq.alt <- object$meq.alt
  iact <- object$iact
  mm <- model.matrix(object$object.rlm)[,,drop=FALSE]
  df.residual <- nrow(mm)-ncol(mm) 
  cov <- object$cov
  df.bar <- NULL
  wt.bar <- NULL
  p.value <- pf.value <- pchi.value <- pt.value <- NULL
  
  if (type != "C" && mix.weights == "mvtnorm") {
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
    if (type == "A") {
      if (mix.weights == "mvtnorm") {
        # compute df
        df.bar <- 0:(nrow(Amat) - meq)
        # p value based on F-distribution
        pf.value <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual, 
                              wt = rev(wt.bar))
        # p value based on chi^2-distribution
        pchi.value <- 1 - pchibar(Ts.org, df1 = df.bar, wt = rev(wt.bar))
      } 
      else if (mix.weights == "none") {
        p.value <- icrlm.pboot(object, constraints = constraints, type = type, 
                               R = R,
                               parallel = parallel, ncpus = ncpus, 
                               cl = cl, seed = seed, verbose = verbose)
        names(p.value) <- Ts.name
      }
      else { stop("method unknown for computing p-value.") }
    }
    else if (type == "B") {          
      if (mix.weights == "mvtnorm") {
        # compute df
        df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
        # p value based on F-distribution
        pf.value <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual, 
                              wt = wt.bar)
        # p value based on chi^2-distribution
        pchi.value <- 1 - pchibar(Ts.org, df1 = df.bar, wt = wt.bar)
      }
      else if (mix.weights == "none") {
        p.value <- icrlm.pboot(object, constraints = constraints, type = type, 
                               R = R,
                               parallel = parallel, ncpus = ncpus, 
                               cl = cl, seed = seed, verbose = verbose)
        names(p.value) <- Ts.name
      }
      else { stop("method unknown for computing p-value.") }
    }
  }
  else if (type == "C") {
    # t-distribution
    pt.value <- 1 - pt(Ts.org, df.residual)
  }
  
  OUT <- list(Ts.org = Ts.org, df.residual = df.residual)
  
  if (type %in% c("A","B")) {
    OUT$df.bar <- df.bar
    OUT$iact <- iact
    OUT$wt.bar <- wt.bar
    OUT$p.value <- p.value
    OUT$pf.value <- pf.value 
    OUT$pchi.value <- pchi.value
  } else if (type == "C") {
    OUT$pt.value <- pt.value
  }
  
  return(OUT)
}
