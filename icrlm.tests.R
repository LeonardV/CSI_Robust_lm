icrlm.tests <- function(object, type="A", meq.alt = 0, tol = 1e-04, test="Fm") {
  
  type <- toupper(type)
  test <- tolower(test)
  if (!"icrlm" %in% class(object)) {
    stop("object must be of class icrlm.")
  }
  if(!(type %in% c("A","B","C"))) {
    stop("type must be \"A\", \"B\" or \"C\"")
  }
  if(!(test %in% c("wald1", "wald2", "score", "fm")) && type %in% c("A","B")) {
    stop("test must be \"wald1\", \"wald2\", \"score\" or \"Fm\"")
  }
  
  # original rlm object fitted by user
  object.rlm <- object$object.rlm
  if(is.null(object.rlm$model)) {
    Y <- object.rlm$residuals + object.rlm$fitted
  } else {
    mfit <- object.rlm$model
    Y <- model.response(mfit)  
  }  
  # data
  X <- model.matrix(object.rlm)[,,drop=FALSE]
  # sample size
  n <- dim(X)[1]
  # unconstraint estimates
  beta.unconstr <- coef(object.rlm)
  # constraint estimates
  beta.constr <- object$beta.constr
  # variable names
  var.names <- names(beta.unconstr)
  #scale parameter unconstrained model
  sigma <- object.rlm$s
  #standard deviation unconstrained model
  tau.hat <- MASS:::summary.rlm(object.rlm)$stddev  
  # inequality constraints matrix
  Amat <- object$Amat
  # rhs constraint matrix
  bvec <- object$bvec
  # number of equality constraints
  meq  <- object$meq
    
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality constraints only.")
  }  
  
  call.rlm <- as.list(object.rlm$call)
  call.rlm <- call.rlm[-1]
  if (call.rlm[["method"]] == "M") {
    if (is.null(call.rlm[["psi"]])) {
      stop("test only applicable with \"psi=psi.bisquare\".")
    } 
    else if (!(call.rlm[["psi"]] == "psi.bisquare")) {
      stop("test only applicable with \"psi=psi.bisquare\".")  
    }
  }
    
  # Silvapulle and Sen, 2005
  if (type == "A") {
    if (is.null(call.rlm[["formula"]])) {
      call.rlm[["data"]] <- NULL
      call.rlm[["x"]] <- NULL
      call.rlm[["y"]] <- NULL
      call.my <- list(x = X, y = Y, Amat = Amat, meq = nrow(Amat), bvec = bvec)      
      CALL <- c(call.my, call.rlm)
      if(any(duplicated(CALL))) {
        stop("duplicated elements in CALL.list")
      }
      rfit <- do.call("icrlm.fit", CALL)
    }
    else {
      call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)  
      CALL <- c(call.rlm, call.my)
      if(any(duplicated(CALL))) {
        stop("duplicated elements in CALL.list")
      }
      Data <- object.rlm$model
      call.rlm[["data"]] <- as.name("Data")
      rfit <- do.call("icrlm.formula", CALL)
    }
    # null model - equality constraints
    beta.eqconstr <- coef(rfit)
      names(beta.eqconstr) <- var.names
    # vector with the indices of the active constraints
    iact <- rfit$iact
    
    if (test %in% c("wald2","score")) {
      # robust Wald (W1, W2) and score test-statistics, Silvapulle, 1996
      wald2Score <- robustWaldScores(x = X, y = Y, beta0 = beta.eqconstr,       #FIXME for no intercept model
                                     beta2 = beta.unconstr, sigma = sigma, 
                                     Amat = Amat, bvec = bvec, meq = meq)
      cov <- wald2Score$V
      #V22 <- wald1Score$V22
      Ts <- unlist(wald2Score)[c("RW2.A", "score.A")]
        names(Ts) <- c("Wald2", "Score")
      if (test == "wald2") {
        Ts <- Ts[1]
      }
      else if(test == "score") {
        Ts <- Ts[2]
      }      
    }
    else if (test == "wald1") {
      # robust Wald, Silvapulle, 1992
      Ts <- robustWald.xx(x = X, beta0 = beta.eqconstr, beta1 = beta.constr, 
                             beta2 = beta.unconstr, tau.hat = tau.hat)
    }
    else if (test == "fm") {
      # robust F-bar, Silvapulle, 1992
      Ts <- robustFm(x = X, y = Y, beta0 = beta.eqconstr, 
                     beta1 = beta.constr, sigma = sigma, tau.hat = tau.hat, 
                     c=4.685061)
    }
  } 
  else if (type == "B") {
    if (meq.alt == 0L) {
      if (test %in% c("wald2","score")) {
        # robust Wald (W1, W2) and score test-statistics, Silvapulle, 1996
        wald2Score <- robustWaldScores(x = X, y = Y, beta0 = beta.constr, 
                                       beta2 = beta.unconstr, sigma = sigma, 
                                       Amat = Amat, bvec = bvec, meq = meq)
        cov <- wald2Score$V
        #V22 <- wald1Score$V22
        Ts <- unlist(wald2Score)[c("RW2.B", "score.B")]
          names(Ts) <- c("Wald2", "Score")
        if (test == "wald2") {
          Ts <- Ts[1]
        }
        else if(test == "score") {
          Ts <- Ts[2]
        }      
      }
      else if (test == "wald1") {
        # robust Wald, Silvapulle, 1992
        Ts <- robustWald.xx(x = X, beta0 = beta.constr, beta1 = beta.unconstr, 
                            beta2 = beta.unconstr, tau.hat = tau.hat)
      }
      else if (test == "fm") {
        # robust F-bar, Silvapulle, 1992
        Ts <- robustFm(x = X, y = Y, beta0 = beta.constr, 
                       beta1 = beta.unconstr, sigma = sigma, tau.hat = tau.hat, 
                       c=4.685061)
      }      
    } 
    else {
      # some equality may be preserved in the alternative hypothesis.           # FIXME for robustWaldScores
      if(meq.alt <= meq) {
        if(is.null(call.rlm[["formula"]])) {
          call.rlm[["data"]] <- NULL
          call.rlm[["x"]] <- NULL
          call.rlm[["y"]] <- NULL
          call.my <- list(x = X, y = Y, Amat = Amat[1:meq.alt,,drop=FALSE], 
                          meq = meq.alt, bvec = bvec[1:meq.alt])          
          CALL <- c(call.my, call.rlm)
          if(any(duplicated(CALL))) {
            stop("duplicated elements in CALL.list")
          }
          rfit <- do.call("icrlm.fit", CALL)
        }
        else {
          call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
                          bvec = bvec[1:meq.alt])          
          CALL <- c(call.rlm, call.my)
          if (any(duplicated(CALL))) {
            stop("duplicated elements in CALL.list")
          }
          rfit <- do.call("icrlm.formula", CALL)
        }        
        beta.constr.alt <- coef(rfit)
          names(beta.constr.alt) <- var.names
        iact <- rfit$iact
        
        if (test %in% c("wald2","score")) {
          # robust Wald (W1, W2) and score test-statistics, Silvapulle, 1996
          wald2Score <- robustWaldScores(x = X, y = Y, beta0 = beta.constr, 
                                         beta2 = beta.unconstr, sigma = sigma, 
                                         Amat = Amat[1:meq.alt,,drop=FALSE], 
                                         bvec = bvec[1:meq.alt], meq = meq.alt)
          cov <- wald2Score$V
          #V22 <- wald1Score$V22
          Ts <- unlist(wald2Score)[c("RW2.B", "score.B")]
            names(Ts) <- c("Wald2", "Score")
          if (test == "wald2") {
            Ts <- Ts[1]
          }
          else if (test == "score") {
            Ts <- Ts[2]
          }      
        }
        else if (test == "wald1") {
          # robust Wald, Silvapulle, 1992
          Ts <- robustWald.xx(x = X, beta0 = beta.constr, 
                              beta1 = beta.constr.alt, beta2 = beta.unconstr, 
                              tau.hat = tau.hat)
        }
        else if (test == "fm") {
          # robust F-bar, Silvapulle, 1992
          Ts <- robustFm(x = X, y = Y, beta0 = beta.constr, 
                         beta1 = beta.constr.alt, sigma = sigma, tau.hat = tau.hat, 
                         c=4.685061)
        }
      } 
      else {
        stop("meq.alt must not be larger than meq.")
      }  
    }
  }
  # Silvapulle and Sen, 2005
  else if (type == "C") {                                                       # FIXME
    if (meq == 0L) {
      # intersection-union test (Sasabuchi, 1980)
      cov <- (t(X)%*%X)                                                         #probably incorrect
      #V22 <- Amat%*%cov%*%t(Amat)
      Ts <- as.vector(min((Amat %*% beta.unconstr - bvec) / 
                            sqrt(diag(Amat%*%cov%*%t(Amat)))))
        names(Ts) <- "Tbar"
    }
    else {
      stop("test not applicable with equality constraints.")
    }
  }
  
  Ts[abs(Ts) < tol] <- 0L
  if (test %in% c("wald1", "fm") || type == "C") {
    cov <- t(X)%*%X
  }
  
  OUT <- list(type = type, object.rlm = object$object.rlm, 
              beta.constr = beta.constr,
              beta.constr.alt = if (type == "B" && meq.alt > 0L) {beta.constr.alt},
              beta.eqconstr = if (type == "A") { beta.eqconstr }, 
              cov = cov, Amat = Amat, bvec = bvec, meq = meq, 
              meq.alt = meq.alt, iact = object$iact, Ts = Ts)
  
    class(OUT) <- "icrlm.test"
  
  OUT

}
