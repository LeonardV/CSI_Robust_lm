wt <- function(cov, meq) {
  if (meq == 0L) {
    wt.bar <- ic.infer:::ic.weights(cov)
  }
  else if (meq > 0) {
    wt.bar <- ic.infer:::ic.weights(solve(solve(cov)[-(1:meq), -(1:meq)]))
  }
  wt.bar
}

