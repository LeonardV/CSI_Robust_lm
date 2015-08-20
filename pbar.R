#acknowlegment: taken from package ibdreg
pfbar <- function(x, df1, df2, wt) {
  if (x <= 0) {
    return(0)
  }
  zed <- df1 == 0
  cdf <- ifelse(any(zed), wt[zed], 0)
  cdf <- cdf + sum(pf(x/df1[!zed], df1[!zed], df2) * wt[!zed])
  return(cdf)
}

pchibar <- function (x, df1, wt) 
{
#acknowlegment: taken from package ibdreg
  if (x <= 0) {
    return(0)
  }
  zed <- df1 == 0
  cdf <- ifelse(any(zed), wt[zed], 0)
  cdf <- cdf + sum(pchisq(x, df1[!zed]) * wt[!zed])
  return(cdf)
}



qchibar <- function (x, df1, wt, crit) 
{
  #acknowlegment: taken from package ibdreg
  if (x <= 0) {
    return(0)
  }
  zed <- df1 == 0
  cdf <- ifelse(any(zed), wt[zed], 0)
  cdf <- cdf + sum(pchisq(x, df1[!zed]) * wt[!zed])
  cdf <- (cdf-(1-crit))^2
  return(cdf)
}

# provides critical value given alpha=.05
#for (i in 1:length(wt)) {
#  startv=0+wt[i]*df[i]
#}

#optim(startv, qchibar, method="BFGS", df1=df, wt=wt, crit=0.05)
