#acknowlegment: taken from package ibdreg
pfbarA <- function(x, df1a, df2a, wt) {
  if (x <= 0) {
    return(0)
  }
  zed <- df1a == 0
  cdf <- ifelse(any(zed), wt[zed], 0)
  cdf <- cdf + sum(pf(x/df1a[!zed], df1a[!zed], df2a) * wt[!zed])
  return(cdf)
}

pchibarA <- function (x, df1a, wt) 
{
#acknowlegment: taken from package ibdreg
  if (x <= 0) {
    return(0)
  }
  zed <- df1a == 0
  cdf <- ifelse(any(zed), wt[zed], 0)
  cdf <- cdf + sum(pchisq(x, df1a[!zed]) * wt[!zed])
  return(cdf)
}
