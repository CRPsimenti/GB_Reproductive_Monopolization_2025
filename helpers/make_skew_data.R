# a helper function to convert tabular skew data into list form for M index
# can also create the null data sets (equal reproductive share, monopoly)

make_skew_data <- function(xdata, party = c("5", "6", "9"), what = c("raw", "monop", "equal")) {
  if (missing(what)) what <- "raw"
  if (missing(party)) stop("need party 5, 6, or 9")
  dat <- list(N = sum(xdata$party == party),
              t = xdata$duration[xdata$party == party]/365.25,
              r = xdata$offspring[xdata$party == party]
  )
  dat$t0 <- rep(0, dat$N)
  names(dat$r) <- xdata$ID[xdata$party == party]
  names(dat$t) <- xdata$ID[xdata$party == party]
  names(dat$t0) <- xdata$ID[xdata$party == party]
  
  if (what == "raw") return(dat)
  
  out <- list(N = dat$N, t = rep(max(dat$t), dat$N), t0 = rep(0, dat$N))
  
  if (what == "monop") out$r <- c(sum(dat$r), rep(0, dat$N - 1))
  if (what == "equal") out$r <- rep(1, dat$N)
  out
}
