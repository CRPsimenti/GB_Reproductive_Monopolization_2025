source("helpers/make_skew_data.R")
# requires SkewCalc package further down for supplem. analyses
# (from https://github.com/ctross/SkewCalc)


# data prep ----
# read raw data
xdata <- read.csv("data/skewdata.csv")

# prep data group-wise for model fitting
dat5 <- make_skew_data(xdata = xdata, party = "5", what = "raw")
dat6 <- make_skew_data(xdata = xdata, party = "6", what = "raw")
dat9 <- make_skew_data(xdata = xdata, party = "9", what = "raw")

# primary analyses -------

# compile model
skewmod <- cmdstanr::cmdstan_model("stan_models/mindex_unconstrained.stan")

# model fits
r5 <- skewmod$sample(data = dat5, parallel_chains = 4, seed = 1, refresh = 0, show_exceptions = FALSE)
r6 <- skewmod$sample(data = dat6, parallel_chains = 4, seed = 1, refresh = 0, show_exceptions = FALSE)
r9 <- skewmod$sample(data = dat9, parallel_chains = 4, seed = 1, refresh = 0, show_exceptions = FALSE)

# table 1 in main text ----
# monopoly and equality values along descriptive data for table
skewdata <- aggregate(xdata$offspring, by = list(xdata$party), function(x)c(males = length(x)))
colnames(skewdata) <- c("party", "males")
skewdata$sires <- aggregate(xdata$offspring, by = list(xdata$party), function(x)c(males = length(x[x>0])))$x
skewdata$offspring <- aggregate(xdata$offspring, by = list(xdata$party), function(x)c(males = sum(x)))$x

skewdata$postmedian <- NA
skewdata$ci <- NA

skewdata$pointestimate <- NA
skewdata$monop <- NA
skewdata$equal <- NA

skewdata$pointestimate <- sapply(c("5", "6", "9"), function(x) {
  aux <- make_skew_data(xdata = xdata, party = x)
  sprintf(SkewCalc::M_index(r = aux$r, t = aux$t), fmt = "%.2f")
})

skewdata$monop <- sapply(c("5", "6", "9"), function(x) {
  aux <- make_skew_data(xdata = xdata, party = x, what = "monop")
  sprintf(SkewCalc::M_index(r = aux$r, t = aux$t), fmt = "%.2f")
})

skewdata$equal <- sapply(c("5", "6", "9"), function(x) {
  aux <- make_skew_data(xdata = xdata, party = x, what = "equal")
  sprintf(SkewCalc::M_index(r = aux$r, t = aux$t), fmt = "%.2f")
})


for (i in 1:3) {
  if (i == 1) p <- c(r5$draws("M"))
  if (i == 2) p <- c(r6$draws("M"))
  if (i == 3) p <- c(r9$draws("M"))
  
  skewdata$postmedian[i] <- sprintf(median(p), fmt = "%.2f")
  ci1 <- sprintf(quantile(p, 0.055), fmt = "%.2f")
  ci2 <- sprintf(quantile(p, 0.945), fmt = "%.2f")
  skewdata$ci[i] <- paste(ci1, ci2, sep = " -- ")
}

skewdata


# figure 2 ----
# M posteriors
p9 <- c(r9$draws("M", format = "draws_matrix"))
p9d <- density(p9)
p6 <- c(r6$draws("M", format = "draws_matrix"))
p6d <- density(p6)
p5 <- c(r5$draws("M", format = "draws_matrix"))
p5d <- density(p5)

par(family = "serif", las = 1, mar = c(2.5, 2.5, 1, 1), mgp = c(1.5, 0.6, 0), tcl = -0.2)
plot(0, 0, xlim = c(-2, 5), ylim = c(0, 0.8), yaxs = "i", type = "n", axes = FALSE, xlab = "M", ylab = "density")
abline(v = 0, lty = 3, col = "grey")
axis(1)
box(bty = "l")
polygon(p5d, border = grey(0), lwd = 0.5, col = adjustcolor("red", 0.4))
polygon(p6d, border = grey(0), lwd = 0.5, col = adjustcolor("blue", 0.4))
polygon(p9d, border = grey(0), lwd = 0.5, col = adjustcolor("black", 0.4))
points(c(0.4, 1.5), c(0.65, 0.7), type = "l")
text(1.6, 0.7, "party 5", adj = 0)
points(c(1.3, 2.0), c(0.23, 0.3), type = "l")
text(2.1, 0.3, "party 6", adj = 0)
points(c(-0.8, -1.2), c(0.15, 0.2), type = "l")
text(-1.2, 0.23, "party 9", adj = 0.5)
segments(median(p5), 0, median(p5), p5d$y[which.min(abs(p5d$x - median(p5)))], col = "red", lwd = 3)
segments(median(p6), 0, median(p6), p6d$y[which.min(abs(p6d$x - median(p6)))], col = "blue", lwd = 3)
segments(median(p9), 0, median(p9), p9d$y[which.min(abs(p9d$x - median(p9)))], col = "black", lwd = 3)

# legend("topright", col = c("red", "blue", "black"), lty = 1, legend = c(5, 6, 9), bty = "n")



# supplements -----

## PP check ----

par(mfrow = c(1, 3), family = "serif", las = 1)
set.seed(1)
for (i in 1:3) {
  if (i == 1) {
    p <- r5$draws("r_rep", format = "draws_matrix")
    dat <- dat5
  }
  if (i == 2) {
    p <- r6$draws("r_rep", format = "draws_matrix")
    dat <- dat6
  }
  if (i == 3) {
    p <- r9$draws("r_rep", format = "draws_matrix")
    dat <- dat9
  }
  
  plot(seq_along(dat$r), dat$r, ylim = c(0, 15), type = "n", xlab = "", ylab = "offspring", las = 1, axes = FALSE)
  axis(1, at = seq_along(dat$r), labels = names(dat$r), tcl = 0, cex.axis = 0.8)
  axis(2)
  box(bty = "l")
  sel <- sample(nrow(p), 50)
  pdata <- apply(p[sel, ], 2, table)
  for (k in seq_len(length(pdata))) {
    points(rep(k, length(pdata[[k]])), as.numeric(names(pdata[[k]])), cex = sqrt(pdata[[k]] * 0.5), 
           pch = 21, bg = "white", col = "black", xpd = TRUE, lwd = 0.5) # bg = grey(0.5), 
  }
  
  # for (i in sel) points(seq_along(p[1, ]), p[i, ] + runif(length(p[i, ]), 0, 1), pch = 16, col = grey(0.2, 0.5))
  points(seq_along(dat$r), dat$r, pch = 4, col = "red", cex = 2, lwd = 3)
  # segments(seq_along(dat$r) - 0.4, dat$r, seq_along(dat$r) + 0.4, dat$r, lwd = 2, col = "red", cex = 2, xpd = TRUE)
  
  if (i == 1) {
    legend("topleft", pch = 21, pt.lwd = 0.5, pt.cex = sqrt(c(1, 5, 10) * 0.5), legend = c(1, 5, 10), bty = "n")
  }
}


## comparison of original stan model with our adaptation using SkewCalc package data ----

# we need the original implementation for validation purposes 
skewmod_ori <- cmdstanr::cmdstan_model("stan_models/mindex_original.stan")

data(ColombiaRS, package = "SkewCalc")
data(KipsigisFemales, package = "SkewCalc")
data(KipsigisMales, package = "SkewCalc")

xdata <- ColombiaRS[ColombiaRS$group == "AFROCOLOMBIAN" & ColombiaRS$sex == "F", ]
s <- skewmod$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_colombia_1 <- c(s$draws("M"))
s <- skewmod_ori$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_orig_colombia_1 <- c(s$draws("M"))


xdata <- ColombiaRS[ColombiaRS$group == "AFROCOLOMBIAN" & ColombiaRS$sex == "M", ]
s <- skewmod$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_colombia_2 <- c(s$draws("M"))
s <- skewmod_ori$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_orig_colombia_2 <- c(s$draws("M"))


xdata <- ColombiaRS[ColombiaRS$group == "EMBERA" & ColombiaRS$sex == "F", ]
s <- skewmod$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_colombia_3 <- c(s$draws("M"))
s <- skewmod_ori$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_orig_colombia_3 <- c(s$draws("M"))


xdata <- ColombiaRS[ColombiaRS$group == "EMBERA" & ColombiaRS$sex == "M", ]
s <- skewmod$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_colombia_4 <- c(s$draws("M"))
s <- skewmod_ori$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_orig_colombia_4 <- c(s$draws("M"))


xdata <- KipsigisFemales
s <- skewmod$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_colombia_5 <- c(s$draws("M"))
s <- skewmod_ori$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_orig_colombia_5 <- c(s$draws("M"))

xdata <- KipsigisMales
s <- skewmod$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_colombia_6 <- c(s$draws("M"))
s <- skewmod_ori$sample(data = list(N = nrow(xdata), r = xdata$rs, t = xdata$age, t0 = rep(0, nrow(xdata))), parallel_chains = 4, seed = 1)
supp_m_orig_colombia_6 <- c(s$draws("M"))


x <- list(supp_m_colombia_1 = supp_m_colombia_1, supp_m_orig_colombia_1 = supp_m_orig_colombia_1,
          supp_m_colombia_2 = supp_m_colombia_2, supp_m_orig_colombia_2 = supp_m_orig_colombia_2,
          supp_m_colombia_3 = supp_m_colombia_3, supp_m_orig_colombia_3 = supp_m_orig_colombia_3,
          supp_m_colombia_4 = supp_m_colombia_4, supp_m_orig_colombia_4 = supp_m_orig_colombia_4,
          supp_m_colombia_5 = supp_m_colombia_5, supp_m_orig_colombia_5 = supp_m_orig_colombia_5,
          supp_m_colombia_6 = supp_m_colombia_6, supp_m_orig_colombia_6 = supp_m_orig_colombia_6)

# visualize posteriors

par(mfcol = c(2, 3))
plot(0, 0, xlim = c(-0.25, 0.8), ylim = c(0, 20), type = "n", yaxs = "i", ylab = "density", xlab = "M", axes = FALSE)
polygon(density(x$supp_m_colombia_1))
polygon(density(x$supp_m_orig_colombia_1), border = "red")
axis(1)
box(bty = "l")
legend("topleft", col = c(NA, "red", "black"), lty = c(NA, 1, 1), legend = c("Afrocolombian females", "original", "modified"), bty = "n", cex = 0.6)

plot(0, 0, xlim = c(-0.25, 0.8), ylim = c(0, 10), type = "n", yaxs = "i", ylab = "density", xlab = "M", axes = FALSE)
polygon(density(x$supp_m_colombia_2))
polygon(density(x$supp_m_orig_colombia_2), border = "blue")
axis(1)
box(bty = "l")
legend("topleft", col = c(NA, "blue", "black"), lty = c(NA, 1, 1), legend = c("Afrocolombian males", "original", "modified"), bty = "n", cex = 0.6)


plot(0, 0, xlim = c(-0.25, 0.8), ylim = c(0, 20), type = "n", yaxs = "i", ylab = "density", xlab = "M", axes = FALSE)
polygon(density(x$supp_m_colombia_3))
polygon(density(x$supp_m_orig_colombia_3), border = "red")
axis(1)
box(bty = "l")
legend("topleft", col = c(NA, "red", "black"), lty = c(NA, 1, 1), legend = c("Embera females", "original", "modified"), bty = "n", cex = 0.6)

plot(0, 0, xlim = c(-0.25, 0.8), ylim = c(0, 10), type = "n", yaxs = "i", ylab = "density", xlab = "M", axes = FALSE)
polygon(density(x$supp_m_colombia_4))
polygon(density(x$supp_m_orig_colombia_4), border = "blue")
axis(1)
box(bty = "l")
legend("topleft", col = c(NA, "blue", "black"), lty = c(NA, 1, 1), legend = c("Embera males", "original", "modified"), bty = "n", cex = 0.6)

plot(0, 0, xlim = c(-0.25, 0.8), ylim = c(0, 20), type = "n", yaxs = "i", ylab = "density", xlab = "M", axes = FALSE)
polygon(density(x$supp_m_colombia_5))
polygon(density(x$supp_m_orig_colombia_5), border = "red")
axis(1)
box(bty = "l")
legend("topleft", col = c(NA, "red", "black"), lty = c(NA, 1, 1), legend = c("Kipsigis females", "original", "modified"), bty = "n", cex = 0.6)

plot(0, 0, xlim = c(-0.25, 0.8), ylim = c(0, 10), type = "n", yaxs = "i", ylab = "density", xlab = "M", axes = FALSE)
polygon(density(x$supp_m_colombia_6))
polygon(density(x$supp_m_orig_colombia_6), border = "blue")
axis(1)
box(bty = "l")
legend("topleft", col = c(NA, "blue", "black"), lty = c(NA, 1, 1), legend = c("Kipsigis males", "original", "modified"), bty = "n", cex = 0.6)

