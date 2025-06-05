# helper function and scripts
source("00_dataprep.R")
source("helpers/pp_foo.R")
source("helpers/pp_foo_stats.R")
source("helpers/rs_plot.R")
source("helpers/elo_illu.R")
source("helpers/domplot.R")

# compile Stan model
library(cmdstanr)
m <- cmdstanr::cmdstan_model("stan_models/rs_model.stan")
m$check_syntax(pedantic = TRUE)

# RUN THIS TO GET A FRESH SET OF MODEL RESULTS
# r <- m$sample(data = standat, parallel_chains = 4, refresh = 100, show_exceptions = FALSE, seed = 1)
# optional directory for storing results:
# if (!dir.exists("results")) dir.create("results")
# r$save_object("results/finalmodel.rds")

# read results instead of refitting
r <- readRDS("results/finalmodel.rds")



# table A5 (numeric model results) ----
xvars <- c("intercept", "elo_slope", "ageclass_slope", "ageclass_elo_interaction")
res <- data.frame(r$summary(variables = xvars))[, c("mean", "median", "rhat")]
res <- data.frame(res, t(apply(r$draws(xvars, format = "draws_matrix"), 2, quantile, probs = c(0.055, 0.945))))
res <- apply(res, 2, round, 2)
res <- res[, c("mean", "median", "X5.5.", "X94.5.", "rhat")] # 
colnames(res)[colnames(res) == "X5.5."] <- "5.5\\%"
colnames(res)[colnames(res) == "X94.5."] <- "94.5\\%"
res


# figure 1 - Elo-illustration ----
# party 6 - 2014-2017
lmat <- matrix(c(1, 6, 1, 6, 2, 6, 2, 6, 3, 6, 3, 5, 4, 5, 4, 5), nrow = 8, byrow = TRUE)
layout(lmat, widths = c(1.5, 1))

# panel A
par(family = "serif", mgp = c(1.1, 0.3, 0), tcl = -0.2, mar = c(2.5, 1, 0.3, 1))
elo_illu(model_env = r, standat = standat, party = "six", year = "2014", xlims = c(-3, 3))
text(-2.9, 2.4, "A", xpd = TRUE, cex = 1.5, font = 2)
text(-1.2, 1, "non-prime-aged", adj = 1, cex = 0.9, col = "#443A83FF", font = 2)
text(2, 0.7, "prime-aged", adj = 0, cex = 0.9, col = "#35B779FF", font = 2)

elo_illu(model_env = r, standat = standat, party = "six", year = "2015", xlims = c(-3, 3))
elo_illu(model_env = r, standat = standat, party = "six", year = "2016", xlims = c(-3, 3))
elo_illu(model_env = r, standat = standat, party = "six", year = "2017", xlims = c(-3, 3))

# panel C
irdata <- data.frame(widths = as.numeric(apply(r$draws("elo_predictor", format = "draws_matrix"), 2,
                                               \(x)diff(quantile(x, probs = c(0.045, 0.955))))))
irdata$cumint <- standat$cumint
irdata$sel <- names(standat$model_party_index_per_obs) == "six" & names(standat$model_year_index_per_obs) %in% (2014:2017)
irdata$nonprime <- standat$nonprime
irdata <- irdata[irdata$sel, ]
irdata$col <- adjustcolor(c(prime = "#35B779FF", nonprime = "#443A83FF")[irdata$nonprime + 1], 0.8)
par(family = "serif", mgp = c(1.1, 0.3, 0), tcl = -0.2, mar = c(2.5, 3, 0.8, 1), las = 1)
plot(irdata$cumint, irdata$widths, xlim = c(0, 100), ylim = c(0, 4), cex.axis = 0.9, 
     xaxs = "i", yaxs = "i", xpd = TRUE, pch = 21, bty = "l",
     lwd = 0.6, bg = irdata$col, cex = 1.3,
     # bg = grey(0.8, 0.4)
     xlab = "cumulative interactions observed", ylab = "width of 89% CI")
text(10, 4, "C", xpd = TRUE, cex = 1.5, font = 2)


# panel B
par(mar = c(3, 4, 1, 1))
dom_plot(model_env = r, standat = standat, party = "six", year = 2014:2017, xlim = c(-7, 7))
axis(1, at = seq(-6, 6, by = 2), labels = abs(seq(-6, 6, 2)))
par(mfrow = c(1, 1))


# figure 3 - model results ----
set.seed(42)
rs_plot(model_env = r, standat = standat, target_party = NULL, 
        n_samples = 50, ylims = c(0, 6), xlims = c(-3, 3), 
        plot_legend = FALSE)
text(2.65, 2.5, "prime-\naged", col = "#35B779FF", adj = 0, font = 2, cex = 0.9)
text(2.37, 0.45, "non-prime-\naged", col = "#443A83FF", adj = 0, font = 2, cex = 0.9)



# supplements ----



## figure SUPP X ----
# as figure 3, but split by party

set.seed(42)
par(mfrow = c(2, 3), family = "serif")
rs_plot(model_env = r, standat = standat, target_party = NULL, n_samples = 50, ylims = c(0, 6), xlims = c(-3, 3))
rs_plot(model_env = r, standat = standat, target_party = "five", n_samples = 50, ylims = c(0, 6), xlims = c(-3, 3))
rs_plot(model_env = r, standat = standat, target_party = "nine", n_samples = 50, ylims = c(0, 6), xlims = c(-3, 3))
rs_plot(model_env = r, standat = standat, target_party = "six", n_samples = 50, ylims = c(0, 6), xlims = c(-3, 3))
rs_plot(model_env = r, standat = standat, target_party = "sixi", n_samples = 50, ylims = c(0, 6), xlims = c(-3, 3))
rs_plot(model_env = r, standat = standat, target_party = "sixw", n_samples = 50, ylims = c(0, 6), xlims = c(-3, 3))
par(mfrow = c(1, 1))

## PP checks ----
# A4
set.seed(42)
par(mfrow = c(1, 3))
pp_foo(model_env = r, standat = standat, split_by_age = FALSE, ytop = 100, ignore_setpar = TRUE)
pp_foo(model_env = r, standat = standat, split_by_age = TRUE, ytop = 70, ignore_setpar = TRUE)
par(mfrow = c(1, 1))

# A5
set.seed(42)
par(family = "serif", mgp = c(1.5, 0.5, 0), las = 1, mar = c(3, 3, 1, 1), tcl = 0)
par(mfrow = c(1, 2))

pp_foo_stats(model_env = r, standat = standat, stat = "max", breaks = seq(0, 26, 2))
title(xlab = "maximum number of females")
pp_foo_stats(model_env = r, standat = standat, stat = "var", breaks = seq(0, 6, 0.5))
title(xlab = "variance in number of females")


