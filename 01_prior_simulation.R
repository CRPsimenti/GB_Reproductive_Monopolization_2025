source("00_dataprep.R")
source("helpers/rs_plot.R")

# simulate 'posteriors' when using just priors (and not conditioning on data)
# the elo part is actually properly modelled from the data,
# it's just the count part (number of females) that is simulated

m <- cmdstanr::cmdstan_model("stan_models/rs_model.stan")
m$check_syntax(pedantic = TRUE)

standat$do_model <- 0 # indicator for prior simulation
r <- m$sample(data = standat, parallel_chains = 4, refresh = 100, show_exceptions = FALSE, seed = 1)
# optional directory for storing results:
# if (!dir.exists("results")) dir.create("results")
# store results
# r$save_object("results/finalmodel_priors.rds")

# plot in the supplement (figure A3)
set.seed(42)
rs_plot(model_env = r, standat = standat, target_party = NULL, n_samples = 50, 
        ylims = c(0, 6), xlims = c(-3, 3), prior_only = TRUE)


