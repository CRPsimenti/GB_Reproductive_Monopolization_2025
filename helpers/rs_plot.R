#' plot elo vs num females
#' 
#' possibly separated by group
#'
#' @param model_env cmdstan model environment
#' @param standat stan data list
#' @param target_party if \code{NULL}: 'average' group, else: one of \code{"six"},
#'                     \code{"sixi"}, \code{"sixw"}, \code{"five"} or \code{"nine"}
#' @param n_samples number of posterior samples to add (default is \code{50})
#' @param ylims,xlims axis limits (default is \code{c(0, 6)}, and \code{c(-3, 3)})
#' @param prior_only logical, if \code{TRUE} plot only one model prediction (not by age)
#' @param plot_legend logical, add the legend
#' @param cols character of length 2, colors for prime-aged and non-prime-aged
#'             ids, default is \code{viridisLite::viridis(7)[c(2, 5)]}
#' @return
#' @export
#'

rs_plot <- function(model_env, 
                    standat, 
                    target_party = NULL, 
                    n_samples = 50,
                    ylims = c(0, 6), 
                    xlims = c(-3, 3), 
                    prior_only = FALSE, 
                    plot_legend = TRUE,
                    cols = c(prime = "#35B779FF", 
                             nonprime = "#443A83FF")) {
  
  obs_sel <- seq_len(standat$n_model_obs)
  n <- standat$n_model_obs
  # character markers for the selected group
  # to find the correct row in the random effects matrix
  grmarker <- NULL
  if (!is.null(target_party)) {
    obs_sel <- which(names(standat$model_party_index_per_obs) %in% target_party)
    n <- length(obs_sel)
    grmarker <- paste0("[", standat$model_party_index_per_obs[target_party][1], ",")
  }
  
  yfem <- standat$model_fems[obs_sel]
  ynonprime <- standat$nonprime[obs_sel]
  
  # plot 'raw' data ('raw' elo ratings)
  edata <- apply(model_env$draws(c("elo_pred_rep")), 3, median)[obs_sel]
  xcols <- rep(cols[2], length(edata))
  xcols[ynonprime == 0] <- cols[1] # prime males color
  plot(edata, yfem , lwd = 1.5, cex = 1.5, pch = 1, # + runif(n, -0.1, 0.1)
       xlab = "Elo-rating (posterior median)", ylab = "number of females", 
       las = 1, col = xcols, xlim = xlims, ylim = ylims)
  
  # extract plotting data for draws
  re <- model_env$draws("blups_mat_party_actual", format = "draws_matrix")
  fe <- model_env$draws(c("intercept", "elo_slope", "ageclass_slope", "ageclass_elo_interaction"), format = "draws_matrix")
  pdata <- fe
  if (!is.null(grmarker)) {
    re <- re[, grepl(grmarker, colnames(re), fixed = TRUE)]
    colnames(re) <- c("intercept", "elo_slope", "ageclass_slope", "ageclass_elo_interaction")
    pdata <- pdata + re
  } else {
    re <- fe - fe # make everything 0
  }
  
  pdata <- pdata[sample(nrow(pdata), n_samples), ]
  
  # prime males ----
  xvals1 <- seq(min(edata[ynonprime == 0]), max(edata[ynonprime == 0]), length.out = 101)
  
  apply(pdata, 1, function(x) {
    points(xvals1, exp(c(x["intercept"]) + c(x["elo_slope"]) * xvals1), 
           type = "l", lwd = 0.5, col = adjustcolor(cols[1], 0.2))
  })
  
  # non prime males ----
  if (!prior_only) {
    xvals2 <- seq(min(edata[ynonprime == 1]), max(edata[ynonprime == 1]), length.out = 101)
    if (length(unique(xvals2)) > 1) {
      apply(pdata, 1, function(x) {
        points(xvals2, exp(c(x["intercept"]) + c(x["elo_slope"]) * xvals2 + 
                             c(x["ageclass_slope"]) * 1 + c(x["ageclass_elo_interaction"]) * xvals2 * 1), 
               type = "l", lwd = 0.5, col = adjustcolor(cols[2], 0.2))
      })
    }
  }
  
  # population effect for prime males
  points(xvals1, exp(colMeans(fe + re)["intercept"] + colMeans(fe + re)["elo_slope"] * xvals1), 
         type = "l", lwd = 3, col = cols[1])
  
  # population effect for non-prime males
  if (!prior_only) {
    if (length(unique(xvals2)) > 1) {
      points(xvals2, exp(colMeans(fe + re)["intercept"] + 
                           (colMeans(fe + re)["elo_slope"]) * xvals2 + 
                           (colMeans(fe + re)["ageclass_slope"]) * 1 + 
                           (colMeans(fe + re)["ageclass_elo_interaction"]) * xvals2 * 1), 
             type = "l", lwd = 3, col = cols[2])
    }
  }
  
  if (plot_legend) {
    partystring <- paste("party", target_party)
    if (isTRUE(target_party == "sixi")) partystring <- paste("party [6I]")
    if (isTRUE(target_party == "sixw")) partystring <- paste("party [6W]")
    if (isTRUE(target_party == "six")) partystring <- paste("party [6]")
    if (isTRUE(target_party == "five")) partystring <- paste("party [5]")
    if (isTRUE(target_party == "nine")) partystring <- paste("party [9]")
    if (is.null(target_party)) partystring <- NA
    legend("topleft", legend = c("prime aged", "non-prime aged"), pch = c(1, 1), lwd = c(1, 1), bty = "n", 
           title = partystring, col = c(cols[1], cols[2]), text.col = c(cols[1], cols[2]), title.col = "black ", cex = 0.8)
  }
  
}
