#' posterior predictions for number of females
#'
#' @param model_env cmdstan model environment
#' @param standat stan data list
#' @param split_by_age logical, do age classes separated or combined
#' @param ytop numeric, upper y-axis limit
#' @param n_samples numeric, the number of posterior samples (default is \code{50})
#' @param target_party character, the party to plot. At its default \code{NULL}: combine all parties.
#' @param target_year character, the year to plot. At its default \code{NULL}: combine all years.
#' @param target_male character, the male to plot. At its default \code{NULL}: combine all males.
#' @param ignore_setpar logical, ignore the setting of \code{par(mfrow = c(1, 2))} that
#'        is executed inside the function
#' @return
#' @export
#'

pp_foo <- function(model_env, 
                   standat, 
                   split_by_age = FALSE, 
                   ytop = 100,
                   n_samples = 50,
                   target_party = NULL,
                   target_year = NULL, 
                   target_male = NULL,
                   ignore_setpar = FALSE) {
  
  if (!ignore_setpar) {
    op <- par(no.readonly = TRUE)
  }
  
  obs_sel <- seq_len(standat$n_model_obs)
  n <- standat$n_model_obs
  
  if (is.null(target_party)) {
    target_party <- unique(names(standat$model_party_index_per_obs))
  }
  if (is.null(target_year)) {
    target_year <- unique(names(standat$model_year_index_per_obs))
  }
  if (is.null(target_male)) {
    target_male <- unique(names(standat$model_male_index_per_obs))
  }
  
  obs_sel <- which(names(standat$model_party_index_per_obs) %in% target_party &
                     names(standat$model_year_index_per_obs) %in% target_year &
                     names(standat$model_male_index_per_obs) %in% target_male)
  
  n <- length(obs_sel)
  if (n == 0) {
    stop("didn't find any values for the specified combination of party, year and male")
  }
  
  
  yfem <- standat$model_fems[obs_sel]
  ynonprime <- standat$nonprime[obs_sel]
  
  
  if (!split_by_age) {
    hist(yfem, xlim = c(-0.5, 10), breaks = seq(-0.5, 50.5, 1), ylim = c(0, ytop), yaxs = "i", xlab = "number of females", las = 1, main = "")
    pdata <- r$draws("femcount_rep", format = "draws_matrix")[, obs_sel]
    pdata <- lapply(apply(pdata[sample(1:nrow(pdata), n_samples), ], 1, hist, plot = FALSE, breaks = seq(-0.5, 20.5, 1)), function(x)x$counts)
    pdata <- do.call("rbind", pdata)
    
    points(seq(0, 20, by = 1), apply(pdata, 2, median), pch = 16, col = "red", cex = 1.5)
    segments(seq(0, 20, by = 1), apply(pdata, 2, quantile, prob = 0.1), seq(0, 20, by = 1), apply(pdata, 2, quantile, prob = 0.9), col = "red", lwd = 3)
    segments(seq(0, 20, by = 1), apply(pdata, 2, min), seq(0, 20, by = 1), apply(pdata, 2, max), col = "red", xpd = TRUE)
    box(bty = "l")
  }
  
  if (split_by_age) {
    if (!ignore_setpar) {
      par(mfrow = c(1, 2))
    }
    
    for (i in c(0, 1)) {
      hist(yfem[ynonprime == i], xlim = c(-0.5, 10), breaks = seq(-0.5, 50.5, 1), 
           ylim = c(0, ytop), yaxs = "i", xlab = "number of females", las = 1, main = "")
      pdata <- r$draws("femcount_rep", format = "draws_matrix")[, obs_sel][, ynonprime == i]
      pdata <- lapply(apply(pdata[sample(1:nrow(pdata), n_samples), ], 1, hist, plot = FALSE, breaks = seq(-0.5, 20.5, 1)), function(x)x$counts)
      pdata <- do.call("rbind", pdata)
      
      points(seq(0, 20, by = 1), apply(pdata, 2, median), pch = 16, col = "red", cex = 1.5)
      segments(seq(0, 20, by = 1), apply(pdata, 2, quantile, prob = 0.1), 
               seq(0, 20, by = 1), apply(pdata, 2, quantile, prob = 0.9), 
               col = "red", lwd = 3)
      segments(seq(0, 20, by = 1), apply(pdata, 2, min), 
               seq(0, 20, by = 1), apply(pdata, 2, max),
               col = "red", xpd = TRUE)
      legend("topright", legend = NA, title = c("prime", "non-prime")[i + 1], bty = "n")
      box(bty = "l")
    }
  }
  
  if (!ignore_setpar) par(op)
}
