#' plot posteriors of ratings for a chosen party and year
#'
#' @param model_env cmdstan model environment
#' @param standat stan data list
#' @param party character: party to plot
#' @param year character: year to plot
#' @param xlims numeric of length 2: fix x-range (if not set: take info from data)
#' @param ymax numeric of length: the max value for the y-axis (if not set: take info from data) 
#' @param adjust adjust value for \code{density()}
#' @param do_plot logical, if FALSE: return the densities as a list without plotting anything
#' @param label_males logical, if TRUE: print labels of individuals
#'
#' @return a plot
#' @export

elo_illu <- function(model_env, standat, party, year, xlims = NULL, ymax = NULL, colorseed = 123, adjust = 1, do_plot = TRUE, label_males = FALSE) {
  
  xx <- standat$model_male_index_per_obs[names(standat$model_party_index_per_obs) == party]
  partymales <- unique(names(xx))

  # viridis colors
  mcols <- viridisLite::viridis(length(partymales) + 7)
  # cut the most yellow again
  mcols <- rev(rev(mcols)[-c(1:3)])

  # do the actual plot
  # all draws
  edraws <- model_env$draws("elo_pred_rep", format = "draws_matrix")
  
  # get x and y limits for target party across all years
  aux <- edraws[, which(names(standat$model_party_index_per_obs) == party)]
  d <- apply(aux, 2, density, adjust = adjust)
  xlims_temp <- range(unlist(lapply(d, function(x)range(x$x))))
  ymax_temp <- max(unlist(lapply(d, function(x)range(x$y))))
  if (interactive()) {
    cat("x-range across all years of party", party, "is: ", round(xlims_temp, 3), "\n")
    cat("max y values across all years of party", party, "is: ", round(ymax_temp, 3), "\n")
  }
  if (is.null(xlims)) {
    xlims <- xlims_temp
  }
  
  # relevant subset for the selected party-year
  index <- which(names(standat$model_year_index_per_obs) == year & names(standat$model_party_index_per_obs) == party)
  plotmales <- names(standat$model_male_index_per_obs)[index]
  edraws <- edraws[, index]
  colnames(edraws) <- plotmales
  
  # age classes
  is_nonprime <- standat$nonprime[index]
  names(is_nonprime) <- plotmales
  
  # assign colors based on age cat
  mm <- c()
  i=1
  for (i in seq_along(is_nonprime)) {
    if (is_nonprime[i] == 1) {
      mm <- c(mm, mcols[1])
      names(mm)[length(mm)] <- names(is_nonprime)[i]
      mcols <- mcols[-1]
    }
    if (is_nonprime[i] == 0) {
      mcols <- rev(mcols)
      mm <- c(mm, mcols[1])
      names(mm)[length(mm)] <- names(is_nonprime)[i]
      mcols <- mcols[-1]
      mcols <- rev(mcols)
    }
  }
  
  mcols <- mm
  
  # sort by column mean (so that labels can be easier arranged)
  edraws <- edraws[, order(colMeans(edraws))]
  
  if (!do_plot) {
    return(edraws)
  }
  
  # compute densities
  pdata <- apply(edraws, 2, density, adjust = adjust)
  if (is.null(ymax)) {
    ymax <- max(unlist(lapply(pdata, function(x)range(x$y))))
  }
  
  plot(0, 0, type = "n", xlim = xlims, ylim = c(0, ymax * 1.25), 
       yaxs = "i", xaxs = "i",
       xlab = paste("Elo-rating in", year),
       axes = FALSE, ylab = "")
  i=names(pdata)[1]
  for (i in names(pdata)) {
    polygon(pdata[[i]], border = NA, col = adjustcolor(mcols[i], 0.5))
    
    if (label_males) {
      text(x = pdata[[i]]$x[which.max(pdata[[i]]$y)], 
           y = ymax * c(1.05, 1.1, 1.15)[(which(names(pdata) == i) - 1) %% 3 + 1],
           labels = i, xpd = TRUE, col = mcols[i], font = c(2, 3)[is_nonprime[i] + 1])
    }
  }
  axis(1)
  abline(h = 0)
}
