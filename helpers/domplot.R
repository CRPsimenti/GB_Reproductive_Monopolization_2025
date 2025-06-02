dom_plot <- function(model_env, 
                     standat,
                     party = "six",
                     year = 2014:2017,
                     xlim = c(-5, 5),
                     cols = c(prime = "#35B779FF", 
                              nonprime = "#443A83FF")
                     ) {
  
  sel <- which(names(standat$model_year_index_per_obs) %in% year & 
                 names(standat$model_party_index_per_obs) %in% party)
  d <- model_env$draws("elo_pred_rep", format = "draws_matrix")[, sel]
  nonprime <- standat$nonprime[sel]
  p <- apply(d, 2, density)
  post_med <- apply(d, 2, median)
  vert_lim <- range(unlist(lapply(p, function(x)range(x$x))))
  
  plot(0, 0, type = "n", ylim = vert_lim, xlim = xlim, axes = FALSE,
       xlab = "density", ylab = "Elo-rating", yaxs = "i")
  abline(v = 0)
  text(-0.7, vert_lim[1]*0.85, "non-prime aged", adj = 1, xpd = TRUE)
  text(0.7, vert_lim[1]*0.85, "prime aged", adj = 0, xpd = TRUE)
  
  for (i in seq_along(p)) {
    
    pdata <- cbind(p[[i]]$x, p[[i]]$y)
    xcol <- cols[1]
    if (nonprime[i] == 1) {
      pdata[, 2] <- pdata[, 2] * (-1)
      # p[[i]]$y <- p[[i]]$y * (-1)
      xcol <- cols[2]
    }
    
    xvals <- c(pdata[, 2], rep(0, length(pdata[, 2])))
    yvals <- c(pdata[, 1], rev(pdata[, 1]))
    
    polygon(x = xvals, y = yvals, col = adjustcolor(xcol, 0.5), lwd = 0.5, border = NA)
  }

  for (i in seq_along(p)) {
    if (nonprime[i] == 1) {
      segments(x0 = 0, y0 = post_med[i], x1 = -1, y1 = post_med[i], col = "black")
    }
    if (nonprime[i] == 0) {
      segments(x0 = 0, y0 = post_med[i], x1 = 1, y1 = post_med[i], col = "black")
    }
  }
  
  abline(v = 0)
  axis(2, las = 1)
}
