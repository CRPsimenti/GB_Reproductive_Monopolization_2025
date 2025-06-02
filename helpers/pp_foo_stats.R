pp_foo_stats <- function(model_env, standat, stat = "max", breaks) {
  if (missing(breaks)) breaks <- "Sturges"
  
  pdata <- model_env$draws("femcount_rep", format = "draws_matrix")
  pd <- apply(pdata, 2, stat)
  hdata <- hist(pd, plot = FALSE, breaks = breaks)
  
  hist(pd, main = "", ylim = c(0, max(hdata$counts) * 1.05), 
       yaxs = "i", xlab = "", ylab = "frequency", breaks = breaks)
  hist_range <- par()$usr[1:2]
  emp_val <- get(stat)(standat$model_fems)
  if (emp_val > hist_range[2]) {
    cat("empirical value is larger than plotting range")
    arrows(x0 = par()$usr[2] - mean(par()$usr[1:2])/2,
           x1 = par()$usr[2],
           y0 = par()$usr[4],
           y1 = par()$usr[4],
           col = "red",
           xpd = TRUE)
    text(x = par()$usr[2] - mean(par()$usr[1:2])/4, y = par()$usr[4] * 0.99, labels = emp_val, adj = c(0.5, 1))
  }
  abline(v = emp_val, col = "red")
  box(bty = "l")
}
