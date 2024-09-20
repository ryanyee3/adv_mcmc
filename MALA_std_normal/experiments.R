# Author: Ryan Yee
# Date: September 19, 2024
# Purpose: standard normal MALA sampler experiments
# Details: 
# Dependencies: mcmcse, dplyr, plyr

source("functions.R")

# settings
lambda = c(0.25, 0.5, 1, 2)
dim = c(1, 10, 100)
iter = seq(1000, 10000, by = 1000) # c(1000, 2000, 3000, 4000, 5000)
settings = expand.grid(
  lambda = lambda,
  dim = dim,
  iter = iter
)

# experiments
results = list()
for (i in 1:nrow(settings)){
  results[[i]] = mala_std_mvn_sample(
    dim = settings$dim[i], 
    lambda = settings$lambda[i], 
    n_burn = 1000, 
    n_post = settings$iter[i], 
    save_samples = FALSE
    )
}

# results
tmp = dplyr::bind_rows(results)
tmp$ess[is.nan(tmp$ess)] = 1
results = dplyr::bind_cols(settings, tmp)

# visualizations
png(filename = "mala_std_normal_diagnostics.png", height = 5, width = 5.5, units = "in", res = 300)

col_list = c("#BB5566", "#DDAA33", "#228833", "#004488")
par(oma = c(3, 1, 1, 1), mar = c(3, 3, 1, 1), mgp = c(1.8, 0.5, 0), mfrow = c(3, 3))

# effective sample size
for (d in dim){
  tmp = subset(results, dim == d)
  lin = split(x = tmp, f = tmp$lambda)
  plot(0, type = "n", xlim = c(min(iter), max(iter)), ylim = c(plyr::round_any(min(tmp$ess), 10), plyr::round_any(max(tmp$ess), 500)), xlab = "Posterior Samples", ylab = "ESS")
  for (i in 1:length(lin)){
    points(lin[[i]]$iter, lin[[i]]$ess, col = col_list[i], pch = i)
    lines(lin[[i]]$iter, lin[[i]]$ess, col = col_list[i])
  }
}

# mean bias
for (d in dim){
  tmp = subset(results, dim == d)
  tmp$mean_bias[is.nan(tmp$mean_bias)] = 1
  lin = split(x = tmp, f = tmp$lambda)
  plot(0, type = "n", xlim = c(min(iter), max(iter)), ylim = c(floor(min(tmp$mean_bias)), ceiling(max(tmp$mean_bias))), xlab = "Posterior Samples", ylab = "Sample mean bias")
  for (i in 1:length(lin)){
    points(lin[[i]]$iter, lin[[i]]$mean_bias, col = col_list[i], pch = i)
    lines(lin[[i]]$iter, lin[[i]]$mean_bias, col = col_list[i])
  }
}

# variance bias
for (d in dim){
  tmp = subset(results, dim == d)
  tmp$var_bias[is.nan(tmp$var_bias)] = 1
  lin = split(x = tmp, f = tmp$lambda)
  plot(0, type = "n", xlim = c(min(iter), max(iter)), ylim = c(floor(min(tmp$var_bias)), ceiling(max(tmp$var_bias))), xlab = "Posterior Samples", ylab = "Sample variance bias")
  for (i in 1:length(lin)){
    points(lin[[i]]$iter, lin[[i]]$var_bias, col = col_list[i], pch = i)
    lines(lin[[i]]$iter, lin[[i]]$var_bias, col = col_list[i])
  }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = lambda, col = col_list, lwd = 1, pch = c(1, 2, 3, 4), xpd = TRUE, horiz = TRUE, cex = 1.1, seg.len=1, bty = 'n', title = expression(lambda))
text(-0.65, 1, substitute(paste(bold("dim = 1"))))
text(0.05, 1, substitute(paste(bold("dim = 10"))))
text(0.75, 1, substitute(paste(bold("dim = 100"))))

dev.off()

# acceptance rate
for (d in dim){
  tmp = subset(results, dim == d)
  tmp$acceptance_rate[is.nan(tmp$acceptance_rate)] = 1
  lin = split(x = tmp, f = tmp$lambda)
  plot(0, type = "n", xlim = c(min(iter), max(iter)), ylim = c(floor(min(tmp$acceptance_rate)), ceiling(max(tmp$acceptance_rate))), xlab = "Posterior Samples", ylab = "Acceptance rate", main = paste0("dim = ", d))
  for (i in 1:length(lin)){
    points(lin[[i]]$iter, lin[[i]]$acceptance_rate, col = col_list[i], pch = i)
    lines(lin[[i]]$iter, lin[[i]]$acceptance_rate, col = col_list[i])
  }
}

