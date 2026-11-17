# Fits missBART1 and missBART2 to the global Amax data of Maire et al. (2015),
# with the settings of Section 7 of Goh, Soh, Parnell and Murphy (2024):
# 5000 burn-in and 5000 kept iterations, 100 regression trees and 50
# missingness trees. The data are not distributed with the package.
#
# Each model is run from four different seeds, because the missBART1 chains
# can settle in different modes on these data.
#
# Saves posterior summaries (no raw data) to
# vignettes/articles/results/global_amax.rds and two figures (first seed of
# each model) to vignettes/articles/figures/, which the "Real data" article
# shows.
#
# Usage, from the root of the repository:
#   Rscript Examples/rerun_global_amax.R path/to/Maire_dataset_for_Andrew.xlsx [burn] [iters]

library(missBART)
library(ggplot2)

args <- commandArgs(TRUE)
burn <- if(length(args) > 1) as.integer(args[2]) else 5000
iters <- if(length(args) > 2) as.integer(args[3]) else 5000

data <- as.data.frame(readxl::read_excel(args[1], sheet = "Data"))
y <- as.matrix(log(data[, c("SLA", "Aarea", "Narea", "Parea", "Gs")]))
x <- as.matrix(data[, which(names(data) == "PPTmean"):ncol(data)])
obs <- !is.na(y)

fit_model <- function(model, seed) {
  set.seed(seed)
  time <- system.time(
    fit <- if(model == "missBART1") {
      missBART1(x, y, n_trees = 100, burn = burn, iters = iters, thin = 1, show_progress = FALSE, save_z = FALSE)
    } else {
      missBART2(x, y, n_reg_trees = 100, n_class_trees = 50, burn = burn, iters = iters, thin = 1,
                show_progress = FALSE, save_z = FALSE)
    })[["elapsed"]]
  fitted_mean <- fitted(fit)
  pred <- simplify2array(fit$y_pred)
  lower <- apply(pred, c(1, 2), quantile, 0.025)
  upper <- apply(pred, c(1, 2), quantile, 0.975)
  imp <- imputed_values(fit)
  imp$column <- factor(imp$column, levels = colnames(y))
  per_response <- data.frame(
    response = colnames(y),
    n_missing = colSums(!obs),
    rmse_observed = sapply(seq_len(ncol(y)), function(j) sqrt(mean((y[obs[, j], j] - fitted_mean[obs[, j], j])^2))),
    coverage_95 = sapply(seq_len(ncol(y)), function(j) mean(y[obs[, j], j] > lower[obs[, j], j] & y[obs[, j], j] < upper[obs[, j], j])),
    mean_observed = colMeans(y, na.rm = TRUE),
    mean_imputed = as.numeric(tapply(imp$mean, imp$column, mean)),
    sd_observed = apply(y, 2, sd, na.rm = TRUE),
    row.names = NULL)
  plot_data <- do.call(rbind, lapply(seq_len(ncol(y)), function(j) {
    data.frame(response = colnames(y)[j], observed = y[obs[, j], j], fitted = fitted_mean[obs[, j], j],
               lower = lower[obs[, j], j], upper = upper[obs[, j], j])
  }))
  plot_data$response <- factor(plot_data$response, levels = colnames(y))
  rug <- data.frame(response = imp$column, imputed = imp$mean)
  fig <- ggplot(plot_data, aes(fitted, observed)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey75", linewidth = 0.2) +
    geom_point(size = 0.4, colour = "#1A9993") +
    geom_abline(linewidth = 0.3) +
    geom_rug(data = rug, aes(y = imputed), inherit.aes = FALSE, sides = "l", colour = "#FD7446", alpha = 0.3) +
    facet_wrap(~response, scales = "free", nrow = 1) +
    labs(x = "Posterior mean prediction", y = "Observed (log scale)", title = model) +
    theme_bw()
  if(seed == 1) ggsave(file.path("vignettes/articles/figures", paste0("amax_", model, ".png")), fig, width = 12, height = 3, dpi = 110)
  list(model = model, seed = seed, seconds = time, per_response = per_response, tree_accept = fit$tree_accept,
       sigma = Reduce("+", fit$sigma_post) / length(fit$sigma_post),
       coef = if(model == "missBART1") coef(fit),
       varimp_missingness = if(model == "missBART2") varimp(fit),
       varimp_regression = varimp(fit, "regression"),
       y_miss_accept = if(model == "missBART2") mean(fit$y_miss_accept))
}

dir.create("vignettes/articles/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("vignettes/articles/results", recursive = TRUE, showWarnings = FALSE)
runs <- expand.grid(seed = 1:4, model = c("missBART1", "missBART2"), stringsAsFactors = FALSE)
res <- parallel::mclapply(seq_len(nrow(runs)), function(i) fit_model(runs$model[i], runs$seed[i]), mc.cores = nrow(runs))
saveRDS(list(runs = res, settings = list(n = nrow(y), q = ncol(x), burn = burn, iters = iters,
                                         version = as.character(packageVersion("missBART")))),
        "vignettes/articles/results/global_amax.rds")
for(r in res) cat(sprintf("%s seed %d: %.0f s, mean imputed: %s\n", r$model, r$seed, r$seconds,
                          paste(sprintf("%.2f", r$per_response$mean_imputed), collapse = " ")))
