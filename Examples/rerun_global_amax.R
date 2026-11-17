# Fits the models of Section 7 of Goh, Soh, Parnell and Murphy (2024) to the
# global Amax data of Maire et al. (2015): 5000 burn-in and 5000 kept
# iterations, 100 regression trees and 50 missingness trees. The data are not
# distributed with the package.
#
# Five analyses, each with several chains because the missBART1 chains can
# settle in different modes on these data:
#   missBART1            probit regression missingness model (4 chains)
#   missBART1_marginal   the same with the update of the missing responses
#                        described in the paper (4 chains)
#   missBART2            probit BART missingness model (4 chains)
#   MAR                  the same response model with the missingness assumed
#                        ignorable (include_y = FALSE; 2 chains)
#   complete_cases       multivariate BART on the complete rows only
#
# Saves posterior summaries (no raw data) to
# vignettes/articles/results/global_amax.rds and figures to
# vignettes/articles/figures/, which the "Real data" article shows.
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
responses <- colnames(y)
miss_col <- factor(responses[col(y)[is.na(y)]], levels = responses)
dir.create("vignettes/articles/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("vignettes/articles/results", recursive = TRUE, showWarnings = FALSE)

# The samples of one chain as a fit of their own.
one_chain <- function(fit, k) {
  idx <- which(fit$chain == k)
  for(nm in c("y_post", "y_pred", "y_impute", "sigma_post", "B_post", "reg_trees", "class_trees")) {
    if(!is.null(fit[[nm]])) fit[[nm]] <- fit[[nm]][idx]
  }
  fit$chain <- fit$chain[idx]
  fit$n_chains <- 1L
  fit
}
pdp_summary <- function(fit, variable, type) {
  pd <- partial_dependence(fit, variable, type = type, n_grid = 25, n_ice = 300, draws = 200)
  do.call(rbind, lapply(responses, function(r) {
    data.frame(response = r, grid = pd$grid, mean = pd$pdp[, r, "mean"], lower = pd$pdp[, r, "lower"], upper = pd$pdp[, r, "upper"])
  }))
}

analyse <- function(name) {
  set.seed(match(name, c("missBART1", "missBART1_marginal", "missBART2", "MAR", "complete_cases")))
  common <- list(burn = burn, iters = iters, thin = 1, show_progress = FALSE)
  time <- system.time(fit <- switch(name,
    missBART1 = do.call(missBART1, c(list(x, y, n_trees = 100, save_z = FALSE, keep_trees = TRUE, n_chains = 4, cores = 4), common)),
    missBART1_marginal = do.call(missBART1, c(list(x, y, n_trees = 100, save_z = FALSE, save_y_pred = FALSE,
                                                    ymis_update = "marginal", n_chains = 4, cores = 4), common)),
    missBART2 = do.call(missBART2, c(list(x, y, n_reg_trees = 100, n_class_trees = 50, save_z = FALSE, keep_trees = TRUE,
                                          n_chains = 4, cores = 4), common)),
    MAR = do.call(missBART1, c(list(x, y, n_trees = 100, include_y = FALSE, save_z = FALSE, keep_trees = TRUE,
                                    n_chains = 2, cores = 2), common)),
    complete_cases = do.call(mvBART, c(list(x[complete.cases(y), ], y[complete.cases(y), ], x_predict = x, n_trees = 100), common))
  ))[["elapsed"]]
  out <- list(name = name, seconds = time, tree_accept = fit$tree_accept)

  if(name == "complete_cases") {
    pred <- predict(fit)
    out$imputed <- data.frame(chain = 1, response = miss_col, value = pred[is.na(y)])
    out$rmse_observed <- sapply(seq_along(responses), function(j) sqrt(mean((y[obs[, j], j] - pred[obs[, j], j])^2)))
    return(out)
  }

  # fit to the observed responses
  fitted_mean <- fitted(fit)
  out$per_response <- data.frame(
    response = responses, n_missing = colSums(!obs),
    rmse_observed = sapply(seq_along(responses), function(j) sqrt(mean((y[obs[, j], j] - fitted_mean[obs[, j], j])^2))),
    r_squared = sapply(seq_along(responses), function(j) cor(y[obs[, j], j], fitted_mean[obs[, j], j])^2),
    mean_observed = colMeans(y, na.rm = TRUE), sd_observed = apply(y, 2, sd, na.rm = TRUE), row.names = NULL)
  if(!is.null(fit$y_pred)) {
    pred <- simplify2array(one_chain(fit, 1)$y_pred)
    lower <- apply(pred, c(1, 2), quantile, 0.025)
    upper <- apply(pred, c(1, 2), quantile, 0.975)
    out$per_response$coverage_95 <- sapply(seq_along(responses), function(j) mean(y[obs[, j], j] > lower[obs[, j], j] & y[obs[, j], j] < upper[obs[, j], j]))
    if(name %in% c("missBART1", "missBART2")) {
      c1 <- one_chain(fit, 1)
      f1 <- fitted(c1)
      plot_data <- do.call(rbind, lapply(seq_along(responses), function(j) {
        data.frame(response = sprintf("%s (R squared = %.2f)", responses[j], cor(y[obs[, j], j], f1[obs[, j], j])^2),
                   observed = y[obs[, j], j], fitted = f1[obs[, j], j], lower = lower[obs[, j], j], upper = upper[obs[, j], j])
      }))
      plot_data$response <- factor(plot_data$response, levels = unique(plot_data$response))
      rug <- data.frame(response = levels(plot_data$response)[as.integer(miss_col)], imputed = imputed_values(c1)$mean)
      rug$response <- factor(rug$response, levels = levels(plot_data$response))
      # predictions (with their intervals) on the vertical axis, data on the horizontal axis
      fig <- ggplot(plot_data, aes(observed, fitted)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey75", linewidth = 0.2) +
        geom_point(size = 0.4, colour = "#1A9993") +
        geom_abline(linewidth = 0.3) +
        geom_rug(data = rug, aes(y = imputed), inherit.aes = FALSE, sides = "l", colour = "#FD7446", alpha = 0.3) +
        facet_wrap(~response, scales = "free", nrow = 1) +
        labs(x = "Observed (log scale)", y = "Predicted (log scale)", title = name) +
        theme_bw()
      ggsave(file.path("vignettes/articles/figures", paste0("amax_", name, ".png")), fig, width = 12, height = 3, dpi = 110)
    }
  }
  out$sigma <- Reduce("+", fit$sigma_post) / length(fit$sigma_post)

  # by chain: imputations, probit coefficients, variable importance
  chains <- lapply(seq_len(fit$n_chains), function(k) one_chain(fit, k))
  out$imputed <- do.call(rbind, lapply(seq_along(chains), function(k) {
    data.frame(chain = k, response = miss_col, value = imputed_values(chains[[k]])$mean)
  }))
  out$diagnostics <- chain_diagnostics(fit)
  if(!is.null(fit$B_post) && fit$include_y) {
    out$own_coef <- t(sapply(chains, function(f) diag(coef(f)[responses, responses, "mean"])))
  }
  keep <- function(k) fit$chain == k
  out$varimp_regression <- t(sapply(seq_along(chains), function(k) {
    v <- fit$var_imp_reg[keep(k), , drop = FALSE]
    colMeans(v / rowSums(v))
  }))
  if(name == "missBART2") {
    out$varimp_missingness <- t(sapply(seq_along(chains), function(k) {
      v <- fit$var_imp[keep(k), , drop = FALSE]
      colMeans(v / rowSums(v))
    }))
    out$interactions_regression <- interactions(chains[[1]])
    out$interactions_missingness <- interactions(chains[[1]], "missingness")
  }

  # partial dependence: responses against soil pH, and missingness against Parea
  if(!is.null(fit$reg_trees)) out$pdp_pH <- pdp_summary(fit, "pH", "response")
  if(name %in% c("missBART1", "missBART2")) {
    out$pdp_missing_Parea <- do.call(rbind, lapply(seq_along(chains), function(k) {
      cbind(chain = k, pdp_summary(chains[[k]], "Parea", "missingness"))
    }))
  }
  out
}

names_ <- c("missBART1", "missBART1_marginal", "missBART2", "MAR", "complete_cases")
res <- parallel::mclapply(names_, analyse, mc.cores = length(names_))
names(res) <- names_
bad <- sapply(res, inherits, what = "try-error")
if(any(bad)) stop("Failed: ", paste(names_[bad], collapse = ", "), "\n", res[[which(bad)[1]]])

# Distribution of the imputed values (posterior means) against the observed values
imp <- do.call(rbind, lapply(names_[names_ != "missBART1_marginal"], function(nm) cbind(model = nm, res[[nm]]$imputed)))
imp$model <- factor(imp$model, levels = c("complete_cases", "MAR", "missBART1", "missBART2"),
                    labels = c("Complete cases", "MAR", "missBART1", "missBART2"))
observed <- data.frame(response = factor(responses[col(y)[obs]], levels = responses), value = y[obs])
fig <- ggplot(imp, aes(value)) +
  geom_density(data = observed, fill = "grey80", colour = NA) +
  geom_density(aes(group = chain, colour = factor(chain)), linewidth = 0.5) +
  facet_grid(model ~ response, scales = "free") +
  scale_colour_brewer(palette = "Dark2", name = "Chain") +
  labs(x = "Value (log scale)", y = "Density") +
  theme_bw()
ggsave("vignettes/articles/figures/amax_imputed.png", fig, width = 12, height = 6.5, dpi = 110)

saveRDS(c(res, list(settings = list(n = nrow(y), q = ncol(x), burn = burn, iters = iters, n_complete = sum(complete.cases(y)),
                                    version = as.character(packageVersion("missBART"))))),
        "vignettes/articles/results/global_amax.rds")
for(r in res) cat(sprintf("%-20s %5.0f s\n", r$name, r$seconds))
