# Re-runs the simulation studies of Goh, Soh, Parnell and Murphy (2024) with the
# C++ samplers, on the simulated data sets and cross-validation folds stored in
# this folder:
#   MAR_1, MAR_2, MNAR_1, MNAR_2          Section 6.1 (bivariate, 4 folds)
#   MNAR_amp0, MNAR_amp1, MNAR_amp2       Section 6.2 (five responses, 8 folds)
#   nshape, ushape                        Appendix D (univariate, 4 folds)
#
# Every (scenario, fold, model) combination is one task; tasks run in parallel.
# The out-of-sample predictions are saved to
# vignettes/articles/results/paper_simulations.rds, which the "Simulation
# checking" article reads.
#
# Usage, from the root of the repository:  Rscript Examples/rerun_paper_simulations.R [cores] [burn] [iters]

library(missBART)
library(parallel)

cores <- as.integer(commandArgs(TRUE)[1])
if(is.na(cores)) cores <- max(1, detectCores() - 2)
# optional second and third arguments shorten the chains, for a quick test
burn <- as.integer(commandArgs(TRUE)[2])
if(is.na(burn)) burn <- 5000
iters <- as.integer(commandArgs(TRUE)[3])
if(is.na(iters)) iters <- 5000
n_reg_trees <- 100

scenarios <- list(MAR_1 = 20, MAR_2 = 20, MNAR_1 = 20, MNAR_2 = 20,
                  MNAR_amp0 = 50, MNAR_amp1 = 50, MNAR_amp2 = 50, nshape = 20, ushape = 20)

load_scenario <- function(name) {
  e <- new.env()
  load(file.path("Examples", paste0(name, ".RData")), e)
  x_imp <- if(!is.null(e$imputed_x_mf) && anyNA(e$x)) as.matrix(e$imputed_x_mf) else e$x
  list(x = as.matrix(e$x), x_imp = x_imp, y = as.matrix(e$y), y_true = as.matrix(e$y_original),
       y_imp = as.matrix(e$imputed_y_mf), folds = e$batch_id)
}

models_for <- function(name, d) {
  p <- ncol(d$y)
  joint <- c("missBART1", "missBART1_marginal", "missBART2", if(anyNA(d$x)) "missBART2_impX")
  if(p == 1) return(c("missBART1", "missBART2", "uniBART_cc", "uniBART_imp"))
  c(joint, "mvBART_cc", "mvBART_imp", "uni_missBART1", "uni_missBART2", "uniBART_cc", "uniBART_imp")
}

# Fits one model to the training folds and predicts the test fold.
run_task <- function(task) {
  d <- load_scenario(task$scenario)
  test <- d$folds[[task$fold]]
  n_class <- scenarios[[task$scenario]]
  set.seed(1000 * task$fold + match(task$model, c("missBART1", "missBART1_marginal", "missBART2", "missBART2_impX",
                                                    "mvBART_cc", "mvBART_imp", "uni_missBART1", "uni_missBART2",
                                                    "uniBART_cc", "uniBART_imp")))
  common <- list(burn = burn, iters = iters, thin = 1, show_progress = FALSE, save_y_pred = FALSE)
  fit_one <- function(type, x, y, x_new) {
    args <- c(list(x = x, y = y, x_predict = x_new), common)
    switch(type,
           missBART1 = do.call(missBART1, c(args, list(n_trees = n_reg_trees, save_z = FALSE))),
           missBART1_marginal = do.call(missBART1, c(args, list(n_trees = n_reg_trees, save_z = FALSE, ymis_update = "marginal"))),
           missBART2 = do.call(missBART2, c(args, list(n_reg_trees = n_reg_trees, n_class_trees = n_class, save_z = FALSE))),
           cc = {
             ok <- stats::complete.cases(y)
             do.call(mvBART, c(list(x = x[ok, , drop = FALSE], y = y[ok, , drop = FALSE], x_predict = x_new, n_trees = n_reg_trees), common))
           },
           imp = do.call(mvBART, c(args, list(n_trees = n_reg_trees))))
  }
  x_na <- d$x[-test, , drop = FALSE]
  x_na_new <- d$x[test, , drop = FALSE]
  x_imp <- d$x_imp[-test, , drop = FALSE]
  x_imp_new <- d$x_imp[test, , drop = FALSE]
  y <- d$y[-test, , drop = FALSE]
  y_imp <- d$y_imp[-test, , drop = FALSE]
  p <- ncol(y)
  extra <- NULL
  time <- system.time({
    pred <- switch(task$model,
      missBART1 = , missBART1_marginal = {
        fit <- fit_one(task$model, x_imp, y, x_imp_new)
        if(task$fold == 1) extra <- list(coef = coef(fit))
        predict(fit)
      },
      missBART2 = {
        fit <- fit_one("missBART2", x_na, y, x_na_new)
        if(task$fold == 1) extra <- list(varimp = varimp(fit), y_miss_accept = mean(fit$y_miss_accept))
        predict(fit)
      },
      missBART2_impX = predict(fit_one("missBART2", x_imp, y, x_imp_new)),
      mvBART_cc = predict(fit_one("cc", x_na, y, x_na_new)),
      mvBART_imp = predict(fit_one("imp", x_na, y_imp, x_na_new)),
      uni_missBART1 = sapply(seq_len(p), function(j) predict(fit_one("missBART1", x_imp, y[, j, drop = FALSE], x_imp_new))),
      uni_missBART2 = sapply(seq_len(p), function(j) predict(fit_one("missBART2", x_na, y[, j, drop = FALSE], x_na_new))),
      uniBART_cc = sapply(seq_len(p), function(j) predict(fit_one("cc", x_na, y[, j, drop = FALSE], x_na_new))),
      uniBART_imp = sapply(seq_len(p), function(j) predict(fit_one("imp", x_na, y_imp[, j, drop = FALSE], x_na_new))))
  })[["elapsed"]]
  c(task, list(pred = unname(as.matrix(pred)), truth = unname(d$y_true[test, , drop = FALSE]),
               observed = unname(!is.na(d$y[test, , drop = FALSE])), seconds = time, extra = extra))
}

tasks <- list()
for(name in names(scenarios)) {
  d <- load_scenario(name)
  for(fold in seq_along(d$folds)) for(model in models_for(name, d)) {
    tasks[[length(tasks) + 1]] <- list(scenario = name, fold = fold, model = model)
  }
}
cat(length(tasks), "tasks on", cores, "cores\n")
# longest tasks first
ord <- order(-sapply(tasks, function(t) (grepl("amp", t$scenario) * 4 + 1) * (1 + 4 * grepl("^uni", t$model))))
start <- Sys.time()
results <- mclapply(tasks[ord], function(t) tryCatch(run_task(t), error = function(e) c(t, list(error = conditionMessage(e)))),
                    mc.cores = cores, mc.preschedule = FALSE)
minutes <- as.numeric(difftime(Sys.time(), start, units = "mins"))
cat("Finished in", round(minutes, 1), "minutes\n")
failed <- Filter(function(r) !is.null(r$error) || is.null(r$scenario), results)
if(length(failed)) { cat(length(failed), "tasks failed:\n"); print(unique(sapply(failed, function(r) paste(r$scenario, r$model, r$error)))) }
dir.create("vignettes/articles/results", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(results = results, settings = list(burn = burn, iters = iters, n_reg_trees = n_reg_trees, n_class_trees = scenarios, cores = cores, minutes = minutes,
                                                version = as.character(packageVersion("missBART")))),
        "vignettes/articles/results/paper_simulations.rds")
