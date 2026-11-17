# Timing and accuracy of the samplers on simulated MNAR data.
#
# The same script runs against version 0.1.0 (pure R) and against the current
# C++ version, so the two can be compared on identical data:
#
#   # version 0.1.0 installed into a separate library
#   remotes::install_github("yongchengoh/missBART@802fb80", lib = "~/R/benchlib")
#   Rscript benchmark.R old 2000 5 25 25 old.rds
#   Rscript benchmark.R new 2000 5 2500 2500 new.rds
#
# Arguments: old|new, n, p (at most 5), burn-in, kept iterations, output file.
# Reported: seconds, RMSE of the posterior mean imputations against the deleted
# values, and RMSE of the fitted mean function against the truth.
a <- commandArgs(TRUE); which <- a[1]; n <- as.integer(a[2]); p <- as.integer(a[3]); burn <- as.integer(a[4]); iters <- as.integer(a[5])
if (which == "old") suppressMessages(library(missBART, lib.loc = "~/R/benchlib")) else suppressMessages(library(missBART))
set.seed(42)
x <- matrix(runif(n*5), ncol=5)
f <- sapply(seq_len(p), function(j) c(10,5,8,6,9)[j]*sin(pi*x[,1]*x[,2]) + c(20,-10,5,15,0)[j]*(x[,3]-.5)^2 + c(10,3,-6,2,7)[j]*x[,4] + c(5,-4,2,0,3)[j]*x[,5])
S <- 0.5^abs(outer(seq_len(p), seq_len(p), "-"))
ytrue <- f + matrix(rnorm(n*p), ncol=p) %*% chol(S)
y <- ytrue
for (j in seq_len(p)) y[runif(n) < pnorm(-0.8 + 0.6*scale(ytrue[,j])), j] <- NA
rm <- function(a,b) sqrt(mean((a-b)^2))
res <- list()
run <- function(name, expr) {
  set.seed(1); t <- system.time(fit <- eval.parent(substitute(expr)))[3]
  mi <- which(is.na(y)); imp <- Reduce("+", fit$y_impute)/length(fit$y_impute); fh <- Reduce("+", fit$y_post)/length(fit$y_post)
  r <- c(time = unname(t), imp_rmse = rm(imp, ytrue[mi]), fit_rmse = rm(fh, f))
  cat(which, name, "n", n, "p", p, "iters", burn+iters, ":", round(r,3), "\n"); r
}
sink(tempfile()); 
r1 <- run("missBART1", missBARTprobit(x, y, n_trees=100, burn=burn, iters=iters, thin=1, show_progress=FALSE, mice_impute=FALSE))
r2 <- run("missBART2", missBART2(x, y, x_predict=NA, predict=FALSE, n_reg_trees=100, n_class_trees=20, burn=burn, iters=iters, thin=1, show_progress=FALSE, mice_impute=FALSE))
r0 <- run("mvBART", {ok <- complete.cases(y); y0 <- y; fit <- mvBART(x, ytrue, n_trees=100, burn=burn, iters=iters, thin=1, show_progress=FALSE, predict=FALSE); fit$y_impute <- list(ytrue[is.na(y)]); fit})
sink()
out <- rbind(missBART1=r1, missBART2=r2, mvBART=r0)
print(round(out,3))
saveRDS(list(which=which,n=n,p=p,total=burn+iters,res=out), a[6])
