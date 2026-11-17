#' missBART: joint models for non-ignorable missing responses using BART
#'
#' The main functions are [missBART1()], [missBART2()] and [mvBART()]. See
#' Goh, Soh, Parnell and Murphy (2024) <doi:10.48550/arXiv.2412.14946>.
#'
#' @useDynLib missBART, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats fitted predict coef
#' @importFrom graphics plot
#' @keywords internal
"_PACKAGE"
