# Helpers kept in R: they are only used to simulate data and are not part of the samplers.


rMVNmu0 <- function(Q) {
  p  <- NCOL(Q)
  z  <- stats::rnorm(p)
  if(p == 1) z/sqrt(Q) else backsolve(PD_chol(Q), z, k=p)
}

rbernoulli <- function(n, p = 0.5) {
  stats::runif(n) > (1 - p)
}

PD_chol <- function(x, ...) tryCatch(chol(x, ...), error = function(e) {
  chol(x + diag(1e-8 * max(abs(diag(x))), nrow(x)), ...)
})

# Get vector for sorting observations into terminal nodes
get_change_points <- function(df, x) {
  n <- nrow(x)
  if(nrow(df) == 1) {
    change_points <- rep(1, n)
  } else {
    change_points <- c()
    terminal_nodes <- as.numeric(setdiff(row.names(df), df$parent))
    terminal_obs <- list()
    for(j in seq_along(terminal_nodes)) {
      y_index <- list()
      parent <- terminal_nodes[j]
      count <- 1
      while(parent != 1) {
        #Add in if/else for type of variable (numeric, integer, categorical)
        #Add in split if missing (not if m == 0/1)
        type <- class(x[,df$split_variable[parent]])
        # switch(EXPR = type,
        #        numeric = {
        #          y_index[[count]] <- if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] <= df$split_value[parent]) else which(x[,df$split_variable[parent]] > df$split_value[parent])
        #          if(any(is.na(x[,df$split_variable[parent]]))){
        #            if(df$NA_direction[parent] == 1) y_index[[count]] <- sort(c(y_index[[count]], which(is.na(x[,df$split_variable[parent]]))))
        #          }
        #        }, factor = {
        #          y_index[[count]] <- if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] == df$split_value[parent]) else which(x[,df$split_variable[parent]] != df$split_value[parent])
        #        }) #End switch
        y_index[[count]] <- if(df$direction[parent] == 0) which(x[,df$split_variable[parent]] <= df$split_value[parent]) else which(x[,df$split_variable[parent]] > df$split_value[parent])
        if(any(is.na(x[,df$split_variable[parent]]))){
          if(df$NA_direction[parent] == 1) y_index[[count]] <- sort(c(y_index[[count]], which(is.na(x[,df$split_variable[parent]]))))
        }
        parent <- df$parent[parent]
        count <- count + 1
      }
      terminal_obs[[j]] <- Reduce(intersect, y_index)
      if(length(terminal_obs[[j]]) == 0) terminal_obs[[j]] <- NA
    }
    terminal_obs <- terminal_obs[!is.na(terminal_obs)]
    change_points <- vector(length=n)

    for(k in seq_along(terminal_obs)) {
      change_points[terminal_obs[[k]]] <- k
    }
    # max_node <- max(unique(change_points))
    # change_points[which(change_points == 0)] <- sample(c(max_node, max_node - 1), 1)
  }
  return(change_points)
}
