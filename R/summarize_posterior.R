#' Compute the posterior similarity matrix (probability that i and j are in same cluster).
#'
#' This function is to evaluate the recovery of clusters
#' @param n total number of subjejcts
#' @param z a matrix of posterior samples, rows for subejcts, columns for MCMC samples kept across
#' multiple iterations
#' @param n_keep the number of columns in z
#'
#' @return a matrix of co-clustering frequencies
#' @export
coclust_mat <- function(n,z,n_keep){
  res     <- matrix(0,n,n)
  count <- 0
  for (index in 1:n_keep){
    for (i in 1:n){
      for (j in 1:n){
        if (z[i,index]==z[j,index]){
          res[i,j] <- res[i,j]+1
        }
      }
    }
    count <- count+1
  }
  res/count
}
