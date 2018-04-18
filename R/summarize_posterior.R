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

#' Frobenius norm: ||AA'-BB'||_F
#'
#'
#' @param A,B two matrixes of the same number of rows
#'
#' @return a positive number
#' @export
#'
distAB <- function(A,B){
  norm(A%*%t(A) - B%*%t(B),"F")
}

#' Convert the cluster indicators to the co-activation patterns.
#'
#' @param z a vector of integers indicating cluster assignments for all subjects
#'
#' @return an N by N binary; an 1 represents that two subjects are in the same
#' cluster; 0 otherwise.
#' @export
#'
#' @examples
#' z2comat(c(1,1,2,2,3,4,5,6,5,7))
z2comat <- function(z){ # for Dahl (2006): simple least square approach to estimate clusters.
  n <- length(z)
  res <- matrix(NA,nrow=n,ncol=n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      res[i,j] <- res[j,i] <- (z[i]==z[j])+0
    }
  }
  diag(res) <- 1
  res
}
