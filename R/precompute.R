#' computing coefficients for mixture of finite mixture model
#'
#' Function to compute log(V_N(T)) under given MFM parameters for t = 1:upto.
#' Pre-computed prior to posterior sampling.
#'
#' @param log_pk log of prior probability: pr(K=k)
#' @param gamma the hyperparameter for symmetric Dirichlet prior for pi_K
#' @param n number of observations
#' @param upto the upper limit (t<= upto) to compute the coefficient
#'
#' @return a vector of log coefficients in mixture of finite mixtures
#' @export
mfm_coefficients <- function(log_pk,gamma,n,upto){
  tol   <- 1e-12
  log_v <- rep(0,upto)
  for (t in 1:upto){
    if (t>n){
      log_v[t] <- -1e99
    } else {
      a <- 0
      c <- -1e99
      k <- 1
      p <- 0
      while( abs(a-c) > tol || p < 1-tol){
        if (k >= t){
          a <- c
          b <- lgamma(k+1)-lgamma(k-t+1)-lgamma(k*gamma+n)+lgamma(k*gamma)+log_pk(k)
          c <- matrixStats::logSumExp(c(a,b))
        }
        p <- p+exp(log_pk(k))
        k <- k+1
      }
      log_v[t] <- c
    }
  }
  log_v
}
