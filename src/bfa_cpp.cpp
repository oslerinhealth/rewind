// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double logsumexp(arma::vec logv_arma)
{
  //int n = logv.size();
  //if(n<=1)
  //	cout<<"Warning in logsumexp"<<endl;
  double max = logv_arma.max();
  double answer = 0;
  // log(sum(exp(logf)) 	= log(sum(exp(logf - max(logf) + max(logf)))
  //			= max(logf) + log(sum(exp(logf - max(logf)))
  answer = max + log(sum(exp(logv_arma-max)));
  return answer;
}

//' R Function to compute the cluster-specific marginal likelihood
//'
//' This R function computes the marginal likelihood by integrating over
//' the distribution of component specific parameter (e.g., machine usage profiles).
//' This function conditions upon a few model parameters: the true and false positive
//' rates (theta and psi), the Q matrix and {p}-the prevalence parameter for each machines.
//'
//' @param Y the data for the current cluster (a subset of observations.)
//' @param Q Q-matrix
//' @param p prevalence parameter for each machine; should be a vector of dimension M.
//' @param theta true positive rates
//' @param psi true positive rates
//'
//' @examples
//' # simulate data:
//' L0 <- 100
//' options_sim0  <- list(N = 200,  # sample size.
//'                      M = 3,   # true number of machines.
//'                      L = L0,   # number of antibody landmarks.
//'                      K = 8,    # number of true components.,
//'                      theta = rep(0.8,L0), # true positive rates
//'                      psi   = rep(0.01,L0), # false positive rates
//'                      alpha1 = 1 # half of the people have the first machine.
//')
//'
//'  simu     <- simulate_data(options_sim0, SETSEED=TRUE)
//'  simu_dat <- simu$datmat
//'  Y <- simu_dat
//'  Q <- simu$Q
//'  p <- c(0.5,0.25,0.1,0.02,0.05)
//'  theta <- options_sim0$theta
//'  psi   <- options_sim0$psi
//'
//' #log_marginal0(Y, Q, p, theta, psi)
//' log_marginal(Y, Q, p, theta, psi) # <-- this is the Rcpp implementation.
//'
//' @return log of marginal likelihood given other model parameters.
//' @export
// [[Rcpp::export]]
double log_marginal(NumericMatrix Y,NumericMatrix Q,
                    NumericVector p,NumericVector theta,NumericVector psi){
  int L = Y.ncol(), n = Y.nrow(), M = Q.nrow();
  NumericVector n1(L), n0(L), p_xi(L), v(2);
  double res = 0;

  for (int l=0; l<L; l++){
    for (int i=0; i<n; i++){
      n1[l] += Y(i,l);
    }
    n0[l] = n-n1[l];
    for (int m=0; m<M; m++){
      p_xi[l] += log(1-p[m])*Q(m,l);
    }
    p_xi[l] = 1-exp(p_xi[l]);
  }

  for (int l=0; l<L; l++){
    v[0] = n1[l]*log(psi[l])+n0[l]*log(1-psi[l])+log(1-p_xi[l]);
    v[1] = n1[l]*log(theta[l])+n0[l]*log(1-theta[l])+log(p_xi[l]);
    res += logsumexp(v);
  }
  return(res);
}

// [[Rcpp::export]]
arma::mat mat_times_vec_by_col(arma::mat m,arma::vec v){
  return(m * arma::diagmat(v));
}

//' Function to compute the full cluster-specific likelihood given latent variables
//'
//' This function computes the likelihood WITHOUT integrating over
//' the distribution of component specific parameter (e.g., machine usage profiles).
//' This function conditions upon a few model parameters: the true and false positive
//' rates (theta and psi), the Q matrix and {p}-the prevalence parameter for each machines.
//'
//' @param Y the data for the current cluster (a subset of observations.)
//' @param eta_star A matrix of M columns (of machines). Multivariate binary indicators of presence or absence of
//' protein landmarks (could be a matrix with rows for multiple )
//' @param Q Q-matrix
//' @param p prevalence parameter for each machine; should be a vector of dimension M.
//' @param theta true positive rates
//' @param psi true positive rates
//'
//' @return a vector of likelihood values (one per cluster) given other model parameters.
//' @export
// [[Rcpp::export]]
arma::vec log_full(arma::mat Y,arma::mat eta_star,arma::mat Q,
                   arma::vec p,arma::vec theta,arma::vec psi){
  int L = Y.n_cols, n = Y.n_rows, M = Q.n_rows, J = eta_star.n_rows;
  arma::vec n1(L), n0(L), res_prior(J),res_lkd(J);n1.zeros();n0.zeros();res_prior.zeros();res_lkd.zeros();
  arma::vec res(J); res.zeros();
  arma::mat xi(J,L);
  arma::mat xi0; xi0.zeros(J,L);
  xi = xi0+(eta_star*Q>0.5);
  for (int l=0; l<L; l++){
    for (int i=0; i<n; i++){
      n1[l] += Y(i,l);
    }
    n0[l] = n-n1[l];
  }

  for (int j=0; j<J; j++){
    for (int m=0; m<M; m++){
      res_prior[j]+= eta_star(j,m)*log(p[m]) + (1-eta_star(j,m))*log(1-p[m]);
    }
    for (int l=0;l<L;l++){
      res_lkd[j] += (n1[l]*log(psi[l])+n0[l]*log(1-psi[l]))*(1-xi(j,l))+
        (n1[l]*log(theta[l])+n0[l]*log(1-theta[l]))*xi(j,l);
    }
    res[j] = res_prior[j]+res_lkd[j];
  }
  return(res);
  //return(exp(res-logsumexp(res))); // exponentiated.
}

// [[Rcpp::export]]
bool equal_unit_vec(arma::vec v,int k){
  int M = v.n_elem;
  arma::vec ek(M); ek.zeros();
  ek[k-1] = 1; // note the index needs to minus one.!!
  return( sum((v-ek)%(v-ek))<0.001);
}

// [[Rcpp::export]]
double compute_Q_condpr(arma::mat Q, arma::mat H, arma::vec Yl,
                        int k, int l, double theta, double psi){
  int n = Yl.n_elem;
  int L = Q.n_cols;
  arma::vec log_lkd(2);log_lkd.zeros();
  arma::mat xi(n,L);
  arma::mat xi0; xi0.zeros(n,L);
  arma::mat Q_cand = Q;
  //set to 0:
  Q_cand(k-1,l-1) = 0;
  xi = xi0+(H*Q_cand>0.5);
  log_lkd[0] = sum(xi.col(l-1)%Yl)*log(theta)+sum(xi.col(l-1)%(1-Yl))*log(1-theta)+
    sum((1-xi.col(l-1))%Yl)*log(psi)+sum((1-xi.col(l-1))%(1-Yl))*log(1-psi);
  //set to 1:
  Q_cand(k-1,l-1) = 1;
  xi = xi0+(H*Q_cand>0.5);
  log_lkd[1] = sum(xi.col(l-1)%Yl)*log(theta)+sum(xi.col(l-1)%(1-Yl))*log(1-theta)+
    sum((1-xi.col(l-1))%Yl)*log(psi)+sum((1-xi.col(l-1))%(1-Yl))*log(1-psi);

  return(exp(log_lkd[1]-logsumexp(log_lkd)));
}

// do_update_Q <- function(Q){
//   M <- nrow(Q)
//   L <- ncol(Q)
//   test <- array(NA,c(M,L,3))
//   ind_unit_Q <- colSums(Q)==1
//   for (k in 1:M){
//     for (l in 1:L){
//       test[k,l,1] <- equal_unit_vec(Q[,l],k)
//       test[k,l,2] <- (sum(Q[k,])==3) && (Q[k,l]==1)
//       Q_cand      <- Q[,ind_unit_Q,drop=FALSE]
//       test[k,l,3] <- (Q[k,l]==0) && sum(Q[,l])==1 &&
//         (sum(sapply(1:ncol(Q_cand), function(c) equal_unit_vec(Q_cand[,c],which(Q[,l]==1))))==2)
//     }
//   }
//   apply(test,c(1,2),function(v) !any(v))
// }
// if ( abs(sum(Q.col(l)-1)) < 0.001){}

// // [[Rcpp::export]]
// arma::mat do_update_Q(arma::mat Q){
//   int M = Q.n_rows, L=Q.n_cols;
//   arma::cube test(M,L,3); test.zeros();
//   arma::mat check_unit(M,L); check_unit.zeros();
//   arma::mat res(M,L); res.zeros();
//   for (int l1=0; l1<L; l1++){
//     for (int m1=0; m1<M; m1++){
//       check_unit(m1,l1) = equal_unit_vec(Q.col(l1),m1);
//     }
//   }
//
//   for (int l=0; l<L; l++){
//     for (int m=0; m<M; m++){
//       test(m,l,0) = check_unit(m,l);
//       test(m,l,1) = (sum(Q.row(m))==3) && (Q(m,l)==1);
//       test(m,l,2) = (Q(m,l)==0) && (sum(Q.col(l))==1);
//       if (test(m,l,2)){
//         test(m,l,2) = (sum(check_unit.t()*check_unit.col(l)) == 2);
//       }
//       //res(m,l)  = 0+(test(m,l,0)+test(m,l,1)+test(m,l,2) < 0.5);
//       res(m,l)  = 0+(test(m,l,0)+test(m,l,1)+test(m,l,2) < 0.5);
//     }
//   }
//   return(res);
// }
//
//
// update_Q <- function(Y,H,Q_old,theta,psi){
//   M  <- nrow(Q_old)
//   N  <- nrow(Y)
//   L  <- ncol(Y)
//   for (k in 1:M){
//     for (l in 1:L){ # begin iteration over elements.
//       if(do_update_Q(Q_old)[k,l]){ # begin an update if needed.
//       Q_cand <- Q_old
//         Q_cand[k,l] <- 1 # set to 1:
//         xi <- (H%*%Q_cand>0.5)+1
//       n11 <- sum(xi[,l]*Y[,l])
//         n10 <- sum(xi[,l]*(1-Y[,l]))
//         n01 <- sum((1-xi[,l])*Y[,l])
//         n00 <- sum((1-xi[,l])*(1-Y[,l]))
//         log_lkd1 <- n11*log(theta[l])+n10*log(1-theta[l])+n01*log(psi[l])+n00*log(1-psi[l])
//
//         Q_cand[k,l] <- 0 # set to 0:
//         xi <- (H%*%Q_cand>0.5)+1
//       n11 <- sum(xi[,l]*Y[,l])
//         n10 <- sum(xi[,l]*(1-Y[,l]))
//         n01 <- sum((1-xi[,l])*Y[,l])
//         n00 <- sum((1-xi[,l])*(1-Y[,l]))
//         log_lkd0 <- n11*log(theta[l])+n10*log(1-theta[l])+n01*log(psi[l])+n00*log(1-psi[l])
//
//
//         curr_prob       <- exp(log_lkd1-matrixStats::logSumExp(c(log_lkd1,log_lkd0)))
//         if (Q_old[k,l]==1){
//           if (curr_prob <= 0.5) {
//             alpha <- 1
//           } else{
//             alpha <- exp(log(1-curr_prob) - log(curr_prob))
//           }
//           Q_old[k,l] <- rbinom(1,1,prob=1-alpha)
//         } else{
//           if (curr_prob >= 0.5){
//             alpha <- 1
//           } else{
//             alpha <- exp(log(curr_prob) - log(1-curr_prob))
//           }
//           Q_old[k,l] <- rbinom(1,1,prob=alpha)
//         } #end flipping.
//       } # end an update if needed.
//     }
//   }  #end iteration over elements.
//   order_mat_byrow(Q_old)
// }
//
// // // [[Rcpp::export]]
// arma::mat update_Q(arma::mat Y, arma::mat H, arma::mat Q_old, arma::vec, arma::vec psi){
//   int M = Q_old.n_rows, N = Y.n_rows, L = Y.n_cols;
//   for (int m = 0; m<M; m++){
//     for (int l=0; l<L; l++){// begin iteration over elements
//
//     }
//   }
// }
//
//
//
//
//

