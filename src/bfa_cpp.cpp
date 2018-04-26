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


//' R Function to compute the cluster-specific marginal likelihood
//'
//' This R function computes the marginal likelihood by integrating over
//' the distribution of component specific parameter (e.g., machine usage profiles).
//' This function conditions upon a few model parameters: the true and false positive
//' rates (theta and psi), the Q matrix and {p}-the prevalence parameter for each machines.
//'
//' @param Y the data for the current cluster (a subset of observations.)
//' @param eta_star_enumerate fixed binary matrix of 2^M rows and M columns. Need to be prespecified.
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
//'  H_enumerate <- as.matrix(expand.grid(rep(list(0:1), options_sim0$M)),ncol=options_sim0$M)
//'
//' #log_marginal0(Y, Q, p, theta, psi)
//' log_marginal(Y,H_enumerate, Q, p, theta, psi) # <-- this is the Rcpp implementation.
//'
//' @return log of marginal likelihood given other model parameters.
//' @export
// [[Rcpp::export]]
double log_marginal(arma::mat Y, arma::mat eta_star_enumerate,
                    arma::mat Q,
                    arma::vec p,arma::vec theta,arma::vec psi){
  int J = eta_star_enumerate.n_rows;
  arma::vec res_enumerate(J); res_enumerate.zeros();
  res_enumerate = log_full(Y,eta_star_enumerate,Q,p,theta,psi);
  return(logsumexp(res_enumerate));
}

//' R Function to compute the cluster-specific marginal likelihood (just for Q=I)
//'
//' This R function computes the marginal likelihood by integrating over
//' the distribution of component specific parameter (e.g., machine usage profiles).
//' This function conditions upon a few model parameters: the true and false positive
//' rates (theta and psi), the Q matrix and {p}-the prevalence parameter for each machines.
//'
//' @param Y the data for the current cluster (a subset of observations.)
//' @param p prevalence parameter for each machine; should be a vector of dimension
//' M=L.
//' @param theta true positive rates
//' @param psi true positive rates
//'
//' @examples
//' # simulate data:
//' L0 <- 100
//' options_sim0  <- list(N = 200,  # sample size.
//'                       M = 3,   # true number of machines.
//'                       L = L0,   # number of antibody landmarks.
//'                       K = 8,    # number of true components.,
//'                      theta = rep(0.8,L0), # true positive rates
//'                      psi   = rep(0.01,L0), # false positive rates
//'                      alpha1 = 1 # half of the people have the first machine.
//')
//'
//'  simu     <- simulate_data(options_sim0, SETSEED=TRUE)
//'  simu_dat <- simu$datmat
//'  Y <- simu_dat
//'  Q <- simu$Q
//'  p <- rep(0.5,L0) #<----- M must equal L.
//'  theta <- options_sim0$theta
//'  psi   <- options_sim0$psi
//'
//' log_marginal_Q_identity(Y, p, theta, psi) # <-- this is the Rcpp implementation.
//'
//' @return log of marginal likelihood given other model parameters.
//' @export
// [[Rcpp::export]]
double log_marginal_Q_identity(arma::mat Y,
                               arma::vec p,
                               arma::vec theta,
                               arma::vec psi){
  int L = Y.n_cols, n = Y.n_rows;
  arma::vec n1(L), n0(L), v(2); n1.zeros(); n0.zeros();v.zeros();
  double res = 0;
  for (int l=0; l<L; l++){
    for (int i=0; i<n; i++){
      n1[l] += Y(i,l);
    }
    n0[l] = n-n1[l];
  }
  // Rcout << "The n1 value is " << n1 << std::endl;
  // Rcout << "The n0 value is " << n0 << std::endl;
  // Rcout << "The p value is " << p << std::endl;
  for (int l=0; l<L; l++){
    v[0] = n1[l]*log(psi[l])+n0[l]*log(1-psi[l])+log(1-p[l]);
    v[1] = n1[l]*log(theta[l])+n0[l]*log(1-theta[l])+log(p[l]);
    res += logsumexp(v);
    // Rcout << "The v value is " << v << std::endl;
  }
  return(res);
}

//' check whether a vector is equal to a unit vector with the one at a particular
//' position
//'
//' This function is used in updating Q matrix if we constrain the updates within
//' the identifiability assumption
//'
//' @param v the vector (a binary vector)
//' @param k the index that is being checked if \code{v[k]} is the only one in
//' vector \code{v}. \code{k} must be smaller than or equal to the length of k
//' @return TRUE if \eqn{v = \mathbf{e}_k}; FALSE otherwise.
//' @examples
//'  equal_unit_vec(c(1,0,0,0,0,0),1)
//'  equal_unit_vec(c(1,0,0,0,0,0),2)
//'  equal_unit_vec(c(0,0,2,0,0,0),3)
//'  equal_unit_vec(c(0,0,1,0,0,0),3)
//' @export
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

// // // [[Rcpp::export]]
// arma::mat update_Q(arma::mat Y, arma::mat H, arma::mat Q_old, arma::vec, arma::vec psi){
//   int M = Q_old.n_rows, N = Y.n_rows, L = Y.n_cols;
//   for (int m = 0; m<M; m++){
//     for (int l=0; l<L; l++){// begin iteration over elements
//
//     }
//   }
// }


