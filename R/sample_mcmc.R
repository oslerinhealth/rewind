# initializations:

#' Randomly generate a Q matrix within the identifying conditions
#'
#' This function is used to initialize the MCMC sampling chain within identifiability
#' cosntraints
#'
#' @param M the number of factors (machines)
#' @param L the dimension of the binary data
#' @param p the Bernoulli probability of an one in the Q matrix
#' @examples
#'
#' simulate_Q(3,100)
#'
#' @return a binary matrix of dimension M by L
#' @export
simulate_Q <- function(M,L,p=0.1){
  Q_col_ordered <-cbind(diag(1,M),diag(1,M),matrix(rbinom(M*(L-2*M),1,p),nrow=M))
  Q_col_ordered[which(rowSums(Q_col_ordered[,-(1:(2*M))]) == 0), sample((2*M+1):(L),1)] <- 1
  Q_col_ordered[,sample(1:L)]
}

#' Randomly generate a Q matrix within the identifying conditions driven by data
#'
#' This function is used to initialize the MCMC sampling chain within identifiability
#' cosntraints; but this is smart - because it will not assign a one to a dimension with
#' few ones in the data
#'
#' @param M the number of factors (machines)
#' @param dat a matrix of binary data (rows for observations, columns for dimensions)
#' @param p the Bernoulli probability of an one in the Q matrix
#' @param frac only initialize those dimension with at least \code{frac}*100%
#' of subjects having a positive response.
#' @examples
#'
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' simu_dat <- simu$datmat
#' simulate_Q_dat(5,simu_dat)
#'
#' @return a binary matrix of dimension M by L
#' @export
simulate_Q_dat <- function(M,dat,p=0.1,frac=1/5){
  #dat <- simu_dat
  ind_all <- which(colSums(dat)>nrow(dat)*frac) # this initial value is very important!. The extra challenge is to known which column is truly all zeros.
  res <- matrix(0,nrow=M,ncol=ncol(dat))
  for (m in 1:M){
    res[m,sample(ind_all,ceiling(p*length(ind_all)))] <- 1
  }
  res
}

# compute likelihood


#' Order a binary matrix by row
#'
#' The order is determined by magnitude of the rows as in a binary system
#'
#' @param mat A binary matrix to be ordered by row
#' @return a list of res - the ordered matrix from large to small, ord - the order
#' @export
order_mat_byrow <- function(mat){
  L <- ncol(mat)
  M <- nrow(mat)
  v <- log(as.matrix(2^{(L-1):0},ncol=1))
  diag_v <- matrix(0,nrow=L,ncol=L)
  diag(diag_v) <- v
  tmp_prod <- mat%*%diag_v
  permute_M_vec <- sapply(1:M,function(i){logSumExp(tmp_prod[i,mat[i,]!=0])})
  list(res = mat[rev(order(permute_M_vec)),,drop=FALSE],
       ord = rev(order(permute_M_vec)))
}

# bin2dec_vec <- function(mat){# This is the log version:
#   #mat <- Q_col_ordered[sample(1:M),]
#   L <- ncol(mat)
#   M <- nrow(mat)
#   v <- log(as.matrix(2^{(L-1):0},ncol=1))
#   diag_v <- matrix(0,nrow=L,ncol=L)
#   diag(diag_v) <- v
#   tmp_prod <- mat%*%diag_v
#   permute_M_vec <- sapply(1:M,function(i){matrixStats::logSumExp(tmp_prod[i,mat[i,]!=0])})
#   permute_M_vec[rev(order(permute_M_vec))]
# }

#
# END of helper functions.
#------------------------------------------------------------------------------#

#' R Function to compute the cluster-specific marginal likelihood
#'
#' This R function computes the marginal likelihood by integrating over
#' the distribution of component specific parameter (e.g., machine usage profiles).
#' This function conditions upon a few model parameters: the true and false positive
#' rates (theta and psi), the Q matrix and {p}-the prevalence parameter for each machines.
#'
#' @param Y the data for the current cluster (a subset of observations.)
#' @param Q Q-matrix
#' @param p prevalence parameter for each machine; should be a vector of dimension M.
#' @param theta true positive rates
#' @param psi true positive rates
#'
#' @examples
#' # simulate data:
#' L0 <- 100
#' options_sim0  <- list(N = 200,  # sample size.
#'                      M = 3,   # true number of machines.
#'                      L = L0,   # number of antibody landmarks.
#'                      K = 8,    # number of true components.,
#'                      theta = rep(0.8,L0), # true positive rates
#'                      psi   = rep(0.01,L0), # false positive rates
#'                      alpha1 = 1 # half of the people have the first machine.
#')
#'
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' simu_dat <- simu$datmat
#' Y <- simu_dat
#'  Q <- simu$Q
#'  p <- c(0.5,0.25,0.1,0.02,0.05)
#'  theta <- options_sim0$theta
#'  psi   <- options_sim0$psi
#'
#' log_marginal0(Y, Q, p, theta, psi)
#' # log_marginal(Y, Q, p, theta, psi) # <-- this is the Rcpp implementation.
#' @return log of marginal likelihood given other model parameters.
#'
# log_marginal0 <- function(Y, Q, p, theta, psi){
#   n1 <- apply(Y,2,sum,na.rm=T)
#   n0 <- apply(1-Y,2,sum,na.rm=T)
#   p_xi <- 1-exp(matrix(log(1-p),nrow=1)%*%Q) # length L.
#
#   mat <- rbind(n1*log(psi)+n0*log(1-psi),
#                n1*log(theta)+n0*log(1-theta))+rbind(log(1-p_xi),log(p_xi))
#   sum(matrixStats::colLogSumExps(mat)) # can remove stuff including and before ::.
# }


#' R Function to compute the cluster-specific marginal likelihood
#'
#' This R function computes the marginal likelihood by integrating over
#' the distribution of component specific parameter (e.g., machine usage profiles).
#' This function conditions upon a few model parameters: the true and false positive
#' rates (theta and psi), the Q matrix and {p}-the prevalence parameter for each machines.
#'
#' @param Y the data for the current cluster (a subset of observations.)
#' @param Q Q-matrix
#' @param p prevalence parameter for each machine; should be a vector of dimension M.
#' @param theta true positive rates
#' @param psi true positive rates
#'
#'
#' # simulate data:
#' L0 <- 100
#' options_sim0  <- list(N = 200,  # sample size.
#'                      M = 3,   # true number of machines.
#'                      L = L0,   # number of antibody landmarks.
#'                      K = 8,    # number of true components.,
#'                      theta = rep(0.8,L0), # true positive rates
#'                      psi   = rep(0.01,L0), # false positive rates
#'                      alpha1 = 1 # half of the people have the first machine.
#')
#'
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' simu_dat <- simu$datmat
#' Y <- simu_dat
#  Q <- simu$Q
#  p <- c(0.5,0.25,0.1,0.02,0.05)
#  theta <- options_sim0$theta
#  psi   <- options_sim0$psi
#
#' l <- 20
#' log_marginal_one_column(Y[,l], Q[,l], p, theta[l], psi[l])
#' @return log of marginal likelihood given other model parameters.
#' @export
log_marginal_one_column <- function(Yl, Ql, p, thetal, psil){
  n1 <- sum(Yl)
  n0 <- sum(1-Yl)
  p_xil <- 1-exp(matrix(log(1-p),nrow=1)%*%matrix(Ql,ncol=1)) # length L.

  mat <- c(n1*log(psil)+n0*log(1-psil),
           n1*log(thetal)+n0*log(1-thetal))+
    c(log(1-p_xil),
      log(p_xil))
  logSumExp(mat) # can remove stuff including and before ::.
}

# mat_times_vec_by_col0 <- function(m,v){m * rep(v, rep(nrow(m),length(v)))}

#' Function to compute the full cluster-specific likelihood given latent variables
#'
#' This function computes the likelihood WITHOUT integrating over
#' the distribution of component specific parameter (e.g., machine usage profiles).
#' This function conditions upon a few model parameters: the true and false positive
#' rates (theta and psi), the Q matrix and {p}-the prevalence parameter for each machines.
#'
#' @param Y the data for the current cluster (a subset of observations.)
#' @param eta_star A matrix of M columns (of machines). Multivariate binary indicators of presence or absence of
#' protein landmarks (could be a matrix with rows for multiple )
#' @param Q Q-matrix
#' @param p prevalence parameter for each machine; should be a vector of dimension M.
#' @param theta true positive rates
#' @param psi true positive rates
#'
#' @importFrom matrixStats colLogSumExps
#' @return log of marginal likelihood given other model parameters.
log_full0 <- function(Y,eta_star,Q, p, theta, psi){# eta_star must be a matrix of M columns.
  n1 <- apply(Y,2,sum,na.rm=T)
  n0 <- apply(1-Y,2,sum,na.rm=T)
  xi <- matrix((eta_star%*%Q > 0.5)+0,ncol=ncol(Q))
  prior <- rowSums(mat_times_vec_by_col(eta_star,log(p))+
                     mat_times_vec_by_col(1-eta_star,log(1-p)))
  PR_mat <- rbind(n1*log(psi)+n0*log(1-psi),
                  n1*log(theta)+n0*log(1-theta))
  likelihood_mat <-mat_times_vec_by_col(1-xi,PR_mat[1,])+
    mat_times_vec_by_col(xi,PR_mat[2,])
  prior+rowSums(likelihood_mat)
}



# eta_star <- as.matrix(expand.grid(rep(list(0:1), options_sim0$M)),ncol=options_sim0$M)
# sum(abs(log_full0(Y,eta_star,Q, p, theta, psi)-
# log_full(Y,eta_star,Q, p, theta, psi)))

#' update true and false positive rates WITHOUT constraints
#'
#' This function samples the theta and psi parameters in bmf
#'
#' @param Y data
#' @param H a matrix whose rows are individual specific machine indicators (all subjects)
#' @param Q Q matrix
#' @param a_theta,a_psi hyperparameters for priors on theta and psi, respectively
#'
#' @return a list of true positive rates theta and false positive rates psi, each of length L.
#' @export
update_positive_rate <- function(Y,H,Q,a_theta,a_psi){
  xi       <- (H%*%Q>0.5)+0
  psi_a1   <- colSums((1-xi)*Y)+a_psi[1]
  psi_a2   <- colSums((1-xi)*(1-Y))+a_psi[2]
  theta_a1 <- colSums(xi*Y)+a_theta[1]
  theta_a2 <- colSums(xi*(1-Y))+a_theta[2]
  theta_samp <- sapply(1:ncol(Q),function(i) {rbeta(1,theta_a1[i],theta_a2[i])})
  psi_samp   <- sapply(1:ncol(Q),function(i) {rbeta(1,psi_a1[i],psi_a2[i])})
  #list(theta=theta_samp,psi=rep(0.1,ncol(Q)))
  list(theta=theta_samp,psi=psi_samp)
}

#' sample alpha - hyperparameter for prevalences p_m
#'
#'  Function to sample the parameters from a grid
#'
#'  @param H_star the matrix of machine usage profiles across clusters
#'  @param t number of clusters
#'  @param M number of machines
#'  @param a,b hyperparameter for beta distribution over reparameterized alpha.
#'  @param show_density Default to FALSE - not showing the full conditional density of alpha
#'  given other unknown parameters; TRUE otherwise.
update_alpha<-function(H_star,t,M,a=1,b=1,show_density=FALSE) {
  th    <- seq(0.001,.999,length=5000)
  sm    <- apply(H_star,2,sum,na.rm=T)
  part1 <- 0
  for (m in 1:M){
    part1 <- part1+lgamma(th/(1-th)/M+sm[m])
  }
  tmp    <- dbeta(th,a,b,log=T)+M*log(th/(1-th))+part1-M*lgamma(th/(1-th)/M+t+1)
  tmp    <- exp(tmp-max(tmp))
  rtheta <- sample(th,1,prob=tmp)
  if (show_density){
    plot(th,tmp,type="l")
    abline(v=rtheta,col="red",lty=2)
  }
  (rtheta)/(1-rtheta)
}

#log_full(simu_dat,as.matrix(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))),
#         Q,p,theta,psi)

#' update prevalence
#'
#' This function updates the prevelance paramters
#'
#' @param H_star machine usage profiles for all clusters
#' @param alpha hyperparameter
#' @param M total number of machines (number of columns for H_star)
#'
#' @return vector of prevalences of length M.
update_prevalence <- function(H_star,alpha,M){
  n1_star <- apply(H_star,2,sum,na.rm=T)
  n0_star <- apply(1-H_star,2,sum,na.rm=T)
  sapply(1:ncol(H_star),function(i) {rbeta(1,n1_star[i]+alpha/M, n0_star[i]+1)})
}

#' check whether a vector is equal to a unit vector with the one at a particular
#' position
#'
#' @param v the vector (a binary vector)
#' @param k the index that is being checked if \code{v[k]} is the only one in
#' vector \code{v}. \code{k} must be smaller than or equal to the length of k
#' @return true for \eqn{v = \mathbf{e}_k}
#' @examples
# equal_unit_vec0(c(1,0,0,0,0,0),1)
# equal_unit_vec0 <- function(v,k){
#   if(k>length(v)){stop("==[rewind]length of vector v is shorter than k!==")}
#   e <- rep(0,length(v))
#   e[k] <- 1
#   (sum(abs(e-v)) == 0)
# }

#' determine to update the current element of Q_ml or not
#'
#'
#' Function to test whether we need to update the current element of Q_ml. This
#' is needed in the constrained Gibbs sampler.
#'
#' @param Q a matrix with rows being machines and columns being protein landmarks (dimension)
#'
#' @return a matrix filled with logical values of dimensions identical to Q. TRUE for updating
#' in constrained Gibbs sampler, FALSE for skipping the updating.
# do_update_Q <- function(Q){
#   M <- nrow(Q)
#   L <- ncol(Q)
#   test <- array(NA,c(M,L,3))
#   ind_unit_Q <- colSums(Q)==1
#   for (k in 1:M){
#     for (l in 1:L){
#       test[k,l,1] <- equal_unit_vec(Q[,l],k) # THIS MEANS A 1 IN E_M NEVER GETS UPDATED.
#       test[k,l,2] <- (sum(Q[k,])==3) && (Q[k,l]==1)
#       Q_cand      <- Q[,ind_unit_Q,drop=FALSE]
#       test[k,l,3] <- (Q[k,l]==0) && sum(Q[,l])==1 &&
#         (sum(sapply(1:ncol(Q_cand), function(c) equal_unit_vec(Q_cand[,c],which(Q[,l]==1))))==2)
#     }
#   }
#   apply(test,c(1,2),function(v) !any(v))
# }

#' determine to update the current element of Q_ml or not
#'
#' Function to test whether we need to update the current element of Q_ml. This
#' is needed in the constrained Gibbs sampler.
#'
#' @param Q a matrix with rows being machines and columns being protein landmarks (dimension)
#' @param k,l row and column indicators for checking wether to update Q_kl
#' @return logical value. TRUE for updating in constrained Gibbs sampler, FALSE
#' for skipping the updating.
#' @export
do_update_Q_one <- function(Q,k,l){ # this saves lots of speed:
  test <- rep(NA,3)
  ind_unit_Q <- colSums(Q)==1
  test[1] <- equal_unit_vec(Q[,l],k)
  test[2] <- (sum(Q[k,])==3) && (Q[k,l]==1)
  Q_cand  <- Q[,ind_unit_Q,drop=FALSE]
  test[3] <- (Q[k,l]==0) && sum(Q[,l])==1 &&
    (sum(sapply(1:ncol(Q_cand), function(c) equal_unit_vec(Q_cand[,c],which(Q[,l]==1))))==2)
  !any(test)
}

#' update the Q matrix element by element
#'
#' this function updates Q matrix by constrained Gibbs sampler
#' NB: - do we need M here? do we modify M after collapsing partner machines.
#'     - need to speed up.
#'
#' @param Y data
#' @param H matrix of machine usage indicators, rows for N subjects, columns for M machines
#' @param Q_old the Q matrix from the last scan in Gibbs sampler (of dimension M by L)
#' @param p prevalence of machines; a vector of length identical to the columns of H.
#' @param theta,psi True and false positive rates. Both are vectors of length L
#'
#' @return \itemize{
#' \item \code{datmat} a matrix for multivariate binary observations from the assumed boolean
#' matrix factorization,
#' \item \code{Q} A Q-matrix,
#' \item \code{H_star} A binary matrix of dimension K by M, for K clusters, and M factors,
#' \item \code{Z} A vector of length N for individual cluster indicators,
#' \item \code{xi} a binary matrix of dimension N by L for true presence or absence of proteins
#' \item \code{Eta} A binary matrix of dimension N by M for presence or absence of factors among
#' all subejcts.
#' }
update_Q_no_H <- function(Y,Q_old,z,t,mylist,p,theta,psi){
  M  <- nrow(Q_old)
  N  <- nrow(Y)
  L  <- ncol(Y)
  for (k in sample(1:M,replace=FALSE)){
    for (l in sample(1:L,replace=FALSE)){ # begin iteration over elements.
      #if(do_update_Q_one(Q_old,k,l)){ # begin an update if needed.
      L0 <- L1 <- 0
      Q_old[k,l] <- 0
      for (j in 1:t){
        L0 <- L0 + log_marginal_one_column(Y[z==mylist[j],l,drop=FALSE],Q_old[,l,drop=FALSE],p,theta[l],psi[l])
      }
      Q_old[k,l] <- 1
      for (j in 1:t){
        L1 <- L1 + log_marginal_one_column(Y[z==mylist[j],l,drop=FALSE],Q_old[,l,drop=FALSE],p,theta[l],psi[l])
      }
      curr_prob <- exp(L1- logSumExp(c(L0,L1)))
      #print(curr_prob)
      Q_old[k,l] <- rbinom(1,1,prob = curr_prob)
      #} # end an update if needed.
    }
  }  #end iteration over elements.
  Q_old#list(res = Q_old, ord = NA)#order_mat_byrow(Q_old)
}

metrop_flip <- function(x,curr_prob){
  if (x==1){
    if (curr_prob <= 0.5) {
      alpha <- 1
    } else{
      alpha <- exp(log(1-curr_prob) - log(curr_prob))
    }
    res <- rbinom(1,1,prob=1-alpha)
  } else{
    if (curr_prob >= 0.5){
      alpha <- 1
    } else{
      alpha <- exp(log(curr_prob) - log(1-curr_prob))
    }
    res <- rbinom(1,1,prob=alpha)
  } #end flipping.
  res
}

#' compute the full conditional probability of Q_ml given others
#'
#' Function to compute this log full conditional density, which will be used
#' in \link{update_Q}
#'
#' @param Y data
#' @param eta_star a matrix of machine usage indicators, rows for clusters, columns for M machines
#' @param Ql l-th column of the Q matrix
#' @param p prevalence of machines; a vector of length identical to the columns of H.
#' @param thetal,psil True and false positive rate (both are numbers between 0 and 1)
#'
#' @return log conditional probability of Q_ml given other unknown parameters
#' @export
log_pr_Qml_cond <- function(Yl,eta_star,Ql,thetal,psil){
  n1 <- sum(Yl)
  n0 <- sum(1-Yl)
  xil <- (eta_star%*%matrix(Ql,ncol=1) > 0.5)+0
  PR_mat <- (1-xil)*(n1*log(psil)+n0*log(1-psil))+
    xil*(n1*log(thetal)+n0*log(1-thetal))
}

#' Element-wise update the Q matrix
#'
#' This function updates Q matrix by Gibbs sampler (with option to do constarined
#' updates in identifiability constrained set or not)
#' NB: - do we need M here? do we modify M after collapsing partner machines.
#'     - need to speed up.
#' @param Y data
#' @param Q_old the Q matrix from the last scan in Gibbs sampler (of dimension M by L)
#' @param H matrix of machine usage indicators, rows for N subjects, columns for M machines
#' @param z a vector of individual cluster indicators
#' @param t the number of clusters at an iteration
#' @param mylist the ordered list that keeps the cluster ids in increasing order
#' @param p prevalence of machines; a vector of length identical to the columns of H.
#' @param theta,psi True and false positive rates. Both are vectors of length L
#' @param constrained Default to FALSE; Set to true if doing constrained Gibbs sampling
#' with identifiability constraints.
#'
#' @return A Q matrix with all its elements updated.
#' @export
update_Q <- function(Y,Q_old,H,z,t,mylist,p,theta,psi,constrained=FALSE){
  M  <- nrow(Q_old)
  N  <- nrow(Y)
  L  <- ncol(Y)
  if (constrained){ # constrained Gibbs sampling:
    for (k in sample(1:M,replace=FALSE)){
      for (l in sample(1:L,replace=FALSE)){ # begin iteration over elements.
        if (do_update_Q_one(Q_old,k,l)){ # begin an update if needed.
          L0 <- L1 <- 0
          Q_old[k,l] <- 0
          for (j in 1:t){ #Yl,eta_star,Ql,thetal,psil
            L0 <- L0 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                       H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
          }
          Q_old[k,l] <- 1
          for (j in 1:t){
            L1 <- L1 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                       H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
          }
          curr_prob <- exp(L1- logSumExp(c(L0,L1)))
          #print(curr_prob)
          #Q_old[k,l] <- metrop_flip(Q_old[k,l],curr_prob) # <-- if doing metroplized flipping.
          Q_old[k,l] <- rbinom(1,1,prob = curr_prob)
        }# end an update if needed.
      }
    }  #end iteration over elements.
  } else {
    for (k in sample(1:M,replace=FALSE)){
      for (l in sample(1:L,replace=FALSE)){ # begin iteration over elements.
        if (do_update_Q_one(Q_old,k,l)){ # begin an update if needed.
          L0 <- L1 <- 0
          Q_old[k,l] <- 0
          for (j in 1:t){ #Yl,eta_star,Ql,thetal,psil
            L0 <- L0 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                       H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
          }
          Q_old[k,l] <- 1
          for (j in 1:t){
            L1 <- L1 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                       H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
          }
          curr_prob <- exp(L1- logSumExp(c(L0,L1)))
          #print(curr_prob)

          #Q_old[k,l] <- metrop_flip(Q_old[k,l],curr_prob)
          Q_old[k,l] <- rbinom(1,1,prob = curr_prob)
        }# end an update if needed.
      }
    }  #end iteration over elements.
  }
  Q_old
}

#' Compute the full conditional probability of column Q_l given other unknown parameters
#'
#' Function to compute the log full conditional density for all patterns of the
#' l-th column of Q. This function is used in \link{update_Q_col_block}.
#'
#' @param Y data
#' @param eta_star a matrix of machine usage indicators, rows for clusters, columns for M machines
#' @param Ql_enumerate l-th column of the Q matrix
#' @param p prevalence of machines; a vector of length identical to the columns of H.
#' @param thetal,psil True and false positive rates. Both are vectors of length L
#'
#' @return a vector of log conditional probability of column Q_l taking each of
#' \code{2^nrow(Ql_enumerate)} given other unknown parameters,
#' the dimension will be identical to \code{2^nrow(Ql_enumerate)}.
#' @export
log_pr_Qml_cond_enum <- function(Yl,eta_star,Ql_enumerate,thetal,psil){
  n1 <- sum(Yl)
  n0 <- sum(1-Yl)
  xil <- (eta_star%*%matrix(Ql_enumerate,nrow=ncol(eta_star)) > 0.5)+0
  PR_mat <- (1-xil)*(n1*log(psil)+n0*log(1-psil))+
    xil*(n1*log(thetal)+n0*log(1-thetal))
}

#' Block update the Q matrix column by column
#'
#' This function updates Q matrix by Gibbs sampler (without identifiability
#' constraints)
#'
#' @param Y data
#' @param Q_old the Q matrix from the last scan in Gibbs sampler (of dimension M by L)
#' @param H matrix of machine usage indicators, rows for N subjects, columns for M machines
#' @param z a vector of individual cluster indicators
#' @param t the number of clusters at an iteration
#' @param mylist the ordered list that keeps the cluster ids in increasing order
#' @param p prevalence of machines; a vector of length identical to the columns of H.
#' @param theta,psi True and false positive rates. Both are vectors of length L
#'
#' @return A Q matrix with all elements updated.
#' @export
update_Q_col_block <- function(Y,Q_old,H,z,t,mylist,p,theta,psi){
  M  <- nrow(Q_old)
  N  <- nrow(Y)
  L  <- ncol(Y)
  Ql_enum  <- t(as.matrix(expand.grid(rep(list(0:1), M)),ncol=M))
  for (l in sample(1:L,replace=FALSE)){ # begin iteration over columns:
    L_enum  <- rep(0,2^M)
    for (j in 1:t){
      L_enum <- L_enum + log_pr_Qml_cond_enum(Y[z==mylist[j],l,drop=FALSE],
                                              H[j,,drop=FALSE],Ql_enum,theta[l],psi[l])
    }
    curr_prob <- exp(L_enum-logSumExp(L_enum))
    Q_old[,l] <- Ql_enum[,sample(1:2^M,1,prob = curr_prob)]
  }
  Q_old
}

#' MCMC sampling designed for binary factor analysis (pre-specified number of factors)
#'
#' This function performs MCMC sampling with user-specified options.
#' NB: 1) add flexibility to specify other parameters as fixed. 2) sample component-specific
#' parameters. 3) sample other model parameters. 4) add timing and printing functionality.
#' 5) add posterior summary functions.
#'
#' @param dat multivariate binary data (row for observations and column for dimensions4)
#' @param model_options Specifying assumed model options:
#' \itemize{
#' \item \code{n} The number of subjects.
#' \item \code{t_max} The maximum (guessed) number of clusters in the data during
#' posterior inference
#' \item \code{m_max} For a model with pre-specified number of factors, \code{m_max};
#' In an infinite factor model with unknown number of factors, \code{m_max} is
#' the maximum (guessed) number of factors during the posterior inference (see \code{\link{slice_sampler}});
#' one can increase this number if this pacakge recommends so;
#' \item \code{a_theta, a_psi} hyperparameters for true and false positive rates;
#' a_theta and a_psi are both a vector of length two.
#' \item \code{a_alpha, b_alpha} Only for infinite factor model  -
#' Gamma hyperparameter for the hyperprior on \code{alpha}. (See \code{\link{slice_sampler}})
#' \item \code{log_v} The charaster string representing the prior
#' distribution for the number of true clusters, e.g.,
#' \code{"function(k) {log(0.1) + (k-1)*log(0.9)}"}. It is usually pre-computed log of the coefficients in Mixture of Finite Mixtures
#' (Miller and Harrison, 2017, JASA). Use this code:
#' \code{coefficients(eval(parse(text=model_options0$log_pk)),
#' model_options0$gamma,
#' model_options0$n,
#' model_options0$t_max+1)}
#' }
#' The following are used if one needs to prefix a few unnown parameters to
#' their respective true or other values
#' \itemize{
#' \item \code{Q} Q matrix
#' \item \code{theta} a vector of true positive rates
#' \item \code{psi} a vector of false positive rates
#' \item \code{p} a vector of machine prevalences
#' \item \code{alpha} For pre-specified number of factors, the hyperparameter
#' for \code{Beta(alpha/m_max,1)} (can set to \code{m_max}); For infinite factor model, the hyperparameter
#' for IBP (can set to 1).
#' }
#' options for specifying data, sample size, max cluster number,
#' coefficients in MFM (Miller and Harrison 2017 JASA), gamma parameter in the MFM
#' Dirchlet prior, number of intermediate Gibbs scan to arrive at the launch state,
#' and other hyperparamter specification if needed, n_total for total number of
#' mcmc iterations and n_keep for the number of samples kept for posterior inference.
#' Note that the options involve other parameters for sampling hyperparameters such as
#' alpha in the Indian Buffet Process.
#' @param mcmc_options Options for MCMC sampling:
#' \itemize{
#' \item \code{n_total} total number of MCMC iterations
#' \item \code{n_keep} number of iterations kept
#' \item \code{n_split} the number of restricted Gibbs scan to arrive at a launch state;
#' see \link{restricted_gibbs}
#' \item \code{print_mod} print machine usage profiles for all clusters and plot
#' the current Q matrix.
#' \item \code{block_update_H} update row of H (if \code{TRUE}) or not
#' (if \code{NULL} or FALSE - must be so for \code{\link{slice_sampler}}). Then no constraints
#' about identifiability is imposed upon Q at any iterations.
#' \item \code{block_update_Q} update columns of Q (if \code{TRUE}) or not
#' (if \code{NULL} or FALSE - must be so for \code{\link{slice_sampler}}). Then no constraints
#' about identifiability is imposed upon Q at any iterations.
#' }
#'
#' @example /inst/example/simulation_fixed_M.R
#'
#' @return posterior samples for quantities of interest. It is a list comprised of the following elements
#' \itemize{
#' \item \code{t_samp}
#' \item \code{z_samp}
#' \item \code{N_samp}
#' \item \code{keepers} indices of MCMC samples kept for inference;
#' \item \code{H_star_samp}
#' \item \code{alpha_samp}
#' } The following are recorded if they are not fixed in a priori:
#' \itemize{
#' \item \code{Q_samp}
#' \item \code{theta_samp}
#' \item \code{psi_samp}
#' \item \code{p_samp}
#' }
#' @export
sampler <- function(dat,model_options,mcmc_options){
  # dat = simu_dat
  # model_options = model_options0
  # mcmc_options = mcmc_options0
  n <- nrow(dat) # number of observations
  L <- ncol(dat) # number of dimensions (protein landmarks).

  n_total <- mcmc_options$n_total # total number of mcmc iterations.
  n_keep  <- mcmc_options$n_keep  # toral number of samples kept for posterior inference.
  keepers <- seq(ceiling(n_total/n_keep),n_total,len=n_keep)
  keep_index <- 0                 # the index to keep during MCMC inference.
  n_split <- mcmc_options$n_split # number of intermediate GS scan to arrive at launch state.
  block_update_Q <- !is.null(mcmc_options$block_update_Q) && mcmc_options$block_update_Q # TRUE for updating columns of Q. FALSE otherwise.
  block_update_H <- !is.null(mcmc_options$block_update_H) && mcmc_options$block_update_H # TRUE for updating rows of H. FALSE otherwise.

  # set options (need to fill in as specific paramters are needed):
  t_max <- model_options$t_max # maximum allowable number of clusters.
  m_max <- model_options$m_max # maximum number of machines (truncated to m_max).
  b     <- model_options$b # gamma parameter in MFM - hyperparameter for Dirichlet distn. needs to be sampled?????
  log_v <- model_options$log_v #coefficients for MFM.

  if (is.null(model_options$Q)){
    Q_samp <- array(0,c(m_max,ncol(dat),n_keep))
    #Q      <- simulate_Q(m_max,L) # random initialization, could be impproved by a smarter starting Q matrix.
    Q      <- simulate_Q_dat(m_max,dat,0.5) # random initialization, could be impproved by a smarter starting Q matrix.
  }else{
    Q     <- model_options$Q # suppose we are given Q.
  }

  if (block_update_H){
    H_star_enumerate <- as.matrix(expand.grid(rep(list(0:1), m_max)),ncol=m_max) # all binary patterns for machine usage profiles. 2^m_max of them.
  }
  # initialize the sampling chain:
  t <- 1        # number of clusters.
  z <- rep(1,n) # z[i] is the cluster ID for observation i.
  mylist <- rep(0,t_max+3); mylist[1] <- 1  # mylist[1:t] is the list of active cluster IDs.
  # mylist is maintained in increasing order for 1:t, and
  # is 0 after that.
  c_next <- 2                    # an available cluster ID to be used next.
  N <- rep(0,t_max+3); N[1] <- n # N[c] is the size of cluster c. Note that N is different from n.

  # H_star <- matrix(0,nrow=t_max+3,ncol=m_max) # parameter for clusters. In our application, the machine usage profiles.
  # NOTE: need to consider SMARTER initial values for this. How to sample H_star? Not needed in conjugate algorithm.

  log_p <- rep(0,n+1) # probability for sampling z; n+1 because at most there could be n clusters;
  # and sometimes one needs to assign an observation to a new cluster - by current bookkeeping,
  # we are not efficient in memory and create a new cluster ID - hence +1.
  zs <- rep(1,n)  # temporary cluster indicators for split-merge assignments.
  S  <- rep(0,n)  # temporary variable for indicies used during split-merge step.

  log_Nb <- log(1:n)+b # the multipliers needed when assigning an observation to an existing cluster. Restaurant process stuff.

  # variables for samples kept:
  t_samp <- rep(0,n_total)
  N_samp <- matrix(0,nrow=t_max+3,ncol=n_total)
  z_samp <- matrix(0,nrow=n,ncol=n_keep) # posterior samples of cluster indicators.
  H_star_samp <- array(0,c(t_max+3,m_max,n_total))# component-specific machine profiles.

  if (is.null(model_options$theta)){
    theta_samp <- matrix(0,nrow=ncol(dat),ncol=n_keep)
    a_theta <- model_options$a_theta
    theta <- sapply(1:L,function(i){rbeta(1,a_theta[1],a_theta[2])}) # initialization.
  } else{
    theta   <- model_options$theta # use specified true positive rates.
  }
  if (is.null(model_options$psi)){
    psi_samp   <- matrix(0,nrow=ncol(dat),ncol=n_keep)
    a_psi      <- model_options$a_psi
    psi <- sapply(1:L,function(i){rbeta(1,a_psi[1],a_psi[2])}) # initialization.
  } else{
    psi     <- model_options$psi # use specified false positive rates.
  }

  if (is.null(model_options$p0)){
    p_samp <- matrix(0,nrow=m_max,ncol=n_total)
    p      <- rep(0.5,m_max) # initialization.
  } else{
    p      <- model_options$p0 # fix prevalence.
  }

  if (is.null(model_options$alpha)){
    alpha_samp <- rep(0,n_total)  # hyperparameter for machine prevalence.
    alpha      <- m_max           # initialization.
  }

  cat("==[rewind] Start MCMC: ==\n")
  for (iter in 1:n_total){
    if (iter%%mcmc_options$print_mod==0){
      cat("==[rewind] iteration: ", iter, "==\n");
      cat("==[rewind] factor indicators (machines usages) for t=",t," clusters:==\n")
      print(H_star[,colSums(H_star)!=0,drop=FALSE])
      image(Q[colSums(H_star)!=0,drop=FALSE,])
    }
    # update cluster indicators z for all subjects - one complete Gibbs scan to refine clusters:
    for (i in 1:n){ # iterate over subjects:
      # remove obs i from its current cluster:
      c    <- z[i]
      N[c] <- N[c]-1
      if (N[c]>0){ # if there are observations left after removal of i, use c_next:
        c_prop <- c_next
      } else{# if i was by itself, use its cluster ID again.
        c_prop <- c
        mylist <- ordered_remove(c,mylist,t)
        t <- t-1 # total number of clusters is now reduced by one.
      }
      # compute probability for Gibbs updating - the probability of assigning a subject
      # to a cluster.
      for (j in 1:t){
        cc       <- mylist[j]
        log_p[j] <- log_Nb[N[cc]]+log_marginal(rbind(dat[(z==cc)[-i],,drop=FALSE],dat[i,]),Q,p,theta,psi)-
          log_marginal(dat[(z==cc)[-i],,drop=FALSE],Q,p, theta,psi) # existing cluster.
      }
      log_p[t+1] <- log_v[t+1]-log_v[t] + log(b) + log_marginal(dat[i,,drop=FALSE],Q,p,theta,psi) # new cluster.

      j <- sample(t+1,1,prob = exp(log_p[1:(t+1)]))

      # add obs i to its new cluster:
      if (j<=t){
        c <- mylist[j]
      } else{
        c      <- c_prop
        mylist <- ordered_insert(c,mylist,t)
        t      <- t+1
        c_next <- ordered_next(mylist)
        if(t>t_max) {stop("==[rewind] Sampled t has exceeded t_max. Increast t_max and retry.==")}
      }
      z[i] <- c
      N[c] <- N[c] + 1

    }# END iteration over subjects.

    # if (use_splitmerge){}
    res_split_merge <- split_merge(dat,z,zs,S,mylist,N,t,b,log_v,n_split,Q,p,theta,psi)
    t <- res_split_merge$t
    z <- res_split_merge$z
    N <- res_split_merge$N
    mylist <- res_split_merge$mylist
    c_next <- ordered_next(mylist)
    if(t>t_max) {stop("==[rewind]Sampled t has exceeded t_max. Increast t_max and retry.==")}

    cat("==[rewind]#subjects for t=",t," clusters:", N[N!=0],"==\n")

    ## Simulate xi for posterior predictive checking:
    # #if (is.null(model_options$Q)){
    # xi_samp <- array(0,c(t_max+3,L,n_total)) # add this to the results.
    # for (j in 1:t){
    #   xi_samp[mylist[j],,iter]  <- update_xi(dat[z==mylist[j],,drop=FALSE],p,theta,psi)
    # }
    # xi_star <- xi_samp[mylist[1:t],,iter]
    # #}

    # update machine usage profile given clusters: (check here for sampling H)
    H_star_redun <- matrix(0,nrow=t_max+3,ncol=m_max)
    for (j in 1:t){
      if (block_update_H){
        curr_pattern_log_p   <- log_full(dat[z==mylist[j],,drop=FALSE],H_star_enumerate,Q,p,theta,psi)
        curr_pattern_p       <- exp(curr_pattern_log_p-logSumExp(curr_pattern_log_p))
        curr_ind             <- sample(1:(2^m_max),1,prob=curr_pattern_p)
        H_star_redun[mylist[j],] <- H_star_enumerate[curr_ind,]
      } else{
        for (m in 1:m_max){
          H_star_redun[mylist[j],m] <- 0
          L0            <- log_full(dat[z==mylist[j],,drop=FALSE],
                                    H_star_redun[mylist[j],,drop=FALSE],Q,p,theta,psi)
          H_star_samp[mylist[j],m,iter] <- 1
          L1            <- log_full(dat[z==mylist[j],,drop=FALSE],
                                    H_star_redun[mylist[j],,drop=FALSE],Q,p,theta,psi)
          curr_eta_p    <- exp(L1-matrixStats::logSumExp(c(L0,L1)))
          H_star_redun[mylist[j],m] <-rbinom(1,1,curr_eta_p)
          #H_star_samp[mylist[j],m,iter] <- metrop_flip(H_star_samp[mylist[j],m,iter],curr_eta_p)
        }
      }
    }

    H_star_samp[,,iter] <- H_star_redun
    H_star <- H_star_redun[mylist[1:t],,drop=FALSE]

    # update Q:
    if (is.null(model_options$Q)){
      #if(iter%%5==0){
      if (block_update_Q){
        Q <- update_Q_col_block(dat,Q,H_star,z,t,mylist,p,theta,psi)
      } else{
        #}else{
        Q <- update_Q(dat,Q,H_star,z,t,mylist,p,theta,psi)
      }
      #}
      #Q <- update_Q_no_H(dat,Q,z,t,mylist,p,theta,psi)
      # if (sum(colSums(H_star)==0)>0){
      #   Q[colSums(H_star)==0,] <- simulate_Q_dat(sum(colSums(H_star)==0),dat)# # reinitialize if there are unused machines vs sample from prior?. Now we can see
      # }
    }

    # update true/false positive rates - theta/psi:
    if (is.null(model_options$theta) && is.null(model_options$psi)){
      res_update_pr <- update_positive_rate(dat,H_star_samp[z,,iter],Q,a_theta,a_psi)
      theta <- res_update_pr$theta
      psi   <- res_update_pr$psi
    }

    # update hyperparameter for {p_m}:
    if (is.null(model_options$alpha)){
      alpha <- update_alpha(H_star,t,m_max)
      alpha_samp[iter] <- alpha
    }

    # update prevalence vector {p_m}:
    if (is.null(model_options$p0)){
      p    <- update_prevalence(H_star,alpha,m_max)
      p_samp[,iter]    <- p
    }

    #
    # Record MCMC results:
    #
    t_samp[iter] <- t
    for (j in 1:t){
      N_samp[mylist[j],iter] <- N[mylist[j]]
    }

    if (iter==keepers[keep_index+1]){
      keep_index <- keep_index +1
      for (i in 1:n){
        z_samp[i,keep_index] <- z[i]
      }
      if (is.null(model_options$theta)){
        theta_samp[,keep_index] <- theta # of length L.
      }
      if (is.null(model_options$psi)){
        psi_samp[,keep_index]   <- psi
      }
      if (is.null(model_options$Q)){
        Q_samp[,,keep_index] <- Q
      }
    }
  }# END mcmc iterations.

  res <- list(t_samp=t_samp,N_samp=N_samp,z_samp=z_samp,keepers=keepers,
              H_star_samp = H_star_samp,
              alpha_samp=alpha_samp)
  if (is.null(model_options$theta)){
    res$theta_samp <- theta_samp # of length L.
  }
  if (is.null(model_options$psi)){
    res$psi_samp   <- psi_samp
  }
  if (is.null(model_options$Q)){
    res$Q_samp <- Q_samp
    #res$xi_samp <- xi_samp
  }
  if (is.null(model_options$p0)){
    res$p_samp <- p_samp
  }
  res
}

#  slice sampler for a unknown number of machines:
#' compute the log of density of inactive probability given
#' the previous one (Teh et al., 07')
#'
#' This function is used in sampling infinite number of columns in H
#'
#' @param log_p0 log probability (the argument of this function)
#' @param alpha hyperparameter for Indian Buffet Process (Infinite version)
#' @param t the number of rows in the IBP
#'
#' @return a value corresponding to the log density f(log_p0)
#' @examples
#'
#' n_grid  <- 10000
#' p0_grid <- seq(log(0.001),0,len=n_grid)
#' y_den  <- log_f_logden(p0_grid,5,8)
#'
#' y <- ars(n_grid,log_f_logden,log_f_logprima,
#' c(-6,-4,-2,-1),m=4,
#' ub=TRUE,xub=0,alpha=5,t=8)
#' hist(y,breaks="Scott",main='Adaptive Rjection Sampling',freq=FALSE)
#' rug(y)
#' points(density(y),type="l",col="blue")
#' points(p0_grid,exp(y_den-rewind:::logsumexp(y_den)-log(diff(p0_grid)[2])),
#'       type="l",col="red") # <-- true density; log_f_logden is only correct upto a proportionality constant.
#'
#' legend("topleft",c("Sample Density","True Density"),lty=c(1,1),col=c("blue","red"),
#'       bty="n")
#'
#' @export
log_f_logden <- function(log_p0,alpha,t){
  #if (p0 > p0_minus_one || p0 < 0){return(-Inf)}
  res <- 0
  for (j in 1:t){
    res  <- res+(1/j)*(1-exp(log_p0))^j
  }
  alpha*res+(alpha-1)*log_p0+
    t*log(1-exp(log_p0))+log_p0
}


#' compute the deriviative of the log of density of inactive probability given
#' the previous one (Teh et al., 07')
#'
#' This function is used in sampling infinite number of columns in H
#'
#' @inheritParams log_f_logden
#'
#' @return a value corresponding to the derivative of the log density d/du f(log_p0)
#' @seealso \code{\link{log_f_logden}} provides an example
#' @export
log_f_logprima <- function(log_p0,alpha,t){
  #if (p0 > p0_minus_one || p0 < 0){return(-Inf)}
  res <- 0
  for (j in 1:t){
    res  <- res - exp(log_p0)*(1-exp(log_p0))^(j-1)
  }
  alpha-t*exp(log_p0)/(1-exp(log_p0))+
    alpha*res
}

#' MCMC sampling designed for binary factor analysis (unknown number of factors)
#'
#' This function performs MCMC sampling with user-specified options.
#' NB: 1) add flexibility to specify other parameters as fixed.
#' 2) sample component-specific
#' parameters.
#' 3) sample other model parameters.
#' 4) add timing and printing functionality.
#' 5) add posterior summary functions.
#'
#' @inheritParams sampler
#'
#' @example /inst/example/simulation_unknown_M.R
#'
#' @return posterior samples for quantities of interest. It is a list comprised of the following elements
#' \itemize{
#' \item \code{t_samp}
#' \item \code{z_samp}
#' \item \code{N_samp}
#' \item \code{keepers} indices of MCMC samples kept for inference;
#' \item \code{H_star_samp}
#' \item \code{alpha_samp}
#' \item \code{m_plus_samp} This is unique to slice_sampler
#' \item \code{m0_samp} This is unique to slice_sampler
#' } The following are recorded if they are not fixed in a priori:
#' \itemize{
#' \item \code{Q_samp}
#' \item \code{theta_samp}
#' \item \code{psi_samp}
#' \item \code{p_samp}
#' }
#' @export
slice_sampler <- function(dat,model_options,mcmc_options){
  # test:
  dat = simu_dat
  model_options = model_options0
  mcmc_options = mcmc_options0

  n <- nrow(dat) # number of observations
  L <- ncol(dat) # number of dimensions (protein landmarks).

  n_total <- mcmc_options$n_total # total number of mcmc iterations.
  n_keep  <- mcmc_options$n_keep  # toral number of samples kept for posterior inference.
  keepers <- seq(ceiling(n_total/n_keep),n_total,len=n_keep)
  keep_index <- 0                 # the index to keep during MCMC inference.
  n_split <- mcmc_options$n_split # number of intermediate GS scan to arrive at launch state.
  block_update_Q <- mcmc_options$block_update_Q # TRUE for updating columns of Q. FALSE otherwise.

  # set options (need to fill in as specific paramters are needed):
  t_max <- model_options$t_max # maximum allowable number of clusters.
  m_max <- model_options$m_max # maximum number of machines (truncated to m_max).
  b     <- model_options$b # gamma parameter in MFM - hyperparameter for Dirichlet distn. needs to be sampled?????
  log_v <- model_options$log_v #coefficients for MFM.

  if (is.null(model_options$m_plus)){
    m_plus_samp <- m0_samp <- rep(NA,n_total) # what does m_plus mean?
    m_plus <- 3 # initialization.
    m0     <- 0 # initialization.
    m_both <- m_plus+m0
  } else{
    m_plus <- model_options$m_plus
    m0     <- 0
    m_both <- m_plus
  }

  if (is.null(model_options$Q)){
    Q_samp <- array(0,c(m_max,ncol(dat),n_keep))
    #Q      <- simulate_Q(m_max,L) # random initialization, could be impproved by a smarter starting Q matrix.
    Q      <- simulate_Q_dat(m_both,dat,0.5) # random initialization, could be impproved by a smarter starting Q matrix.
  }else{
    Q     <- model_options$Q # suppose we are given Q.
  }

  # initialize the sampling chain:
  t <- 1        # number of clusters.
  z <- rep(1,n) # z[i] is the cluster ID for observation i.
  mylist <- rep(0,t_max+3); mylist[1] <- 1  # mylist[1:t] is the list of active cluster IDs.
  # mylist is maintained in increasing order for 1:t, and
  # is 0 after that.
  c_next <- 2                    # an available cluster ID to be used next.
  N <- rep(0,t_max+3); N[1] <- n # N[c] is the size of cluster c. Note that N is different from n.

  # H_star <- matrix(0,nrow=t_max+3,ncol=m_max) # parameter for clusters. In our application, the machine usage profiles.
  # NOTE: need to consider SMARTER initial values for this. How to sample H_star? Not needed in conjugate algorithm.

  log_p <- rep(0,n+1) # probability for sampling z; n+1 because at most there could be n clusters;
  # and sometimes one needs to assign an observation to a new cluster - by current bookkeeping,
  # we are not efficient in memory and create a new cluster ID - hence +1.
  zs <- rep(1,n)  # temporary cluster indicators for split-merge assignments.
  S  <- rep(0,n)  # temporary variable for indicies used during split-merge step.

  log_Nb <- log(1:n)+b # the multipliers needed when assigning an observation to an existing cluster. Restaurant process stuff.

  # variables for samples kept:
  t_samp <- rep(0,n_total)
  N_samp <- matrix(0,nrow=t_max+3,ncol=n_total)
  z_samp <- matrix(0,nrow=n,ncol=n_keep) # posterior samples of cluster indicators.
  H_star_samp <- array(0,c(t_max+3,m_max,n_total))# component-specific machine profiles.

  if (is.null(model_options$theta)){
    theta_samp <- matrix(0,nrow=ncol(dat),ncol=n_keep)
    a_theta <- model_options$a_theta
    theta <- sapply(1:L,function(i){rbeta(1,a_theta[1],a_theta[2])}) # initialization.
  } else{
    theta   <- model_options$theta # use specified true positive rates.
  }
  if (is.null(model_options$psi)){
    psi_samp   <- matrix(0,nrow=ncol(dat),ncol=n_keep)
    a_psi      <- model_options$a_psi
    psi <- sapply(1:L,function(i){rbeta(1,a_psi[1],a_psi[2])}) # initialization.
  } else{
    psi     <- model_options$psi # use specified false positive rates.
  }

  if (is.null(model_options$p_both)){
    p_samp <- matrix(0,nrow=m_max,ncol=n_total)
    p      <- rep(0.5,m_both) # initialization.
  } else{
    p      <- model_options$p_both # fix prevalence.
  }

  if (is.null(model_options$alpha)){
    alpha_samp <- rep(0,n_total)  # hyperparameter for machine prevalence.
    alpha      <- 5           # initialization.
  }

  cat("==[rewind] Start MCMC: ==\n")
  for (iter in 1:n_total){
    if (iter%%mcmc_options$print_mod==0){
      cat("==[rewind] iteration: ", iter, "==\n");
      cat("==[rewind] factor indicators (machines usages) for t=",t," clusters:==\n")
      print(H_star[,colSums(H_star)!=0,drop=FALSE])
      image(Q[colSums(H_star)!=0,drop=FALSE,])
    }
    # update cluster indicators z for all subjects - one complete Gibbs scan to refine clusters:
    for (i in 1:n){ # iterate over subjects:
      # remove obs i from its current cluster:
      c    <- z[i]
      N[c] <- N[c]-1
      if (N[c]>0){ # if there are observations left after removal of i, use c_next:
        c_prop <- c_next
      } else{# if i was by itself, use its cluster ID again.
        c_prop <- c
        mylist <- ordered_remove(c,mylist,t)
        t <- t-1 # total number of clusters is now reduced by one.
      }
      # compute probability for Gibbs updating - the probability of assigning a subject
      # to a cluster.
      for (j in 1:t){
        cc       <- mylist[j]
        log_p[j] <- log_Nb[N[cc]]+log_marginal(rbind(dat[(z==cc)[-i],,drop=FALSE],dat[i,]),Q,p,theta,psi)-
          log_marginal(dat[(z==cc)[-i],,drop=FALSE],Q,p, theta,psi) # existing cluster.
      }
      log_p[t+1] <- log_v[t+1]-log_v[t] + log(b) + log_marginal(dat[i,,drop=FALSE],Q,p,theta,psi) # new cluster.

      j <- sample(t+1,1,prob = exp(log_p[1:(t+1)]))

      # add obs i to its new cluster:
      if (j<=t){
        c <- mylist[j]
      } else{
        c      <- c_prop
        mylist <- ordered_insert(c,mylist,t)
        t      <- t+1
        c_next <- ordered_next(mylist)
        if(t>t_max) {stop("==[rewind] Sampled t has exceeded t_max. Increast t_max and retry.==")}
      }
      z[i] <- c
      N[c] <- N[c] + 1

    }# END iteration over subjects.

    # if (use_splitmerge){}
    res_split_merge <- split_merge(dat,z,zs,S,mylist,N,t,b,log_v,n_split,Q,p,theta,psi)
    t <- res_split_merge$t
    z <- res_split_merge$z
    N <- res_split_merge$N
    mylist <- res_split_merge$mylist
    c_next <- ordered_next(mylist)
    if(t>t_max) {stop("==[rewind]Sampled t has exceeded t_max. Increast t_max and retry.==")}

    cat("==[rewind]#subjects for t=",t," clusters:", N[N!=0],"==\n")

    # update machine usage profile given clusters: (check here for sampling H)
    # H_star_enumerate <- as.matrix(expand.grid(rep(list(0:1), m_both)),ncol=m_both) # all binary patterns for machine usage profiles. 2^m_max of them.

    if (iter ==1){
      H_star_redun <- matrix(0,nrow=t_max+3,ncol=m_max)
      for (m in 1:m_both){
        for (j in 1:t){
          H_star_redun[mylist[j],m] <- 0
          L0            <- log_full(dat[z==mylist[j],,drop=FALSE],
                                    H_star_redun[mylist[j],1:m_both,drop=FALSE],Q,p,theta,psi)
          H_star_redun[mylist[j],m] <- 1
          L1            <- log_full(dat[z==mylist[j],,drop=FALSE],
                                    H_star_redun[mylist[j],1:m_both,drop=FALSE],Q,p,theta,psi)
          curr_eta_p    <- exp(L1-matrixStats::logSumExp(c(L0,L1)))
          H_star_redun[mylist[j],m] <- rbinom(1,1,curr_eta_p)
          #H_star_samp[mylist[j],m,iter] <- metrop_flip(H_star_samp[mylist[j],m,iter],curr_eta_p)
        }
      }
      # good code:
      # for (m in 1:m_both){
      #   for (j in 1:t){
      #     H_star_samp[mylist[j],m,iter] <- 0
      #     L0            <- log_full(dat[z==mylist[j],,drop=FALSE],
      #                               matrix(H_star_samp[mylist[j],1:m_both,iter],nrow=1),Q,p,theta,psi)
      #     H_star_samp[mylist[j],m,iter] <- 1
      #     L1            <- log_full(dat[z==mylist[j],,drop=FALSE],
      #                               matrix(H_star_samp[mylist[j],1:m_both,iter],nrow=1),Q,p,theta,psi)
      #     curr_eta_p    <- exp(L1-matrixStats::logSumExp(c(L0,L1)))
      #     H_star_samp[mylist[j],m,iter] <- rbinom(1,1,curr_eta_p)
      #     #H_star_samp[mylist[j],m,iter] <- metrop_flip(H_star_samp[mylist[j],m,iter],curr_eta_p)
      #   }
      #   # curr_pattern_log_p   <- log_full(dat[z==mylist[j],,drop=FALSE],H_star_enumerate,Q,p,theta,psi)
      #   # curr_pattern_p       <- exp(curr_pattern_log_p-logSumExp(curr_pattern_log_p))
      #   # curr_ind             <- sample(1:(2^m_both),1,prob=curr_pattern_p)
      #   # H_star_samp[mylist[j],1:m_both,iter] <- H_star_enumerate[curr_ind,]
      # }
    } else{
      for (m in 1:m_both){
        for (j in 1:t){
          H_star_redun[mylist[j],m] <- 0
          ind_prop_mplus <- which(colSums(H_star_redun[,1:m_plus,drop=FALSE])!=0)
          L0            <- -log(min(p[ind_prop_mplus]))+log_full(dat[z==mylist[j],,drop=FALSE],
                                                                 H_star_redun[mylist[j],1:m_both,drop=FALSE],Q,p,theta,psi)
          H_star_redun[mylist[j],m] <- 1
          ind_prop_mplus <- which(colSums(H_star_redun[,1:m_plus,drop=FALSE])!=0)
          L1            <- -log(min(p[ind_prop_mplus]))+log_full(dat[z==mylist[j],,drop=FALSE],
                                                                 H_star_redun[mylist[j],1:m_both,drop=FALSE],Q,p,theta,psi)
          curr_eta_p    <- exp(L1-matrixStats::logSumExp(c(L0,L1)))
          H_star_redun[mylist[j],m] <- rbinom(1,1,curr_eta_p)
          #H_star_samp[mylist[j],m,iter] <- metrop_flip(H_star_samp[mylist[j],m,iter],curr_eta_p)
        }
        # curr_pattern_log_p   <- log_full(dat[z==mylist[j],,drop=FALSE],H_star_enumerate,Q,p,theta,psi)
        # curr_pattern_p       <- exp(curr_pattern_log_p-logSumExp(curr_pattern_log_p))
        # curr_ind             <- sample(1:(2^m_both),1,prob=curr_pattern_p)
        # H_star_samp[mylist[j],1:m_both,iter] <- H_star_enumerate[curr_ind,]
      }
    }

    H_star_samp[,,iter] <- H_star_redun
    H_star <- H_star_redun[mylist[1:t],1:m_both,drop=FALSE]

    # update Q:
    if (is.null(model_options$Q)){
      #if(iter%%5==0){
      if (block_update_Q){
        Q <- update_Q_col_block(dat,Q,H_star,z,t,mylist,p,theta,psi)
      } else{
        #}else{
        Q <- update_Q(dat,Q,H_star,z,t,mylist,p,theta,psi)
      }
      #}
      #Q <- update_Q_no_H(dat,Q,z,t,mylist,p,theta,psi)
      # if (sum(colSums(H_star)==0)>0){
      #   Q[colSums(H_star)==0,] <- simulate_Q_dat(sum(colSums(H_star)==0),dat)# # reinitialize if there are unused machines vs sample from prior?. Now we can see
      # }
    }

    # drop inactive ones and add newly active ones:
    ind_plus <- which(colSums(H_star)!=0)
    m_plus   <- length(ind_plus)  # <-- what if ind_plus is empty?
    if (m_plus==0){stop("==[rewind] no active factors (machines), initialize Q with more
                        plausible values and retry.==")}
    H_star <- H_star[,ind_plus,drop=FALSE]

    # deal with inactive parts to be added to columns of H_star:
    # sample slice variable:
    p <- rep(NA,m_plus)
    for (m in 1:(m_plus)){
      sm <- sum(H_star[,m]) # number of clusters using machine m.
      p[m] <- rbeta(1,sm,1+t-sm)
    }
    s      <- runif(1,0,min(p)) # <---- p for active ones. we have already just focused on active ones.
    cat("s = ",s," ==\n")
    # sample inactive p0 from 1 to m0 (adaptive rejection sampling):
    p0    <- 1
    count <- 1
    curr_alpha <- alpha
    curr_t     <- t
    cat("alpha = ", alpha,"==\n")
    while (p0[count] >= s){
      curr_log_p0_samp <- ars(1,log_f_logden,log_f_logprima,
                              log(c(p0[count]/10,p0[count]/2,p0[count]/10*9)),m=3,
                              lb=TRUE,xlb=log(0.00001),
                              ub=TRUE,xub=log(p0[count]),alpha=curr_alpha,t=curr_t)
      p0 <- c(p0, exp(curr_log_p0_samp))
      count <- count + 1
    }
    if (length(p0)>1){
      p0 <- p0[-1]
      m0 <- length(p0)
      p  <- c(p,p0)
    } else{
      p0 <- m0 <- 0
    }
    m_both <- m_plus + m0 # <- what if both are zero? A: will show an error meesage above.
    print(m_plus)
    print(m0)
    if (m_both > m_max){stop("==[rewind] total #factors (machines) exceeds m_max=",m_max,", increase m_max and retry.==")}

    Q <- Q[ind_plus,,drop=FALSE]
    if (m0>0){
      H_star <- cbind(H_star,matrix(0,nrow=t,ncol=m0))
      Q <- rbind(Q,simulate_Q_dat(m0,dat,p = 0.5)) # pad Q with extra rows.
    }

    # sample p given the newly updated active feature matrix H_star_samp:
    # 1. The p vector corresponds to those probabilities greater than s (given the slice).
    # 2. When integrating over xi_il, we only need to consider these p.

    #
    # how to change H_star here?????
    #



    # update true/false positive rates - theta/psi:
    if (is.null(model_options$theta) && is.null(model_options$psi)){
      res_update_pr <- update_positive_rate(dat,H_star_samp[z,1:m_both,iter],Q,a_theta,a_psi)
      theta <- res_update_pr$theta
      psi   <- res_update_pr$psi
    }

    # update hyperparameter for {p_m}:
    if (is.null(model_options$alpha)){ # <--- change this as well, use gamma prior, things will be easier!!!
      alpha <- rgamma(1,shape=model_options$a_alpha+m_plus,
                      rate= model_options$b_alpha+sum(1/(1:t)))
      alpha_samp[iter] <- alpha
    }

    # update prevalence vector {p_m}:
    if (is.null(model_options$p_both)){
      #p    <- update_prevalence(H_star,alpha,m_max)
      p_samp[1:m_both,iter]    <- p
    }

    #
    # Record MCMC results:
    #
    t_samp[iter] <- t
    for (j in 1:t){
      N_samp[mylist[j],iter] <- N[mylist[j]]
    }

    if (iter==keepers[keep_index+1]){
      keep_index <- keep_index +1
      for (i in 1:n){
        z_samp[i,keep_index] <- z[i]
      }
      if (is.null(model_options$theta)){
        theta_samp[,keep_index] <- theta # of length L.
      }
      if (is.null(model_options$psi)){
        psi_samp[,keep_index]   <- psi
      }
      if (is.null(model_options$Q)){
        Q_samp[1:m_both,,keep_index] <- Q
      }
    }
    }# END mcmc iterations.

  res <- list(t_samp=t_samp,N_samp=N_samp,z_samp=z_samp,keepers=keepers,
              H_star_samp = H_star_samp,
              alpha_samp=alpha_samp,
              m_plus_samp = m_plus,
              m0_samp = m0)
  if (is.null(model_options$theta)){
    res$theta_samp <- theta_samp # of length L.
  }
  if (is.null(model_options$psi)){
    res$psi_samp   <- psi_samp
  }
  if (is.null(model_options$Q)){
    res$Q_samp <- Q_samp
    #res$xi_samp <- xi_samp
  }
  if (is.null(model_options$p0)){
    res$p_samp <- p_samp
  }
  res
}

