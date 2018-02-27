# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# Define helper functions:

#' insert a cluster ID prior to "index"; mylist keeps a list of cluster IDs
#'
#' This function is called when adding a new cluster with more than one subjects
#'
#' @param index an integer of cluster ID to be inserted
#' @param mylist a vector of integer cluster IDs
#' @param t the effective number of clusters
#'
#' @examples
#' ordered_insert(4,c(1:3,5:10),7)
#' ordered_insert(4,c(1:10),7)
#'
#' @return an vector of cluster IDs with an inserted ID \code{index}.
#' @export
ordered_insert <- function(index,mylist,t){
  j <- t
  while (j>0 && mylist[j]>index){
    mylist[j+1] <- mylist[j]
    j <- j-1
  }
  mylist[j+1] <- index # all elements greater than index have been moved one position to the right,
  # leaving space for index.
  mylist
}


#' remove an index from a vector of cluster ID
#'
#' This function is called when removing a cluster.
#'
#' @inheritParams ordered_insert
#'
#' @examples
#' ordered_remove(3,c(1:3,5:10),7)
#' ordered_remove(3,c(1:10),7)
#'
#'
#'# some quick examples:
#'
#' list0 <- c(1:7,0,0,0)
#' k <- 1
#' while (list0[k]==k){ k=k+ 1}; cm = k
#' while (list0[k]==k+1){ k =k+ 1}; ci = k+1
#' while (list0[k]==k+2) {k =k+ 1};  cj = k+2
#'
#' a1 <- ordered_remove(3,list0,7)
#' a2 <- ordered_remove(5,a1,6)
#' a3 <- ordered_insert(cm,a2,5)
#'
#' # end of examples.
#'
#' @export
ordered_remove <- function(index,mylist,t){
  for (j in 1:t){ # it stopped at t, which might be smaller than length of mylist.
    if (mylist[j] >= index){
      mylist[j] <- mylist[j+1]
    }
  }
  mylist # only t-1 matters now.
}


#' Find the next available cluster ID to use:
#'
#' @param mylist a vector of integer cluster IDs
#' @examples
#' ordered_next(c(1:3,5:10))
#' ordered_next(c(1:3,5:6,9:10))
#' @return an index to use next for an observation during complete Gibbs scan
#' @export
ordered_next <- function(mylist){
  j <- 1
  while (mylist[j]==j) {
    j <- j+1
  }
  j # identify the first gap.
}

#' Order a binary matrix by row
#'
#' The order is determined by magnitude of the rows as in a binary system
#'
#' @param mat A binary matrix to be ordered by row
#' @return a list of res - the ordered matrix from large to small, ord - the order
#' @export

order_mat_byrow <- function(mat){
  #mat <- Q_col_ordered[sample(1:M),]
  L <- ncol(mat)
  M <- nrow(mat)
  v <- log(as.matrix(2^{(L-1):0},ncol=1))
  diag_v <- matrix(0,nrow=L,ncol=L)
  diag(diag_v) <- v
  tmp_prod <- mat%*%diag_v
  permute_M_vec <- sapply(1:M,function(i){matrixStats::logSumExp(tmp_prod[i,mat[i,]!=0])})
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
  matrixStats::logSumExp(mat) # can remove stuff including and before ::.
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


equal_unit_vec0 <- function(v,k){
  if(k>length(v)){stop("==[rewind]length of vector v is shorter than k!==")}
  e <- rep(0,length(v))
  e[k] <- 1
  (sum(abs(e-v)) == 0)
}

# equal_unit_vec0(c(1,0,0,0,0,0),1)
# equal_unit_vec(c(1,0,0,0,0,0),1)



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
do_update_Q <- function(Q){
  M <- nrow(Q)
  L <- ncol(Q)
  test <- array(NA,c(M,L,3))
  ind_unit_Q <- colSums(Q)==1
  for (k in 1:M){
    for (l in 1:L){
      test[k,l,1] <- equal_unit_vec(Q[,l],k) # THIS MEANS A 1 IN E_M NEVER GETS UPDATED.
      test[k,l,2] <- (sum(Q[k,])==3) && (Q[k,l]==1)
      Q_cand      <- Q[,ind_unit_Q,drop=FALSE]
      test[k,l,3] <- (Q[k,l]==0) && sum(Q[,l])==1 &&
        (sum(sapply(1:ncol(Q_cand), function(c) equal_unit_vec(Q_cand[,c],which(Q[,l]==1))))==2)
    }
  }
  apply(test,c(1,2),function(v) !any(v))
}

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
      curr_prob <- exp(L1- matrixStats::logSumExp(c(L0,L1)))
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

#' update the Q matrix element by element
#'
#' This function updates Q matrix by Gibbs sampler (with option to do constarined
#' updates within identifiability constraint or not)
#' NB: - do we need M here? do we modify M after collapsing partner machines.
#'     - need to speed up.
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
          curr_prob <- exp(L1- matrixStats::logSumExp(c(L0,L1)))
          #print(curr_prob)

          #Q_old[k,l] <- metrop_flip(Q_old[k,l],curr_prob)
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
          curr_prob <- exp(L1- matrixStats::logSumExp(c(L0,L1)))
          #print(curr_prob)

          #Q_old[k,l] <- metrop_flip(Q_old[k,l],curr_prob)
          Q_old[k,l] <- rbinom(1,1,prob = curr_prob)
        }# end an update if needed.
      }
    }  #end iteration over elements.
  }
  Q_old
}

#' compute the full conditional probability of column Q_l given others
#'
#' Function to compute this log full conditional density, which will be used
#' in \link{update_Q}
#'
#' @param Y data
#' @param eta_star a matrix of machine usage indicators, rows for clusters, columns for M machines
#' @param Ql_enumerate l-th column of the Q matrix
#' @param p prevalence of machines; a vector of length identical to the columns of H.
#' @param thetal,psil True and false positive rates. Both are vectors of length L
#'
#' @return a vector of log conditional probability of column Q_l given other unknown parameters,
#' the dimension will be identical to \code{2^(length(Ql))}.s
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
#' This function updates Q matrix by Gibbs sampler (with option to do constarined
#' updates within identifiability constraint or not)
#' NB: - do we need M here? do we modify M after collapsing partner machines.
#'     - need to speed up.
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
  for (l in sample(1:L,replace=FALSE)){ # begin iteration over elements.
    L_enum  <- rep(0,2^M)
    for (j in 1:t){
      L_enum <- L_enum + log_pr_Qml_cond_enum(Y[z==mylist[j],l,drop=FALSE],
                                              H[j,,drop=FALSE],Ql_enum,theta[l],psi[l])
    }
    curr_prob <- exp(L_enum- matrixStats::logSumExp(L_enum))
    #print(round(curr_prob,2))
    Q_old[,l] <- Ql_enum[,sample(1:2^M,1,prob = curr_prob)]
  }
  Q_old
}

#' Restricted Gibbs scan to update assignment indicators in S.
#'
#' Function to perform Gibbs scan only for the observations associated
#' with either i or j (randomly chosen), but not including the cluster indicator
#' for i or j. This function does not change the cluster observation i or j currently is in,
#' because the restricted Gibbs scan only updates other observations in S by randomly
#' assigning them to the cluster with i or j.
#'
#' @param Y multivariate binary data (row for subjects, column for dimensions)
#' @param zsa cluster assignment indicators prior to updating (transition from in MH)
#' @param zsb cluster assignment indicators after updating (transition to in MH; same below.)
#'            both zsa and zsb will be set to zs when active is true, because we want the zs to be
#'            modified/updated.
#' @param cia,cja cluster indicator for subject i,j prior to updating
#' @param cib,cjb cluster indicator for subject i,j after updating
#' @param ni,nj the number of subects clustered with i (j) - including i (j).
#' @param i,j two distinct subejcts' indices
#' @param S the set of indices clustered with either i or j (including i and j)
#' @param ns size of S
#' @param b gamma in the dirichlet distribution in MFM.
#' @param active TRUE for performing a split update (both updating cluster indicators and computing the transition probabilities), FALSE otherwise (only for calculating
#' proposal density; used in q(Z|Z_merge) = q(Z|Z_launch)).
#' @param Q a Q matrix of dimension M by L
#' @param p a vector of machine (factor) prevalences of length M
#' @param theta a vector of true positive rates of length L
#' @param psi a vector of false positive rates of length L
#'
#'
#' @return \itemize{
#' \item \code{log_p} the log probability of transitioning from zsa to zsb (product of all subjects in S);
#' \item \code{ni,nj} the number of subjects assigned to cluster for i,j
#' \item \code{zsb} updated assignment indicators for subjects in S
#' }
#'
#' @export
restricted_gibbs <- function(Y,zsa,zsb,cia,cib,cja,cjb,ni,nj,i,j,S,ns,b,active,Q,p,theta,psi){
  # ??need to update the observations to be in sync with zs.
  # ?? i must be different from j?
  log_p <- 0 # variable for log transition probabilities.
  for (ks in 1:ns){ # beging iteration over observations in S.
    k <- S[ks]
    if (k!=i && k!=j){
      if (zsa[k]==cia){# <--mod # if the current subject to be updated beldongs to cluster of i:
        ni <- ni-1 # <--mod
      } else{
        nj <- nj-1 # <--mod
      }
      # compute probability to assign observation k to the cluster to which obs i belongs:
      Li <- log_marginal(rbind(Y[(zsa==cia)[-k],,drop=FALSE],Y[k,]),Q,p,theta,psi) - log_marginal(Y[(zsa==cia)[-k],,drop=FALSE],Q,p,theta,psi) # from restaurant process.
      Lj <- log_marginal(rbind(Y[(zsa==cja)[-k],,drop=FALSE],Y[k,]),Q,p,theta,psi) - log_marginal(Y[(zsa==cja)[-k],,drop=FALSE],Q,p,theta,psi)
      Pi <- exp(log(ni+b)+Li- matrixStats::logSumExp(c(log(ni+b)+Li,log(nj+b)+Lj))) # the (ni+b) also comes from restaurant process.

      # if we need to update the assignment indicators:
      if (active){
        zsb[k] <- ifelse(runif(1) < Pi, cib, cjb)
      }
      if (zsb[k]==cib){
        ni <- ni+1
        log_p <-log_p+log(Pi)
      } else{
        nj <- nj+1
        log_p <- log_p + log(1-Pi)
      }
    }
  } #end iteration over observations in S.
  return(list(log_p=log_p,ni=ni,nj=nj,zsb=zsb)) # zsb will be updated and passed to the next intermediate Gibbs scan. zsb indicates which observations are used for calculating log marginal likelihood.
}

#' Split-merge type update of clusterings
#'
#' This function uses restricted_gibbs to arrive at a launch state and then propose
#' a split-merge move. Based on the computed acceptance probability, this function accepts
#' or rejects the proposal arriving at an old or new state for the cluster indicators.
#'
#' @param Y multivariate binary data (row for subjects, column for dimensions)
#' @param z original cluster indicators
#' @param zs temporary variable used for split-merge assignments
#' @param S temporary set of indicies for restricted Gibbs
#' @param mylist the cluster ID list
#' @param N a vector of the numbers of subjects in the list of clusters
#' @param t total number of non-empty clusters
#' @param b gamma parameter in the mixture of finite mixture formulation
#' @param log_v coefficients in MFM (check Miller and Harrison 2016 JASA)
#' @param n_split the number of intermediate Gibbs scan to arrive at the launch split state (usually 5)
#' @param Q a Q matrix of dimension M by L
#' @param p a vector of machine (factor) prevalences of length M
#' @param theta a vector of true positive rates of length L
#' @param psi a vector of false positive rates of length L
#' @return Returns the values at the end of the current iteration \itemize{
#' \item \code{t} the number of clusters;
#' \item \code{z} the cluster assignment indicators for all subjects (taking values from \code{mylist});
#' \item \code{N} the number of subjects for each cluster; please refer to mylist to get the matching cluster IDs;
#' \item \code{mylist} the non-empty and empty cluster IDs.
#' }
#'
#' @export
split_merge <- function(Y,z,zs,S,mylist,N,t,b,log_v,n_split,Q,p,theta,psi){
  # for non-conjugate sampler, there needs to be n_merge
  n <- nrow(Y)
  # randomly choose a pair of indices:
  rand_pair <- sample(n,2,replace=FALSE)
  i <- rand_pair[1]
  j <- rand_pair[2]
  ci0 <- z[i]
  cj0 <- z[j] # original states.

  # set S[1],...,S[ns] to the indices of the observations in clusters ci0 and cj0:
  ns <- 0
  for (k in 1:n){
    if (z[k]==ci0 ||  z[k]==cj0){
      ns <- ns+1
      S[ns] <- k
    }
  }

  # find available cluster IDs for split and merge parameters:
  k <- 1
  while (mylist[k]==k){ k=k+ 1};     cm = k
  while (mylist[k]==k+1){ k =k+ 1};  ci = k+1
  while (mylist[k]==k+2) {k =k+ 1};  cj = k+2

  # merge state:
  # for (ks in 1:ns){
  #   #??
  # }

  # randomly choose the launch split state:
  zs[i] <- ci; ni <- 1
  zs[j] <- cj; nj <- 1

  for (ks in 1:ns){# start with a uniformly chosen split:
    k <- S[ks]
    if (k !=i && k!= j){
      if (runif(1)<0.5){
        zs[k] <- ci
        ni<- ni+1 # need adjoin function to collect observations in the same cluster.
      } else{
        zs[k] <- cj
        nj <- nj+1
      }
    }
  }

  for (rep in 1:n_split){#make several restricte Gibb scan moves:
    res_intermediate_GS <- restricted_gibbs(Y,zs,zs,ci,ci,cj,cj,ni,nj,i,j,S,ns,b,TRUE,Q,p,theta,psi) # make sure arguments are updated.
    log_p <- res_intermediate_GS$log_p
    ni <- res_intermediate_GS$ni
    nj <- res_intermediate_GS$nj
    zs <- res_intermediate_GS$zsb
  }

  # make MH proposal:
  if (ci0==cj0){#propose a split:
    # make one final sweep and compute the probability density:
    res_final_GS <- restricted_gibbs(Y,zs,zs,ci,ci,cj,cj,ni,nj,i,j,S,ns,b,TRUE,Q,p,theta,psi) # make sure arguments are updated.
    log_prop_ab <- res_intermediate_GS$log_p
    ni <- res_intermediate_GS$ni
    nj <- res_intermediate_GS$nj
    zs <- res_intermediate_GS$zsb

    # probability of transitioning from the merge state to original state:
    log_prop_ba <- 0 #log(1). Because given current split state, with probability 1 we get the merged state.

    # compute MH acceptance probability
    log_prior_b <- log_v[t+1]+lgamma(ni+b)+lgamma(nj+b)-2*lgamma(b)
    log_prior_a <- log_v[t] + lgamma(ns+b)-lgamma(b)
    log_lik_ratio <- log_marginal(Y[zs==ci,,drop=FALSE],Q,p,theta,psi)+log_marginal(Y[zs==cj,,drop=FALSE],Q,p,theta,psi)-
      log_marginal(Y[z==ci0,,drop=FALSE],Q,p,theta,psi)
    p_accept <- min(1,exp(log_prop_ba-log_prop_ab+log_prior_b-log_prior_a+log_lik_ratio))

    # accept or reject:
    if (runif(1)<p_accept){#accept split:
      for (ks in 1:ns){
        z[S[ks]] <- zs[S[ks]]
      }
      mylist <- ordered_remove(ci0,mylist,t)
      mylist <- ordered_insert(ci,mylist,t-1)
      mylist <- ordered_insert(cj,mylist,t)
      N[ci0] <- 0
      N[ci]  <- ni
      N[cj]  <- nj
      t <- t+1
    }
  } else{ # propose a merge:
    # probability of transitioning to merge state:
    log_prop_ab <- 0  # log(1)

    # compute probability density of going from split launch state to original state:
    res_imaginary_GS <- restricted_gibbs(Y,zs,z,ci,ci0,cj,cj0,ni,nj,i,j,S,ns,b,FALSE,Q,p,theta,psi)
    log_prop_ba <- res_imaginary_GS$log_p
    ni <- res_imaginary_GS$ni
    nj <- res_imaginary_GS$nj

    # compute acceptance probability:
    log_prior_b <- log_v[t-1]+lgamma(ns+b)-lgamma(b)
    log_prior_a <- log_v[t]+lgamma(ni+b)+lgamma(nj+b)-2*lgamma(b)
    log_lik_ratio <- log_marginal(Y[S,,drop=FALSE],Q,p,theta,psi) -
      log_marginal(Y[z==ci0,,drop=FALSE],Q,p,theta,psi)-log_marginal(Y[z==cj0,,drop=FALSE],Q,p,theta,psi) # computed for original (not launch state) and proposed states.
    p_accept <- min(1,log_prop_ba-log_prop_ab+log_prior_b-log_prior_a+log_lik_ratio)

    #accept or reject:
    if (runif(1)<p_accept){
      for (ks in 1:ns){
        z[S[ks]] <- cm
      }
      mylist <- ordered_remove(ci0,mylist,t)
      mylist <- ordered_remove(cj0,mylist,t-1)
      mylist <- ordered_insert(cm,mylist,t-2)
      N[cm] <- ns
      N[ci0] <- 0
      N[cj0] <- 0
      t <- t-1
    }
  }

  return(list(t=t,z=z,N=N,mylist=mylist)) # be extremely clear about what variables got updated.
}















# ------------------
# not used
#

# update_xi <- function(Y,Q,p,theta,psi){
#   n1 <- apply(Y,2,sum,na.rm=T)
#   n0 <- apply(1-Y,2,sum,na.rm=T)
#   #p_xi <- 1-exp(matrix(log(1-p),nrow=1)%*%Q) # length L.
#   tmp <- compute_marg_xi_prob(length(p))
#
#   mat <- rbind(n1*log(psi)+n0*log(1-psi),
#                n1*log(theta)+n0*log(1-theta))+
#     rbind(log(1-p_xi),
#           log(p_xi))
#   pr_gs_xi <- exp(mat[2,] - matrixStats::colLogSumExps(mat))
#   sapply(1:ncol(Y), function(i){rbinom(1,1,prob=pr_gs_xi[i])})
# }

# # not used:
# compute_marg_xi_prob <- function(M){
#   res <- 0
#   for (j in 1:M){
#     res <- res+choose(M,j)*(1-1/2^j)
#   }
#   res/2^M
# }
#
# # not used:
# update_xi <- function(Y,p,theta,psi){
#   n1 <- apply(Y,2,sum,na.rm=T)
#   n0 <- apply(1-Y,2,sum,na.rm=T)
#   tmp <- compute_marg_xi_prob(length(p))
#   p_xi <- rep(tmp,ncol(Y))
#
#   mat <- rbind(n1*log(psi)+n0*log(1-psi),
#                n1*log(theta)+n0*log(1-theta))+
#     rbind(log(1-p_xi),
#           log(p_xi))
#   pr_gs_xi <- exp(mat[2,] - matrixStats::colLogSumExps(mat))
#   sapply(1:ncol(Y), function(i){rbinom(1,1,prob=pr_gs_xi[i])})
# }
